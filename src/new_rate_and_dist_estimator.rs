use std::{
    collections::{HashMap, HashSet},
    path::PathBuf,
    time::Instant,
};

use mhgl::{HGraph, HyperGraph};
use serde::{ser::SerializeStruct, Deserialize, Serialize};

use crate::{
    math::{
        coset_complex_bfs::{bfs, find_complete_nodes, NodeData},
        finite_field::{FFRep, FiniteField},
        polynomial::FFPolynomial,
    },
    matrices::{
        mat_trait::RankMatrix,
        parallel_matrix::ParallelFFMatrix,
        sparse_ffmatrix::{MemoryLayout, SparseFFMatrix},
        sparse_sparse_ffmatrix::SparseSparseFFMatrix,
        sparse_vec::SparseVector,
    },
    reed_solomon::ReedSolomon,
};

const HGRAPH_CACHE_FILE_NAME: &str = "hgraph.hg";
const PARALLEL_MATRIX_CACHE: &str = "parallel_matrix_cache";
const BORDER_CACHE: &str = "border_cache";
const RATE_AND_DIST_CACHE: &str = "rate_and_dist_cache";
const INTERIOR_CACHE: &str = "interior_cache";

fn entries_to_rows(
    entries: Vec<(usize, usize, FFRep)>,
    field_mod: FFRep,
) -> Vec<(usize, SparseVector)> {
    let mut hm = HashMap::new();
    for (row_ix, col_ix, entry) in entries {
        let e = hm.entry(row_ix).or_insert(SparseVector::new_empty());
        e.add_entry(col_ix, entry, field_mod);
    }
    hm.into_iter().collect()
}

#[derive(Debug)]
pub struct NewRateAndDistConfig {
    directory: PathBuf,
    quotient: FFPolynomial,
    dim: usize,

    data_id_to_col_ix: HashMap<u64, usize>,
    col_ix_to_data_id: HashMap<usize, u64>,
    next_col_ix: usize,
    /// Stored in increasing order.
    completed_nodes: Vec<u32>,
    /// Stored in decreasing order so that way nodes can just be popped
    /// off to process.
    remaining_nodes: Vec<u32>,
    /// How many nodes to process before computing an (n,k,d) tuple.
    num_nodes_per_checkpoint: usize,

    /// How many times to compute the distance of the code
    num_distance_calcs: usize,

    /// Map from a checkpoint node to the value of n, k, and d
    /// after processing all nodes up to and including the node.
    /// a d is only Some() if it is computed at that checkpoint.
    code_checkpoints: Vec<(u32, (usize, usize, Option<usize>))>,

    row_ix_to_check_id: HashMap<usize, u64>,
    check_id_to_row_ixs: HashMap<u64, Vec<usize>>,
    next_row_ix: usize,

    pivot_row_to_col_ix: HashMap<usize, usize>,
    pivot_col_to_row_ix: HashMap<usize, usize>,
    matrix: SparseSparseFFMatrix,
}

impl NewRateAndDistConfig {
    pub fn new(
        quotient: FFPolynomial,
        dim: usize,
        truncation: Option<usize>,
        hgraph_log_rate: Option<usize>,
        num_distance_calcs: Option<usize>,
        directory: PathBuf,
        num_checkpoints: usize,
    ) -> Self {
        if directory.is_dir() == false {
            panic!("Provided directory is not a directory.")
        }
        let mut hg_output_file = directory.clone();
        hg_output_file.push(HGRAPH_CACHE_FILE_NAME);
        let hg = bfs(
            quotient.clone(),
            dim,
            truncation,
            Some(hg_output_file),
            hgraph_log_rate,
        );
        let mut remaining_nodes = find_complete_nodes(&hg, &quotient, dim);
        remaining_nodes.sort();
        remaining_nodes.reverse();
        let num_nodes_per_checkpoint = remaining_nodes.len() / num_checkpoints;
        let field_mod = quotient.field_mod;
        // let num_threads = 1;
        NewRateAndDistConfig {
            directory,
            quotient,
            dim,
            data_id_to_col_ix: HashMap::new(),
            col_ix_to_data_id: HashMap::new(),
            next_col_ix: 0,
            completed_nodes: Vec::new(),
            remaining_nodes,
            num_nodes_per_checkpoint,
            num_distance_calcs: num_distance_calcs.unwrap_or(1),
            code_checkpoints: Vec::new(),
            row_ix_to_check_id: HashMap::new(),
            check_id_to_row_ixs: HashMap::new(),
            next_row_ix: 0,
            pivot_col_to_row_ix: HashMap::new(),
            pivot_row_to_col_ix: HashMap::new(),
            matrix: SparseSparseFFMatrix::new(field_mod),
        }
    }

    fn process_node_batch(
        &mut self,
        hg: &HGraph<NodeData, ()>,
        local_code: &ReedSolomon,
        compute_distance: bool,
    ) {
        let mut next_node_batch = Vec::new();
        while self.remaining_nodes.is_empty() == false
            && next_node_batch.len() < self.num_nodes_per_checkpoint
        {
            next_node_batch.push(self.remaining_nodes.pop().unwrap());
        }
        let start_interior = Instant::now();
        for node in next_node_batch.iter() {
            self.process_node(hg, local_code, *node)
        }

        let elapsed_interior = start_interior.elapsed().as_secs_f64();

        let last_node = *next_node_batch.last().unwrap();
        self.completed_nodes.append(&mut next_node_batch);

        let n = self.col_ix_to_data_id.len();
        let tot_num_pivots = self.pivot_row_to_col_ix.len();
        let k = n - tot_num_pivots;
        if compute_distance {
            let start_distance = Instant::now();
            let mut col_weights = Vec::new();
            for col_ix in self.col_ix_to_data_id.keys() {
                if self.pivot_col_to_row_ix.contains_key(col_ix) {
                    continue;
                }
                let col_weight = self.matrix.get_col_weight(*col_ix);
                col_weights.push(col_weight);
            }
            let d = col_weights
                .iter()
                // .filter(|w| **w > 0)
                .min()
                .cloned()
                .unwrap_or(0);
            let elapsed_distance = start_distance.elapsed().as_secs_f64();

            println!(
                "[n, k, d] = [{:}, {:}, {:}], k/n = {:.6}, d/n = {:.6}",
                n,
                k,
                d,
                (k as f64 / n as f64),
                (d as f64) / (n as f64)
            );
            let tot_time = start_interior.elapsed().as_secs_f64();
            println!("tot time for batch: {:}", tot_time);
            println!("% time spent interior: {:}", elapsed_interior / tot_time);
            println!("% time spent distance: {:}", elapsed_distance / tot_time);

            self.code_checkpoints.push((last_node, (n, k, Some(d))));
        } else {
            let tot_time = start_interior.elapsed().as_secs_f64();
            println!(
                "[n, k, d] = [{:}, {:}, ??], k/n = {:.6}, d/n = ??",
                n,
                k,
                (k as f64 / n as f64),
            );
            println!("tot time for batch: {:}", tot_time);
            println!("% time spent interior: {:}", elapsed_interior / tot_time);
            self.code_checkpoints.push((last_node, (n, k, None)));
        }
    }

    /// Returns pivots added to the matrix.
    fn process_node(&mut self, hg: &HGraph<NodeData, ()>, local_code: &ReedSolomon, node: u32) {
        let containing_edges = hg.containing_edges_of_nodes([node]);
        let data_ids = containing_edges
            .iter()
            .filter(|id| {
                if let Some(nodes) = hg.query_edge(id) {
                    nodes.len() == self.dim
                } else {
                    false
                }
            })
            .cloned()
            .collect::<Vec<u64>>();
        let mut checks = Vec::new();
        let mut interior_checks = containing_edges
            .iter()
            .filter(|id| {
                if let Some(nodes) = hg.query_edge(id) {
                    nodes.len() == self.dim - 1
                } else {
                    false
                }
            })
            .cloned()
            .collect::<Vec<u64>>();
        checks.append(&mut interior_checks);
        let mut border_checks = data_ids
            .iter()
            .filter_map(|id| {
                let data_nodes = hg.query_edge(id).unwrap();
                let border_nodes = data_nodes
                    .into_iter()
                    .filter(|node| *hg.get_node(node).unwrap() != 0)
                    .collect::<Vec<u32>>();
                hg.find_id(border_nodes)
            })
            .filter(|id| {
                let border_check_view = hg.link(id);
                if border_check_view.len() != self.quotient.field_mod as usize {
                    return false;
                }
                let max_node_of_border = border_check_view
                    .iter()
                    .max_by_key(|(_id, nodes)| nodes[0])
                    .unwrap()
                    .1[0];
                max_node_of_border == node
            })
            .collect::<Vec<u64>>();
        checks.append(&mut border_checks);
        for data_id in data_ids.iter() {
            if self.data_id_to_col_ix.contains_key(data_id) {
                continue;
            }
            self.data_id_to_col_ix.insert(*data_id, self.next_col_ix);
            self.col_ix_to_data_id.insert(self.next_col_ix, *data_id);
            self.next_col_ix += 1;
        }

        let local_parity_check = local_code.parity_check_matrix();
        let mut row_ixs_added = HashSet::new();
        for check in checks {
            let mut data_ids_visible = hg.maximal_edges(&check);
            data_ids_visible.sort();
            for data_id in data_ids_visible.iter() {
                if self.data_id_to_col_ix.contains_key(data_id) {
                    continue;
                }
                self.data_id_to_col_ix.insert(*data_id, self.next_col_ix);
                self.col_ix_to_data_id.insert(self.next_col_ix, *data_id);
                self.next_col_ix += 1;
            }
            let mut new_entries = Vec::new();
            assert_eq!(local_parity_check.n_cols, data_ids_visible.len());
            for ix in 0..data_ids_visible.len() {
                let col = local_parity_check.get_col(ix);
                for offset in 0..col.len() {
                    let new_row_ix = self.next_row_ix + offset;
                    row_ixs_added.insert(new_row_ix);
                    new_entries.push((
                        new_row_ix,
                        self.data_id_to_col_ix[&data_ids_visible[ix]],
                        col[offset],
                    ));
                }
            }
            self.next_row_ix += local_parity_check.n_rows;
            let mut new_rows = HashMap::new();
            for (row, col, entry) in new_entries {
                let vec = new_rows.entry(row).or_insert(SparseVector::new_empty());
                vec.add_entry(col, entry, self.matrix.field_mod);
            }
            for (row_ix, row) in new_rows.into_iter() {
                if let Some((pivot_row, pivot_col)) = self.matrix.add_row_and_pivot(row_ix, row) {
                    self.pivot_col_to_row_ix.insert(pivot_col, pivot_row);
                    self.pivot_row_to_col_ix.insert(pivot_row, pivot_col);
                }
            }
        }
    }

    /// returns time taken
    pub fn run(&mut self) -> f64 {
        let mut hg_file = self.directory.clone();
        hg_file.push(HGRAPH_CACHE_FILE_NAME);
        let hg = HGraph::from_file(&hg_file).unwrap();
        let local_code = ReedSolomon::new(
            self.quotient.field_mod,
            (self.quotient.field_mod - 1) as usize,
        );
        let pc_mat = local_code.parity_check_matrix();
        println!("pc_mat:\n{:}", pc_mat);
        let num_batches = self
            .remaining_nodes
            .len()
            .div_ceil(self.num_nodes_per_checkpoint);
        let num_batches_per_dist = num_batches / self.num_distance_calcs;
        let mut num_batches_done = 0;
        let run_start = Instant::now();
        while self.remaining_nodes.len() > 0 {
            num_batches_done += 1;
            let batch_start = Instant::now();
            if num_batches_done % num_batches_per_dist == 0 {
                println!("{:}", "=".repeat(77));
                self.process_node_batch(&hg, &local_code, true);
            } else {
                println!("{:}", "-".repeat(77));
                self.process_node_batch(&hg, &local_code, false);
            }
            let batch_time_taken = batch_start.elapsed().as_secs_f64();
            let num_batches_remaining = num_batches - num_batches_done;
            let time_remaining = num_batches_remaining as f64 * batch_time_taken;
            println!("Estimated time remaining: {:}", time_remaining);
        }
        println!("Complete!");
        let time_taken = run_start.elapsed().as_secs_f64();
        println!("Total time taken: {:}", time_taken);
        println!(
            "Time per node: {:}",
            time_taken / self.completed_nodes.len() as f64
        );
        time_taken
        // let mut pivots = Vec::new();
        // for (row_ix, col_ix) in self.border_manager.pivot_row_to_col_ix.iter() {
        //     pivots.push((*row_ix, *col_ix));
        // }
        // for (row_ix, col_ix) in self.interior_manager.pivot_row_to_col_ix.iter() {
        //     pivots.push((*row_ix, *col_ix));
        // }
        // println!("Checking pivots");
        // for pivot in pivots {
        //     let good_interior = self.interior_manager.matrix.assert_pivot(pivot).is_ok();
        //     let good_border = self.border_manager.matrix.assert_pivot(pivot).is_ok();
        //     assert!(good_interior);
        //     assert!(good_border);
        // }
    }
}

#[cfg(test)]
mod tests {
    use std::{path::PathBuf, str::FromStr, time::Instant};

    use simple_logger::SimpleLogger;

    use crate::{
        math::polynomial::FFPolynomial, new_rate_and_dist_estimator::NewRateAndDistConfig,
        rate_and_dist_estimator::RateAndDistConfig,
    };

    #[test]
    fn single_node_example_sparse_sparse() {
        let q = FFPolynomial::from_str("x^2 + 2x + 2 % 3").unwrap();
        let dim = 3;
        let dir = PathBuf::from_str("/Users/matt/repos/qec/tmp/single_node_test").unwrap();
        let t = Some(718000);
        let l = Some(100);
        let l = None;
        let c = Some(1);
        let checks = 10;
        // let _logger = SimpleLogger::new().init().unwrap();
        let mut rate_estimator =
            NewRateAndDistConfig::new(q.clone(), dim, t, l, c, dir.clone(), checks);

        let new_taken = rate_estimator.run();

        let mut old_rate_estimator = RateAndDistConfig::new(q, dim, t, l, c, dir, checks);
        let old_taken = old_rate_estimator.run();

        println!("new / old: {:}", new_taken / old_taken);

        let (mut interior, mut border) = old_rate_estimator.quit();
        interior.append(&mut border);
        let new_mat = rate_estimator.matrix.to_sparse();
        println!("new == old? {:}", new_mat == interior);
        // println!("mat: {:}", new_mat.to_dense());
        // println!("interior: {:}", i.to_dense());
        // println!("border: {:}", b.to_dense());
        // dbg!(&rate_estimator);
    }
}
