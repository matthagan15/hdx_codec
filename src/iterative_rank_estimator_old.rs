use fxhash::{FxHashMap, FxHashSet};
use indexmap::IndexMap;
use mhgl::{HGraph, HyperGraph};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::io::Write;
use std::{collections::HashSet, fs::File, path::PathBuf, str::FromStr, time::Instant};

use crate::{
    code::Code,
    math::{finite_field::FiniteField, iterative_bfs_new::GroupBFS, polynomial::FFPolynomial},
    matrices::sparse_ffmatrix::{MemoryLayout, SparseFFMatrix},
    reed_solomon::ReedSolomon,
};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RankEstimatorConfig {
    quotient_poly: String,
    dim: usize,
    reed_solomon_degree: usize,
    output_dir: String,
    hgraph_file: PathBuf,
}

impl RankEstimatorConfig {
    pub fn new(
        quotient_poly: String,
        dim: usize,
        reed_solomon_degree: usize,
        output_dir: String,
        hgraph_file: PathBuf,
    ) -> Self {
        Self {
            quotient_poly,
            dim,
            reed_solomon_degree,
            output_dir,
            hgraph_file,
        }
    }
}
pub struct IterativeRankEstimator {
    hgraph: HGraph<u16, ()>,
    message_id_to_col_ix: IndexMap<u64, usize>,
    parity_check_to_row_ixs: IndexMap<u64, Vec<usize>>,
    pivotized_local_checks: HashSet<u32>,
    column_ix_to_pivot_row: BTreeMap<usize, usize>,
    parity_check_to_col_ixs: FxHashMap<u64, Vec<usize>>,
    pub parity_check_matrix: SparseFFMatrix,
    local_code: ReedSolomon,
    dim: usize,
}

impl IterativeRankEstimator {
    pub fn new(conf: RankEstimatorConfig) -> Self {
        let quotient =
            FFPolynomial::from_str(&conf.quotient_poly).expect("Could not parse polynomial.");
        log::trace!("Loaded polynomial: {:}", quotient);
        let local_code = ReedSolomon::new(quotient.field_mod, conf.reed_solomon_degree);
        local_code.print();
        let parity_check_matrix =
            SparseFFMatrix::new(0, 0, quotient.field_mod, MemoryLayout::RowMajor);
        let pathbuf = PathBuf::from(conf.output_dir);
        log::trace!(
            "Attempting to read hgraph from: {:}",
            &conf.hgraph_file.to_str().unwrap()
        );
        let hgraph =
            HGraph::from_file(&conf.hgraph_file).expect("Could not read hgraph from file.");

        // TODO: Need to assert that the local_code message_len is the same as what we
        // expect from the coset complex parameters.

        Self {
            hgraph,
            message_id_to_col_ix: IndexMap::new(),
            parity_check_to_row_ixs: IndexMap::new(),
            pivotized_local_checks: HashSet::new(),
            parity_check_matrix,
            column_ix_to_pivot_row: BTreeMap::new(),
            parity_check_to_col_ixs: FxHashMap::default(),
            local_code,
            dim: conf.dim,
        }
    }

    fn get_complete_border_checks(&self, local_check: u32) -> Vec<u64> {
        let mut ret = Vec::new();
        let maximal_faces = self.hgraph.maximal_edges_of_nodes([local_check]);
        for max_face in maximal_faces {
            let border_nodes: Vec<u32> = self
                .hgraph
                .query_edge(&max_face)
                .unwrap()
                .into_iter()
                .filter(|node| *node != local_check)
                .collect();
            let border_check = self
                .hgraph
                .find_id(border_nodes)
                .expect("Border nodes should have edge id.");
            let star: Vec<u32> = self
                .hgraph
                .link(&border_check)
                .into_iter()
                .map(|(_, rest)| rest[0])
                .collect();
            if star.into_iter().fold(true, |acc, e| {
                acc & self.pivotized_local_checks.contains(&e)
            }) {
                ret.push(border_check);
            }
        }
        ret
    }

    pub fn compute_rate(&mut self) {
        let mut counter = 0;
        let check_rate = 1000;
        let output_filename = String::from("/Users/matt/repos/qec/tmp/rate.txt");
        let mut output_file = File::create(output_filename).expect("could not create rate file.");
        let rate_start = Instant::now();
        let mut rate_upper_bound: Option<f64> = None;
        let nodes = self.hgraph.nodes();
        log::trace!("Found {:} nodes to check.", nodes.len());
        for node in nodes {
            log::trace!("Checking node: {node}");
            self.node_step(node);
            log::trace!(
                "rate: {:}",
                self.column_ix_to_pivot_row.len() as f64 / self.message_id_to_col_ix.len() as f64
            );
        }
        println!("{:}", "#".repeat(80));
        println!("{:}BORDER CHECKS{:}", " ".repeat(34), " ".repeat(34));
        println!("{:}", "#".repeat(80));
        self.pivotize_border_checks();

        log::trace!(
            "Parity check dims: {:}, {:}",
            self.parity_check_matrix.n_rows,
            self.parity_check_matrix.n_cols
        );
        log::trace!(
            "rate: {:}",
            self.column_ix_to_pivot_row.len() as f64 / self.message_id_to_col_ix.len() as f64
        )
    }

    fn old_compute_rate(&mut self) {

        // let mut counter = 0;
        // let check_rate = 1000;
        // let output_filename = String::from("/Users/matt/repos/qec/tmp/rate.txt");
        // let mut output_file = File::create(output_filename).expect("could not create rate file.");
        // let rate_start = Instant::now();
        // let mut rate_upper_bound: Option<f64> = None;
        // let nodes = self.hgraph.nodes();
        // if discovered.len() == 0 {
        //     let rate = 1.0
        //         - (self.column_ix_to_pivot_row.len() as f64
        //             / self.message_id_to_col_ix.len() as f64);
        //     log::debug!(
        //         "FINAL pivots, triangles, rate: {:}, {:}, {rate}",
        //         self.column_ix_to_pivot_row.len(),
        //         self.message_id_to_col_ix.len()
        //     );
        //     log::info!(
        //         "time taken (seconds): {:}",
        //         rate_start.elapsed().as_secs_f64()
        //     );
        //     let r = write!(
        //         output_file,
        //         "{:},{:}\n",
        //         self.column_ix_to_pivot_row.len(),
        //         self.message_id_to_col_ix.len()
        //     );
        //     if r.is_err() {
        //         log::error!("Could not write final output to disk!");
        //     }
        //     break;
        // }
        // local_checks_to_pivotize.insert(discovered[0]);
        // counter += 1;
        // if counter % check_rate == 0 {
        //     println!("{:}", "#".repeat(75));
        //     log::warn!(
        //         "Counter: {counter}, total_time: {:}",
        //         rate_start.elapsed().as_secs_f64()
        //     );
        //     let mut local_checks_not_ready = Vec::new();
        //     let mut total_time_interior = 0.0;
        //     let mut num_interior = 0;
        //     let mut time_border = 0.0;
        //     let mut num_border = 0;
        //     for local_check in local_checks_to_pivotize.drain() {
        //         let local_check_complete = self.is_local_view_complete(local_check);
        //         if local_check_complete {
        //             let instant_interior = Instant::now();
        //             let pivots_before = self.column_ix_to_pivot_row.len();
        //             let upper_bound = self.pivotize_interior_checks(local_check);
        //             if rate_upper_bound.is_none() {
        //                 rate_upper_bound = Some(upper_bound);
        //             }
        //             total_time_interior += instant_interior.elapsed().as_secs_f64();
        //             num_interior += self.column_ix_to_pivot_row.len() - pivots_before;

        //             self.pivotized_local_checks.insert(local_check);
        //             let border_checks = self.get_complete_border_checks(local_check);
        //             if border_checks.is_empty() == false {
        //                 for border_check in border_checks.into_iter() {
        //                     let temp = Instant::now();
        //                     self.pivotize_border_check(border_check);
        //                     time_border += temp.elapsed().as_secs_f64();
        //                     num_border += 1;
        //                 }
        //             }
        //         } else {
        //             local_checks_not_ready.push(local_check);
        //         }
        //     }
        //     for check in local_checks_not_ready.into_iter() {
        //         local_checks_to_pivotize.insert(check);
        //     }
        //     let rate = 1.0
        //         - (self.column_ix_to_pivot_row.len() as f64
        //             / self.message_id_to_col_ix.len() as f64);
        //     log::info!(
        //         "time interior: {:}, num_interior: {num_interior}, avg: {:}",
        //         total_time_interior,
        //         total_time_interior / (num_interior as f64)
        //     );
        //     log::info!(
        //         "time border: {:}, num_border: {:}, avg: {:}",
        //         time_border,
        //         num_border,
        //         time_border / (num_border as f64)
        //     );
        //     log::debug!(
        //         "pivots, triangles, rate: {:}, {:}, {rate}",
        //         self.column_ix_to_pivot_row.len(),
        //         self.message_id_to_col_ix.len()
        //     );
        //     if rate_upper_bound.is_some() {
        //         log::trace!("rate / upper_bound: {:}", rate / rate_upper_bound.unwrap());
        //     }
        //     let r = write!(
        //         output_file,
        //         "{:},{:}\n",
        //         self.column_ix_to_pivot_row.len(),
        //         self.message_id_to_col_ix.len()
        //     );
        //     if r.is_err() {
        //         log::error!("Could not write output to disk!");
        //     }
        // }
    }

    fn is_parity_check_complete(&self, parity_check: u64) -> bool {
        let maximal_edges = self.hgraph.maximal_edges(&parity_check);
        maximal_edges.len() == (self.parity_check_matrix.field_mod as usize)
    }

    fn is_local_view_complete(&self, local_check: u32) -> bool {
        let containing_edges = self.hgraph.containing_edges_of_nodes([local_check]);
        let mut is_complete = true;
        for containing_edge in containing_edges {
            let containing_edge_nodes = self.hgraph.query_edge(&containing_edge).unwrap();
            if containing_edge_nodes.len() == self.dim {
                let boundary_check: Vec<u32> = containing_edge_nodes
                    .iter()
                    .filter(|node| **node != local_check)
                    .cloned()
                    .collect();
                let boundary_check_id = self.hgraph.find_id(&boundary_check[..]).expect("Boundary check should still be a valid edge in the hgraph. If the traingle contains this local_check node then it must contain it's link within the triangle.");
                if self.is_parity_check_complete(boundary_check_id) == false {
                    is_complete = false;
                    break;
                }
            } else {
                if self.is_parity_check_complete(containing_edge) == false {
                    is_complete = false;
                    break;
                }
            }
        }
        is_complete
    }

    /// Gets the containing data edges for this parity check first.
    fn add_parity_check_to_matrix(&mut self, parity_check: u64) {
        let mut maximal_edges = self.hgraph.maximal_edges(&parity_check);
        maximal_edges.sort();
        if maximal_edges.len() != self.local_code.encoded_len() {
            return;
        }

        let row_indices_of_outputs = self
            .parity_check_to_row_ixs
            .get(&parity_check)
            .expect("Parity check needs rows before adding to matrix");
        for message_id in maximal_edges.iter() {
            let basis_vec: Vec<FiniteField> = maximal_edges
                .clone()
                .into_iter()
                .map(|temp_id: u64| {
                    if temp_id == *message_id {
                        FiniteField::new(1, self.parity_check_matrix.field_mod)
                    } else {
                        FiniteField::new(0, self.parity_check_matrix.field_mod)
                    }
                })
                .collect();
            let parity_check_output = self.local_code.parity_check(&basis_vec);
            if row_indices_of_outputs.len() != parity_check_output.len() {
                panic!("Number of rows assigned to parity check does not match output length of parity check.")
            }
            for ix in 0..row_indices_of_outputs.len() {
                self.parity_check_matrix.insert(
                    row_indices_of_outputs[ix] as usize,
                    *self.message_id_to_col_ix.get(message_id).unwrap() as usize,
                    parity_check_output[ix].0,
                );
            }
        }
    }

    fn parity_check_handler(&mut self, parity_check: &u64) {
        if self.parity_check_to_row_ixs.contains_key(parity_check) {
            return;
        }
        let current_ix_start = if self.parity_check_to_row_ixs.is_empty() {
            0
        } else {
            self.parity_check_to_row_ixs
                .get_index(self.parity_check_to_row_ixs.len() - 1)
                .unwrap()
                .1
                .last()
                .unwrap()
                + 1
        };
        let new_row_indices: Vec<usize> =
            (current_ix_start..(current_ix_start + self.local_code.parity_check_len())).collect();
        self.parity_check_to_row_ixs
            .insert(*parity_check, new_row_indices);

        self.add_parity_check_to_matrix(*parity_check);
    }

    pub fn step(&mut self) -> Vec<u32> {
        // let triangle = self.bfs.step();
        // if triangle.is_empty() {
        //     return Vec::new();
        // }
        // let hg = &self.hgraph;
        // let e1 = hg
        //     .find_id([triangle[0], triangle[1]])
        //     .expect("should have been allocated in bfs.");
        // let e2 = hg
        //     .find_id([triangle[0], triangle[2]])
        //     .expect("should have been allocated in bfs.");
        // let e3 = hg
        //     .find_id([triangle[1], triangle[2]])
        //     .expect("should have been allocated in bfs.");
        // let triangle_id = hg
        //     .find_id(&triangle[0..3])
        //     .expect("should have been allocated in bfs.");
        // let message_ix = if self.message_id_to_col_ix.is_empty() {
        //     0
        // } else {
        //     self.message_id_to_col_ix.last().unwrap().1 + 1
        // };
        // self.message_id_to_col_ix.insert(triangle_id, message_ix);
        // self.parity_check_handler(&e1);
        // self.parity_check_handler(&e2);
        // self.parity_check_handler(&e3);
        // triangle
        todo!()
    }

    fn node_step(&mut self, node: u32) {
        // first need to make sure that each row in the nodes view
        // has matrix rows filled out.
        let type_ix = self.hgraph.get_node(&node).unwrap();
        if *type_ix != 0 {
            return;
        }
        let star = self.hgraph.containing_edges_of_nodes([node]);
        let mut interior_checks = Vec::new();
        let mut message_ixs = Vec::new();
        let mut border_checks = HashMap::new();
        for edge in star.iter() {
            let edge_nodes = self.hgraph.query_edge(&edge).unwrap();
            if edge_nodes.len() == 2 {
                interior_checks.push(*edge);
            } else if edge_nodes.len() == 3 {
                // make sure the message_id has a column index
                let message_ix = if self.message_id_to_col_ix.contains_key(edge) == false {
                    let message_ix = if self.message_id_to_col_ix.is_empty() {
                        0
                    } else {
                        self.message_id_to_col_ix.last().unwrap().1 + 1
                    };
                    self.message_id_to_col_ix.insert(*edge, message_ix);
                    message_ixs.push(message_ix);
                    message_ix
                } else {
                    let ix = self.message_id_to_col_ix.get(edge).unwrap();
                    message_ixs.push(*ix);
                    *ix
                };
                // get the border check
                let border_nodes: Vec<u32> = edge_nodes
                    .into_iter()
                    .filter(|node| *self.hgraph.get_node(&node).unwrap() != 0)
                    .collect();
                if border_nodes.len() != 2 {
                    panic!("Border nodes should have length 2");
                }
                let border_check = self.hgraph.find_id(&border_nodes[..]).unwrap();
                border_checks.insert(border_check, false);
                let e = self
                    .parity_check_to_col_ixs
                    .entry(border_check)
                    .or_default();
                e.push(message_ix);
                if e.len() == self.local_code.message_len() {
                    let border_check_complete = border_checks.get_mut(&border_check).unwrap();
                    *border_check_complete = true;
                }
            }
        }

        if message_ixs.len() != 2 * self.parity_check_matrix.field_mod.pow(2) as usize {
            println!("local check does not have complete triangles.");
            return;
        } else {
            log::trace!("Local node has enough triangles.");
        }

        // Now its time to pivotize the interior checks
        for interior_check in interior_checks.iter() {
            log::trace!("interior_check getting handled: {:}", interior_check);
            self.parity_check_handler(interior_check);
        }
        log::trace!("interior checks handled.");
        for border_check in border_checks
            .iter()
            .filter_map(|(k, v)| if *v { Some(k) } else { None })
        {
            log::trace!("border check getting handled: {:}", border_check);
            self.parity_check_handler(border_check);
        }
        log::trace!("border checks handled.");
        let mut interior_ixs = Vec::new();
        for check in interior_checks {
            let mut ixs = self.parity_check_to_row_ixs.get(&check).unwrap().clone();
            interior_ixs.append(&mut ixs);
        }
        for ix in interior_ixs.iter() {
            if let Some(col) = self
                .parity_check_matrix
                .pivotize_row_within_range(*ix, &interior_ixs[..])
            {
                log::trace!("Adding Pivot: {:}, {:}", col, *ix);
                self.column_ix_to_pivot_row.insert(col, *ix);
            } else {
                log::trace!("Found nonzero row {:}", *ix);
            }
        }
        // let mut border_ixs = Vec::new();
        for check in border_checks
            .iter()
            .filter_map(|(k, v)| if *v { Some(k) } else { None })
        {
            let border_ixs = self.parity_check_to_row_ixs.get(check).unwrap();
            let col_ixs = self.parity_check_to_col_ixs.get(check).unwrap();
            for col_ix in col_ixs.iter() {
                if let Some(pivot_row) = self.column_ix_to_pivot_row.get(col_ix) {
                    self.parity_check_matrix
                        .eliminate_rows((*pivot_row, *col_ix), border_ixs.clone());
                }
            }
        }
    }

    fn pivotize_border_checks(&mut self) {
        for (border_check, message_ixs) in self.parity_check_to_col_ixs.clone() {
            self.parity_check_handler(&border_check);
            for border_ix in self.parity_check_to_row_ixs.get(&border_check).unwrap() {
                let col = self.parity_check_matrix.pivotize_row(*border_ix);
                if let Some(col) = col {
                    log::trace!("Adding border pivot: {:}, {:}", *border_ix, col);
                    self.column_ix_to_pivot_row.insert(col, *border_ix);
                }
            }
        }
    }

    // pub fn relative_rate_upper_bound(&mut self) -> f64 {
    //     let mut first_vertex = None;
    //     for step in 0..usize::MAX {
    //         let t = self.step();
    //         if first_vertex.is_none() {
    //             first_vertex = Some(t[0]);
    //         }
    //         if self.is_local_view_complete(first_vertex.unwrap()) {
    //             println!("First check is done! step number {:}", step);
    //             return self.pivotize_interior_checks(first_vertex.unwrap());
    //         }
    //     }
    //     panic!("Reached the end of the complex without finding a completed local check?")
    // }

    fn compute_interior_and_border_ixs(&self, local_check: u32) -> (Vec<usize>, Vec<usize>) {
        todo!()
    }
    /// Computes a *lower* bound on the (normalized) rate of the error correcting code
    /// from the resulting complex, given a completed local_check.
    pub fn pivotize_interior_checks(
        &mut self,
        local_check: u32,
        // containing_edges: Vec<u64>,
        // interior_ixs: Vec<usize>,
        // border_ixs: Vec<usize>,
    ) -> f64 {
        let mut interior_checks = Vec::new();
        let mut border_checks = Vec::new();
        let mut message_ids = Vec::new();
        let containing_edges = self.hgraph.containing_edges_of_nodes([local_check]);
        for edge in containing_edges {
            let nodes = self.hgraph.query_edge(&edge).unwrap();
            if nodes.len() == 2 {
                interior_checks.push(edge);
            } else if nodes.len() == 3 {
                message_ids.push(edge);
                let border_nodes: Vec<u32> = nodes
                    .into_iter()
                    .filter(|node| *node != local_check)
                    .collect();
                let border_id = self.hgraph.find_id(border_nodes).unwrap();
                border_checks.push(border_id);
            }
        }
        let mut border_ixs = Vec::new();
        for border_check in border_checks.iter() {
            border_ixs.append(
                &mut self
                    .parity_check_to_row_ixs
                    .get(border_check)
                    .unwrap()
                    .clone(),
            );
        }
        let mut interior_ixs = Vec::new();
        for interior_check in interior_checks.iter() {
            interior_ixs.append(
                &mut self
                    .parity_check_to_row_ixs
                    .get(interior_check)
                    .unwrap()
                    .clone(),
            );
        }
        let mut total_ixs = interior_ixs.clone();
        total_ixs.append(&mut border_ixs.clone());

        // iterate through each triangle, attempt to find a pivot on the interior checks.
        let mut num_pivots_found = 0;
        let mut pivots = HashSet::new();
        for message_id in message_ids.iter() {
            let message_ix = self.message_id_to_col_ix.get(message_id).unwrap();
            let possible_ixs: HashSet<_> = interior_ixs
                .iter()
                .filter(|x| pivots.contains(*x) == false)
                .cloned()
                .collect();
            if let Some(pivot) = self
                .parity_check_matrix
                .find_nonzero_entry_among_rows(*message_ix, possible_ixs.into_iter().collect())
            {
                // log::trace!("pivot! {:}", pivot);
                self.parity_check_matrix
                    .eliminate_rows((pivot, *message_ix), total_ixs.clone());
                pivots.insert(pivot);
                self.column_ix_to_pivot_row.insert(*message_ix, pivot);
                num_pivots_found += 1;
            }
        }
        let ret = 1.0 - (num_pivots_found as f64) / (message_ids.len() as f64);
        ret
    }

    fn print_border_ixs(&self, border_ixs: Vec<usize>) {
        let mut cols = HashSet::new();
        for row_ix in border_ixs.iter() {
            let row = self.parity_check_matrix.row(*row_ix);
            for (c, e) in row.to_vec() {
                cols.insert(c);
            }
        }
        let mut mat = SparseFFMatrix::new(
            border_ixs.len(),
            cols.len(),
            self.parity_check_matrix.field_mod,
            MemoryLayout::RowMajor,
        );
        for row_ix in border_ixs.iter() {
            let row = self.parity_check_matrix.row(*row_ix);
            for (c, e) in row.to_vec() {
                mat.insert(*row_ix, c, e);
            }
        }
        mat.dense_print();
    }

    fn pivotize_border_check(&mut self, border_check: u64) {
        let mut rows = self
            .parity_check_to_row_ixs
            .get(&border_check)
            .expect("border check does not have rows?")
            .clone();
        rows.sort();
        for row in rows {
            let r = self.parity_check_matrix.row(row);
            if let Some((ix, _)) = r.first_nonzero() {
                self.parity_check_matrix.eliminate_all_rows((row, ix));

                self.column_ix_to_pivot_row.insert(ix, row);
            }
        }
    }

    fn print_sub_matrix(
        &self,
        border_ixs: Vec<usize>,
        interior_ixs: Vec<usize>,
        message_ids: Vec<u64>,
    ) {
        let mut entries = Vec::new();
        let mut new_row_ix: usize = 0;
        for row_ix in interior_ixs.iter() {
            let mut new_col_ix: usize = 0;
            for message_id in &message_ids {
                let message_ix = self.message_id_to_col_ix.get(message_id).unwrap();
                let entry = self.parity_check_matrix.get(*row_ix, *message_ix);
                entries.push((new_row_ix, new_col_ix, entry.0));
                new_col_ix += 1;
            }
            new_row_ix += 1;
        }
        let mut new_col_ix: usize = 0;
        for _ in &message_ids {
            entries.push((
                new_row_ix,
                new_col_ix,
                self.parity_check_matrix.field_mod - 1,
            ));
            new_col_ix += 1;
        }
        new_row_ix += 1;
        for row_ix in border_ixs {
            let mut new_col_ix: usize = 0;
            for message_id in &message_ids {
                let message_ix = self.message_id_to_col_ix.get(message_id).unwrap();
                let entry = self.parity_check_matrix.get(row_ix, *message_ix);
                entries.push((new_row_ix, new_col_ix, entry.0));
                new_col_ix += 1;
            }
            new_row_ix += 1;
        }
        let sub_matrix = SparseFFMatrix::new_with_entries(
            new_row_ix - 1,
            message_ids.len(),
            self.parity_check_matrix.field_mod,
            MemoryLayout::RowMajor,
            entries,
        );
        println!("dense");
        sub_matrix.dense_print();
    }
    fn print_hgraph(&self) {
        println!("{:?}", self.hgraph);
    }
}
