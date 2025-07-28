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
struct InteriorManager {
    matrix: SparseFFMatrix,
    col_ix_to_pivot_row: HashMap<usize, usize>,
    pivot_row_to_col_ix: HashMap<usize, usize>,
}

impl InteriorManager {
    pub fn new(field_mod: FFRep) -> Self {
        InteriorManager {
            matrix: SparseFFMatrix::new(0, 0, field_mod, MemoryLayout::RowMajor),
            col_ix_to_pivot_row: HashMap::new(),
            pivot_row_to_col_ix: HashMap::new(),
        }
    }

    pub fn add_entries(&mut self, entries: Vec<(usize, usize, FFRep)>) {
        self.matrix.insert_entries(entries);
    }

    pub fn reduce_external_row(&self, row: &mut SparseVector) {
        let mut next_reducable_col = None;
        for (col_ix, _entry) in row.0.iter() {
            if self.col_ix_to_pivot_row.contains_key(col_ix) {
                next_reducable_col = Some(*col_ix);
                break;
            }
        }
        'outer_loop: while next_reducable_col.is_some() {
            let pivot_row = self.matrix.get_row(
                *self
                    .col_ix_to_pivot_row
                    .get(next_reducable_col.as_ref().unwrap())
                    .unwrap(),
            );
            let mut scalar = FiniteField::new(
                row.query(next_reducable_col.as_ref().unwrap()),
                self.matrix.field_mod,
            );
            scalar = -1 * scalar;
            row.add_scaled_row_to_self(scalar, &pivot_row);
            for (col_ix, _entry) in row.0.iter() {
                if self.col_ix_to_pivot_row.contains_key(col_ix) {
                    next_reducable_col = Some(*col_ix);
                    continue 'outer_loop;
                }
            }
            next_reducable_col = None;
        }
    }

    pub fn add_pivots(&mut self, pivots: Vec<(usize, usize)>) {
        for (row_ix, col_ix) in pivots {
            self.pivot_row_to_col_ix.insert(row_ix, col_ix);
            self.col_ix_to_pivot_row.insert(col_ix, row_ix);
        }
    }
}

impl Serialize for InteriorManager {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let col_ix_to_pivot_row: Vec<(usize, usize)> = self
            .col_ix_to_pivot_row
            .iter()
            .map(|(k, v)| (*k, *v))
            .collect();
        let pivot_row_to_col_ix: Vec<_> = self
            .pivot_row_to_col_ix
            .iter()
            .map(|(k, v)| (*k, *v))
            .collect();
        let mut s = serializer.serialize_struct("InteriorManager", 3)?;
        s.serialize_field("matrix", &self.matrix)?;
        s.serialize_field("col_ix_to_pivot_row", &col_ix_to_pivot_row)?;
        s.serialize_field("pivot_row_to_col_ix", &pivot_row_to_col_ix)?;
        s.end()
    }
}

impl<'de> Deserialize<'de> for InteriorManager {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let mut data = <serde_json::Value>::deserialize(deserializer)?;
        let matrix: Result<SparseFFMatrix, _> = serde_json::from_value(data["matrix"].take());
        let col_ix_to_pivot_row: Result<Vec<(usize, usize)>, _> =
            serde_json::from_value(data["col_ix_to_pivot_row"].take());
        let pivot_row_to_col_ix: Result<Vec<(usize, usize)>, _> =
            serde_json::from_value(data["pivot_row_to_col_ix"].take());
        Ok(InteriorManager {
            matrix: matrix.unwrap(),
            col_ix_to_pivot_row: col_ix_to_pivot_row.unwrap().into_iter().collect(),
            pivot_row_to_col_ix: pivot_row_to_col_ix.unwrap().into_iter().collect(),
        })
    }
}

#[derive(Debug, Serialize, Deserialize)]
struct BorderCache {
    col_ix_to_pivot_row: Vec<(usize, usize)>,
    pivot_row_to_col_ix: Vec<(usize, usize)>,
}

#[derive(Debug)]
struct BorderManager {
    col_ix_to_pivot_row: HashMap<usize, usize>,
    pivot_row_to_col_ix: HashMap<usize, usize>,
    matrix: ParallelFFMatrix,
}

impl BorderManager {
    pub fn new(field_mod: FFRep, num_threads: usize) -> Self {
        Self {
            col_ix_to_pivot_row: HashMap::new(),
            pivot_row_to_col_ix: HashMap::new(),
            matrix: ParallelFFMatrix::new(
                (0..num_threads)
                    .map(|_ix| SparseFFMatrix::new(0, 0, field_mod, MemoryLayout::RowMajor))
                    .collect(),
            ),
        }
    }

    pub fn add_entries(&mut self, entries: Vec<(usize, usize, FFRep)>) {
        self.matrix.add_entries(entries);
    }

    pub fn add_row_and_pivot(
        &mut self,
        row_ix: usize,
        row: SparseVector,
    ) -> Option<(usize, usize)> {
        self.matrix.add_row_and_pivot(row_ix, row)
    }

    pub fn save_matrix_to_disk(&self, directory: PathBuf) {
        self.matrix.cache(&directory);
    }

    pub fn cache_data(&self) -> String {
        let bc = BorderCache {
            col_ix_to_pivot_row: self
                .col_ix_to_pivot_row
                .iter()
                .map(|(k, v)| (*k, *v))
                .collect(),
            pivot_row_to_col_ix: self
                .pivot_row_to_col_ix
                .iter()
                .map(|(k, v)| (*k, *v))
                .collect(),
        };
        serde_json::to_string(&bc).unwrap()
    }

    pub fn from_cache(directory: PathBuf, num_threads: usize) -> Option<Self> {
        let mut border_cache_file = directory.clone();
        border_cache_file.push(BORDER_CACHE);
        if let Ok(border_cache_data) = std::fs::read_to_string(border_cache_file.as_path()) {
            if let Ok(border_cache) = serde_json::from_str::<BorderCache>(&border_cache_data[..]) {
                let mut parallel_matrix_cache = directory.clone();
                parallel_matrix_cache.push(PARALLEL_MATRIX_CACHE);
                if let Some(parallel_matrix) =
                    ParallelFFMatrix::from_disk(parallel_matrix_cache, num_threads)
                {
                    return Some(BorderManager {
                        col_ix_to_pivot_row: border_cache.col_ix_to_pivot_row.into_iter().collect(),
                        pivot_row_to_col_ix: border_cache.pivot_row_to_col_ix.into_iter().collect(),
                        matrix: parallel_matrix,
                    });
                } else {
                    return None;
                }
            } else {
                return None;
            }
        } else {
            return None;
        }
    }

    pub fn add_pivots(&mut self, pivots: Vec<(usize, usize)>) {
        for (row_ix, col_ix) in pivots {
            self.pivot_row_to_col_ix.insert(row_ix, col_ix);
            self.col_ix_to_pivot_row.insert(col_ix, row_ix);
        }
    }
}

#[derive(Debug, Serialize, Deserialize)]
struct RateDistCache {
    directory: PathBuf,
    quotient: FFPolynomial,
    dim: usize,

    data_id_to_col_ix: Vec<(u64, usize)>,
    col_ix_to_data_id: Vec<(usize, u64)>,
    next_col_ix: usize,
    /// Stored in increasing order.
    completed_nodes: Vec<u32>,
    /// Stored in decreasing order so that way nodes can just be popped
    /// off to process.
    remaining_nodes: Vec<u32>,
    /// How many nodes to process before computing an (n,k,d) tuple.
    num_nodes_per_checkpoint: usize,
    /// Map from a checkpoint node to the value of n, k, and d
    /// after processing all nodes up to and including the node.
    code_checkpoints: Vec<(u32, (usize, usize, Option<usize>))>,
}

#[derive(Debug)]
pub struct RateAndDistConfig {
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

    interior_manager: InteriorManager,
    border_manager: BorderManager,
}

impl RateAndDistConfig {
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
        let num_threads = std::thread::available_parallelism().unwrap().into();
        // let num_threads = 1;
        RateAndDistConfig {
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
            interior_manager: InteriorManager::new(field_mod),
            border_manager: BorderManager::new(field_mod, num_threads),
        }
    }

    pub fn create_cache(&self) -> RateDistCache {
        RateDistCache {
            directory: self.directory.clone(),
            quotient: self.quotient.clone(),
            dim: self.dim,
            data_id_to_col_ix: self
                .data_id_to_col_ix
                .iter()
                .map(|(k, v)| (*k, *v))
                .collect(),
            col_ix_to_data_id: self
                .col_ix_to_data_id
                .iter()
                .map(|(k, v)| (*k, *v))
                .collect(),
            next_col_ix: self.next_col_ix,
            completed_nodes: self.completed_nodes.clone(),
            remaining_nodes: self.remaining_nodes.clone(),
            num_nodes_per_checkpoint: self.num_nodes_per_checkpoint,
            code_checkpoints: self
                .code_checkpoints
                .iter()
                .map(|(k, v)| (*k, *v))
                .collect(),
        }
    }

    pub fn cache(&self) {
        let cache = self.create_cache();
        let mut cache_file = self.directory.clone();
        cache_file.push(RATE_AND_DIST_CACHE);
        std::fs::write(cache_file, serde_json::to_string(&cache).unwrap()).unwrap();

        let interior_cache = serde_json::to_string(&self.interior_manager).unwrap();
        let mut interior_cache_file = self.directory.clone();
        interior_cache_file.push(INTERIOR_CACHE);
        std::fs::write(interior_cache_file, &interior_cache[..]).unwrap();

        let border_cache_data = self.border_manager.cache_data();
        let mut border_cache_file = self.directory.clone();
        border_cache_file.push(BORDER_CACHE);
        std::fs::write(border_cache_file, &border_cache_data[..]).unwrap();

        let mut parallel_matrix_cache = self.directory.clone();
        parallel_matrix_cache.push(PARALLEL_MATRIX_CACHE);
        self.border_manager
            .save_matrix_to_disk(parallel_matrix_cache);
    }

    pub fn from_cache(directory: PathBuf, num_threads: usize) -> Option<Self> {
        let mut cache_file = directory.clone();
        cache_file.push(RATE_AND_DIST_CACHE);
        let cache_data = if let Ok(data) = std::fs::read_to_string(cache_file) {
            data
        } else {
            return None;
        };
        let cache: RateDistCache = if let Ok(cache) = serde_json::from_str(&cache_data[..]) {
            cache
        } else {
            return None;
        };

        let mut interior_cache_file = directory.clone();
        interior_cache_file.push(INTERIOR_CACHE);
        let interior_data = if let Ok(data) = std::fs::read_to_string(interior_cache_file) {
            data
        } else {
            return None;
        };
        let interior_manager: InteriorManager =
            if let Ok(interior) = serde_json::from_str(&interior_data[..]) {
                interior
            } else {
                return None;
            };

        let border_manager = BorderManager::from_cache(directory.clone(), num_threads)?;
        todo!()
        // Some(RateAndDistConfig {
        //     directory,
        //     quotient: cache.quotient,
        //     dim: cache.dim,
        //     data_id_to_col_ix: cache.data_id_to_col_ix.into_iter().collect(),
        //     col_ix_to_data_id: cache.col_ix_to_data_id.into_iter().collect(),
        //     next_col_ix: cache.next_col_ix,
        //     completed_nodes: cache.completed_nodes,
        //     remaining_nodes: cache.remaining_nodes,
        //     num_nodes_per_checkpoint: cache.num_nodes_per_checkpoint,
        //     code_checkpoints: cache.code_checkpoints,
        //     row_ix_to_check_id: HashMap::new(),
        //     check_id_to_row_ixs: HashMap::new(),
        //     next_row_ix: 0,
        //     interior_manager,
        //     border_manager,
        // })
    }
    fn eliminate_interior_from_border_pivots(&mut self, border_pivots: Vec<(usize, usize)>) {
        for (row_ix, col_ix) in border_pivots {
            let pivot_row = self.border_manager.matrix.get_row(row_ix);
            self.interior_manager
                .matrix
                .eliminate_col_with_pivot(&pivot_row, col_ix);
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
        let mut interior_pivots_added = Vec::new();
        for node in next_node_batch.iter() {
            interior_pivots_added.append(&mut self.process_node_interior(hg, local_code, *node));
        }
        self.interior_manager.add_pivots(interior_pivots_added);
        let elapsed_interior = start_interior.elapsed().as_secs_f64();

        let start_border = Instant::now();
        let mut border_pivots_added = Vec::new();
        for node in next_node_batch.iter() {
            border_pivots_added.append(&mut self.process_node_border(hg, local_code, *node));
        }

        self.border_manager.add_pivots(border_pivots_added.clone());
        self.eliminate_interior_from_border_pivots(border_pivots_added);
        let elapsed_border = start_border.elapsed().as_secs_f64();

        let last_node = *next_node_batch.last().unwrap();
        self.completed_nodes.append(&mut next_node_batch);

        let n = self.col_ix_to_data_id.len();
        let tot_num_pivots = self.interior_manager.pivot_row_to_col_ix.len()
            + self.border_manager.pivot_row_to_col_ix.len();
        let k = n - tot_num_pivots;
        if compute_distance {
            let start_distance = Instant::now();
            let mut col_weights = Vec::new();
            for col_ix in self.col_ix_to_data_id.keys() {
                if self
                    .interior_manager
                    .col_ix_to_pivot_row
                    .contains_key(col_ix)
                    || self.border_manager.col_ix_to_pivot_row.contains_key(col_ix)
                {
                    continue;
                }
                let col_weight = self.interior_manager.matrix.get_col_weight(*col_ix)
                    + self.border_manager.matrix.get_col_weight(*col_ix);
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
            println!("% time spent border: {:}", elapsed_border / tot_time);
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
            println!("% time spent border: {:}", elapsed_border / tot_time);
            self.code_checkpoints.push((last_node, (n, k, None)));
        }
    }

    fn process_node_border(
        &mut self,
        hg: &HGraph<NodeData, ()>,
        local_code: &ReedSolomon,
        node: u32,
    ) -> Vec<(usize, usize)> {
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
        let border_checks = data_ids
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
        let mut new_entries = Vec::new();
        let local_parity_check = local_code.parity_check_matrix();
        let mut row_ixs_added = HashSet::new();
        for border_check in border_checks {
            let mut data_ids_visible = hg.maximal_edges(&border_check);
            data_ids_visible.sort();
            for data_id in data_ids_visible.iter() {
                if self.data_id_to_col_ix.contains_key(data_id) {
                    continue;
                }
                self.data_id_to_col_ix.insert(*data_id, self.next_col_ix);
                self.col_ix_to_data_id.insert(self.next_col_ix, *data_id);
                self.next_col_ix += 1;
            }

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
        }

        entries_to_rows(new_entries, self.quotient.field_mod)
            .into_iter()
            .filter_map(|(row_ix, mut row)| {
                self.interior_manager.reduce_external_row(&mut row);
                self.border_manager.add_row_and_pivot(row_ix, row)
            })
            .collect()
    }

    /// Returns pivots added to the interior matrix.
    fn process_node_interior(
        &mut self,
        hg: &HGraph<NodeData, ()>,
        local_code: &ReedSolomon,
        node: u32,
    ) -> Vec<(usize, usize)> {
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
        let interior_checks = containing_edges
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
        for data_id in data_ids.iter() {
            if self.data_id_to_col_ix.contains_key(data_id) {
                continue;
            }
            self.data_id_to_col_ix.insert(*data_id, self.next_col_ix);
            self.col_ix_to_data_id.insert(self.next_col_ix, *data_id);
            self.next_col_ix += 1;
        }

        let mut new_entries = Vec::new();
        let local_parity_check = local_code.parity_check_matrix();
        let mut row_ixs_added = HashSet::new();
        for interior_check in interior_checks {
            let mut data_ids_visible = hg.maximal_edges(&interior_check);
            data_ids_visible.sort();
            for data_id in data_ids_visible.iter() {
                if self.data_id_to_col_ix.contains_key(data_id) {
                    continue;
                }
                self.data_id_to_col_ix.insert(*data_id, self.next_col_ix);
                self.col_ix_to_data_id.insert(self.next_col_ix, *data_id);
                self.next_col_ix += 1;
            }

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
        }
        self.interior_manager.add_entries(new_entries);
        let row_ixs_added = row_ixs_added.into_iter().collect::<Vec<usize>>();
        let mut new_row_ixs = row_ixs_added.clone();
        new_row_ixs.sort();
        let mut pivots_added = Vec::new();
        for row_ix in new_row_ixs {
            if let Some(col_ix) = self
                .interior_manager
                .matrix
                .pivotize_row_within_range(row_ix, &row_ixs_added[..])
            {
                pivots_added.push((row_ix, col_ix));
            }
        }
        pivots_added
    }

    pub fn run(&mut self) {
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
        while self.remaining_nodes.len() > 0 {
            num_batches_done += 1;
            if num_batches_done % num_batches_per_dist == 0 {
                self.process_node_batch(&hg, &local_code, true);
            } else {
                self.process_node_batch(&hg, &local_code, false);
            }
        }
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
    /// returns the resulting interior and border matrices as `(interior, border)`
    pub fn quit(self) -> (SparseFFMatrix, SparseFFMatrix) {
        (
            self.interior_manager.matrix,
            self.border_manager.matrix.quit(),
        )
    }
}

#[cfg(test)]
mod tests {
    use std::{path::PathBuf, str::FromStr};

    use simple_logger::SimpleLogger;

    use crate::{math::polynomial::FFPolynomial, rate_and_dist_estimator::RateAndDistConfig};

    #[test]
    fn single_node_example() {
        let q = FFPolynomial::from_str("x^2 + 2x + 2 % 3").unwrap();
        let dim = 3;
        let dir = PathBuf::from_str("/Users/matt/repos/qec/tmp/single_node_test").unwrap();
        let _logger = SimpleLogger::new().init().unwrap();
        let mut rate_estimator =
            RateAndDistConfig::new(q, dim, Some(5_000), Some(1_000), Some(2), dir, 10);
        dbg!(&rate_estimator.remaining_nodes.len());
        rate_estimator.run();
        let (i, b) = rate_estimator.quit();
        // println!("interior: {:}", i.to_dense());
        // println!("border: {:}", b.to_dense());
        // dbg!(&rate_estimator);
    }
}
