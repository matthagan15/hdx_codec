use indexmap::IndexMap;
use mhgl::{HGraph, HyperGraph};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::fs::{self, read, read_to_string};
use std::io::Read;
use std::path::Path;
use std::thread::{available_parallelism, panicking};
use std::{collections::HashSet, fs::File, path::PathBuf, str::FromStr};

use crate::hdx_code;
use crate::math::coset_complex_bfs::bfs;
use crate::math::finite_field::FFRep;
use crate::matrices::mat_trait::RankMatrix;
use crate::matrices::parallel_matrix::ParallelFFMatrix;
use crate::{
    code::Code,
    math::{finite_field::FiniteField, polynomial::FFPolynomial},
    matrices::sparse_ffmatrix::{MemoryLayout, SparseFFMatrix},
    reed_solomon::ReedSolomon,
};

const RANK_CONFIG_FILE_NAME: &str = "rank_config.json";
const HGRAPH_CACHE_FILE_NAME: &str = "hgraph_cache";
const INTERIOR_MATRIX_CACHE_NAME: &str = "interior_matrix_cache";
const BORDER_MATRIX_CACHE_NAME: &str = "border_matrix_cache";

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
enum ComputationState {
    Start,
    BFS,
    ComputeMatrices,
    BorderRank,
    Done,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct RankConfig {
    quotient_poly: FFPolynomial,
    dim: usize,
    rs_degree: usize,
    truncation: usize,
    rate_upper_bound: Option<f64>,
    computation_state: ComputationState,
    step_size: usize,
    truncation_to_rate: Vec<(usize, f64)>,
}

impl RankConfig {
    pub fn new(
        quotient_poly: FFPolynomial,
        dim: usize,
        rs_degree: usize,
        truncation: usize,
        step_size: usize,
    ) -> Self {
        Self {
            quotient_poly,
            dim,
            rs_degree,
            truncation,
            rate_upper_bound: None,
            computation_state: ComputationState::Start,
            step_size,
            truncation_to_rate: Vec::new(),
        }
    }
    pub fn from_disk(config_dir: impl AsRef<Path>) -> Option<Self> {
        let mut config_path = PathBuf::from(config_dir.as_ref());
        if !config_path.is_dir() {
            panic!("Provided configuration path is not a directory.")
        }
        config_path.push(RANK_CONFIG_FILE_NAME);
        let config_file_data = std::fs::read_to_string(config_path.as_path()).ok();
        let config = config_file_data
            .map(|config_data| serde_json::from_str::<RankConfig>(&config_data[..]).ok());
        config.flatten()
    }

    pub fn save_to_disk(&self, config_dir: impl AsRef<Path>) {
        let mut config_path = PathBuf::from(config_dir.as_ref());
        if !config_path.is_dir() {
            panic!("Provided configuration path is not a directory.")
        }
        config_path.push(RANK_CONFIG_FILE_NAME);
        let config_serialized =
            serde_json::to_string_pretty(self).expect("Could not serialize RankConfig.");
        std::fs::write(config_path, &config_serialized[..])
            .expect("Could not write RankConfig to disk.");
    }

    fn next_truncation(&self) -> Option<usize> {
        if self.truncation_to_rate.len() > 0 {
            let (prev_truncation, _) = self.truncation_to_rate.last().unwrap();
            let next = self.truncation.min(prev_truncation + self.step_size);
            if next == self.truncation {
                None
            } else {
                Some(next)
            }
        } else {
            Some(self.step_size)
        }
    }

    pub fn run(&mut self, config_dir: impl AsRef<Path>, num_threads: usize) {
        let mut current_truncation = 0;
        let mut hgraph = None;
        let mut border_matrix = None;
        let mut current_interior_pivots = None;
        let mut num_cols = None;
        // let mut num_cols = None;
        while self.computation_state != ComputationState::Done {
            match self.computation_state {
                ComputationState::Start => {
                    // log::trace!("START");
                    println!("START");
                    let next_truncation = self.next_truncation();
                    if next_truncation.is_none() {
                        self.computation_state = ComputationState::Done;
                    } else {
                        current_truncation = next_truncation.unwrap();
                        self.computation_state = ComputationState::BFS;
                    }
                    log::trace!("Next truncation: {:?}", next_truncation);
                    self.save_to_disk(config_dir.as_ref());
                    // self.computation_state = ComputationState::BFS;
                }
                ComputationState::BFS => {
                    log::trace!("BFS");
                    let mut hgraph_cache = PathBuf::from(config_dir.as_ref());
                    hgraph_cache.push(HGRAPH_CACHE_FILE_NAME);
                    let (hg, _new_edges) = bfs(
                        self.quotient_poly.clone(),
                        self.dim,
                        Some(current_truncation),
                        Some(hgraph_cache),
                        None,
                    );
                    hgraph = Some(hg);
                    self.computation_state = ComputationState::ComputeMatrices;
                    self.save_to_disk(config_dir.as_ref());
                }
                ComputationState::ComputeMatrices => {
                    println!("COMPUTE MATRICES");
                    if hgraph.is_none() {
                        println!("No hgraph?");
                        self.computation_state = ComputationState::BFS;
                        continue;
                    }
                    let hgraph = hgraph.as_ref().unwrap();
                    let local_rs_code =
                        ReedSolomon::new(self.quotient_poly.field_mod, self.rs_degree);
                    #[derive(Debug, Serialize, Deserialize)]
                    struct MatrixCache {
                        matrix: SparseFFMatrix,
                        pivots: Vec<(usize, usize)>,
                        /// Maps a parity check id to the start row ix in the matrix
                        used_edges: Vec<(u64, usize)>,
                    }
                    let mut message_ids = hgraph.edges_of_size(self.dim);
                    message_ids.sort();
                    let message_ids_len = message_ids.len();
                    let message_id_to_ix: HashMap<u64, usize> = message_ids
                        .into_iter()
                        .zip((0..message_ids_len).into_iter())
                        .collect();
                    let mut interior_matrix_cache_file = PathBuf::from(config_dir.as_ref());
                    interior_matrix_cache_file.push(INTERIOR_MATRIX_CACHE_NAME);
                    let interior_matrix_data = read_to_string(interior_matrix_cache_file.as_path());
                    let interior_matrix_parsed: Option<MatrixCache> =
                        if let Ok(data) = interior_matrix_data {
                            serde_json::from_str(&data[..]).ok()
                        } else {
                            None
                        };
                    let mut interior_matrix_cache = if let Some(cache) = interior_matrix_parsed {
                        cache
                    } else {
                        MatrixCache {
                            matrix: SparseFFMatrix::new(
                                0,
                                0,
                                self.quotient_poly.field_mod,
                                MemoryLayout::RowMajor,
                            ),
                            pivots: Vec::new(),
                            used_edges: Vec::new(),
                        }
                    };

                    let mut border_matrix_cache_file = PathBuf::from(config_dir.as_ref());
                    border_matrix_cache_file.push(BORDER_MATRIX_CACHE_NAME);
                    let border_matrix_data = read_to_string(border_matrix_cache_file.as_path());
                    let border_matrix_parsed: Option<MatrixCache> =
                        if let Ok(data) = border_matrix_data {
                            serde_json::from_str(&data[..]).ok()
                        } else {
                            None
                        };
                    let mut border_matrix_cache = if let Some(cache) = border_matrix_parsed {
                        cache
                    } else {
                        MatrixCache {
                            matrix: SparseFFMatrix::new(
                                0,
                                0,
                                self.quotient_poly.field_mod,
                                MemoryLayout::RowMajor,
                            ),
                            pivots: Vec::new(),
                            used_edges: Vec::new(),
                        }
                    };
                    let border_edge_to_index: HashMap<u64, usize> =
                        border_matrix_cache.used_edges.iter().cloned().collect();

                    let hgraph_edges = hgraph.edges_of_size(self.dim - 1);
                    let mut interior_edge_to_index: HashMap<u64, usize> =
                        interior_matrix_cache.used_edges.iter().cloned().collect();

                    let mut new_interior_checks = Vec::new();
                    let mut new_border_checks = Vec::new();
                    let mut local_checks = HashMap::new();
                    for edge in hgraph_edges.into_iter().filter(|e_id| {
                        hgraph.maximal_edges(e_id).len() == local_rs_code.encoded_len()
                    }) {
                        if interior_edge_to_index.contains_key(&edge)
                            || border_edge_to_index.contains_key(&edge)
                        {
                            continue;
                        }
                        let nodes = hgraph.query_edge(&edge).unwrap();
                        let type_zero_node: Vec<u32> = nodes
                            .into_iter()
                            .filter(|node| *hgraph.get_node(node).unwrap() == 0)
                            .collect();
                        if type_zero_node.len() == 1 {
                            new_interior_checks.push(edge);
                            let local_check_edges =
                                local_checks.entry(type_zero_node[0]).or_insert(Vec::new());
                            local_check_edges.push(edge);
                        } else {
                            new_border_checks.push(edge);
                        }
                    }

                    let mut new_row_indices = Vec::new();
                    let pivot_col_to_pivot_row: HashMap<usize, usize> = interior_matrix_cache
                        .pivots
                        .iter()
                        .map(|(row, col)| (*col, *row))
                        .collect();
                    let mut pivots_to_use = HashSet::new();
                    for new_interior_check in new_interior_checks.into_iter() {
                        let row_start = interior_edge_to_index
                            .entry(new_interior_check)
                            .or_insert(interior_matrix_cache.matrix.n_rows);
                        let mut message_ids = hgraph.maximal_edges(&new_interior_check);
                        message_ids.sort();
                        let message_ixs: Vec<usize> = message_ids
                            .into_iter()
                            .map(|id| *message_id_to_ix.get(&id).unwrap())
                            .collect();
                        let local_parity_check = local_rs_code.parity_check_matrix();
                        let mut new_entries = Vec::new();
                        let mut new_row_ix = *row_start;
                        for row_ix in 0..local_parity_check.n_rows {
                            new_row_indices.push(new_row_ix);
                            for col_ix in 0..local_parity_check.n_cols {
                                // dbg!(&local_parity_check);
                                // dbg!(col_ix);
                                // dbg!(&message_ixs);
                                let new_col_ix = message_ixs[col_ix];
                                if pivot_col_to_pivot_row.contains_key(&new_col_ix) {
                                    pivots_to_use.insert(new_col_ix);
                                }
                                new_entries.push((
                                    new_row_ix,
                                    new_col_ix,
                                    local_parity_check[[row_ix, col_ix]].0,
                                ));
                            }
                            new_row_ix += 1;
                        }
                        interior_matrix_cache.matrix.insert_entries(new_entries);
                    }
                    for pivot_col in pivots_to_use {
                        let pivot_row = pivot_col_to_pivot_row.get(&pivot_col).unwrap();
                        interior_matrix_cache
                            .matrix
                            .eliminate_rows((*pivot_row, pivot_col), &new_row_indices[..]);
                    }

                    let mut new_border_rows = Vec::new();
                    let pivot_col_to_pivot_row: HashMap<usize, usize> = border_matrix_cache
                        .pivots
                        .iter()
                        .map(|(row, col)| (*col, *row))
                        .collect();
                    let mut pivots_to_use = HashSet::new();
                    for new_border_check in new_border_checks.into_iter() {
                        let row_start =
                            if let Some(row_ix) = border_edge_to_index.get(&new_border_check) {
                                *row_ix
                            } else {
                                border_matrix_cache.matrix.n_rows
                            };
                        let mut message_ids = hgraph.maximal_edges(&new_border_check);
                        message_ids.sort();
                        let message_ixs: Vec<usize> = message_ids
                            .into_iter()
                            .map(|id| *message_id_to_ix.get(&id).unwrap())
                            .collect();
                        let local_parity_check = local_rs_code.parity_check_matrix();
                        let mut new_entries = Vec::new();
                        let mut new_row_ix = row_start;
                        for row_ix in 0..local_parity_check.n_rows {
                            new_border_rows.push(new_row_ix);
                            for col_ix in 0..local_parity_check.n_cols {
                                let new_col_ix = message_ixs[col_ix];
                                if pivot_col_to_pivot_row.contains_key(&new_col_ix) {
                                    pivots_to_use.insert(new_col_ix);
                                }
                                new_entries.push((
                                    new_row_ix,
                                    new_col_ix,
                                    local_parity_check[[row_ix, col_ix]].0,
                                ));
                            }
                            new_row_ix += 1;
                        }
                        border_matrix_cache.matrix.insert_entries(new_entries);
                    }
                    for pivot_col in pivots_to_use {
                        let pivot_row = pivot_col_to_pivot_row.get(&pivot_col).unwrap();
                        border_matrix_cache
                            .matrix
                            .eliminate_rows((*pivot_row, pivot_col), &new_row_indices[..]);
                    }

                    // Now need to reduce the new interior checks.
                    for (local_check, new_interior_edges) in local_checks {
                        let mut rows_to_pivot = Vec::new();
                        for interior_edge in new_interior_edges.iter() {
                            let row_start = interior_edge_to_index.get(interior_edge).unwrap();
                            for new_row in *row_start..*row_start + local_rs_code.parity_check_len()
                            {
                                rows_to_pivot.push(new_row);
                            }
                        }
                        let mut new_pivots = Vec::new();
                        for new_row in rows_to_pivot.iter() {
                            if let Some(new_pivot) = interior_matrix_cache
                                .matrix
                                .pivotize_row_within_range(*new_row, &rows_to_pivot[..])
                            {
                                new_pivots.push((*new_row, new_pivot));
                            }
                        }
                        interior_matrix_cache.pivots.append(&mut new_pivots);
                    }
                    let interior_pivot_col_to_row: HashMap<usize, usize> = interior_matrix_cache
                        .pivots
                        .iter()
                        .map(|(row, col)| (*col, *row))
                        .collect();
                    for row in new_border_rows {
                        if let Some(mut old_row) = border_matrix_cache.matrix.remove_row(row) {
                            let mut first_nnz = old_row.first_nonzero();
                            while first_nnz.is_some() {
                                let (pivot_col, to_pivot_entry) = first_nnz.unwrap();
                                if let Some(pivot_row) = interior_pivot_col_to_row.get(&pivot_col) {
                                    let interior_row =
                                        interior_matrix_cache.matrix.get_row(*pivot_row);
                                    let scalar = -1
                                        * FiniteField::new(
                                            to_pivot_entry,
                                            self.quotient_poly.field_mod,
                                        );
                                    old_row.add_scaled_row_to_self(scalar, &interior_row);
                                    first_nnz = old_row.first_nonzero();
                                } else {
                                    first_nnz = None;
                                }
                            }
                        }
                    }
                    // let (local_checks, interior_checks, border_checks, message_ids) =
                    //     split_hg_into_checks(hgraph, &local_rs_code);
                    // let (interior, border, num_interior_pivots) = build_matrices_from_checks(
                    //     local_checks,
                    //     interior_checks,
                    //     border_checks,
                    //     message_ids,
                    //     hgraph,
                    //     &local_rs_code,
                    // );

                    border_matrix = Some(border_matrix_cache.matrix);
                    current_interior_pivots = Some(interior_matrix_cache.pivots.len());
                    // current_interior_pivots = Some(num_interior_pivots);
                    num_cols = Some(interior_matrix_cache.matrix.n_cols);
                    // let new_interior_rate = 1.0
                    // - (interior_matrix_cache.pivots.len() as f64 / num_cols.unwrap() as f64);
                    // if self.rate_upper_bound.is_none() {
                    //     self.rate_upper_bound = Some(new_interior_rate);
                    // } else {
                    //     self.rate_upper_bound =
                    //         Some(new_interior_rate.max(self.rate_upper_bound.unwrap()));
                    // }

                    // TODO update the caches, don't forget about the used_edges field
                    // println!("INTERIOR MATRIX:\n");
                    // interior_matrix_cache.matrix.dense_print();
                    // println!("BORDER MATRIX:\n");
                    // border_matrix_cache.matrix.dense_print();
                    self.computation_state = ComputationState::BorderRank;
                    // println!("{:?}", config_dir.as_ref());
                    // self.save_to_disk(config_dir.as_ref());
                    // self.truncation_to_rate.push();
                }
                ComputationState::BorderRank => {
                    println!("BORDER RANK");
                    let mut parallel_border_mat = match ParallelFFMatrix::from_disk(
                        PathBuf::from(config_dir.as_ref()),
                        num_threads,
                    ) {
                        Some(mut par_mat) => {
                            // although the border matrix was just read from disk we will instead use the parallelized
                            // cached version to avoid keeping two copies on hand. Keep border around until this point
                            // in case the reading from disk fails though.
                            if par_mat.num_threads() != num_threads {
                                log::error!("Number of threads provided does not match number of cache points found.");
                                let mut big_mat = par_mat.quit();
                                par_mat = big_mat.split_into_parallel(
                                    big_mat.row_ixs().into_iter().collect(),
                                    num_threads,
                                );
                            }
                            par_mat
                        }
                        None => {
                            if border_matrix.is_none() {
                                self.computation_state = ComputationState::Start;
                                continue;
                            }
                            let mut border: SparseFFMatrix = border_matrix.unwrap();
                            let row_ixs = border.row_ixs();
                            border.split_into_parallel(row_ixs.into_iter().collect(), num_threads)
                        }
                    };
                    border_matrix = None;
                    println!("............ Reducing Border Matrix ............");
                    // moving this into a new thread so that we can increase the stack size and name the thread
                    // to try and debug this stack overflow issue.
                    let stack_size = if current_truncation > 10_000_000 {
                        64 * 1024
                    } else {
                        32 * 1024
                    };
                    let thread_builder = std::thread::Builder::new()
                        .stack_size(stack_size)
                        .name("Coordinator".into());
                    let cache_dir = PathBuf::from(config_dir.as_ref());

                    let handle = thread_builder
                        .spawn(move || {
                            let pivots = parallel_border_mat.row_echelon_form(
                                Some(&cache_dir),
                                Some(10_000),
                                Some(1_000),
                            );
                            parallel_border_mat.quit();
                            pivots
                        })
                        .expect("Could not create coordinator thread");
                    let pivots = handle.join().expect("Parallel matrix solver had error.");
                    log::trace!("Found {:} pivots for the border matrix.", pivots.len());
                    let rate = 1.0
                        - ((current_interior_pivots.unwrap() + pivots.len()) as f64
                            / num_cols.unwrap() as f64);
                    println!("truncation: {:}, rate: {:}", current_truncation, rate);
                    self.truncation_to_rate.push((current_truncation, rate));
                    self.save_to_disk(config_dir.as_ref());
                    hgraph = None;
                    border_matrix = None;
                    current_interior_pivots = None;
                    num_cols = None;
                    self.computation_state = ComputationState::Start;
                    self.save_to_disk(config_dir.as_ref());
                }
                ComputationState::Done => {
                    break;
                }
            }
        }
        println!("DONE");
        println!("results: {:?}", self.truncation_to_rate);
    }
}

/// Computes the upper and lower bounds of the rate for a given quotient polynomia, matrix
/// dimension, and reed solomon degree. Truncates the coset complex truncation, if no truncation
/// is given it computes the full thing.
/// Returns - (`num_pivots_lower_bound`, `num_pivots_upper_bound`, `num_cols`)
pub fn compute_rank_bounds(
    quotient_poly: FFPolynomial,
    dim: usize,
    rs_degree: usize,
    cache_dir: PathBuf,
    truncation: Option<usize>,
    cache_rate: Option<usize>,
    log_rate: Option<usize>,
    num_threads: Option<usize>,
) -> f64 {
    println!("--------------- RANK COMPUTATION ---------------");
    // First check if there is a config in the directory:
    // If config present - make sure parameters match. Then try to load cache
    // If config not present - we are starting from scratch.
    let mut config_file = cache_dir.clone();
    let num_threads = num_threads.unwrap_or(available_parallelism().unwrap().into());
    config_file.push("config.json");
    let coset_complex_size =
        crate::math::coset_complex_bfs::size_of_coset_complex(&quotient_poly, dim);
    let input_config = RankConfig {
        quotient_poly: quotient_poly.clone(),
        dim,
        rs_degree,
        rate_upper_bound: None,
        computation_state: ComputationState::Start,
        truncation: truncation.unwrap_or(coset_complex_size + 1),
        truncation_to_rate: Vec::new(),
        step_size: 1,
    };
    let mut bfs_steps_needed = truncation.unwrap_or(coset_complex_size + 1);
    log::trace!("............ Managing Config ............");
    if config_file.is_file() {
        let config_file_data = std::fs::read_to_string(config_file.as_path())
            .expect("Could not read present config file.");
        if let Ok(config) = serde_json::from_str::<RankConfig>(&config_file_data[..]) {
            if config != input_config {
                if (
                    &input_config.quotient_poly,
                    input_config.dim,
                    input_config.rs_degree,
                ) == (&config.quotient_poly, config.dim, config.rs_degree)
                    && input_config.truncation > config.truncation
                {
                    log::trace!("New truncation found, updating config on disk.");
                    bfs_steps_needed -= config.truncation;
                    std::fs::write(
                        config_file,
                        serde_json::to_string_pretty(&input_config)
                            .expect("Could not serialize config."),
                    )
                    .expect("Could not write config to disk.");
                } else {
                    dbg!(("file_config: ", config));
                    dbg!(("input config:", input_config));
                    log::error!("Mismatch between cached config and input config: bailing.");
                    panic!()
                }
            } else {
                bfs_steps_needed = 0;
            }
        } else {
            log::error!("There exists a rank config file, but it cannot be parsed.");
            panic!()
        }
    } else {
        std::fs::write(
            config_file,
            serde_json::to_string_pretty(&input_config).expect("Could not serialize config."),
        )
        .expect("Could not write config to disk.");
    }
    log::trace!(
        "Quotient Polynomial: {:}\ndim = {:}\nReed-Solomon Degree: {:}\ntruncation: {:}\nnum_threads: {:}",
        quotient_poly,
        dim,
        rs_degree,
        truncation.unwrap(),
        num_threads
    );
    let local_rs_code = ReedSolomon::new(quotient_poly.field_mod, rs_degree);
    /// A valid cache has already had the interior matrices pivotized properly and the border
    /// matrix reduced as much as possible.
    #[derive(Serialize, Deserialize)]
    struct MatrixCache {
        interior_mat: SparseFFMatrix,
        border_mat: SparseFFMatrix,
        num_interior_pivots: usize,
    }
    log::trace!("............ Checking Matrix Cache ............");
    let mut matrix_cache_file = cache_dir.clone();
    matrix_cache_file.push("matrix_cache");
    if bfs_steps_needed > 0 {
        if matrix_cache_file.is_file() {
            std::fs::remove_file(matrix_cache_file.as_path())
                .expect("Could not remove matrix cache.");
        }
    }
    let mut matrix_cache = None;
    if matrix_cache_file.is_file() {
        let matrix_cache_data = std::fs::read_to_string(matrix_cache_file.as_path())
            .expect("Could not read existing matrix cache file.");
        if let Ok(mat_cache) = serde_json::from_str::<MatrixCache>(&matrix_cache_data[..]) {
            log::trace!("Matrix cache successfully retrieved!");
            matrix_cache = Some(mat_cache);
        }
    }
    if matrix_cache.is_none() {
        log::trace!("No Matrix Cache found, starting from scratch.");
        let mut check_cache = cache_dir.clone();
        check_cache.push("check_cache");
        #[derive(Deserialize)]
        struct Checks {
            local: Vec<u32>,
            interior: Vec<u64>,
            border: Vec<u64>,
            message_ids: Vec<u64>,
        }
        let mut checks = None;
        if check_cache.is_file() {
            let check_cache_data = std::fs::read_to_string(check_cache.as_path())
                .expect("Could not read check cache even though file is present.");
            if let Ok(cached_checks) = serde_json::from_str::<Checks>(&check_cache_data[..]) {
                log::trace!("Retrieved cached parity check data.");
                checks = Some(cached_checks);
            }
        }
        let mut hgraph_cache_file = cache_dir.clone();
        hgraph_cache_file.push("hgraph_cache");
        log::trace!("Querying BFS manager.");
        let (hg, _new_edges) = bfs(
            quotient_poly.clone(),
            dim,
            truncation,
            Some(hgraph_cache_file),
            None,
        );
        if checks.is_none() {
            log::trace!("............ Preprocessing HGraph for Matrix ............");
            let splitted_checks = split_hg_into_checks(&hg, &local_rs_code);
            log::trace!("Found {:} local checks, {:} interior checks, {:} border checks, and {:} message ids.", splitted_checks.0.len(), splitted_checks.1.len(), splitted_checks.2.len(), splitted_checks.3.len());
            checks = Some(Checks {
                local: splitted_checks.0,
                interior: splitted_checks.1,
                border: splitted_checks.2,
                message_ids: splitted_checks.3,
            });
        };
        let checks = checks.unwrap();
        log::trace!("............ Building Matrices ............");
        let matrices = build_matrices_from_checks(
            checks.local,
            checks.interior,
            checks.border,
            checks.message_ids,
            &hg,
            &local_rs_code,
        );

        matrix_cache = Some(MatrixCache {
            interior_mat: matrices.0,
            border_mat: matrices.1,
            num_interior_pivots: matrices.2,
        });
        serde_json::to_writer(
            std::fs::File::create(matrix_cache_file.as_path()).unwrap(),
            matrix_cache.as_ref().unwrap(),
        )
        .expect("Could not cache matrices");
    }

    let mats = matrix_cache.unwrap();
    let (interior, mut border) = (mats.interior_mat, mats.border_mat);
    let num_interior_pivots = mats.num_interior_pivots;
    log::trace!(
        "Matrices retreived. Interior Matrix has dims {:}x{:} and {:} pivots were found.",
        interior.n_rows,
        interior.n_cols,
        num_interior_pivots
    );
    log::trace!("Border Matrix is {:}x{:}", border.n_rows, border.n_cols);
    let num_cols = interior.n_cols;
    drop(interior);

    let mut parallel_border_mat = match ParallelFFMatrix::from_disk(cache_dir.clone(), num_threads)
    {
        Some(mut par_mat) => {
            // although the border matrix was just read from disk we will instead use the parallelized
            // cached version to avoid keeping two copies on hand. Keep border around until this point
            // in case the reading from disk fails though.
            drop(border);
            if par_mat.num_threads() != num_threads {
                log::error!(
                    "Number of threads provided does not match number of cache points found."
                );
                let mut big_mat = par_mat.quit();
                par_mat = big_mat
                    .split_into_parallel(big_mat.row_ixs().into_iter().collect(), num_threads);
            }
            par_mat
        }
        None => border.split_into_parallel(border.row_ixs().into_iter().collect(), num_threads),
    };
    log::trace!("............ Reducing Border Matrix ............");
    // moving this into a new thread so that we can increase the stack size and name the thread
    // to try and debug this stack overflow issue.
    let thread_builder = std::thread::Builder::new()
        .stack_size(64 * 1024)
        .name("Coordinator".into());
    let handle = thread_builder
        .spawn(move || {
            let pivots = parallel_border_mat.row_echelon_form(
                Some(cache_dir.as_path()),
                cache_rate,
                log_rate,
            );
            parallel_border_mat.quit();
            pivots
        })
        .expect("Could not create coordinarot thread");
    let pivots = handle.join().expect("Parallel matrix solver had error.");
    log::trace!("Found {:} pivots for the border matrix.", pivots.len());
    let rate = 1.0 - (pivots.len() + num_interior_pivots) as f64 / num_cols as f64;
    log::info!("Final computed rate: {:}", rate);
    rate
}

/// Scans a hgraph for local checks, interior checks, and border checks. Needs reference to
/// the polynomial to determine when a check is complete.
/// Returns (local_checks, interior_checks, border_checks, message_ids)
pub fn split_hg_into_checks(
    hgraph: &HGraph<u16, ()>,
    local_code: &ReedSolomon,
) -> (Vec<u32>, Vec<u64>, Vec<u64>, Vec<u64>) {
    let mut local_checks = HashSet::new();
    let mut message_ids = HashSet::new();
    let mut interior_checks = HashSet::new();
    let mut border_checks = HashSet::new();
    let nodes = hgraph.nodes();
    log::trace!("Found {:} nodes to check.", nodes.len());
    for node in nodes {
        let mut one_complete_edge = false;
        // We only want to get local checks of the first type.
        if *hgraph.get_node(&node).unwrap() != 0 {
            continue;
        }
        let star = hgraph.containing_edges_of_nodes([node]);
        for edge in star {
            let edge_nodes = hgraph.query_edge(&edge).unwrap();
            if edge_nodes.len() == 2 {
                if hgraph.link(&edge).len() == local_code.encoded_len() {
                    interior_checks.insert(edge);
                    one_complete_edge = true;
                }
            } else if edge_nodes.len() == 3 {
                message_ids.insert(edge);
                let border_nodes: Vec<u32> = edge_nodes
                    .into_iter()
                    .filter(|node| *hgraph.get_node(&node).unwrap() != 0)
                    .collect();
                if border_nodes.len() != 2 {
                    panic!("Border nodes should have length 2");
                }
                let border_check = hgraph.find_id(&border_nodes[..]).unwrap();

                if hgraph.maximal_edges(&border_check).len() == local_code.encoded_len() {
                    border_checks.insert(border_check);
                    one_complete_edge = true;
                }
            }
        }
        if one_complete_edge {
            local_checks.insert(node);
        }
    }
    let message_ids: Vec<u64> = message_ids
        .into_iter()
        .filter(|id| {
            let mut contains_one_good_edge = false;
            let nodes = hgraph.query_edge(&id).unwrap();
            'outer: for node_1 in nodes.iter() {
                for node_2 in nodes.iter() {
                    if *node_1 == *node_2 {
                        continue;
                    }
                    let check_id = hgraph.find_id([*node_1, *node_2]);
                    if check_id.is_none() {
                        println!(
                            "node_1: {:}, node_2: {:}, message_id: {:}",
                            node_1, node_2, id
                        );
                        dbg!(nodes);
                        println!("PANIC: BAD HGRAPH.\n{:}", hgraph);
                        dbg!(hgraph);
                        panic!()
                    }
                    let check_id = check_id.unwrap();
                    if interior_checks.contains(&check_id) || border_checks.contains(&check_id) {
                        contains_one_good_edge = true;
                        break 'outer;
                    }
                }
            }
            contains_one_good_edge
        })
        .collect();
    (
        local_checks.into_iter().collect(),
        interior_checks.into_iter().collect(),
        border_checks.into_iter().collect(),
        message_ids,
    )
}

/// Constructs the interior and border check matrices from the given parity checks.
/// The interior matrix will be in row-echelon form and the border matrix will be as
/// reduced as possible from the interior check pivots.
///
/// Assumes provided inputs satisfy the following:
/// - each parity check, whether interior or border, has a complete check view for the provided
///     local code.
/// - Each message_id has at least one parity check that sees the given message symbol.
fn build_matrices_from_checks(
    local_checks: Vec<u32>,
    interior_checks: Vec<u64>,
    border_checks: Vec<u64>,
    message_ids: Vec<u64>,
    hgraph: &HGraph<u16, ()>,
    local_code: &ReedSolomon,
) -> (SparseFFMatrix, SparseFFMatrix, usize) {
    let num_messages = message_ids.len();
    let message_id_to_ix: HashMap<u64, usize> =
        message_ids.into_iter().zip(0..num_messages).collect();
    let mut message_id_to_checks: HashMap<u64, HashSet<u64>> =
        HashMap::with_capacity(message_id_to_ix.len());
    let interior_check_to_message_ixs: HashMap<_, _> = interior_checks
        .into_iter()
        .map(|check| {
            let max_edges = hgraph.maximal_edges(&check);
            let mut max_ixs = Vec::new();
            for max_edge in max_edges {
                let max_edge_to_checks = message_id_to_checks.entry(max_edge).or_default();
                max_edge_to_checks.insert(check);
                max_ixs.push(*message_id_to_ix.get(&max_edge).expect("Interior check was collected that did not have a corresponding message_id accounted for."));
            }
            max_ixs.sort();
            max_ixs.dedup();
            (check, max_ixs)
        })
        .collect();
    let border_check_to_message_ixs: HashMap<_, _> = border_checks
    .into_iter()
    .map(|check| {
        let max_edges = hgraph.maximal_edges(&check);
        let mut max_ixs = Vec::new();
        for max_edge in max_edges {
            let max_edge_to_checks = message_id_to_checks.entry(max_edge).or_default();
            max_edge_to_checks.insert(check);
            max_ixs.push(*message_id_to_ix.get(&max_edge).expect("Interior check was collected that did not have a corresponding message_id accounted for."));
        }
        max_ixs.sort();
        max_ixs.dedup();
        (check, max_ixs)
    })
    .collect();
    let mut interior_entries: Vec<(usize, usize, FFRep)> = Vec::new();
    let mut border_entries: Vec<(usize, usize, FFRep)> = Vec::new();
    let num_rows_per_check = local_code.parity_check_len();
    let mut counter = 0;
    let interior_check_to_row_ix_offset: HashMap<u64, usize> = interior_check_to_message_ixs
        .keys()
        .map(|check| {
            let ret = counter * num_rows_per_check;
            counter += 1;
            (*check, ret)
        })
        .collect();
    counter = 0;
    let border_check_to_row_ix_offset: HashMap<u64, usize> = border_check_to_message_ixs
        .keys()
        .map(|check| {
            let ret = counter * num_rows_per_check;
            counter += 1;
            (*check, ret)
        })
        .collect();
    for (message_id, message_ix) in message_id_to_ix {
        // These are the checks that can see the current message_id
        let viewing_checks = message_id_to_checks.get(&message_id).unwrap();
        for viewing_check in viewing_checks.iter() {
            let is_interior = interior_check_to_message_ixs.contains_key(viewing_check);
            if is_interior && border_check_to_message_ixs.contains_key(viewing_check)
                || !is_interior && !border_check_to_message_ixs.contains_key(viewing_check)
            {
                log::error!(
                    "Parity check encountered which is not in border or interior check data."
                );
                panic!()
            }
            let local_view = if is_interior {
                interior_check_to_message_ixs.get(viewing_check).unwrap()
            } else {
                border_check_to_message_ixs.get(viewing_check).unwrap()
            };
            let parity_check_input = local_view
                .into_iter()
                .map(|ix| {
                    FiniteField::new(
                        if *ix == message_ix { 1 } else { 0 },
                        local_code.field_mod(),
                    )
                })
                .collect();
            let parity_check = local_code.parity_check(&parity_check_input);
            for ix in 0..parity_check.len() {
                let row_ix = if is_interior {
                    ix + interior_check_to_row_ix_offset.get(viewing_check).unwrap()
                } else {
                    ix + border_check_to_row_ix_offset.get(viewing_check).unwrap()
                };
                let col_ix = message_ix;
                let entry = parity_check[ix].0;
                if is_interior {
                    interior_entries.push((row_ix, col_ix, entry));
                } else {
                    border_entries.push((row_ix, col_ix, entry));
                }
            }
        }
    }
    let mut interior_matrix = SparseFFMatrix::new_with_entries(
        0,
        0,
        local_code.field_mod(),
        MemoryLayout::RowMajor,
        interior_entries,
    );
    let mut border_matrix = SparseFFMatrix::new_with_entries(
        0,
        0,
        local_code.field_mod(),
        MemoryLayout::RowMajor,
        border_entries,
    );
    let mut num_interior_pivots = 0;
    for node in local_checks {
        let containing_edges = hgraph.containing_edges_of_nodes([node]);
        let mut messages = Vec::new();
        let mut interior_checks = Vec::new();
        let mut rows_to_pivot = Vec::new();
        for containing_edge in containing_edges {
            if message_id_to_checks.contains_key(&containing_edge) {
                messages.push(containing_edge);
            } else if interior_check_to_row_ix_offset.contains_key(&containing_edge) {
                interior_checks.push(containing_edge);
                let mut rows: Vec<usize> = (0..local_code.parity_check_len())
                    .map(|ix| {
                        ix + interior_check_to_row_ix_offset
                            .get(&containing_edge)
                            .unwrap()
                    })
                    .collect();
                rows_to_pivot.append(&mut rows);
            }
        }
        let border_checks: HashSet<u64> = hgraph
            .link_of_nodes([node])
            .into_iter()
            .filter(|(_, nodes)| nodes.len() == 2)
            .map(|(id, _)| id)
            .filter(|id| border_check_to_row_ix_offset.contains_key(id))
            .collect();
        let mut border_rows_to_eliminate = Vec::new();
        for border_check in border_checks.into_iter() {
            let mut rows: Vec<usize> = (0..local_code.parity_check_len())
                .map(|ix| ix + border_check_to_row_ix_offset.get(&border_check).unwrap())
                .collect();
            border_rows_to_eliminate.append(&mut rows);
        }
        for row_ix in rows_to_pivot.clone() {
            if let Some(pivot_col) =
                interior_matrix.pivotize_row_within_range(row_ix, &rows_to_pivot[..])
            {
                num_interior_pivots += 1;
                let pivot_row = interior_matrix.get_row(row_ix);
                border_matrix.pivotize_with_row((row_ix, pivot_col), pivot_row);
            }
        }
    }
    (interior_matrix, border_matrix, num_interior_pivots)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RankEstimatorConfig {
    pub quotient_poly: String,
    pub dim: usize,
    pub reed_solomon_degree: usize,
    pub cache_file: Option<String>,
    pub hgraph_file: PathBuf,
    pub num_threads: usize,
}

impl RankEstimatorConfig {
    pub fn new(
        quotient_poly: String,
        dim: usize,
        reed_solomon_degree: usize,
        cache_file: Option<String>,
        hgraph_file: PathBuf,
        num_threads: usize,
    ) -> Self {
        Self {
            quotient_poly,
            dim,
            reed_solomon_degree,
            cache_file,
            hgraph_file,
            num_threads,
        }
    }

    pub fn from_disk(path: &Path) -> Self {
        let mut s = String::new();
        let mut file = File::open(path).expect("Could not open file for matrix read.");
        file.read_to_string(&mut s)
            .expect("Could not read matrix from disk");
        serde_json::from_str::<RankEstimatorConfig>(&s).expect("Could not deserialize file.")
    }
}

#[cfg(test)]
mod test {
    use std::{path::PathBuf, str::FromStr};

    use serde::Deserialize;
    use simple_logger::SimpleLogger;

    use crate::{
        math::{
            finite_field::{FFRep, FiniteField},
            polynomial::FFPolynomial,
        },
        matrices::{
            ffmatrix::FFMatrix,
            mat_trait::RankMatrix,
            parallel_matrix::ParallelFFMatrix,
            sparse_ffmatrix::{MemoryLayout, SparseFFMatrix},
        },
    };

    use super::{compute_rank_bounds, RankConfig};

    #[test]
    fn small_example() {
        let _ = SimpleLogger::new()
            .with_level(log::LevelFilter::Warn)
            .init()
            .unwrap();
        let q = FFPolynomial::from_str("1*x^2 + 2*x^ 1 + 2*x^0 % 3").unwrap();
        let mut rc = RankConfig::new(q, 3, 2, 30000, 1000);
        println!("here");
        rc.run("/Users/matt/repos/qec/tmp", 1);
    }

    #[test]
    fn load_galois_test() {
        #[derive(Debug, Clone, Deserialize)]
        struct TestMatrices {
            num_samples: usize,
            n_rows: usize,
            n_cols: usize,
            finite_field: FFRep,
            memory_layout: MemoryLayout,
            input: Vec<u32>,
            output: Vec<u32>,
        }
        let filename = PathBuf::from("/Users/matt/repos/qec/scripts/galois_bench_data.json");
        if let Some(test_data) = std::fs::read_to_string(filename.as_path()).ok() {
            let x: serde_json::Value =
                serde_json::from_str(&test_data[..]).expect("cannot deserialize");
            let test_matrices: TestMatrices = serde_json::from_value(x).unwrap();
            dbg!(test_matrices.num_samples);
            dbg!(test_matrices.n_rows);
            dbg!(test_matrices.n_cols);
            dbg!(test_matrices.finite_field);
            dbg!(test_matrices.memory_layout);
            let matrix_len = test_matrices.n_rows * test_matrices.n_cols;
            for sample_ix in 0..test_matrices.num_samples {
                let input_matrix_raw =
                    &test_matrices.input[sample_ix * matrix_len..(sample_ix + 1) * matrix_len];
                let output_matrix_raw =
                    &test_matrices.output[sample_ix * matrix_len..(sample_ix + 1) * matrix_len];
                assert_eq!(input_matrix_raw.len(), matrix_len);
                assert_eq!(output_matrix_raw.len(), matrix_len);
                let input_parsed_entries = input_matrix_raw
                    .iter()
                    .map(|entry| FiniteField::from((*entry, test_matrices.finite_field)))
                    .collect();
                let output_parsed_entries = output_matrix_raw
                    .iter()
                    .map(|entry| FiniteField::from((*entry, test_matrices.finite_field)))
                    .collect();
                let mut dense_rref = FFMatrix::new(
                    input_parsed_entries,
                    test_matrices.n_rows,
                    test_matrices.n_cols,
                );
                let mut sparse = SparseFFMatrix::from(dense_rref.clone());
                let mut parallel = sparse.clone().split_into_parallel(
                    sparse.row_ixs().into_iter().collect(),
                    std::thread::available_parallelism().unwrap().into(),
                );
                parallel.rref(None, None, None);
                let parallel_rref = parallel.quit().to_dense();
                sparse.rref();
                let sparse_rref = sparse.to_dense();
                dense_rref.rref();
                let galois_rref = FFMatrix::new(
                    output_parsed_entries,
                    test_matrices.n_rows,
                    test_matrices.n_cols,
                );
                assert_eq!(galois_rref, dense_rref);
                assert_eq!(galois_rref, sparse_rref);
                assert_eq!(galois_rref, parallel_rref);
            }
        }
    }

    #[test]
    fn functionized() {
        let _ = SimpleLogger::new().init().unwrap();
        let q = FFPolynomial::from_str("1*x^2 + 2*x^ 1 + 2*x^0 % 3").unwrap();
        let cache_dir = PathBuf::from_str("/Users/matt/repos/qec/tmp/rank/").unwrap();
        dbg!(compute_rank_bounds(
            q,
            3,
            2,
            cache_dir,
            Some(1000),
            None,
            None,
            None
        ),);
    }
}
