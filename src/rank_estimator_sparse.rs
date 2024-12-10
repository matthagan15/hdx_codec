use indexmap::IndexMap;
use mhgl::{HGraph, HyperGraph};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::fs;
use std::io::Read;
use std::path::Path;
use std::thread::available_parallelism;
use std::{collections::HashSet, fs::File, path::PathBuf, str::FromStr};

use crate::hdx_code;
use crate::math::coset_complex_bfs::bfs;
use crate::math::finite_field::FFRep;
use crate::matrices::mat_trait::RankMatrix;
use crate::{
    code::Code,
    math::{finite_field::FiniteField, polynomial::FFPolynomial},
    matrices::sparse_ffmatrix::{MemoryLayout, SparseFFMatrix},
    reed_solomon::ReedSolomon,
};

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
struct RankConfig {
    quotient_poly: FFPolynomial,
    dim: usize,
    rs_degree: usize,
    truncation: usize,
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
) {
    // First check if there is a config in the directory:
    // If config present - make sure parameters match. Then try to load cache
    // If config not present - we are starting from scratch.
    let mut config_file = cache_dir.clone();
    config_file.push("config.json");
    let coset_complex_size =
        crate::math::coset_complex_bfs::size_of_coset_complex(&quotient_poly, dim);
    let input_config = RankConfig {
        quotient_poly: quotient_poly.clone(),
        dim,
        rs_degree,
        truncation: truncation.unwrap_or(coset_complex_size + 1),
    };
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
                    std::fs::write(
                        config_file,
                        serde_json::to_string_pretty(&input_config)
                            .expect("Could not serialize config."),
                    )
                    .expect("Could not write config to disk.");
                } else {
                    log::error!("Mismatch between cached config and input config: bailing.");
                    panic!()
                }
            }
        }
    } else {
        std::fs::write(
            config_file,
            serde_json::to_string_pretty(&input_config).expect("Could not serialize config."),
        )
        .expect("Could not write config to disk.");
    }
    let local_rs_code = ReedSolomon::new(quotient_poly.field_mod, rs_degree);
    /// A valid cache has already had the interior matrices pivotiized properly and the border
    /// matrix reduced as much as possible.
    #[derive(Serialize, Deserialize)]
    struct MatrixCache {
        interior_mat: SparseFFMatrix,
        border_mat: SparseFFMatrix,
        num_interior_pivots: usize,
    }
    let mut matrix_cache_file = cache_dir.clone();
    matrix_cache_file.push("matrix_cache");
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
        let (hg, _new_edges) = bfs(
            quotient_poly.clone(),
            dim,
            truncation,
            Some(hgraph_cache_file),
            None,
        );
        if checks.is_none() {
            let splitted_checks = split_hg_into_checks(&hg, &local_rs_code);
            checks = Some(Checks {
                local: splitted_checks.0,
                interior: splitted_checks.1,
                border: splitted_checks.2,
                message_ids: splitted_checks.3,
            });
        };
        let checks = checks.unwrap();
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
        // Cache the newly created matrices.
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
    drop(interior);
    let mut parallel_border_mat = border.split_into_parallel(
        border.row_ixs().into_iter().collect(),
        available_parallelism().unwrap().into(),
    );
    let pivots = parallel_border_mat.row_echelon_form();
    log::trace!("Found {:} pivots for the border matrix.", pivots.len());
    let reduced_border = parallel_border_mat.quit();
    log::info!(
        "Final computed rate: {:}",
        1.0 - (pivots.len() + num_interior_pivots) as f64 / reduced_border.n_cols as f64
    );
}

/// Scans a hgraph for local checks, interior checks, and border checks. Needs reference to
/// the polynomial to determine when a check is complete.
/// Returns (local_checks, interior_checks, border_checks, message_ids)
fn split_hg_into_checks(
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
/// - each parity check, whether interior or border, has a complete local view for the provided
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
            if let Some(_) = interior_matrix.pivotize_row_within_range(row_ix, &rows_to_pivot[..]) {
                num_interior_pivots += 1;
                let pivot_row = interior_matrix.row(row_ix);
                border_matrix.eliminate_col_with_range(&pivot_row, &border_rows_to_eliminate[..]);
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

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IterativeRankEstimator {
    hgraph_file: PathBuf,
    message_id_to_col_ix: IndexMap<u64, usize>,
    parity_check_to_row_ixs: IndexMap<u64, Vec<usize>>,
    column_ix_to_pivot_row: BTreeMap<usize, usize>,
    border_checks: Vec<u64>,
    processed_border_checks: Vec<u64>,
    pub parity_check_matrix: SparseFFMatrix,
    local_code: ReedSolomon,
    rate_upper_bound: Option<f64>,
    num_threads: usize,
    cache_file: Option<PathBuf>,
    field_mod: FFRep,
}

impl IterativeRankEstimator {
    pub fn new(conf: RankEstimatorConfig) -> Self {
        let quotient =
            FFPolynomial::from_str(&conf.quotient_poly).expect("Could not parse polynomial.");
        log::trace!("Loaded polynomial: {:}", quotient);
        let local_code = ReedSolomon::new(quotient.field_mod, conf.reed_solomon_degree);
        local_code.print();
        let cache_file = conf.cache_file.map(|s| PathBuf::from_str(&s[..]).unwrap());
        log::trace!(
            "Attempting to read hgraph from: {:}",
            &conf.hgraph_file.to_str().unwrap()
        );

        // TODO: Need to assert that the local_code message_len is the same as what we
        // expect from the coset complex parameters.

        Self {
            hgraph_file: conf.hgraph_file,
            message_id_to_col_ix: IndexMap::new(),
            parity_check_to_row_ixs: IndexMap::new(),
            parity_check_matrix: SparseFFMatrix::new(
                0,
                0,
                quotient.field_mod,
                MemoryLayout::RowMajor,
            ),
            column_ix_to_pivot_row: BTreeMap::new(),
            border_checks: Vec::new(),
            processed_border_checks: Vec::new(),
            local_code,
            rate_upper_bound: None,
            num_threads: conf.num_threads,
            cache_file,
            field_mod: quotient.field_mod,
        }
    }

    fn cache(&self) {
        if self.cache_file.is_none() {
            log::trace!("No cache file provided so skipping caching.");
            return;
        }
        let s = serde_json::to_string(self).expect("Could not serialize self.");
        let path = self.cache_file.clone().unwrap();
        fs::write(path.as_path(), s).expect("Could not write cache");
        log::trace!("Successfully cached!")
    }

    /// Currently panics if a rank estimator cannot be constructed.
    pub fn load_from_cache(cache_file: PathBuf) -> Option<Self> {
        if cache_file.is_file() {
            // todo: serialize this to JsonObject
            let file_data = fs::read_to_string(&cache_file).expect("Could not read file.");
            log::trace!("Successfully loaded from cache!");
            return Some(serde_json::from_str(&file_data).expect("could not parse serde"));
        }
        log::trace!("Could not load cache.");
        None
    }

    fn initial_pass(&mut self) {
        let mut count = 0;
        let hgraph: HGraph<u16, ()> =
            HGraph::from_file(&self.hgraph_file).expect("Could not read hgraph from file.");
        let nodes = hgraph.nodes();
        log::trace!("Found {:} nodes to check.", nodes.len());
        let mut border_check_to_cols: HashMap<u64, Vec<usize>> = HashMap::new();
        for node in nodes {
            // first need to make sure that each row in the nodes view
            // has matrix rows filled out.
            let type_ix = hgraph.get_node(&node).unwrap();
            if *type_ix != 0 {
                continue;
            }
            let star = hgraph.containing_edges_of_nodes([node]);
            let mut interior_checks = Vec::new();
            let mut message_ixs = Vec::new();
            let mut border_checks = HashMap::new();
            for edge in star {
                let edge_nodes = hgraph.query_edge(&edge).unwrap();
                if edge_nodes.len() == 2 {
                    interior_checks.push(edge);
                } else if edge_nodes.len() == 3 {
                    // make sure the message_id has a column index
                    let message_ix = if self.message_id_to_col_ix.contains_key(&edge) == false {
                        let message_ix = if self.message_id_to_col_ix.is_empty() {
                            0
                        } else {
                            self.message_id_to_col_ix.last().unwrap().1 + 1
                        };
                        self.message_id_to_col_ix.insert(edge, message_ix);
                        message_ixs.push(message_ix);
                        message_ix
                    } else {
                        let ix = self.message_id_to_col_ix.get(&edge).unwrap();
                        message_ixs.push(*ix);
                        *ix
                    };
                    // get the border check
                    let border_nodes: Vec<u32> = edge_nodes
                        .into_iter()
                        .filter(|node| *hgraph.get_node(&node).unwrap() != 0)
                        .collect();
                    if border_nodes.len() != 2 {
                        panic!("Border nodes should have length 2");
                    }
                    let border_check = hgraph.find_id(&border_nodes[..]).unwrap();
                    border_checks.insert(border_check, false);
                    let e = border_check_to_cols.entry(border_check).or_default();
                    e.push(message_ix);
                    if e.len() == self.local_code.encoded_len() {
                        let border_check_complete = border_checks.get_mut(&border_check).unwrap();
                        *border_check_complete = true;
                    }
                }
            }

            if message_ixs.len() != self.field_mod.pow(3) as usize {
                continue;
            }

            // Now its time to pivotize the interior checks
            for interior_check in interior_checks.iter() {
                self.parity_check_handler(interior_check);
                let mut maximal_edges = hgraph.maximal_edges(interior_check);
                maximal_edges.sort();
                self.add_parity_check_to_matrix(*interior_check, maximal_edges);
            }
            for border_check in border_checks
                .iter()
                .filter_map(|(k, v)| if *v { Some(k) } else { None })
            {
                self.parity_check_handler(border_check);
                let mut maximal_edges = hgraph.maximal_edges(border_check);
                maximal_edges.sort();
                self.add_parity_check_to_matrix(*border_check, maximal_edges);
                self.border_checks.push(*border_check);
            }

            let mut interior_ixs = Vec::new();
            for check in interior_checks {
                let mut ixs = self.parity_check_to_row_ixs.get(&check).unwrap().clone();
                interior_ixs.append(&mut ixs);
            }
            let mut pivots_found = 0;
            for ix in interior_ixs.iter() {
                if let Some(col) = self
                    .parity_check_matrix
                    .pivotize_row_within_range(*ix, &interior_ixs[..])
                {
                    self.column_ix_to_pivot_row.insert(col, *ix);
                    pivots_found += 1;
                }
            }

            if self.rate_upper_bound.is_none() {
                self.rate_upper_bound = Some(1.0 - pivots_found as f64 / interior_ixs.len() as f64);
            }

            for check in border_checks
                .iter()
                .filter_map(|(k, v)| if *v { Some(k) } else { None })
            {
                let border_ixs = self.parity_check_to_row_ixs.get(check).unwrap();
                let col_ixs = border_check_to_cols.get(check).unwrap();
                for col_ix in col_ixs.iter() {
                    if let Some(pivot_row) = self.column_ix_to_pivot_row.get(col_ix) {
                        self.parity_check_matrix
                            .eliminate_rows((*pivot_row, *col_ix), border_ixs.clone());
                    }
                }
            }
            count += 1;
            if count % 1000 == 0 {
                log::trace!("Count: {count}");
                log::trace!(
                    "rate: {:.4}",
                    1.0 - self.column_ix_to_pivot_row.len() as f64
                        / self.message_id_to_col_ix.len() as f64,
                );
            }
        }
        log::trace!("Initial pass complete. Processed {:} nodes, found {:} pivots among interior checks, and have {:} border checks to process.", count, self.column_ix_to_pivot_row.len(), self.border_checks.len());
        log::trace!("Caching re-ordered border check matrix.");
        let border_ixs = self
            .border_checks
            .iter()
            .filter_map(|check| self.parity_check_to_row_ixs.get(check))
            .fold(HashSet::<usize>::new(), |mut acc, v| {
                for ix in v {
                    acc.insert(*ix);
                }
                acc
            });
        let mut border_matrix = self.parity_check_matrix.split(border_ixs.clone());
        border_matrix.reorder_row_and_col_indices();
        let mut cache_file = PathBuf::from(
            self.hgraph_file
                .parent()
                .expect("Hgraph file has no parent."),
        );
        cache_file.push("border_matrix.cache");
        log::trace!("writing border matrix to: {:?}", &cache_file);
        border_matrix
            .to_disk(&cache_file)
            .expect("Could not write border check matrix to disk.");
        self.cache_file = Some(cache_file);
    }

    pub fn compute_rate(&mut self) {
        if self.cache_file.is_none() {
            self.initial_pass();
        }
        let mut border_check_matrix = SparseFFMatrix::from_disk(&self.cache_file.as_ref().unwrap());

        println!("{:}", "#".repeat(80));
        println!("{:}BORDER CHECKS{:}", " ".repeat(34), " ".repeat(34));
        println!("{:}", "#".repeat(80));
        let mut count = 0;
        let mut num_pivots_found = 0;
        let num_border_checks = self.border_checks.len();
        let mut tot_time_border_checks = 0.0;
        let pivots_found = border_check_matrix.row_echelon_form();
        let tot_pivots = pivots_found.len() + self.column_ix_to_pivot_row.len();
        log::trace!(
            "final rate: {:}",
            1.0 - tot_pivots as f64 / self.message_id_to_col_ix.len() as f64
        );
    }

    /// Gets the containing data edges for this parity check first.
    ///
    /// **ASSUMES MAXIMAL EDGES ARE SORTED!**
    fn add_parity_check_to_matrix(&mut self, parity_check: u64, maximal_edges: Vec<u64>) {
        // let mut maximal_edges = self.hgraph.maximal_edges(&parity_check);
        // maximal_edges.sort();
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
                        FiniteField::new(1, self.field_mod)
                    } else {
                        FiniteField::new(0, self.field_mod)
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
    }
}

#[cfg(test)]
mod test {
    use std::{path::PathBuf, str::FromStr};

    use simple_logger::SimpleLogger;

    use crate::math::polynomial::FFPolynomial;

    use super::compute_rank_bounds;

    #[test]
    fn functionized() {
        let _ = SimpleLogger::new().init().unwrap();
        let q = FFPolynomial::from_str("1*x^2 + 2*x^ 1 + 2*x^0 % 3").unwrap();
        let cache_dir = PathBuf::from_str("/Users/matt/repos/qec/tmp/rank/").unwrap();
        dbg!(compute_rank_bounds(q, 3, 2, cache_dir, Some(10000)));
    }
}
