use fxhash::{FxHashMap, FxHashSet};
use indexmap::IndexMap;
use mhgl::{HGraph, HyperGraph};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::fs;
use std::io::{Read, Write};
use std::path::Path;
use std::rc::Rc;
use std::{collections::HashSet, fs::File, path::PathBuf, str::FromStr, time::Instant};

use crate::math::finite_field::FFRep;
use crate::matrices::ffmatrix::FFMatrix;
use crate::matrices::mat_trait::RankMatrix;
use crate::matrices::sparse_ffmatrix::ParallelFFMatrix;
use crate::{
    code::Code,
    math::{coset_complex_bfs::GroupBFS, finite_field::FiniteField, polynomial::FFPolynomial},
    matrices::sparse_ffmatrix::{MemoryLayout, SparseFFMatrix},
    reed_solomon::ReedSolomon,
};

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
        let pivots_found = border_check_matrix.rank_via_upper_triangular();
        let tot_pivots = pivots_found + self.column_ix_to_pivot_row.len();
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
