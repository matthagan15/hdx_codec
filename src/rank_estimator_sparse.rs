use fxhash::{FxHashMap, FxHashSet};
use indexmap::IndexMap;
use mhgl::{HGraph, HyperGraph};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::io::Write;
use std::path::Path;
use std::rc::Rc;
use std::{collections::HashSet, fs::File, path::PathBuf, str::FromStr, time::Instant};

use crate::math::finite_field::FFRep;
use crate::matrices::mat_trait::RankMatrix;
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

pub struct IterativeRankEstimator<M: RankMatrix> {
    hgraph: HGraph<u16, ()>,
    message_id_to_col_ix: IndexMap<u64, usize>,
    parity_check_to_row_ixs: IndexMap<u64, Vec<usize>>,
    column_ix_to_pivot_row: BTreeMap<usize, usize>,
    parity_check_to_col_ixs: FxHashMap<u64, Vec<usize>>,
    pub parity_check_matrix: M,
    local_code: ReedSolomon,
    rate_upper_bound: Option<f64>,
    dim: usize,
    field_mod: FFRep,
}

impl<M: RankMatrix> IterativeRankEstimator<M> {
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
            parity_check_matrix: M::new(quotient.field_mod),
            column_ix_to_pivot_row: BTreeMap::new(),
            parity_check_to_col_ixs: FxHashMap::default(),
            local_code,
            rate_upper_bound: None,
            dim: conf.dim,
            field_mod: quotient.field_mod,
        }
    }

    pub fn compute_rate(&mut self) {
        let rate_start = Instant::now();
        let mut tot_node_time = 0.0;
        let mut count = 0;
        let nodes = self.hgraph.nodes();
        let num_nodes = nodes.len();
        log::trace!("Found {:} nodes to check.", nodes.len());
        for node in nodes {
            let node_start = rate_start.elapsed();
            self.node_step(node);
            count += 1;
            let time_for_node = (rate_start.elapsed() - node_start).as_secs_f64();
            tot_node_time += time_for_node;
            let time_remaining = time_for_node * (num_nodes as f64 - count as f64);
            if count % 1000 == 0 {
                log::trace!("Count: {count}");
                log::trace!(
                    "rate: {:.4}, estimated time remaining: {:}",
                    1.0 - self.column_ix_to_pivot_row.len() as f64
                        / self.message_id_to_col_ix.len() as f64,
                    time_remaining,
                );
                log::trace!(
                    "num pivots found: {:}, num triangles: {:}",
                    self.column_ix_to_pivot_row.len(),
                    self.message_id_to_col_ix.len()
                );
            }
        }
        println!("{:}", "#".repeat(80));
        println!("{:}BORDER CHECKS{:}", " ".repeat(34), " ".repeat(34));
        println!("{:}", "#".repeat(80));
        let mut count = 0;
        let mut num_pivots_found = 0;
        let num_border_checks = self.parity_check_to_col_ixs.len();
        let mut tot_time_border_checks = 0.0;
        let border_ixs = self
            .parity_check_to_col_ixs
            .iter()
            .filter_map(|(check, col_ixs)| self.parity_check_to_row_ixs.get(check))
            .fold(HashSet::<usize>::new(), |mut acc, v| {
                for ix in v {
                    acc.insert(*ix);
                }
                acc
            });
        let mut border_parity_check = self
            .parity_check_matrix
            .split_into_parallel(border_ixs.clone(), 2);
        for border_ix in border_ixs.iter() {
            let border_start = Instant::now();
            let col = border_parity_check.pivotize_row(*border_ix);
            if let Some(col) = col {
                self.column_ix_to_pivot_row.insert(col, *border_ix);
                num_pivots_found += 1;
            }
            count += 1;
            let current_border_time_taken = border_start.elapsed().as_secs_f64();
            if count % 400 == 0 {
                let time_per_border_check = tot_time_border_checks / count as f64;
                let time_remaining =
                    time_per_border_check * (num_border_checks as f64 - count as f64);
                let possible_ixs_left = border_ixs.len() - count;
                let most_pivots_possible = self.column_ix_to_pivot_row.len() + possible_ixs_left;
                let worst_case_rate =
                    1.0 - (most_pivots_possible as f64 / self.message_id_to_col_ix.len() as f64);
                let best_case_rate = 1.0
                    - (self.column_ix_to_pivot_row.len() as f64
                        / self.message_id_to_col_ix.len() as f64);
                log::trace!(
                    "avg_border check: {:}, Estimated time remaining: {:}",
                    time_per_border_check,
                    time_remaining
                );
                log::trace!(
                    "time spent on this border check: {:}",
                    current_border_time_taken
                );
                log::trace!("borders processed: {:}", count);
                log::trace!("borders left: {:}", num_border_checks - count);
                log::trace!("rate best case scenario: {:}", best_case_rate);
                log::trace!("rate worst case scenario: {:}", worst_case_rate);
            }

            tot_time_border_checks += current_border_time_taken;
        }
        let out = border_parity_check.quit();

        log::trace!(
            "Parity check dims: {:}, {:}",
            self.parity_check_matrix.n_rows(),
            self.parity_check_matrix.n_cols()
        );
        log::trace!(
            "final rate: {:}",
            1.0 - self.column_ix_to_pivot_row.len() as f64 / self.message_id_to_col_ix.len() as f64
        )
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

        self.add_parity_check_to_matrix(*parity_check);
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
        let node_start = Instant::now();
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
                if e.len() == self.local_code.encoded_len() {
                    let border_check_complete = border_checks.get_mut(&border_check).unwrap();
                    *border_check_complete = true;
                }
            }
        }

        if message_ixs.len() != self.field_mod.pow(3) as usize {
            return;
        }

        // Now its time to pivotize the interior checks
        for interior_check in interior_checks.iter() {
            self.parity_check_handler(interior_check);
        }
        for border_check in border_checks
            .iter()
            .filter_map(|(k, v)| if *v { Some(k) } else { None })
        {
            self.parity_check_handler(border_check);
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
            let col_ixs = self.parity_check_to_col_ixs.get(check).unwrap();
            for col_ix in col_ixs.iter() {
                if let Some(pivot_row) = self.column_ix_to_pivot_row.get(col_ix) {
                    self.parity_check_matrix
                        .eliminate_rows((*pivot_row, *col_ix), border_ixs.clone());
                }
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
                entries.push((new_row_ix, new_col_ix, entry));
                new_col_ix += 1;
            }
            new_row_ix += 1;
        }
        let mut new_col_ix: usize = 0;
        for _ in &message_ids {
            entries.push((new_row_ix, new_col_ix, self.field_mod - 1));
            new_col_ix += 1;
        }
        new_row_ix += 1;
        for row_ix in border_ixs {
            let mut new_col_ix: usize = 0;
            for message_id in &message_ids {
                let message_ix = self.message_id_to_col_ix.get(message_id).unwrap();
                let entry = self.parity_check_matrix.get(row_ix, *message_ix);
                entries.push((new_row_ix, new_col_ix, entry));
                new_col_ix += 1;
            }
            new_row_ix += 1;
        }
        let sub_matrix = SparseFFMatrix::new_with_entries(
            new_row_ix - 1,
            message_ids.len(),
            self.field_mod,
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
