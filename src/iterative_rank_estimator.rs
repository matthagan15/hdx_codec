use std::{
    borrow::BorrowMut,
    collections::{HashMap, HashSet},
    path::PathBuf,
    str::FromStr,
    thread::current,
    usize,
};

use fxhash::FxHashMap;
use indexmap::IndexMap;
use mhgl::HyperGraph;
use serde::{Deserialize, Serialize};

use crate::{
    code::Code,
    math::{
        finite_field::{FFRep, FiniteField},
        iterative_bfs_new::GroupBFS,
        polynomial::FFPolynomial,
        quotient_polynomial,
    },
    matrices::sparse_ffmatrix::{MemoryLayout, SparseFFMatrix, SparseVector},
    reed_solomon::ReedSolomon,
};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RankEstimatorConfig {
    quotient_poly: String,
    dim: usize,
    reed_solomon_degree: usize,
    output_dir: String,
}

impl RankEstimatorConfig {
    pub fn new(
        quotient_poly: String,
        dim: usize,
        reed_solomon_degree: usize,
        output_dir: String,
    ) -> Self {
        Self {
            quotient_poly,
            dim,
            reed_solomon_degree,
            output_dir,
        }
    }
}
pub struct IterativeRankEstimator {
    bfs: GroupBFS,
    message_id_to_col_ix: IndexMap<u64, usize>,
    parity_check_to_row_ixs: IndexMap<u64, Vec<usize>>,
    pivotized_local_checks: HashSet<u32>,
    pivots: Vec<(usize, usize)>,
    pub parity_check_matrix: SparseFFMatrix,
    local_code: ReedSolomon,
    dim: usize,
}

impl IterativeRankEstimator {
    pub fn new(conf: RankEstimatorConfig) -> Self {
        let quotient =
            FFPolynomial::from_str(&conf.quotient_poly).expect("Could not parse polynomial.");
        let local_code = ReedSolomon::new(quotient.field_mod, conf.reed_solomon_degree);
        let parity_check_matrix = SparseFFMatrix::new(
            usize::MAX,
            usize::MAX,
            quotient.field_mod,
            MemoryLayout::RowMajor,
        );
        let pathbuf = PathBuf::from(conf.output_dir);
        let bfs = GroupBFS::new(&pathbuf, String::from("temporary"), &quotient);
        Self {
            bfs,
            message_id_to_col_ix: IndexMap::new(),
            parity_check_to_row_ixs: IndexMap::new(),
            pivotized_local_checks: HashSet::new(),
            parity_check_matrix,
            pivots: Vec::new(),
            local_code,
            dim: conf.dim,
        }
    }

    fn get_complete_border_checks(&self, local_check: u32) -> Vec<u64> {
        let mut ret = Vec::new();
        let maximal_faces = self.bfs.hgraph().maximal_edges_of_nodes([local_check]);
        for max_face in maximal_faces {
            let border_nodes: Vec<u32> = self
                .bfs
                .hgraph()
                .query_edge(&max_face)
                .unwrap()
                .into_iter()
                .filter(|node| *node != local_check)
                .collect();
            let border_check = self
                .bfs
                .hgraph()
                .find_id(border_nodes)
                .expect("Border nodes should have edge id.");
            let star: Vec<u32> = self
                .bfs
                .hgraph()
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
        let mut local_checks_to_pivotize = HashSet::new();
        loop {
            let discovered = self.step();
            local_checks_to_pivotize.insert(discovered[0]);
            counter += 1;
            if counter % check_rate == 0 {
                log::warn!("Counter: {counter}");
                let mut local_checks_not_ready = Vec::new();
                for local_check in local_checks_to_pivotize.drain() {
                    if self.is_local_view_complete(local_check) {
                        self.pivotize_interior_checks(local_check);
                        self.pivotized_local_checks.insert(local_check);
                        let border_checks = self.get_complete_border_checks(local_check);
                        if border_checks.is_empty() == false {
                            for border_check in border_checks.into_iter() {
                                log::info!("pivotizing border check: {border_check}");
                                self.pivotize_border_check(border_check);
                            }
                        }
                    } else {
                        local_checks_not_ready.push(local_check);
                    }
                }
                for check in local_checks_not_ready.into_iter() {
                    local_checks_to_pivotize.insert(check);
                }
                let rate =
                    1.0 - (self.pivots.len() as f64 / self.message_id_to_col_ix.len() as f64);
                log::debug!(
                    "pivots, triangles, rate: {:}, {:}, {rate}",
                    self.pivots.len(),
                    self.message_id_to_col_ix.len()
                );
            }
        }
    }

    fn is_parity_check_complete(&self, parity_check: u64) -> bool {
        let maximal_edges = self.bfs.hgraph().maximal_edges(&parity_check);
        maximal_edges.len() == (self.bfs.field_mod() as usize)
    }

    fn is_local_view_complete(&self, local_check: u32) -> bool {
        let containing_edges = self.bfs.hgraph().containing_edges_of_nodes([local_check]);
        let mut is_complete = true;
        for containing_edge in containing_edges {
            let containing_edge_nodes = self.bfs.hgraph().query_edge(&containing_edge).unwrap();
            if containing_edge_nodes.len() == self.dim {
                let boundary_check: Vec<u32> = containing_edge_nodes
                    .iter()
                    .filter(|node| **node != local_check)
                    .cloned()
                    .collect();
                let boundary_check_id = self.bfs.hgraph().find_id(&boundary_check[..]).expect("Boundary check should still be a valid edge in the hgraph. If the traingle contains this local_check node then it must contain it's link within the triangle.");
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

    fn add_parity_check_to_matrix(&mut self, parity_check: u64) {
        let mut maximal_edges = self.bfs.hgraph().maximal_edges(&parity_check);
        maximal_edges.sort();
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
                        FiniteField::new(1, self.bfs.field_mod())
                    } else {
                        FiniteField::new(0, self.bfs.field_mod())
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
        if self.is_parity_check_complete(*parity_check) == false {
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
        let triangle = self.bfs.step();
        let hg = self.bfs.hgraph();
        let e1 = hg
            .find_id([triangle[0], triangle[1]])
            .expect("should have been allocated in bfs.");
        let e2 = hg
            .find_id([triangle[0], triangle[2]])
            .expect("should have been allocated in bfs.");
        let e3 = hg
            .find_id([triangle[1], triangle[2]])
            .expect("should have been allocated in bfs.");
        let triangle_id = hg
            .find_id(&triangle[0..3])
            .expect("should have been allocated in bfs.");
        let message_ix = if self.message_id_to_col_ix.is_empty() {
            0
        } else {
            self.message_id_to_col_ix.last().unwrap().1 + 1
        };
        self.message_id_to_col_ix.insert(triangle_id, message_ix);
        self.parity_check_handler(&e1);
        self.parity_check_handler(&e2);
        self.parity_check_handler(&e3);
        triangle
    }

    pub fn relative_rate_upper_bound(&mut self) -> f64 {
        let mut first_vertex = None;
        for step in 0..usize::MAX {
            let t = self.step();
            if first_vertex.is_none() {
                first_vertex = Some(t[0]);
            }
            if self.is_local_view_complete(first_vertex.unwrap()) {
                println!("First check is done! step number {:}", step);
                return self.pivotize_interior_checks(first_vertex.unwrap());
            }
        }
        panic!("Reached the end of the complex without finding a completed local check?")
    }

    /// Computes a *lower* bound on the (normalized) rate of the error correcting code
    /// from the resulting complex, given a completed local_check.
    pub fn pivotize_interior_checks(&mut self, local_check: u32) -> f64 {
        log::trace!("pivotizing local_check: {:}", local_check);
        let containing_edges = self.bfs.hgraph().containing_edges_of_nodes([local_check]);
        let mut interior_checks = Vec::new();
        let mut border_checks = Vec::new();
        let mut message_ids = Vec::new();
        for edge in containing_edges {
            let nodes = self.bfs.hgraph().query_edge(&edge).unwrap();
            if nodes.len() == 2 {
                interior_checks.push(edge);
            } else if nodes.len() == 3 {
                message_ids.push(edge);
                let border_nodes: Vec<u32> = nodes
                    .into_iter()
                    .filter(|node| *node != local_check)
                    .collect();
                let border_id = self.bfs.hgraph().find_id(border_nodes).unwrap();
                border_checks.push(border_id);
                if self.is_parity_check_complete(border_id) == false {
                    panic!("Border check failed: {:}", border_id);
                }
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
                self.pivots.push((pivot, *message_ix));
                num_pivots_found += 1;
            } else {
                // log::trace!(
                //     "could not find pivot among interior rows for message_ix: {:}",
                //     message_ix
                // );
            }
        }
        log::trace!("Done!");
        let ret = 1.0 - (num_pivots_found as f64) / (message_ids.len() as f64);
        // self.print_border_ixs(border_ixs.clone());
        // self.print_sub_matrix(border_ixs, interior_ixs, message_ids);

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
            self.bfs.field_mod(),
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

                self.pivots.push((row, ix));
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
                let entry = self.parity_check_matrix.query(*row_ix, *message_ix);
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
                let entry = self.parity_check_matrix.query(row_ix, *message_ix);
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
        println!("{:?}", self.bfs.hgraph());
    }
}

mod tests {

    use super::{IterativeRankEstimator, RankEstimatorConfig};

    #[test]
    fn printing() {
        let conf = RankEstimatorConfig {
            quotient_poly: String::from("1*x^2 + 2*x^1 + 2*x^0 % 3"),
            dim: 3,
            reed_solomon_degree: 2,
            output_dir: String::from("/Users/matt/repos/qec/tmp"),
        };
        let mut iterator = IterativeRankEstimator::new(conf);
        let mut first_vertex = None;
        for step in 0..400000 {
            let t = iterator.step();
            if first_vertex.is_none() {
                first_vertex = Some(t[0]);
            }
            if iterator.is_local_view_complete(first_vertex.unwrap()) {
                println!("First check is done! step number {:}", step);
                iterator.pivotize_interior_checks(first_vertex.unwrap());
                return;
            }
        }
        println!("Finished steps before completed?");
    }

    #[test]
    fn rate_entry() {
        let conf = RankEstimatorConfig {
            quotient_poly: String::from("1*x^2 + 2*x^1 + 2*x^0 % 3"),
            dim: 3,
            reed_solomon_degree: 2,
            output_dir: String::from("/Users/matt/repos/qec/tmp"),
        };
        let mut iterator = IterativeRankEstimator::new(conf);
        iterator.compute_rate();
    }
}
