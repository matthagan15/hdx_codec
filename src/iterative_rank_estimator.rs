use std::{borrow::BorrowMut, collections::HashMap, path::PathBuf, str::FromStr, thread::current};

use indexmap::IndexMap;
use mhgl::HyperGraph;
use serde::{Deserialize, Serialize};

use crate::{
    code::Code,
    math::{
        finite_field::{FFRep, FiniteField},
        iterative_bfs_new::GroupBFS,
        polynomial::FiniteFieldPolynomial,
        sparse_ffmatrix::{MemoryLayout, SparseFFMatrix, SparseVector},
    },
    reed_solomon::ReedSolomon,
};

#[derive(Debug, Clone, Serialize, Deserialize)]
struct RankEstimatorConfig {
    quotient_poly: String,
    dim: usize,
    reed_solomon_degree: usize,
    output_dir: String,
}

struct IterativeRankEstimator {
    bfs: GroupBFS,
    message_id_to_col_ix: IndexMap<u64, usize>,
    parity_check_to_row_ixs: IndexMap<u64, Vec<usize>>,
    pub parity_check_matrix: SparseFFMatrix,
    local_code: ReedSolomon,
    dim: usize,
}

impl IterativeRankEstimator {
    pub fn new(conf: RankEstimatorConfig) -> Self {
        let quotient = FiniteFieldPolynomial::from_str(&conf.quotient_poly)
            .expect("Could not parse polynomial.");
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
            parity_check_matrix,
            local_code,
            dim: conf.dim,
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
        self.parity_check_handler(&e2);
        triangle
    }

    pub fn attempt_block_diagonalize(&mut self, local_check: u32) {
        println!("local_check: {:}", local_check);
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
        for message_id in message_ids {
            if let Some(pivot) = self.parity_check_matrix.find_nonzero_entry_among_rows(
                *self.message_id_to_col_ix.get(&message_id).unwrap(),
                interior_ixs.clone(),
            ) {}
        }
    }

    fn print_hgraph(&self) {
        println!("{:}", self.bfs.hgraph());
    }
}

mod tests {
    use std::iter;

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
        for step in 0..100000 {
            let t = iterator.step();
            if first_vertex.is_none() {
                first_vertex = Some(t[0]);
            }
            if iterator.is_local_view_complete(first_vertex.unwrap()) {
                println!("First check is done! step number {:}", step);
                iterator.attempt_block_diagonalize(first_vertex.unwrap());
            }
        }
    }
}
