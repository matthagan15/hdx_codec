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
    triangle_id_to_col_ix: IndexMap<u64, u64>,
    parity_check_to_row_ixs: IndexMap<u64, Vec<usize>>,
    pub parity_check_matrix: SparseFFMatrix,
    local_code: ReedSolomon,
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
            triangle_id_to_col_ix: IndexMap::new(),
            parity_check_to_row_ixs: IndexMap::new(),
            parity_check_matrix,
            local_code,
        }
    }

    fn is_parity_check_complete(&self, parity_check: u64) -> bool {
        let maximal_edges = self.bfs.hgraph().maximal_edges(&parity_check);
        maximal_edges.len() == (self.bfs.field_mod() as usize)
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
            println!("basis_vec: {:}", basis_vec.len());
            let parity_check_output = self.local_code.parity_check(&basis_vec);
            if row_indices_of_outputs.len() != parity_check_output.len() {
                panic!("Number of rows assigned to parity check does not match output length of parity check.")
            }
            for ix in 0..row_indices_of_outputs.len() {
                self.parity_check_matrix.insert(
                    row_indices_of_outputs[ix] as usize,
                    *self.triangle_id_to_col_ix.get(message_id).unwrap() as usize,
                    parity_check_output[ix].0,
                );
            }
        }
    }

    pub fn step(&mut self) {
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
        let message_ix = if self.triangle_id_to_col_ix.is_empty() {
            0
        } else {
            self.triangle_id_to_col_ix.last().unwrap().1 + 1
        };
        self.triangle_id_to_col_ix.insert(triangle_id, message_ix);
        let mut handle_parity_check = |&parity_check| {
            if self.is_parity_check_complete(parity_check) {
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
                let new_row_indices: Vec<usize> = (current_ix_start
                    ..(current_ix_start + self.local_code.parity_check_len()))
                    .collect();
                self.parity_check_to_row_ixs
                    .insert(parity_check, new_row_indices);
                self.add_parity_check_to_matrix(parity_check);
            }
        };
        handle_parity_check(&e1);
        handle_parity_check(&e2);
        handle_parity_check(&e2);
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
        for _ in 0..30 {
            iterator.step();
        }
        iterator.print_hgraph();
        println!("{:}", iterator.parity_check_matrix.to_dense());
    }
}
