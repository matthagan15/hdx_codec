use std::{borrow::BorrowMut, collections::HashMap, path::PathBuf, str::FromStr};

use indexmap::IndexMap;
use mhgl::HyperGraph;

use crate::math::{
    iterative_bfs_new::GroupBFS, polynomial::FiniteFieldPolynomial, sparse_ffmatrix::SparseVector,
};

struct IterativeRankEstimator {
    bfs: GroupBFS,
    /// for now each sparse vector is going to have the triangle edge index be the
    /// matrix index
    parity_check_to_row: IndexMap<u64, SparseVector>,
}

impl IterativeRankEstimator {
    pub fn new() -> Self {
        let dir = PathBuf::from("/Users/matt/repos/qec/tmp");
        let q = FiniteFieldPolynomial::from_str("1*x^2 + 2 * x^1 + 2 * x^0 % 3").unwrap();
        let bfs = GroupBFS::new(&dir, String::from("temporary"), &q);
        // let rows = bfs.hgraph().edges_of_size(2);
        // let triangle = bfs.hgraph().edges_of_size(3)[0];
        Self {
            bfs,
            parity_check_to_row: IndexMap::new(),
            // parity_check_to_row: rows
            //     .into_iter()
            //     .map(|row| {
            //         let mut sparse_vec = SparseVector::new_empty();
            //         // sparse_vec.insert(triangle as usize, 1);
            //         (row, sparse_vec)
            //     })
            //     .collect(),
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

        for e in [e1, e2, e3].iter() {
            let v = self
                .parity_check_to_row
                .entry(*e)
                .or_insert(SparseVector::new_empty());
            v.insert(triangle_id as usize, 1);
        }
        // for vertex in triangle {
        //     let vertex_view = hg.maximal_edges_of_nodes([vertex]);
        //     if vertex_view.len() == self.bfs.field_mod().pow(3) as usize {
        //         // do the row reduction of each row in this view
        //         // get each edge that also has this vertex
        //         let mut containing_lines: Vec<_> = hg
        //             .containing_edges_of_nodes([vertex])
        //             .into_iter()
        //             .filter(|e| hg.query_edge(e).unwrap().len() == 2)
        //             .collect();
        //         containing_lines.sort();
        //         self.rref(containing_lines);
        //     }
        // }
    }

    fn rref(&mut self, checks: Vec<u64>) {}

    fn print_hgraph(&self) {
        println!("{:}", self.bfs.hgraph());
    }

    fn print_matrix(&self) {
        let cols = self.bfs.hgraph().edges_of_size(3);
        let mut rows: Vec<u64> = self
            .parity_check_to_row
            .keys()
            .into_iter()
            .cloned()
            .collect();
        rows.sort();
        let mut row_prefix = 0;
        for row in rows.iter() {
            let row_string = row.to_string();
            row_prefix = row_prefix.max(row_string.len());
        }
        let mut s = String::new();
        let mut col_width = 0;
        s.push_str(&" ".repeat(row_prefix)[..]);
        s.push_str(" | ");
        for col in cols.iter() {
            let col_string = col.to_string();
            col_width = col_width.max(col_string.len());
        }
        for col in cols.iter() {
            let mut col_string = col.to_string();
            col_string.push_str(&" ".repeat(col_width - col_string.len())[..]);
            s.push_str(&col_string[..]);
            s.push_str(" | ");
        }
        s.pop();
        let row_len = s.len();
        s.push('\n');
        s.push_str(&"_".repeat(row_len)[..]);
        s.push('\n');
        for row in rows {
            let mut row_string = row.to_string();
            let pad = row_prefix - row_string.len();
            row_string.push_str(&" ".repeat(pad)[..]);
            s.push_str(&row_string[..]);
            s.push_str(" | ");
            for col in cols.clone() {
                let v = self.parity_check_to_row.get(&row).unwrap();
                let entry = v.query(&(col as usize));
                let mut col_string = entry.to_string();
                col_string.push_str(&" ".repeat(col_width - col_string.len())[..]);
                s.push_str(&col_string[..]);
                s.push_str(" | ");
            }
            s.pop();
            s.push('\n');
        }
        s.push_str(&"_".repeat(row_len)[..]);
        println!("{:}", s);
    }
}

mod tests {
    use super::IterativeRankEstimator;

    #[test]
    fn printing() {
        let mut iterator = IterativeRankEstimator::new();
        for _ in 0..30 {
            iterator.step();
        }
        iterator.print_hgraph();
        iterator.print_matrix();
    }
}
