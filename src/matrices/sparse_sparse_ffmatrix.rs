use fxhash::FxHashMap;
use mhgl::{HGraph, HyperGraph};
use serde::{Deserialize, Serialize};

use super::sparse_ffmatrix::{MemoryLayout, SparseFFMatrix};
use crate::{
    math::finite_field::{FFRep, FiniteField as FF},
    matrices::sparse_vec::SparseVector,
};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SparseSparseFFMatrix {
    row_nodes: FxHashMap<usize, u32>,
    col_nodes: FxHashMap<usize, u32>,
    hgraph: HGraph<usize, FFRep>,
    pub field_mod: FFRep,
}

impl SparseSparseFFMatrix {
    pub fn new(field_mod: FFRep) -> Self {
        Self {
            row_nodes: FxHashMap::default(),
            col_nodes: FxHashMap::default(),
            hgraph: HGraph::new(),
            field_mod,
        }
    }
    pub fn new_with_entries(field_mod: FFRep, entries: Vec<(usize, usize, FFRep)>) -> Self {
        let mut new = SparseSparseFFMatrix::new(field_mod);
        new.insert_entries(entries);
        new
    }

    pub fn to_sparse(self) -> SparseFFMatrix {
        let mut entries = Vec::new();
        for (row_ix, row_node) in self.row_nodes {
            let containing_edges = self.hgraph.containing_edges_of_nodes([row_node]);
            for edge in containing_edges {
                let entry = self.hgraph.get_edge(&edge).unwrap();
                let nodes = self.hgraph.query_edge(&edge).unwrap();
                assert_eq!(nodes.len(), 2);
                let col_node = if nodes[0] == row_node {
                    nodes[1]
                } else {
                    nodes[0]
                };
                let col_ix = self.hgraph.get_node(&col_node).unwrap();
                entries.push((row_ix, *col_ix, *entry));
            }
        }
        SparseFFMatrix::new_with_entries(0, 0, self.field_mod, MemoryLayout::RowMajor, entries)
    }

    /// Adds the entry to (row_ix, col_ix) in the matrix
    pub fn insert(&mut self, row_ix: usize, col_ix: usize, entry: FFRep) {
        let r = self
            .row_nodes
            .entry(row_ix)
            .or_insert_with(|| self.hgraph.add_node(row_ix));
        let c = self
            .col_nodes
            .entry(col_ix)
            .or_insert_with(|| self.hgraph.add_node(col_ix));
        let id = self
            .hgraph
            .find_id([*r, *c])
            .or_else(|| Some(self.hgraph.add_edge([*r, *c], 0)))
            .unwrap();

        let edge_mut = self.hgraph.get_edge_mut(&id).unwrap();
        *edge_mut += entry;
        *edge_mut %= self.field_mod;
        if *edge_mut == 0 {
            self.hgraph.remove_edge(id);
            if self.hgraph.containing_edges_of_nodes([*r]).len() == 0 {
                self.hgraph.remove_node(*r);
            }
        }
        if self.hgraph.containing_edges_of_nodes([*c]).len() == 0 {
            self.hgraph.remove_node(*c);
        }
    }

    /// Adds a collection of (`row_ix`, `col_ix`, `entry`) to the matrix.
    pub fn insert_entries(&mut self, entries: Vec<(usize, usize, FFRep)>) {
        for (row_ix, col_ix, entry) in entries.into_iter() {
            self.insert(row_ix, col_ix, entry);
        }
    }

    pub fn get(&self, row_ix: usize, col_ix: usize) -> FF {
        let r = self.row_nodes.get(&row_ix);
        let c = self.col_nodes.get(&col_ix);
        if r.is_none() || c.is_none() {
            return FF::new(0, self.field_mod);
        }
        let id = self.hgraph.find_id([*r.unwrap(), *c.unwrap()]);
        if id.is_none() {
            return FF::new(0, self.field_mod);
        }
        FF::new(*self.hgraph.get_edge(&id.unwrap()).unwrap(), self.field_mod)
    }

    pub fn get_row(&self, row_ix: &usize) -> Vec<(usize, FFRep)> {
        if self.row_nodes.contains_key(row_ix) == false {
            return Vec::new();
        }
        let row_node = self.row_nodes.get(row_ix).unwrap();
        self.hgraph
            .link_of_nodes([*row_node])
            .into_iter()
            .map(|(edge_id, col_node)| {
                let col_ix = *self.hgraph.get_node(&col_node[0]).unwrap();
                let entry = *self.hgraph.get_edge(&edge_id).unwrap();
                (col_ix, entry)
            })
            .collect()
    }

    pub fn get_col(&self, col_ix: usize) -> SparseVector {
        if self.col_nodes.contains_key(&col_ix) == false {
            return SparseVector::new_empty();
        }
        let col_node = self.col_nodes.get(&col_ix).unwrap();
        let entries = self
            .hgraph
            .link_of_nodes([*col_node])
            .into_iter()
            .map(|(edge_id, row_node)| {
                let row_ix = *self.hgraph.get_node(&row_node[0]).unwrap();
                let entry = *self.hgraph.get_edge(&edge_id).unwrap();
                (row_ix, entry)
            })
            .collect();
        SparseVector::new_with_entries(entries)
    }

    /// TODO: Could speed this up by not checking for non-zero entries, but keeping for correctness
    /// for now.
    pub fn get_col_weight(&self, col_ix: usize) -> usize {
        if self.col_nodes.contains_key(&col_ix) == false {
            return 0;
        }
        let col_node = self.col_nodes.get(&col_ix).unwrap();
        self.hgraph
            .link_of_nodes([*col_node])
            .iter()
            .filter(|(edge_id, _row_node)| *self.hgraph.get_edge(edge_id).unwrap() != 0)
            .count()
    }

    /// Finds the smallest row that satisfies the provided filter.
    pub fn find_smallest_nonzero_row(
        &self,
        col_ix: usize,
        filter: impl Fn(usize) -> bool,
    ) -> Option<usize> {
        let col_node = self.col_nodes.get(&col_ix);
        if col_node.is_none() {
            return None;
        }
        let col_node = *col_node.unwrap();
        let link = self.hgraph.link_of_nodes([col_node]);
        link.into_iter().fold(None, |acc: Option<usize>, x| {
            let node_data = self.hgraph.get_node(&x.1[0]).unwrap();
            if filter(*node_data) {
                if acc.is_none() {
                    Some(*node_data)
                } else {
                    Some(acc.unwrap().min(*node_data))
                }
            } else {
                acc
            }
        })
    }

    pub fn add_scaled_row_to_target(
        &mut self,
        source_row: usize,
        target_row: usize,
        scalar: FFRep,
    ) {
        let source_row_node = self.row_nodes.get(&source_row);
        if source_row_node.is_none() {
            return;
        }
        let source_row_node = source_row_node.unwrap();
        let link = self.hgraph.link_of_nodes([*source_row_node]);

        let t_row = self
            .row_nodes
            .entry(target_row)
            .or_insert_with(|| self.hgraph.add_node(target_row));
        let mut edges_to_remove = Vec::new();
        for (edge_id, col_node_vec) in link {
            let source_entry = *self.hgraph.get_edge(&edge_id).unwrap();
            let col_node = col_node_vec[0];
            let target_edge_id = self
                .hgraph
                .find_id([*t_row, col_node])
                .or_else(|| Some(self.hgraph.add_edge([*t_row, col_node], 0)))
                .unwrap();
            let e = self.hgraph.get_edge_mut(&target_edge_id).unwrap();
            *e = (*e + (source_entry * scalar)) % self.field_mod;
            if *e == 0 {
                edges_to_remove.push(target_edge_id);
            }
        }
        for id in edges_to_remove.into_iter() {
            self.hgraph.remove_edge(id);
        }
        if self.hgraph.containing_edges_of_nodes([*t_row]).len() == 0 {
            self.hgraph.remove_node(*t_row);
            self.row_nodes.remove(&target_row);
        }
    }

    pub fn scale_row(&mut self, row_ix: usize, scalar: FFRep) {
        let row_node = self.row_nodes.get(&row_ix);
        if row_node.is_none() {
            return;
        }
        let row_node = row_node.unwrap();
        let edges = self.hgraph.containing_edges_of_nodes([*row_node]);
        for edge_id in edges {
            let e = self.hgraph.get_edge_mut(&edge_id).unwrap();
            *e = (*e * scalar) % self.field_mod;
        }
    }

    pub fn add_row_and_pivot(
        &mut self,
        row_ix: usize,
        row: SparseVector,
    ) -> Option<(usize, usize)> {
        self.insert_entries(row.0.into_iter().map(|x| (row_ix, x.0, x.1)).collect());
        self.pivotize_row(row_ix).map(|col_ix| (row_ix, col_ix))
    }

    pub fn pivotize_row(&mut self, row_ix: usize) -> Option<usize> {
        let row_node = self.row_nodes.get(&row_ix);
        if row_node.is_none() {
            return None;
        }
        let row_node = *row_node.unwrap();
        let link = self.hgraph.link_of_nodes([row_node]);
        if link.is_empty() {
            self.row_nodes.remove(&row_ix);
            return None;
        }

        let pivot_col_node = link.iter().min_by_key(|x| x.1[0]).unwrap().1[0];

        self.pivotize(row_ix, *self.hgraph.get_node(&pivot_col_node).unwrap())
    }
    /// Returns the column of the successfully pivotized row, or None if
    /// the row is empty.
    pub fn pivotize(&mut self, row_ix: usize, col_ix: usize) -> Option<usize> {
        let row_node = self.row_nodes.get(&row_ix);
        if row_node.is_none() {
            return None;
        }
        let row_node = *row_node.unwrap();
        let pivot_col_node = *self.col_nodes.get(&col_ix).unwrap();
        let pivot_id = self.hgraph.find_id([row_node, pivot_col_node]).unwrap();
        let pivot_entry_pre_scaling = self.hgraph.get_edge(&pivot_id).unwrap();
        let scalar = FF::new(*pivot_entry_pre_scaling, self.field_mod)
            .modular_inverse()
            .0;
        self.scale_row(row_ix, scalar);
        let rows_link = self.hgraph.link_of_nodes([pivot_col_node]);
        for (edge_id, row_node_vec) in rows_link {
            let target_row_node = row_node_vec[0];
            if target_row_node == row_node {
                continue;
            }
            let target_entry = self.hgraph.get_edge(&edge_id).unwrap();
            let scalar = -1 * FF::new(*target_entry, self.field_mod);
            let target_row_ix = self.hgraph.get_node(&target_row_node).unwrap();
            self.add_scaled_row_to_target(row_ix, *target_row_ix, scalar.0);
        }
        Some(*self.hgraph.get_node(&pivot_col_node).unwrap())
    }

    pub fn print(&self) {
        let mut n_rows = 0;
        let mut n_cols = 0;
        let entries =
            self.row_nodes
                .keys()
                .fold(Vec::<(usize, usize, FFRep)>::new(), |mut acc, row_ix| {
                    let mut v = self
                        .get_row(row_ix)
                        .into_iter()
                        .map(|(col_ix, entry)| {
                            n_rows = n_rows.max(row_ix + 1);
                            n_cols = n_cols.max(col_ix + 1);
                            (*row_ix, col_ix, entry)
                        })
                        .collect();
                    acc.append(&mut v);
                    acc
                });
        let sparse_mat = SparseFFMatrix::new_with_entries(
            n_rows,
            n_cols,
            self.field_mod,
            MemoryLayout::RowMajor,
            entries,
        );
        println!("{:}", sparse_mat.to_dense());
    }
}

impl From<SparseFFMatrix> for SparseSparseFFMatrix {
    fn from(value: SparseFFMatrix) -> Self {
        let field_mod = value.field_mod;
        let entries = value.into_entries();
        let mut ret = SparseSparseFFMatrix::new(field_mod);
        ret.insert_entries(entries);
        ret
    }
}

#[cfg(test)]
mod tests {
    use crate::matrices::sparse_sparse_ffmatrix::SparseSparseFFMatrix;

    #[test]
    fn basic_row_ops() {
        let mut m = SparseSparseFFMatrix::new(5);
        m.insert_entries(vec![(0, 0, 1), (0, 1, 2), (0, 2, 3), (0, 3, 1)]);
        println!("{:?}", m.get_row(&0));
        m.insert(1, 3, 3);
        m.print();
        m.add_scaled_row_to_target(0, 1, 2);
        m.print();
        m.scale_row(0, 3);
        m.print();
        for ix in 2..7 {
            m.add_scaled_row_to_target(0, ix, ix as u32);
        }
        m.print();
        println!("pivot 0");
        println!("col: {:?}", m.pivotize_row(0));
        m.print();
        println!("pivot 1");
        println!("col: {:?}", m.pivotize_row(1));
        m.print();
    }
}
