use core::num;
use std::{
    collections::HashMap,
    fmt::Display,
    fs::File,
    io::{Read, Write},
    ops::{AddAssign, Index, Mul},
    path::Path,
};

use fxhash::{FxHashMap, FxHashSet};
use mhgl::{HGraph, HyperGraph};
use serde::{Deserialize, Serialize};

use super::{ffmatrix::FFMatrix, sparse_ffmatrix::SparseVector};
use crate::math::finite_field::{FFRep, FiniteField as FF};

// TODO: Should I implement a Compressed Sparse Row version of this?
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SparseSparseFFMatrix {
    row_nodes: FxHashMap<usize, u32>,
    col_nodes: FxHashMap<usize, u32>,
    hgraph: HGraph<usize, FFRep>,
    field_mod: FFRep,
}

impl SparseSparseFFMatrix {
    /// Creates an all zeros sparse matrix with the specified memory layout.
    /// If your matrix is row-sparse then use `MemoryLayout::RowMajor`, if it is
    /// column sparse use `MemoryLayout::ColMajor`.
    pub fn new(field_mod: FFRep) -> Self {
        Self {
            row_nodes: FxHashMap::default(),
            col_nodes: FxHashMap::default(),
            hgraph: HGraph::new(),
            field_mod,
        }
    }
    pub fn new_with_entries(
        n_rows: usize,
        n_cols: usize,
        field_mod: FFRep,
        entries: Vec<(usize, usize, FFRep)>,
    ) -> Self {
        let mut new = SparseSparseFFMatrix::new(field_mod);
        new.insert_entries(entries);
        new
    }

    pub fn row(&self, row_ix: usize) -> SparseVector {
        todo!()
    }

    /// Adds the entry to (row_ix, col_ix) in the matrix
    pub fn insert(&mut self, row_ix: usize, col_ix: usize, entry: FFRep) {
        let r = self
            .row_nodes
            .entry(row_ix)
            .or_insert(self.hgraph.add_node(row_ix));
        let c = self
            .col_nodes
            .entry(col_ix)
            .or_insert(self.hgraph.add_node(col_ix));
        let id = self
            .hgraph
            .find_id([*r, *c])
            .or(Some(self.hgraph.add_edge([*r, *c], 0)))
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

    pub fn query(&self, row_ix: usize, col_ix: usize) -> FF {
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

    /// Finds the smallest row that satisfies the provided filter.
    fn find_smallest_nonzero_row(
        &self,
        col_ix: usize,
        filter: impl Fn(usize) -> bool,
    ) -> Option<usize> {
        self.col_nodes
            .get(&col_ix)
            .map(|node| self.hgraph.link_of_nodes([*node]))
            .map(|v| {
                v.into_iter().fold(usize::MAX, |acc, x| {
                    let node_data = self.hgraph.get_node(&x.1[0]).unwrap();
                    if filter(*node_data) {
                        acc.min(*node_data)
                    } else {
                        acc
                    }
                })
            })
            .filter(|r| *r == usize::MAX)
    }

    pub fn add_row_to_other(&mut self, source_row: usize, target_row: usize) {
        let source_row_node = self.row_nodes.get(&source_row);
        if source_row_node.is_none() {
            return;
        }
        let source_row_node = source_row_node.unwrap();
        let link = self.hgraph.link_of_nodes([*source_row_node]);
        let t_row = self
            .col_nodes
            .entry(target_row)
            .or_insert(self.hgraph.add_node(target_row));
        for (edge_id, col_node_vec) in link {
            let entry = *self.hgraph.get_edge(&edge_id).unwrap();
            let col_node = col_node_vec[0];
            let target_edge_id = self
                .hgraph
                .find_id([*t_row, col_node])
                .or(Some(self.hgraph.add_edge([*t_row, col_node], 0)))
                .unwrap();
            let e = self.hgraph.get_edge_mut(&target_edge_id).unwrap();
            *e = (*e + entry) % self.field_mod;
            if *e == 0 {
                self.hgraph.remove_edge(target_edge_id);
            }
        }
    }

    /// Returns the column of the successfully pivotized row, or None if
    /// the row is empty.
    pub fn pivotize_row(&mut self, row_ix: &usize) -> Option<usize> {
        todo!()
    }
}

mod tests {}
