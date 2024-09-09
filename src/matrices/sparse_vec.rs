use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::math::finite_field::{FFRep, FiniteField};

#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
/// Meant to be used as a row or column of a sparse matrix. Stores the index and
/// value of the nonzero entries in the vector. The dimension of the vector
/// is not stored, although that probably wouldn't be too hard to add.
pub struct SparseVector(pub Vec<(usize, FFRep)>);
impl SparseVector {
    pub fn new_empty() -> Self {
        SparseVector(Vec::new())
    }

    pub fn first_nonzero(&self) -> Option<(usize, FFRep)> {
        for (col_ix, entry) in self.0.iter() {
            if *entry != 0 {
                return Some((*col_ix, *entry));
            }
        }
        None
    }

    pub fn to_vec(self) -> Vec<(usize, FFRep)> {
        self.0
    }

    pub fn new_with_entries(entries: Vec<(usize, FFRep)>) -> Self {
        let mut new_entries = entries;
        new_entries.sort_by(|a, b| a.0.cmp(&b.0));
        SparseVector(new_entries)
    }

    pub fn is_zero(&self) -> bool {
        self.first_nonzero().is_none()
    }

    pub fn max_index(&self) -> usize {
        self.0.iter().fold(0, |acc, x| acc.max(x.0))
    }

    pub fn clone_and_scale(&self, scalar: FiniteField) -> Self {
        let v: Vec<(usize, FFRep)> = self
            .0
            .iter()
            .map(|(ix, entry)| (*ix, (&scalar * *entry).0))
            .collect();
        SparseVector(v)
    }

    /// Inserts the given `entry` at position `ix`, overwriting existing
    /// entries if they exist.
    pub fn insert(&mut self, ix: usize, entry: FFRep) {
        let insertion_point = self.0.partition_point(|&x| x.0 < ix);
        if insertion_point >= self.0.len() {
            self.0.insert(insertion_point, (ix, entry));
        } else if self.0[insertion_point].0 == ix {
            self.0[insertion_point].1 = entry;
        } else {
            self.0.insert(insertion_point, (ix, entry));
        }
    }
    pub fn query(&self, ix: &usize) -> FFRep {
        if let Ok(position) = self.0.binary_search_by_key(ix, |&(ix, _)| ix) {
            self.0[position].1
        } else {
            0
        }
    }

    pub fn sparsity(&self) -> usize {
        self.0.len()
    }

    pub fn scale(&mut self, scalar: FFRep, field_mod: FFRep) {
        for (_, entry) in self.0.iter_mut() {
            *entry *= scalar;
            *entry %= field_mod;
        }
    }

    pub fn add_to_self(&mut self, rhs: &SparseVector, field_mod: FFRep) {
        let mut hm: HashMap<usize, u32> = self.0.clone().into_iter().collect();
        for (ix, entry) in rhs.0.iter() {
            let e = hm.entry(*ix).or_default();
            *e += *entry;
            *e %= field_mod;
        }
        let mut v: Vec<(usize, FFRep)> =
            hm.into_iter().filter(|(_ix, entry)| *entry != 0).collect();
        v.sort_by(|a, b| a.0.cmp(&b.0));
        self.0 = v;
    }

    pub fn add_scaled_row_to_self(&mut self, scalar: FiniteField, rhs: &SparseVector) {
        let mut self_ix = 0;
        let mut rhs_ix = 0;
        let mut new_vec = Vec::new();
        loop {
            if self_ix == self.0.len() && rhs_ix < rhs.0.len() {
                let (rhs_pos, rhs_entry) = rhs.0[rhs_ix];
                new_vec.push((rhs_pos, (&scalar * rhs_entry).0));
                rhs_ix += 1;
            } else if self_ix < self.0.len() && rhs_ix == rhs.0.len() {
                new_vec.push(self.0[self_ix]);
                self_ix += 1;
            } else if self_ix == self.0.len() && rhs_ix == rhs.0.len() {
                break;
            } else {
                let (self_pos, self_entry) = self.0[self_ix];
                let (rhs_pos, rhs_entry) = rhs.0[rhs_ix];
                if self_pos == rhs_pos {
                    let new_entry = &(&scalar * rhs_entry) + self_entry;
                    if new_entry.0 != 0 {
                        new_vec.push((self_pos, new_entry.0));
                    }
                    self_ix += 1;
                    rhs_ix += 1;
                } else if self_pos < rhs_pos {
                    new_vec.push((self_pos, self_entry));
                    self_ix += 1;
                } else {
                    new_vec.push((rhs_pos, (&scalar * rhs_entry).0));
                    rhs_ix += 1;
                }
            }
        }
        self.0 = new_vec;
    }
}
