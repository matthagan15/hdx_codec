use std::{collections::HashMap, ops::Mul};

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
    pub fn nnz(&self) -> usize {
        self.0.len()
    }

    /// Finds the first nonzero entry with an index larger than (or equal to) the provided
    /// `min_ix` cutoff.
    pub fn first_nonzero_larger_than(&self, min_ix: usize) -> Option<(usize, FFRep)> {
        self.0
            .iter()
            .filter(|(ix, _)| *ix >= min_ix)
            .min_by(|a, b| a.0.cmp(&b.0))
            .copied()
    }

    pub fn capacity(&self) -> usize {
        self.0.capacity()
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

    /// Inserts the given `entry` at position `ix`, overwriting existing
    /// entries if they exist.
    pub fn insert(&mut self, ix: usize, entry: FFRep) {
        if entry == 0 {
            return;
        }
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

    pub fn swap(&mut self, ix: usize, jx: usize) {
        let memory_ix = self.0.binary_search_by_key(&ix, |&(ix, _)| ix);
        let memory_jx = self.0.binary_search_by_key(&jx, |&(jx, _)| jx);
        match (memory_ix, memory_jx) {
            (Ok(mem_ix), Ok(mem_jx)) => {
                let tmp = self.0[mem_ix].1;
                self.0[mem_ix].1 = self.0[mem_jx].1;
                self.0[mem_jx].1 = tmp;
            }
            (Ok(mem_ix), Err(_)) => {
                let old_entry = self.0.remove(mem_ix);
                self.insert(jx, old_entry.1);
            }
            (Err(_), Ok(mem_jx)) => {
                let old_entry = self.0.remove(mem_jx);
                self.insert(ix, old_entry.1);
            }
            (Err(_), Err(_)) => {}
        }
    }

    /// Adds scalar * rhs to self. Does nothing if scalar is 0.
    pub fn add_scaled_row_to_self(&mut self, scalar: FiniteField, rhs: &SparseVector) {
        if scalar.0 == 0 {
            return;
        }
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

impl Mul<&SparseVector> for &SparseVector {
    type Output = u32;

    fn mul(self, rhs: &SparseVector) -> Self::Output {
        let mut self_ix = 0;
        let mut rhs_ix = 0;
        let mut tot = 0;
        loop {
            if self_ix >= self.0.len() || rhs_ix >= rhs.0.len() {
                break;
            }
            let (rhs_dim_ix, rhs_entry) = &rhs.0[rhs_ix];
            let (self_dim_ix, self_entry) = &self.0[self_ix];
            if self_dim_ix == rhs_dim_ix {
                tot += rhs_entry * self_entry;
                rhs_ix += 1;
                self_ix += 1;
            } else if self_dim_ix > rhs_dim_ix {
                rhs_ix += 1;
            } else {
                self_ix += 1;
            }
        }
        tot
    }
}

pub trait SparseVec {
    fn empty_set() -> Self;
    /// returns the (logical) index of the first nonzero.
    /// ```rust
    /// let s = SparseVec::from([3, 4, 5]);
    /// assert_eq!(s.first_nonzero(), Some(3));
    /// assert!(SparseVec::empty_set().first_nonzero().is_none());
    /// ```
    fn first_nonzero(&self) -> Option<usize>;
    fn nnz(&self) -> usize;
    /// returns the first logical ix with nonzero entry.
    fn first_nonzero_larger_than(&self, min_ix: usize) -> Option<usize>;
    fn is_zero(&self) -> bool;
    fn max_index(&self) -> Option<usize>;
    fn insert(&mut self, ix: usize, entry: FFRep);
    fn query(&self, ix: usize) -> Option<FFRep>;
    fn swap(&mut self, ix: usize, jx: usize);
    fn scale(&mut self, scaler: FFRep, field_mod: FFRep);
    fn add_scaled_row_to_self(&mut self, scaler: FFRep, field_mod: FFRep, other: &Self);
}

#[cfg(test)]
mod test {
    use crate::matrices::sparse_vec::SparseVector;

    #[test]
    fn test_sparse_section() {
        let v = vec![(0, 10), (200, 4), (10, 30), (1000, 2)];
        let mut ss = SparseVector::new_with_entries(v);
        assert_eq!(ss.query(&7), 0);
        assert_eq!(ss.query(&1000), 2);
        ss.insert(4, 17);
        assert_eq!(ss.query(&4), 17);
    }

    #[test]
    fn sparse_vector_additions() {
        let p = 7_u32;
        let mut v = Vec::new();
        for ix in 0..10 {
            v.push((ix as usize, 1));
        }
        let mut u = Vec::new();
        for ix in 10..20 {
            u.push((ix as usize, 1));
        }

        let adder = SparseVector::new_with_entries(vec![(3, 4), (4, 6), (9, 3)]);
        let mut sparse_v = SparseVector::new_with_entries(v);
        let mut sparse_u = SparseVector::new_with_entries(u);
        sparse_v.add_scaled_row_to_self(crate::math::finite_field::FiniteField::new(1, p), &adder);
        sparse_u
            .add_scaled_row_to_self(crate::math::finite_field::FiniteField::new(1, p), &sparse_v);
        dbg!(sparse_u);
        dbg!(sparse_v);
    }
}
