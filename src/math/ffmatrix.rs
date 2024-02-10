use std::ops::{Index, Mul};

use serde::{Deserialize, Serialize};

use crate::math::finite_field::FiniteField as FF;

fn dot_prod(left: Vec<FF>, right: Vec<FF>) -> FF {
    if left.len() != right.len() {
        panic!("Cannot take dot product of two unequal length vectors.")
    }
    let mut ret = FF::new(0, left[0].1);
    for ix in 0..left.len() {
        ret += left[ix] * right[ix];
    }
    ret
}

pub fn vandermonde(elements: &Vec<FF>, n_rows: usize) -> FFMatrix {
    let mut entries = Vec::new();
    for k in 0..n_rows {
        let mut tmp = elements.clone().into_iter().map(|x| x.pow(k as u32)).collect();
        entries.append(&mut tmp);
    }
    FFMatrix::new(entries, n_rows, elements.len())
}

#[derive(Debug, Clone, Hash, PartialEq, Eq, Serialize, Deserialize)]
pub struct FFMatrix {
    pub entries: Vec<FF>,
    pub n_rows: usize,
    pub n_cols: usize,
    pub field_mod: u32,
}

impl FFMatrix {
    pub fn new(entries: Vec<FF>,
        n_rows: usize,
        n_cols: usize) -> Self {
        if entries.len() == 0 {
            panic!("Why did you try to make an empty matrix?")
        }
        let p = entries[0].1;
        Self {
            entries,
            n_rows,
            n_cols,
            field_mod: p,
        }
    }
    pub fn zero(n_rows: usize, n_cols: usize, field_mod: u32) -> Self {
        let mut entries = Vec::with_capacity(n_rows * n_cols);
        for _ in 0..(n_rows * n_cols) {
            entries.push(FF::new(0, field_mod));
        }
        Self { entries, n_rows, n_cols, field_mod }
    }

    /// ix is the row index (starts at 0) and jx is the col index (also 
    /// starts at 0)
    fn convert_indices(&self, ix: usize, jx: usize) -> usize {
        ((ix % self.n_rows) * self.n_cols) + (jx % self.n_cols)
    }
    /// Clones the specified row of the matrix as `FiniteFieldPolynomials`, does not include information about the quotient polynomial.
    pub fn get_row(&self, row_ix: usize) -> Vec<FF> {
        let mut ret = Vec::with_capacity(self.n_cols);
        for jx in 0..self.n_cols {
            ret.push(self.entries[self.convert_indices(row_ix, jx)].clone());
        }
        ret
    }

    /// Clones the specified column of the matrix as `FiniteFieldPolynomials`, does not include information about the quotient polynomial.
    pub fn get_col(&self, col_ix: usize) -> Vec<FF> {
        let mut ret = Vec::with_capacity(self.n_rows);
        for ix in 0..self.n_rows {
            ret.push(self.entries[self.convert_indices(ix, col_ix)].clone());
        }
        ret
    }
}

impl PartialOrd for FFMatrix {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        if self.n_cols != other.n_cols {
            return None;
        }
        if self.n_rows != other.n_rows {
            return None;
        }

        for ix in 0..self.entries.len() {
            if self.entries[ix] == other.entries[ix] {
                continue;
            } else {
                return self.entries[ix].partial_cmp(&other.entries[ix]);
            }
        }
        Some(std::cmp::Ordering::Equal)
    }
    
}

impl Ord for FFMatrix {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(other).expect("Cannot order non-equal shape matrices")
    }
}

impl Mul<&Vec<FF>> for &FFMatrix {
    type Output = Vec<FF>;

    fn mul(self, rhs: &Vec<FF>) -> Self::Output {
        if rhs.len() != self.n_cols {
            panic!("Mismatched shapes in matrix-vector multiplication.")
        }
        let mut ret = Vec::new();
        let p = self.field_mod;
        for ix in 0..self.n_rows {
            let mut tmp = FF::new(0, p);
            for jx in 0..self.n_cols {
                tmp += self.entries[self.convert_indices(ix, jx)]
            }
            ret.push(tmp);
        }
        ret
    }
}

impl Mul for &FFMatrix {
    type Output = FFMatrix;

    fn mul(self, rhs: &FFMatrix) -> Self::Output {
        if self.n_cols != rhs.n_rows {
            panic!("Tried to multiply incompatible matrices.")
        }
        let mut entries = Vec::with_capacity(self.entries.len());
        let field_mod = self.field_mod;
        for ix in 0..self.n_rows {
            for jx in 0..rhs.n_cols {
                let mut entry = FF::new(0, self.field_mod);
                for kx in 0..self.n_cols {
                    entry += self[[ix, kx]] * rhs[[kx, jx]];
                }
                entries.push(entry);
            }
        }
        FFMatrix {
            entries,
            n_rows: self.n_rows,
            n_cols: rhs.n_cols,
            field_mod,
        }
    }
}

impl Index<[usize; 2]> for FFMatrix {
    type Output = FF;

    fn index(&self, index: [usize; 2]) -> &Self::Output {
        self.entries
            .get(self.convert_indices(index[0], index[1]))
            .expect("Matrix indexing out of bounds.")
    }
}