use std::{cmp::Ordering, hash::Hash, ops::Mul, rc::Rc};

use super::{finite_field::FFRep, galois_field::GaloisField, polynomial::FiniteFieldPolynomial};

#[derive(Debug, Clone)]
pub struct GaloisMatrix {
    pub entries: Vec<FFRep>,
    pub n_rows: usize,
    pub n_cols: usize,
    pub field_mod: u32,
    lookup_table: Rc<GaloisField>,
}

impl GaloisMatrix {
    pub fn new(n_rows: usize, n_cols: usize, lookup: Rc<GaloisField>) -> Self {
        let entries = (0..n_rows * n_cols).into_iter().map(|_| 0).collect();
        GaloisMatrix {
            entries,
            n_rows,
            n_cols,
            field_mod: lookup.field_mod,
            lookup_table: lookup,
        }
    }

    /// ix is the row index (starts at 0) and jx is the col index (also
    /// starts at 0)
    fn convert_indices(&self, row_ix: usize, col_ix: usize) -> usize {
        ((row_ix % self.n_rows) * self.n_cols) + (col_ix % self.n_cols)
    }

    pub fn id(dim: usize, lookup: Rc<GaloisField>) -> Self {
        let mut entries = Vec::new();
        for ix in 0..dim {
            for jx in 0..dim {
                entries.push(if ix == jx { 1 } else { 0 });
            }
        }
        GaloisMatrix {
            entries,
            n_rows: dim,
            n_cols: dim,
            field_mod: lookup.field_mod,
            lookup_table: lookup,
        }
    }

    pub fn set_entry(&mut self, row_ix: usize, col_ix: usize, new_entry: FiniteFieldPolynomial) {
        let rem = &new_entry % &self.lookup_table.quotient;
        let num = rem.to_number();
        let ix = self.convert_indices(row_ix, col_ix);
        self.entries[ix] = num;
    }

    pub fn get_num(&self, row_ix: usize, col_ix: usize) -> FFRep {
        self.entries[self.convert_indices(row_ix, col_ix)]
    }

    pub fn get_poly(&self, row_ix: usize, col_ix: usize) -> FiniteFieldPolynomial {
        let num = self.entries[self.convert_indices(row_ix, col_ix)];
        FiniteFieldPolynomial::from_number(num, self.lookup_table.field_mod)
    }
}

impl Mul<&GaloisMatrix> for &GaloisMatrix {
    type Output = GaloisMatrix;

    fn mul(self, rhs: &GaloisMatrix) -> Self::Output {
        if self.n_cols != rhs.n_rows {
            panic!("Tried to multiply incompatible matrices.")
        }
        let mut entries = Vec::with_capacity(self.n_rows * rhs.n_cols);
        for ix in 0..self.n_rows {
            for jx in 0..rhs.n_cols {
                let mut entry = 0;
                for kx in 0..self.n_cols {
                    let left_entry = self.get_num(ix, kx);
                    let right_entry = rhs.get_num(kx, jx);
                    let additive_amt = self.lookup_table.mul(left_entry, right_entry);
                    entry = self.lookup_table.add(&entry, &additive_amt);
                }
                entries.push(entry);
            }
        }
        GaloisMatrix {
            entries,
            n_rows: self.n_rows,
            n_cols: rhs.n_cols,
            field_mod: self.field_mod,
            lookup_table: self.lookup_table.clone(),
        }
    }
}

impl Hash for GaloisMatrix {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.entries.hash(state);
        self.n_rows.hash(state);
        self.n_cols.hash(state);
        self.field_mod.hash(state);
    }
}

impl PartialEq for GaloisMatrix {
    fn eq(&self, other: &Self) -> bool {
        self.entries == other.entries
            && self.n_rows == other.n_rows
            && self.n_cols == other.n_cols
            && self.field_mod == other.field_mod
    }
}

impl PartialOrd for GaloisMatrix {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        if (self.n_rows, self.n_cols) != (other.n_rows, other.n_cols) {
            return None;
        }

        for ix in 0..self.entries.len() {
            let entry_cmp = self.entries[ix].cmp(&other.entries[ix]);
            if entry_cmp != Ordering::Equal {
                return Some(entry_cmp);
            }
        }
        Some(Ordering::Equal)
    }
}

impl Eq for GaloisMatrix {}

impl Ord for GaloisMatrix {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other)
            .expect("Incorrect shapes for comparing.")
    }
}
