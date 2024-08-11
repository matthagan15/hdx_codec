use std::{
    cmp::Ordering,
    fmt::{Display, Write},
    hash::Hash,
    ops::Mul,
    rc::Rc,
    sync::Arc,
};

use serde::{Deserialize, Serialize};

use crate::math::{finite_field::FFRep, galois_field::GaloisField, polynomial::FFPolynomial};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GaloisMatrix {
    pub entries: Vec<FFRep>,
    pub n_rows: usize,
    pub n_cols: usize,
}

impl GaloisMatrix {
    pub fn new(n_rows: usize, n_cols: usize) -> Self {
        let entries = (0..n_rows * n_cols).into_iter().map(|_| 0).collect();
        GaloisMatrix {
            entries,
            n_rows,
            n_cols,
        }
    }

    /// ix is the row index (starts at 0) and jx is the col index (also
    /// starts at 0)
    fn convert_indices(&self, row_ix: usize, col_ix: usize) -> usize {
        ((row_ix % self.n_rows) * self.n_cols) + (col_ix % self.n_cols)
    }

    pub fn id(dim: usize) -> Self {
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
        }
    }

    pub fn set_entry(
        &mut self,
        row_ix: usize,
        col_ix: usize,
        new_entry: FFPolynomial,
        lookup: Arc<GaloisField>,
    ) {
        let rem = &new_entry % &lookup.quotient;
        let num = rem.to_number();
        let ix = self.convert_indices(row_ix, col_ix);
        self.entries[ix] = num;
    }

    pub fn get_num(&self, row_ix: usize, col_ix: usize) -> FFRep {
        self.entries[self.convert_indices(row_ix, col_ix)]
    }

    pub fn get_poly(&self, row_ix: usize, col_ix: usize, lookup: Arc<GaloisField>) -> FFPolynomial {
        let num = self.entries[self.convert_indices(row_ix, col_ix)];
        FFPolynomial::from_number(num, lookup.field_mod)
    }

    pub fn mul(&self, rhs: &GaloisMatrix, lookup: Arc<GaloisField>) -> GaloisMatrix {
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
                    let additive_amt = lookup.mul(left_entry, right_entry);
                    entry = lookup.add(&entry, &additive_amt);
                }
                entries.push(entry);
            }
        }
        GaloisMatrix {
            entries,
            n_rows: self.n_rows,
            n_cols: rhs.n_cols,
        }
    }

    pub fn pretty_print(&self, lookup: Arc<GaloisField>) -> String {
        let mut max_string_len = 0;
        let strings: Vec<String> = self
            .entries
            .iter()
            .map(|q| {
                let p = FFPolynomial::from_number(*q, lookup.field_mod);
                let s = p.to_string();
                max_string_len = max_string_len.max(s.len());
                s
            })
            .collect();
        let mut rows = Vec::new();
        let mut col_counter = 0;
        let mut row = String::from("| ");
        let mut row_len = 0;
        for mut string in strings {
            let diff_len = max_string_len - string.len();
            string.push_str(&" ".repeat(diff_len));
            col_counter += 1;
            col_counter %= self.n_cols;
            if col_counter == 0 {
                row.push_str(&string);
                row.push_str(" |\n");
                rows.push(row.clone());
                row_len = row.len();
                row.clear();
                row.push_str("| ");
            } else {
                row.push_str(&string);
                row.push_str(" ; ");
            }
        }
        let mut ret = String::new();
        ret.push_str(&"_".repeat(row_len - 1));
        ret.push('\n');
        for row in rows {
            ret.push_str(&row);
        }
        ret.push_str(&"-".repeat(row_len - 1));
        ret.push_str(&format!(" modulo F_{:}", lookup.field_mod));
        ret
    }
}

impl Hash for GaloisMatrix {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.entries.hash(state);
        self.n_rows.hash(state);
        self.n_cols.hash(state);
    }
}

impl PartialEq for GaloisMatrix {
    fn eq(&self, other: &Self) -> bool {
        self.entries == other.entries && self.n_rows == other.n_rows && self.n_cols == other.n_cols
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
