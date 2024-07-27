use core::num;
use std::{
    fmt::{Display, Write},
    ops::{Index, Mul},
};

use serde::{Deserialize, Serialize};

use crate::math::finite_field::FiniteField as FF;
use crate::math::finite_field::FiniteFieldExt as FFX;

use crate::math::finite_field::FFRep;

/// A matrix with entries in a
#[derive(Debug, Clone, Hash, PartialEq, Eq, Serialize, Deserialize)]
pub struct FFXMatrix {
    entries: Vec<FFRep>,
    n_rows: usize,
    n_cols: usize,
    prime_base: FFRep,
    prime_power: FFRep,
}

impl FFXMatrix {
    pub fn new(entries: Vec<FFRep>, dim: usize, prime_base: FFRep, prime_power: FFRep) -> Self {
        if entries.len() != dim * dim {
            panic!("Only allow new square matrices for consistency.")
        }
        Self {
            entries,
            n_rows: dim,
            n_cols: dim,
            prime_base,
            prime_power,
        }
    }
    pub fn id(dim: usize, prime_base: FFRep, prime_power: FFRep) -> Self {
        let mut entries = Vec::with_capacity(dim * dim);
        for ix in 0..dim {
            for jx in 0..dim {
                entries.push(if ix == jx { 1 } else { 0 })
            }
        }
        FFXMatrix {
            entries,
            n_rows: dim,
            n_cols: dim,
            prime_base,
            prime_power,
        }
    }
    pub fn nrows(&self) -> usize {
        self.n_rows
    }
    pub fn ncols(&self) -> usize {
        self.n_cols
    }

    pub fn set_entry(&mut self, row_ix: usize, col_ix: usize, new_entry: FFX) {
        let ix = self.convert_indices(row_ix, col_ix);
        self.entries[ix] = new_entry.0;
    }

    /// ix is the row index (starts at 0) and jx is the col index (also
    /// starts at 0)
    fn convert_indices(&self, ix: usize, jx: usize) -> usize {
        ((ix % self.n_rows) * self.n_cols) + (jx % self.n_cols)
    }

    pub fn k_type_subgroup(dim: usize) {}
}
impl PartialOrd for FFXMatrix {
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

impl Ord for FFXMatrix {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(other)
            .expect("Cannot order non-equal shape matrices")
    }
}

impl From<(Vec<FFRep>, FFRep, FFRep)> for FFXMatrix {
    /// Takes the mod of each entry and returns a proper `FFXMatrix`
    fn from(value: (Vec<FFRep>, FFRep, FFRep)) -> Self {
        let q = (value.1).pow(value.2);
        let entries: Vec<FFRep> = value.0.into_iter().map(|x| x % q).collect();
        let dim = (entries.len() as f64).sqrt().round() as usize;
        FFXMatrix {
            entries,
            n_rows: dim,
            n_cols: dim,
            prime_base: value.1,
            prime_power: value.2,
        }
    }
}

impl Mul for &FFXMatrix {
    type Output = FFXMatrix;

    fn mul(self, rhs: &FFXMatrix) -> Self::Output {
        if self.n_cols != rhs.n_rows {
            panic!("Tried to multiply incompatible matrices.")
        }
        let mut entries = Vec::with_capacity(self.n_rows * rhs.n_cols);
        let modulator = self.prime_base.pow(self.prime_power);
        for ix in 0..self.n_rows {
            for jx in 0..rhs.n_cols {
                let mut entry = 0;
                for kx in 0..self.n_cols {
                    entry += self[[ix, kx]] * rhs[[kx, jx]];
                }
                entries.push(entry % modulator);
            }
        }
        FFXMatrix {
            entries,
            n_rows: self.n_rows,
            n_cols: rhs.n_cols,
            prime_base: self.prime_base,
            prime_power: self.prime_power,
        }
    }
}

impl Index<[usize; 2]> for FFXMatrix {
    type Output = FFRep;

    fn index(&self, index: [usize; 2]) -> &Self::Output {
        let ix = self.convert_indices(index[0], index[1]);
        &self.entries[ix]
    }
}

impl Display for FFXMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut max_string_len = 0;
        let strings: Vec<String> = self
            .entries
            .iter()
            .map(|q| {
                let s = q.to_string();
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
        f.write_char('\n')?;
        f.write_str(&"_".repeat(row_len - 1))?;
        f.write_char('\n')?;
        for row in rows {
            f.write_str(&row)?;
        }
        f.write_str(&"-".repeat(row_len - 1))?;
        f.write_str(&format!(
            " modulo F_( {:} ^ {:})",
            self.prime_base, self.prime_power
        ))
    }
}

mod tests {
    use crate::math::finite_field::FiniteField;

    use super::FFXMatrix;

    #[test]
    fn test_ffx_matrix() {
        let a = vec![1, 2, 3, 4, 5, 6, 7, 8, 9];
        let m = FFXMatrix::new(a, 3, 3, 2);
        println!("m: {:}", m);
    }
}
