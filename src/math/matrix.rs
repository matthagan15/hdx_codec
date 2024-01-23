use std::{
    fmt::{Display, Write},
    ops::{Add, Index, Mul, MulAssign},
};

use deepsize::DeepSizeOf;

use crate::math::polynomial::*;

#[derive(Debug, Clone, Hash, PartialEq, Eq, DeepSizeOf)]
pub struct PolyMatrix {
    // Todo: use a sparse representation.
    // one possible problem, hash maps using up lots of entropy?
    entries: Vec<FiniteFieldPolynomial>,
    n_rows: usize,
    n_cols: usize,
    field_mod: u32,
    quotient: FiniteFieldPolynomial,
}

impl PolyMatrix {
    fn zero(dim: usize, quotient: FiniteFieldPolynomial) -> PolyMatrix {
        let mut entries = Vec::with_capacity(dim * dim);
        let z = FiniteFieldPolynomial::zero(quotient.field_mod);
        for _ in 0..dim * dim {
            entries.push(z.clone());
        }
        PolyMatrix {
            entries,
            n_rows: dim,
            n_cols: dim,
            field_mod: z.field_mod,
            quotient,
        }
    }

    fn swap_rows(&mut self, row_1: usize, row_2: usize) {
        let r1 = row_1 % self.n_rows;
        let r2 = row_2 % self.n_rows;
        let r1_start_ix = self.convert_indices(r1, 0);
        let r2_start_ix = self.convert_indices(r2, 0);
        for k in 0..self.n_cols {
            let tmp = self.entries[r1_start_ix + k].clone();
            self.entries[r1_start_ix + k] = self.entries[r2_start_ix + k].clone();
            self.entries[r2_start_ix + k] = tmp;
        }

    }

    pub fn clean(&mut self) {
        for p in self.entries.iter_mut() {
            p.clean();
        }
    }

    pub fn id(dim: usize, quotient: FiniteFieldPolynomial) -> PolyMatrix {
        let mut entries = Vec::with_capacity(dim * dim);
        let zero = FiniteFieldPolynomial::zero(quotient.field_mod);
        let one = FiniteFieldPolynomial::constant(1, quotient.field_mod);
        for ix in 0..dim {
            for jx in 0..dim {
                entries.push(if ix == jx { one.clone() } else { zero.clone() })
            }
        }
        PolyMatrix {
            entries,
            n_rows: dim,
            n_cols: dim,
            field_mod: quotient.field_mod,
            quotient,
        }
    }

    // fn det(&self) -> () {let a = QuotientPoly::new(self);}

    pub fn basis_state(ix: usize, jx: usize, dim: usize, quotient: FiniteFieldPolynomial) -> Self {
        let mut entries = Vec::with_capacity(dim * dim);
        let zero = FiniteFieldPolynomial::zero(quotient.field_mod);
        let one = FiniteFieldPolynomial::constant(1, quotient.field_mod);
        for ix_2 in 0..dim {
            for jx_2 in 0..dim {
                entries.push(if ix_2 == jx_2 || (ix == ix_2 && jx == jx_2) {
                    one.clone()
                } else {
                    zero.clone()
                })
            }
        }
        PolyMatrix {
            entries,
            n_rows: dim,
            n_cols: dim,
            field_mod: quotient.field_mod,
            quotient,
        }
    }

    fn convert_indices(&self, ix: usize, jx: usize) -> usize {
        ((ix % self.n_rows) * self.n_cols) + (jx % self.n_cols)
    }

    pub fn get_mut(&mut self, ix: usize, jx: usize) -> &mut FiniteFieldPolynomial {
        let idx = self.convert_indices(ix, jx);
        self.entries.get_mut(idx).expect("Could not index entry.")
    }

    pub fn generator_matrix(  
        ix: usize,
        jx: usize,
        dim: usize,
        alpha: u32,
        quotient: FiniteFieldPolynomial,
    ) -> Self {
        if ix == jx {
            panic!("Indices cannot be equal for a generator matrix!")
        }
        let p = FiniteFieldPolynomial::monomial((alpha, quotient.field_mod).into() , 1);
        let mut m = PolyMatrix::id(dim, quotient.clone());
        let e = m.get_mut(ix, jx);
        *e = p;
        m
    }
}

impl PartialOrd for PolyMatrix {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        if self.quotient != other.quotient {
            return None;
        }
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

impl Ord for PolyMatrix {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(other).expect("Cannot order")
    }
}

impl Mul<&PolyMatrix> for &PolyMatrix {
    type Output = PolyMatrix;

    fn mul(self, rhs: &PolyMatrix) -> Self::Output {
        let mut entries = Vec::with_capacity(self.entries.len());
        let field_mod = self.field_mod;
        let q = self.quotient.clone();
        for ix in 0..self.n_rows {
            for jx in 0..rhs.n_cols {
                let mut entry = FiniteFieldPolynomial::zero(self.field_mod);
                for kx in 0..self.n_cols {
                    let left_entry = &self[[ix, kx]];
                    let right_entry = &rhs[[kx, jx]];
                    entry += left_entry * right_entry;
                }
                entry = &entry % &self.quotient;
                entries.push(entry);
            }
        }
        PolyMatrix {
            entries,
            n_rows: self.n_rows,
            n_cols: self.n_cols,
            field_mod,
            quotient: q,
        }
    }
}

impl MulAssign<&PolyMatrix> for PolyMatrix {
    fn mul_assign(&mut self, rhs: &PolyMatrix) {
        let out = (self as &PolyMatrix) * rhs;
        *self = out;
    }
}

impl Add<&PolyMatrix> for PolyMatrix {
    type Output = PolyMatrix;

    fn add(self, rhs: &PolyMatrix) -> Self::Output {
        if (self.n_rows, self.n_cols) != (rhs.n_rows, rhs.n_cols) {
            panic!("Shapes are not equal squares!")
        }
        let mut new_entries = self.entries.clone();
        for ix in 0..self.n_rows {
            for jx in 0..self.n_cols {
                let idx = self.convert_indices(ix, jx);
                new_entries[idx] += rhs[[ix, jx]].clone();
                new_entries[idx] = &new_entries[idx] % &self.quotient;
            }
        }
        PolyMatrix {
            entries: new_entries,
            n_rows: self.n_rows,
            n_cols: self.n_cols,
            field_mod: self.field_mod,
            quotient: self.quotient,
        }
    }
}


impl From<&[QuotientPoly]> for PolyMatrix {
    fn from(value: &[QuotientPoly]) -> Self {
        let l = value.len();
        let mut is_square = false;
        let mut n_cols = 0;
        for ix in 1..=l {
            if ix * ix == l {
                is_square = true;
                n_cols = ix;
                break;
            }
        }
        if is_square == false {
            panic!("Attempted to make matrix that is not square!")
        }
        let new_entries = value.iter().map(|f| {
            &f.poly % &f.quotient
        }).collect();
        PolyMatrix {
            entries: new_entries,
            n_cols,
            n_rows: n_cols,
            field_mod: value[0].field_mod,
            quotient: value[0].quotient.clone(),
        }
    }
}

impl Index<[usize; 2]> for PolyMatrix {
    type Output = FiniteFieldPolynomial;

    fn index(&self, index: [usize; 2]) -> &Self::Output {
        self.entries
            .get(self.convert_indices(index[0], index[1]))
            .expect("Matrix indexing out of bounds.")
    }
}

impl Display for PolyMatrix {
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
        f.write_str(&format!(" modulo [{:}] over F_{:}", self.quotient, self.field_mod))
    }
}

mod tests {
    use crate::math::polynomial::{FiniteFieldPolynomial, QuotientPoly};

    use super::PolyMatrix;

    fn basic_matrix() -> PolyMatrix {
        let q = FiniteFieldPolynomial::monomial((1, 199).into(), 4);
        let p = FiniteFieldPolynomial::monomial((1, 199).into(), 1);
        let r = FiniteFieldPolynomial::monomial((7, 199).into(), 3);
        let data = [
            QuotientPoly::from((&p, &q)),
            QuotientPoly::from((&r, &q)),
            QuotientPoly::from((&r + &p, q.clone())),
            QuotientPoly::from((&p * &r, q.clone())),
        ];
        PolyMatrix::from(&data[..])
    }
    #[test]
    fn test_indexing() {
        let m = basic_matrix();
        println!("m:\n{:}", m);
        for ix in 0..m.n_rows {
            for jx in 0..m.n_cols {
                println!("m[[{:}, {:}]] = {:}", ix, jx, m[[ix, jx]]);
            }
        }
        // println!("m[[0,0]] = {:}", m[[0,0]]);
        // println!("m:\n{:}", m);
    }

    #[test]
    fn test_matrix_mul() {
        let m = basic_matrix();
        let m2 = basic_matrix();
        let out = &m * &m2;
        println!("input:\n{:}", m2);
        println!("out:\n{:}", out);
    }

    #[test]
    fn test_generator_matrices() {
        let q = FiniteFieldPolynomial::monomial((1_u32, 3_u32).into(), 3);
        let g = PolyMatrix::generator_matrix(1, 2, 3, 1, q);
        println!("g = {:}", g);
        let g_2 = &g * &g;
        let g_3 = &g * &g_2;
        let g_4 = &g * &g_3;
        println!("g_2 = {:}", g_2);
        println!("g_3 = {:}", g_3);
        println!("g_4 = {:}", g_4);
    }
}
