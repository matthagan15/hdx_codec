use std::{
    fmt::{Display, Write},
    ops::{Add, Index, Mul, MulAssign},
};

use crate::math::polynomial::*;

#[derive(Debug, Clone)]
struct Matrix {
    // Todo: use a sparse representation.
    // one possible problem, hash maps using up lots of entropy?
    entries: Vec<QuotientPoly>,
    n_rows: usize,
    n_cols: usize,
    field_mod: u32,
    quotient: FiniteFieldPolynomial,
}

impl Matrix {
    fn zero(dim: usize, quotient: FiniteFieldPolynomial) -> Matrix {
        let mut entries = Vec::with_capacity(dim * dim);
        let z = QuotientPoly::zero(quotient.field_mod, quotient);
        for _ in 0..dim * dim {
            entries.push(z.clone());
        }
        Matrix {
            entries,
            n_rows: dim,
            n_cols: dim,
            field_mod: z.field_mod,
            quotient: z.quotient,
        }
    }

    pub fn id(dim: usize, quotient: FiniteFieldPolynomial) -> Matrix {
        let mut entries = Vec::with_capacity(dim * dim);
        let zero = QuotientPoly::zero(quotient.field_mod, quotient.clone());
        let one = QuotientPoly::constant(1, quotient.clone());
        for ix in 0..dim {
            for jx in 0..dim {
                entries.push(if ix == jx { one.clone() } else { zero.clone() })
            }
        }
        Matrix {
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
        let zero = QuotientPoly::zero(quotient.field_mod, quotient.clone());
        let one = QuotientPoly::constant(1, quotient.clone());
        for ix_2 in 0..dim {
            for jx_2 in 0..dim {
                entries.push(if ix_2 == jx_2 || (ix == ix_2 && jx == jx_2) {
                    one.clone()
                } else {
                    zero.clone()
                })
            }
        }
        Matrix {
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

    pub fn get_mut(&mut self, ix: usize, jx: usize) -> &mut QuotientPoly {
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
        let p = QuotientPoly::monomial(alpha , 1, quotient.clone());
        let mut m = Matrix::id(dim, quotient.clone());
        let e = m.get_mut(ix, jx);
        *e = p;
        m
    }
}

impl Mul<&Matrix> for &Matrix {
    type Output = Matrix;

    fn mul(self, rhs: &Matrix) -> Self::Output {
        let mut entries = Vec::with_capacity(self.entries.len());
        let field_mod = self.field_mod;
        let q = self.quotient.clone();
        for ix in 0..self.n_rows {
            for jx in 0..rhs.n_cols {
                let mut entry = QuotientPoly::zero(self.field_mod, self.quotient.clone());
                for kx in 0..self.n_cols {
                    let left_entry = &self[[ix, kx]];
                    let right_entry = &rhs[[kx, jx]];
                    entry += left_entry * right_entry;
                }
                entries.push(entry);
            }
        }
        Matrix {
            entries,
            n_rows: self.n_rows,
            n_cols: self.n_cols,
            field_mod,
            quotient: q,
        }
    }
}

impl MulAssign<&Matrix> for Matrix {
    fn mul_assign(&mut self, rhs: &Matrix) {
        let out = (self as &Matrix) * rhs;
        *self = out;
    }
}

impl Add<&Matrix> for Matrix {
    type Output = Matrix;

    fn add(self, rhs: &Matrix) -> Self::Output {
        todo!()
    }
}

impl From<(QuotientPoly, usize)> for Matrix {
    fn from(value: (QuotientPoly, usize)) -> Self {
        todo!()
    }
}
impl From<&[QuotientPoly]> for Matrix {
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
        Matrix {
            entries: value.clone().to_vec(),
            n_cols,
            n_rows: n_cols,
            field_mod: value[0].field_mod,
            quotient: value[0].quotient.clone(),
        }
    }
}

impl Index<[usize; 2]> for Matrix {
    type Output = QuotientPoly;

    fn index(&self, index: [usize; 2]) -> &Self::Output {
        self.entries
            .get(self.convert_indices(index[0], index[1]))
            .expect("Matrix indexing out of bounds.")
    }
}

impl Display for Matrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut max_string_len = 0;
        let strings: Vec<String> = self
            .entries
            .iter()
            .map(|q| {
                let s = q.poly.to_string();
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
        f.write_str(&format!(" modulo {:}", self.quotient))
    }
}

mod tests {
    use crate::math::polynomial::{FiniteFieldPolynomial, QuotientPoly};

    use super::Matrix;

    fn basic_matrix() -> Matrix {
        let q = FiniteFieldPolynomial::monomial((1, 199).into(), 4);
        let p = FiniteFieldPolynomial::monomial((1, 199).into(), 1);
        let r = FiniteFieldPolynomial::monomial((7, 199).into(), 3);
        let data = [
            QuotientPoly::from((&p, &q)),
            QuotientPoly::from((&r, &q)),
            QuotientPoly::from((&r + &p, q.clone())),
            QuotientPoly::from((&p * &r, q.clone())),
        ];
        Matrix::from(&data[..])
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
        let g = Matrix::generator_matrix(1, 2, 3, 1, q);
        println!("g = {:}", g);
        let g_2 = &g * &g;
        let g_3 = &g * &g_2;
        let g_4 = &g * &g_3;
        println!("g_2 = {:}", g_2);
        println!("g_3 = {:}", g_3);
        println!("g_4 = {:}", g_4);
    }
}
