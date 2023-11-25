use std::{collections::HashMap, ops::{Mul, IndexMut, Index}};

use super::group_ring_field::Ring;

struct Matrix<R: Ring> {
    data: Vec<R>,
    n_rows: usize,
    n_cols: usize,
}

impl<R: Ring> Index<(usize, usize)> for Matrix<R> {
    type Output = R;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        let ix = index.1 + self.n_cols * (index.0 - 1);
        &self.data[ix]
    }
}

impl<R: Ring> IndexMut<(usize, usize)> for Matrix<R> {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        let ix = index.1 + self.n_cols * (index.0 - 1);
        &mut self.data[ix]
    }
}

impl<R: Ring> Matrix<R> {
    pub fn mul(&self, rhs: &Matrix<R>) -> Matrix<R> {
        if self.n_cols != rhs.n_rows {
            panic!("Matrix Dimension mismatch.")
        }
        let mut data = Vec::with_capacity(rhs.n_cols * self.n_rows);
        for _ in 0..self.n_rows * rhs.n_cols {
            data.push(R::zero());
        }
        let mut ret = Matrix {
            data,
            n_rows: self.n_rows,
            n_cols: rhs.n_cols,
        };
        for ix in 0..self.n_rows {
            for jx in 0..rhs.n_cols {
                let mut tot = R::zero();
                for kx in 0.. self.n_cols {
                    tot += self[(ix, kx)] * rhs[(kx, jx)];
                }
                ret[(ix, jx)] = tot;
            }
        }
        ret
    }
}

/// Polynomial in single indeterminate
#[derive(Debug)]
pub struct Polynomial<R> where R: Ring {
    pub coeffs: HashMap<i32, R>,
}

impl<R: Ring> Polynomial<R> {
    pub fn division(&self, rhs: &Polynomial<R>) {
        
    }
}

impl<R: Ring> Mul for Polynomial<R> {
    type Output = Polynomial<R>;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut ret = HashMap::new();
        for (deg_1, coeff_1) in self.coeffs.iter() {
            for (deg_2, coeff_2) in rhs.coeffs.iter() {
                let e = ret.entry(deg_1 + deg_2).or_insert(R::zero());
                *e += *coeff_1 * *coeff_2;
            }
        }
        Polynomial { coeffs: ret }
    }
}

impl Ring for f64 {
    fn zero() -> Self {
        0.0
    }

    fn one() -> Self {
        1.0
    }

    fn additive_inv(&self) -> Self {
        self * (-1.0)
    }
}

mod tests {
    use std::collections::HashMap;

    use super::Polynomial;

    #[test]
    fn test_polynomial_multiplication() {
        let p1 = Polynomial {
            coeffs : HashMap::from([
                (0, 1.0),
                (1, 2.),
                (3, 3.),
            ])
        };
        let p2 = Polynomial {
            coeffs : HashMap::from( [
                (1, 2.)
            ])
        };
        dbg!(p1 * p2);
    }
}