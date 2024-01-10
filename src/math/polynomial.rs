use std::{collections::HashMap, ops::{Mul, IndexMut, Index, Add, Sub, SubAssign}};


use super::group_ring_field::{Ring, Field, self};

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

#[derive(Debug)]
struct QuotientedPoly<F: Field> {
    pub coeffs: Vec<F>,
    pub divisor: Polynomial<F>,
}

/// Polynomial in single indeterminate. Uses dense storage (vec)
#[derive(Debug, Clone, PartialEq)]
pub struct Polynomial<F> where F: Field {
    pub coeffs: Vec<F>,
}

impl<F> Polynomial<F>
where F: Field {
    pub fn deg(&self) -> usize {
        let n = self.coeffs.len();
        for ix in 0..self.coeffs.len() {
            if self.coeffs[n - ix - 1] != <F as group_ring_field::Field>::zero() {
                return n - ix - 1;
            }
        }
        0
    }

    pub fn leading_coeff(&self) -> F {
        for coeff in self.coeffs.iter().rev() {
            if *coeff != <F as group_ring_field::Field>::zero() {
                return *coeff;
            }
        }
        <F as group_ring_field::Field>::zero()
    }

    pub fn monomial(coefficient: F, degree: usize) -> Polynomial<F> {
        let mut ret: Vec<F> = (0..degree).into_iter().map(|_| <F as group_ring_field::Field>::zero()).collect();
        ret.push(coefficient);
        Polynomial {coeffs: ret}
    }
    /// Performs Euclidean division, returing a tuple of (quotient, remainder)
    pub fn division(&self, denominator: &Polynomial<F>) -> (Polynomial<F>, Polynomial<F>) {
        let mut quotient = Polynomial {coeffs: vec![<F as group_ring_field::Ring>::zero()]};
        let mut remainder = self.clone();
        let d = denominator.deg();
        let mut lc = *denominator.coeffs.last().expect("Tried to divide by zero :'( ");
        while remainder.deg() >= d {
            let coeff_s = remainder.leading_coeff() * lc.mul_inv();
            let deg_s = remainder.deg() - d;
            let s = Polynomial::monomial(coeff_s, deg_s);
            let remainder_sub = s.clone() * denominator.clone();
            quotient = quotient + s;
            remainder = remainder -  remainder_sub;
        }
        remove_trailing_zeros(&mut quotient.coeffs);
        remove_trailing_zeros(&mut remainder.coeffs);
        (quotient, remainder)
    }
}

impl<F: Field> Mul for Polynomial<F> {
    type Output = Polynomial<F>;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut ret: Vec<F> = (1..=(self.coeffs.len() + rhs.coeffs.len() - 1)).into_iter().map(|_| <F as group_ring_field::Field>::zero()).collect();
        for deg_1 in 0..self.coeffs.len() {
            for deg_2 in 0..rhs.coeffs.len() {
                let coeff_1 = self.coeffs[deg_1];
                let coeff_2 = rhs.coeffs[deg_2];
                if let Some( e) = ret.get_mut(deg_1 + deg_2) {
                    *e += coeff_1 * coeff_2;
                }
            }
        }
        Polynomial { coeffs: ret }
    }
}

impl<F: Field> Add for Polynomial<F> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let deg = usize::max(self.deg(), rhs.deg());
        let mut ret: Vec<F> = (0..(deg + 1)).into_iter().map(|_| <F as group_ring_field::Field>::zero()).collect();
        for ix in 0..self.coeffs.len() {
            ret[ix] += self.coeffs[ix];
        }
        for ix in 0..rhs.coeffs.len() {
            ret[ix] += rhs.coeffs[ix];
        }
        Polynomial {coeffs: ret}
    }
}

impl<F: Field> Sub for Polynomial<F> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let deg = usize::max(self.deg(), rhs.deg());
        let mut ret: Vec<F> = (0..(deg + 1)).into_iter().map(|_| <F as group_ring_field::Field>::zero()).collect();
        for ix in 0..self.deg() + 1 {
            ret[ix] += self.coeffs[ix];
        }
        for ix in 0..rhs.deg() + 1 {
            ret[ix] += rhs.coeffs[ix].additive_inv();
        }
        Polynomial {coeffs: ret}
    }
}

impl<F: Field> SubAssign for Polynomial<F> {
    fn sub_assign(&mut self, rhs: Self) {
        let deg = usize::max(self.deg(), rhs.deg());
        let mut ret: Vec<F> = (0..(deg + 1)).into_iter().map(|_| <F as group_ring_field::Field>::zero()).collect();
        for ix in 0..self.deg() + 1 {
            ret[ix] += self.coeffs[ix];
        }
        for ix in 0..rhs.deg() + 1 {
            ret[ix] += rhs.coeffs[ix].additive_inv();
        }
        remove_trailing_zeros(&mut ret);
        self.coeffs = ret;
    }
}

fn remove_trailing_zeros<F: Field>(coeffs: &mut Vec<F>) {
    let zero = <F as group_ring_field::Field>::zero();
    let n = coeffs.len();
    let mut is_trailing_zeros = coeffs[n - 1] == zero;
    while is_trailing_zeros {
        if let Some(coeff) = coeffs.pop() {
            if coeff != zero {
                coeffs.push(coeff);
                is_trailing_zeros = false;
            }
        } else {
            break;
        }
    }
    coeffs.shrink_to_fit();
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
impl Field for f64 {
    fn mul_inv(&self) -> Self {
        1. / *self
    }

    fn zero() -> Self {
        0.0_f64
    }

    fn one() -> Self {
        1.0_f64
    }
}

mod tests {
    use std::collections::HashMap;

    use super::{Polynomial, remove_trailing_zeros};

    #[test]
    fn test_polynomial_multiplication() {
        let p1 = Polynomial {
            coeffs : vec![1.0, 2., 0.0, 3.]
        };
        let p2 = Polynomial {
            coeffs : vec![0.0, 5.0, 7.0]
        };
        dbg!(p1 * p2);
    }

    #[test]
    fn test_remove_trailing_zeros() {
        let mut r = vec![1., 2., 3., 4., 5., 0., 0., 2., 0., 0.,];
        remove_trailing_zeros(&mut r);
        println!("r capacity: {:}", r.capacity());
        dbg!(r);
    }

    #[test]
    fn test_constructors() {
        let p = Polynomial::monomial(3.2, 3);
        dbg!(&p);
        println!("p degree: {:}, leading_coeff: {:}", p.deg(), p.leading_coeff());
    }
    #[test]
    fn test_polynomial_division() {
        let b = Polynomial {
            coeffs: vec![1., 3., 21., 22., 23.]
        };
        let q = Polynomial {
            coeffs: vec![17., 9., 13.]
        };
        let r = Polynomial {
            coeffs: vec![5., 11.]
        };
        let a = b.clone() * q.clone() + r.clone();
        let (q_comp, r_comp) = a.division(&b);
        assert_eq!(q, q_comp);
        assert_eq!(r, r_comp);
        dbg!(q_comp);
        dbg!(r_comp);
    }
}