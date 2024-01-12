use std::{collections::HashMap, ops::{Mul, IndexMut, Index, Add, Sub, SubAssign}, fmt::Display};


use super::{group_ring_field::{Ring, Field, self}, finite_field::FiniteField};

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
struct QuotientedPoly {
    pub coeffs: Vec<FiniteField>,
    pub divisor: FiniteFieldPolynomial,
}

/// Polynomial in single indeterminate. Uses dense storage (vec)
#[derive(Debug, Clone, PartialEq)]
pub struct FiniteFieldPolynomial {
    pub coeffs: Vec<FiniteField>,
    pub field_mod: u32,
}

impl FiniteFieldPolynomial {
    pub fn deg(&self) -> usize {
        let n = self.coeffs.len();
        for ix in 0..self.coeffs.len() {
            if self.coeffs[n - ix - 1].0 != 0 {
                return n - ix - 1;
            }
        }
        0
    }

    pub fn leading_coeff(&self) -> FiniteField {
        for coeff in self.coeffs.iter().rev() {
            if coeff.0 != 0 {
                return *coeff;
            }
        }
        (0_u32, self.field_mod).into()
    }

    pub fn monomial(coefficient: FiniteField, degree: usize, field_mod: u32) -> FiniteFieldPolynomial {
        let mut ret: Vec<FiniteField> = (0..degree).into_iter().map(|_| (0, field_mod).into()).collect();
        ret.push(coefficient);
        FiniteFieldPolynomial {coeffs: ret, field_mod}
    }

    /// Performs Euclidean division, returing a tuple of (quotient, remainder)
    pub fn division(&self, denominator: &FiniteFieldPolynomial) -> (FiniteFieldPolynomial, FiniteFieldPolynomial) {
        let mut quotient = FiniteFieldPolynomial {coeffs: vec![(0, self.field_mod).into()], field_mod: self.field_mod};
        let mut remainder = self.clone();
        let d = denominator.deg();
        let mut lc = *denominator.coeffs.last().expect("Tried to divide by zero :'( ");
        while remainder.deg() >= d {
            let coeff_s = remainder.leading_coeff() * lc.mul_inv();
            let deg_s = remainder.deg() - d;
            let s = FiniteFieldPolynomial::monomial(coeff_s, deg_s, self.field_mod);
            let remainder_sub = s.clone() * denominator.clone();
            quotient = quotient + s;
            remainder = remainder -  remainder_sub;
        }
        remove_trailing_zeros(&mut quotient.coeffs);
        remove_trailing_zeros(&mut remainder.coeffs);
        (quotient, remainder)
    }
}

impl Display for FiniteFieldPolynomial {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&format!("{:} + ", self.coeffs[0].0)).expect("couldn't write polynomial?");
        for ix in 1..self.coeffs.len() - 1 {
            if self.coeffs[ix].0 != 0 {
                f.write_str(&format!("{:} x^{:} + ", self.coeffs[ix].0, ix)).expect("couldn't write polynomial?");
            }
        }
        f.write_str(&format!("{:} x^{:} mod {:}", self.coeffs[self.coeffs.len() - 1].0, self.coeffs.len() - 1, self.field_mod))
    }
}

impl Mul for FiniteFieldPolynomial {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut ret: Vec<FiniteField> = (1..=(self.coeffs.len() + rhs.coeffs.len() - 1)).into_iter().map(|_| (0, self.field_mod).into() ).collect();
        for deg_1 in 0..self.coeffs.len() {
            for deg_2 in 0..rhs.coeffs.len() {
                let coeff_1 = self.coeffs[deg_1];
                let coeff_2 = rhs.coeffs[deg_2];
                if let Some( e) = ret.get_mut(deg_1 + deg_2) {
                    *e += coeff_1 * coeff_2;
                }
            }
        }
        FiniteFieldPolynomial { coeffs: ret , field_mod: self.field_mod}
    }
}

impl Add for FiniteFieldPolynomial {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let deg = usize::max(self.deg(), rhs.deg());
        let mut ret: Vec<FiniteField> = (0..(deg + 1)).into_iter().map(|_| FiniteField(0, self.field_mod) ).collect();
        for ix in 0..self.coeffs.len() {
            ret[ix] += self.coeffs[ix];
        }
        for ix in 0..rhs.coeffs.len() {
            ret[ix] += rhs.coeffs[ix];
        }
        FiniteFieldPolynomial {coeffs: ret, field_mod: self.field_mod}
    }
}

impl Sub for FiniteFieldPolynomial {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let deg = usize::max(self.deg(), rhs.deg());
        let mut ret: Vec<FiniteField> = (0..(deg + 1)).into_iter().map(|_| FiniteField(0, self.field_mod)).collect();
        for ix in 0..self.deg() + 1 {
            ret[ix] += self.coeffs[ix];
        }
        for ix in 0..rhs.deg() + 1 {
            ret[ix] += rhs.coeffs[ix].additive_inv();
        }
        FiniteFieldPolynomial {coeffs: ret, field_mod: self.field_mod}
    }
}

impl SubAssign for FiniteFieldPolynomial {
    fn sub_assign(&mut self, rhs: Self) {
        let deg = usize::max(self.deg(), rhs.deg());
        let mut ret: Vec<FiniteField> = (0..(deg + 1)).into_iter().map(|_| FiniteField(0, self.field_mod)).collect();
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

fn remove_trailing_zeros(coeffs: &mut Vec<FiniteField>) {
    let n = coeffs.len();
    let mut is_trailing_zeros = coeffs[n - 1].0 == 0;
    while is_trailing_zeros {
        if let Some(coeff) = coeffs.pop() {
            if coeff.0 != 0 {
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

    use crate::math::finite_field::FiniteField;

    use super::{FiniteFieldPolynomial, remove_trailing_zeros};

    #[test]
    fn test_polynomial_multiplication() {
        let p1 = FiniteFieldPolynomial {
            coeffs : vec![(1, 199).into(), (2, 199).into(), (0, 199).into(), (3, 199).into()],
            field_mod: 199,
        };
        let z = FiniteField(0, 199);
        let p2 = FiniteFieldPolynomial {
            coeffs : vec![z, z + 5, z + 7],
            field_mod: 199,
        };
        dbg!(p1 * p2);
    }

    #[test]
    fn test_remove_trailing_zeros() {
        let z = FiniteField(0, 199);
        let mut r = vec![z + 1, z + 2, z + 3, z + 4, z + 5, z, z, z + 2, z + 0,z + 0];
        remove_trailing_zeros(&mut r);
        println!("r capacity: {:}", r.capacity());
        dbg!(r);
    }

    #[test]
    fn test_constructors() {
        let p = FiniteFieldPolynomial::monomial((0, 199).into(), 3, 199);
        dbg!(&p);
        println!("p degree: {:}, leading_coeff: {:}", p.deg(), p.leading_coeff());
    }
    #[test]
    fn test_polynomial_division() {
        let z = FiniteField(0, 199);
        let b = FiniteFieldPolynomial {
            coeffs: vec![z + 1, z + 3, z + 21, z + 22, z + 23],
            field_mod: 199,
        };
        let q = FiniteFieldPolynomial {
            coeffs: vec![z + 17, z + 9, z + 13],
            field_mod: 199,
        };
        let r = FiniteFieldPolynomial {
            coeffs: vec![z + 5, z + 11],
            field_mod: 199,
        };
        let a = b.clone() * q.clone() + r.clone();
        let (q_comp, r_comp) = a.division(&b);
        assert_eq!(q, q_comp);
        assert_eq!(r, r_comp);
        println!("quotient - {:}", q_comp);
        println!("remainder - {:}", r_comp);
    }
}