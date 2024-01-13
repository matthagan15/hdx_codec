use gcd::binary_u32;
use serde::de;
use std::{
    cmp::Ordering,
    collections::{HashMap, HashSet},
    fmt::Display,
    ops::{Add, Div, Index, IndexMut, Mul, Rem, Sub, SubAssign, AddAssign, MulAssign},
};

use super::{
    finite_field::FiniteField,
    group_ring_field::{self, Field, Ring},
};

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
                for kx in 0..self.n_cols {
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
    pub poly: FiniteFieldPolynomial,
    pub divisor: FiniteFieldPolynomial,
    pub field_mod: u32,
}

impl QuotientedPoly {
    // fn new(poly: FiniteFieldPolynomial, divisor: FiniteFieldPolynomial) -> Self {
    // Currently do not run check to save time
    // if poly.field_mod != divisor.field_mod {
    //     panic!("Polynomials defined over different fields.")
    // }
    // TODO: need to make the divisor primitive or check if it is primitive.
    // let p = poly.field_mod;
    // let (_, r) = poly / divisor.clone();
    // QuotientedPoly {
    //     poly: r,
    //     divisor,
    //     field_mod: p,
    // }
    // }
}

impl Add for QuotientedPoly {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        todo!()
    }
}

impl Sub for QuotientedPoly {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        todo!()
    }
}

impl Mul for QuotientedPoly {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        todo!()
    }
}

impl Div for QuotientedPoly {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        todo!()
    }
}

fn get_smallest_divisor(n: u32) -> Option<u32> {
    for x in 2..(n / 2) {
        if n % x == 0 {
            return Some(x);
        }
    }
    None
}

/// Returns all the prime numbers that divide the provided number `n`. Note will return `n` if itself is prime, but will not return `1`.
fn get_divisors(n: u32) -> Vec<u32> {
    let mut ret = Vec::new();
    if n == 1 {
        return ret;
    }
    if let Some(l) = get_smallest_divisor(n) {
        ret.push(l);
        let mut new_n = n;
        while new_n % l == 0 {
            new_n /= l;
        }
        ret.append(&mut get_divisors(new_n));
    } else {
        ret.push(n);
    }
    ret
}

/// Polynomial in single indeterminate. Uses dense storage (vec)
#[derive(Debug, Clone, PartialEq)]
pub struct FiniteFieldPolynomial {
    // pub coeffs: Vec<FiniteField>,
    pub coeffs: HashMap<usize, FiniteField>,
    degree: usize,
    pub field_mod: u32,
}

impl FiniteFieldPolynomial {
    /// Removes zeros and updates the degree
    fn clean(&mut self) {
        let mut new_degree = 0;
        let mut hm: HashMap<usize, FiniteField> = self
            .coeffs
            .clone()
            .into_iter()
            .filter(|(d, c)| {
                new_degree = new_degree.max(*d);
                c.0 != 0
            })
            .collect();
        if hm.len() == 0 {
            hm.insert(0, (0, self.field_mod).into());
        }
        self.degree = new_degree;
        self.coeffs = hm;
    }
    pub fn scale(&mut self, scalar: &FiniteField) {
        for (_, v) in self.coeffs.iter_mut() {
            *v *= scalar;
        }
    }

    /// Creates a new zero polynomial
    pub fn new(field_mod: u32) -> Self {
        let z = FiniteField::new(0, field_mod);
        let hm = HashMap::from([(0, z)]);
        FiniteFieldPolynomial {
            coeffs: hm,
            degree: 0,
            field_mod,
        }
    }

    //     pub fn gcd(&self, rhs: &FiniteFieldPolynomial) -> FiniteFieldPolynomial {
    //         if rhs.is_zero() {
    //             self.clone()
    //         } else {
    //             let (_, r) = self.division(rhs);
    //             rhs.gcd(&r)
    //         }
    //     }

    //     /// Implementation of the [Rabin Irreducibility Test ](https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields#Rabin's_test_of_irreducibility) to determine if `self` can be factored into two polynomials p(x) q(x) each of degree less than `self.deg()`.
    //     pub fn is_irreducible(&self) -> bool {
    //         let n = self.deg() as u32;
    //         let divisors = get_divisors(n);
    //         let powers: Vec<u32> = divisors
    //             .into_iter()
    //             .filter(|p| *p != n)
    //             .map(|p| n / p)
    //             .collect();
    //         for power in powers {
    //             let tmp_degree = self.field_mod.pow(power);
    //             let tmp_term_1 = FiniteFieldPolynomial::monomial(
    //                 (1, self.field_mod).into(),
    //                 tmp_degree as usize,
    //                 self.field_mod,
    //             );
    //             let tmp_term_2 =
    //                 FiniteFieldPolynomial::monomial((-1, self.field_mod).into(), 1, self.field_mod);
    //             let t = tmp_term_1 + tmp_term_2;
    //             let (_, tmp) = t.division(self);
    //             let g = self.gcd(&tmp);
    //             if g.deg() == 0 && g.coeffs[0].0 == 1 {
    //                 return false;
    //             }
    //         }
    //         let coeff_1: FiniteField = (1, self.field_mod).into();
    //         let deg_1 = self.field_mod.pow(n) as usize;
    //         let tmp_term_1 = FiniteFieldPolynomial::monomial(
    //             coeff_1,
    //             deg_1,
    //             self.field_mod,
    //         );
    //         let tmp_term_2 =
    //             FiniteFieldPolynomial::monomial((-1, self.field_mod).into(), 1, self.field_mod);
    //         let t = tmp_term_1 + tmp_term_2;

    //         let (_, g) = t.division(self);

    //         g.is_zero()
    //     }

    //     pub fn is_primitive(&self) -> bool {
    //         let d = self.deg();
    //         let upper_bound = self.field_mod.pow(d as u32) - 1;
    //         for n in 1..upper_bound {
    //             let t1 = FiniteFieldPolynomial::monomial(
    //                 (1, self.field_mod).into(),
    //                 n as usize,
    //                 self.field_mod,
    //             );
    //             let t2 = FiniteFieldPolynomial::new(vec![(-1, self.field_mod).into()]);
    //             let g = t1 + t2;
    //             let (_, r) = self.division(&g);
    //             if r.is_zero() {
    //                 return false;
    //             }
    //         }
    //         let t1 = FiniteFieldPolynomial::monomial(
    //             (1, self.field_mod).into(),
    //             upper_bound as usize,
    //             self.field_mod,
    //         );
    //         let t2 = FiniteFieldPolynomial::new(vec![(-1, self.field_mod).into()]);
    //         let g = t1 + t2;
    //         let (_, r) = self.division(&g);
    //         r.is_zero()
    //     }

    //     pub fn is_zero(&self) -> bool {
    //         if self.coeffs.len() == 0 {
    //             return true;
    //         }
    //         let mut coeffs_are_zero = self.coeffs[0].0 == 0;
    //         for coeff in self.coeffs.iter() {
    //             coeffs_are_zero |= coeff.0 == 0;
    //         }
    //         coeffs_are_zero
    //     }

    //     pub fn deg(&self) -> usize {
    //         let n = self.coeffs.len();
    //         for ix in 0..self.coeffs.len() {
    //             if self.coeffs[n - ix - 1].0 != 0 {
    //                 return n - ix - 1;
    //             }
    //         }
    //         0
    //     }

    //     pub fn leading_coeff(&self) -> FiniteField {
    //         for coeff in self.coeffs.iter().rev() {
    //             if coeff.0 != 0 {
    //                 return *coeff;
    //             }
    //         }
    //         (0_u32, self.field_mod).into()
    //     }

    //     pub fn monomial(
    //         coefficient: FiniteField,
    //         degree: usize,
    //         field_mod: u32,
    //     ) -> FiniteFieldPolynomial {
    //         println!("degree: {:}", degree);
    //         let mut ret: Vec<FiniteField> = (0..degree)
    //             .into_iter()
    //             .map(|d| {
    //                 (0, field_mod).into()
    //             })
    //             .collect();
    //         ret.push(coefficient);
    //         dbg!("post push?");
    //         FiniteFieldPolynomial {
    //             coeffs: ret,
    //             field_mod,
    //         }
    //     }

    //     /// Performs Euclidean division, returing a tuple of (quotient, remainder)
    //     pub fn division(
    //         &self,
    //         denominator: &FiniteFieldPolynomial,
    //     ) -> (FiniteFieldPolynomial, FiniteFieldPolynomial) {
    //         let mut quotient = FiniteFieldPolynomial {
    //             coeffs: vec![(0, self.field_mod).into()],
    //             field_mod: self.field_mod,
    //         };
    //         let mut remainder = self.clone();
    //         let d = denominator.deg();
    //         let mut lc = *denominator
    //             .coeffs
    //             .last()
    //             .expect("Tried to divide by zero :'( ");
    //         while remainder.deg() >= d {
    //             let coeff_s = remainder.leading_coeff() * lc.modular_inverse();
    //             let deg_s = remainder.deg() - d;
    //             let s = FiniteFieldPolynomial::monomial(coeff_s, deg_s, self.field_mod);
    //             let remainder_sub = s.clone() * denominator.clone();
    //             quotient = quotient + s;
    //             remainder = remainder - remainder_sub;
    //         }
    //         remove_trailing_zeros(&mut quotient.coeffs);
    //         remove_trailing_zeros(&mut remainder.coeffs);
    //         (quotient, remainder)
    //     }
}

impl Display for FiniteFieldPolynomial {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.coeffs.len() == 0 {
            return f.write_str("zero polynomial encountered.");
        }
        let mut coeffs: Vec<(usize, FiniteField)> = self.coeffs.clone().into_iter().collect();
        coeffs.sort_by(|x, y| x.0.cmp(&y.0).reverse());

        if coeffs[0].0 == 0 {
            f.write_str(&format!("{:}", coeffs[0].1 .0))
                .expect("couldn't write polynomial?");
        } else {
            f.write_str(&format!("{:} x^{:}", coeffs[0].1 .0, coeffs[0].0))
                .expect("couldn't write polynomial?");
        }
        if coeffs.len() >= 2 {
            for ix in 1..coeffs.len() {
                if coeffs[ix].1 .0 != 0 {
                    if coeffs[ix].0 == 0 {
                        f.write_str(&format!(" + {:}", coeffs[ix].1 .0))
                            .expect("couldn't write polynomial?");
                    } else {
                        f.write_str(&format!(" + {:} x^{:}", coeffs[ix].1 .0, coeffs[ix].0))
                            .expect("couldn't write polynomial?");
                    }
                }
            }
        }
        f.write_str(&format!(" mod {:}", self.field_mod))
    }
}

// impl Mul for FiniteFieldPolynomial {
//     type Output = Self;

//     fn mul(self, rhs: Self) -> Self::Output {
//         let mut ret: Vec<FiniteField> = (1..=(self.coeffs.len() + rhs.coeffs.len() - 1))
//             .into_iter()
//             .map(|_| (0, self.field_mod).into())
//             .collect();
//         for deg_1 in 0..self.coeffs.len() {
//             for deg_2 in 0..rhs.coeffs.len() {
//                 let coeff_1 = self.coeffs[deg_1];
//                 let coeff_2 = rhs.coeffs[deg_2];
//                 if let Some(e) = ret.get_mut(deg_1 + deg_2) {
//                     *e += coeff_1 * coeff_2;
//                 }
//             }
//         }
//         FiniteFieldPolynomial {
//             coeffs: ret,
//             field_mod: self.field_mod,
//         }
//     }
// }

impl Mul for FiniteFieldPolynomial {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        if self.field_mod != rhs.field_mod {
            panic!("Tried to multiply polynomials of different fields.")
        }
        let mut hm = HashMap::new();
        for (d1, c1) in self.coeffs.into_iter() {
            for (d2, c2) in rhs.coeffs.iter() {
                let e = hm.entry(d1 + d2).or_insert((0_u32, self.field_mod).into());
                *e += c1 * c2;
            }
        }
        let mut p = FiniteFieldPolynomial {
            coeffs: hm,
            degree: self.degree * rhs.degree,
            field_mod: self.field_mod,
        };
        p.clean();
        p
    }
}

impl Mul<&FiniteFieldPolynomial> for &FiniteFieldPolynomial {
    type Output = FiniteFieldPolynomial;

    fn mul(self, rhs: &FiniteFieldPolynomial) -> Self::Output {
        if self.field_mod != rhs.field_mod {
            panic!("Tried to multiply polynomials of different fields.")
        }
        let mut hm = HashMap::new();
        for (d1, c1) in self.coeffs.iter() {
            for (d2, c2) in rhs.coeffs.iter() {
                let e = hm.entry(d1 + d2).or_insert((0_u32, self.field_mod).into());
                *e += *c1 * c2;
            }
        }
        let mut p = FiniteFieldPolynomial {
            coeffs: hm,
            degree: self.degree * rhs.degree,
            field_mod: self.field_mod,
        };
        p.clean();
        p
    }
}

impl MulAssign for FiniteFieldPolynomial {
    fn mul_assign(&mut self, rhs: Self) {
        if self.field_mod != rhs.field_mod {
            panic!("Tried to multiply polynomials of different fields.")
        }
        let mut hm = HashMap::new();
        for (d1, c1) in self.coeffs.iter() {
            for (d2, c2) in rhs.coeffs.iter() {
                let e = hm.entry(d1 + d2).or_insert((0_u32, self.field_mod).into());
                *e += *c1 * c2;
            }
        }
        self.coeffs = hm;
        self.clean();
    }
}

// impl Div for FiniteFieldPolynomial {
//     type Output = (Self, Self);

//     fn div(self, rhs: Self) -> Self::Output {
//         self.division(&rhs)
//     }
// }

impl<'a, 'b> Add<&'a FiniteFieldPolynomial> for &'b FiniteFieldPolynomial {
    type Output = FiniteFieldPolynomial;

    fn add(self, rhs: &'a FiniteFieldPolynomial) -> Self::Output {
        if self.field_mod != rhs.field_mod {
            panic!("Tried to add polynomials of different fields.")
        }
        let mut hm = HashMap::new();
        for (k, v) in self.coeffs.iter() {
            let e = hm.entry(*k).or_insert(FiniteField::new(0, self.field_mod));
            *e += v;
        }
        for (k, v) in rhs.coeffs.iter() {
            let e = hm.entry(*k).or_insert(FiniteField::new(0, self.field_mod));
            *e += v;
        }
        let mut p = FiniteFieldPolynomial {
            coeffs: hm,
            degree: self.degree.max(rhs.degree),
            field_mod: self.field_mod,
        };
        p.clean();
        p
    }
}

impl Add for FiniteFieldPolynomial {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        if self.field_mod != rhs.field_mod {
            panic!("Tried to add polynomials of different fields.")
        }
        let mut hm = rhs.coeffs;
        for (k, v) in self.coeffs {
            let e = hm.entry(k).or_insert(FiniteField::new(0, self.field_mod));
            *e += v;
        }
        let hm = clear_zeros(hm);
        let mut p = FiniteFieldPolynomial {
            coeffs: hm,
            degree: self.degree.max(rhs.degree),
            field_mod: self.field_mod,
        };
        p.clean();
        p
    }
}

impl AddAssign for FiniteFieldPolynomial {
    fn add_assign(&mut self, rhs: Self) {
        if self.field_mod != rhs.field_mod {
            panic!("Tried to add polynomials of different fields.")
        }
        for (k, v) in rhs.coeffs.into_iter() {
            let e = self
                .coeffs
                .entry(k)
                .or_insert((0_u32, self.field_mod).into());
            *e += v;
        }
        self.clean();
    }
}

impl AddAssign<&FiniteFieldPolynomial> for FiniteFieldPolynomial {
    fn add_assign(&mut self, rhs: &FiniteFieldPolynomial) {
        if self.field_mod != rhs.field_mod {
            panic!("Tried to add polynomials of different fields.")
        }
        for (k, v) in rhs.coeffs.iter() {
            let e = self
                .coeffs
                .entry(*k)
                .or_insert((0_u32, self.field_mod).into());
            *e += v;
        }
        self.clean();
    }
}

impl Sub for FiniteFieldPolynomial {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let mut hm = self.coeffs;
        for (k, v) in rhs.coeffs {
            let e = hm.entry(k).or_insert(FiniteField::new(0, self.field_mod));
            *e -= v;
        }
        let hm = clear_zeros(hm);
        let mut p = FiniteFieldPolynomial {
            coeffs: hm,
            degree: self.degree.max(rhs.degree),
            field_mod: self.field_mod,
        };
        p.clean();
        p
    }
}

impl Sub<&FiniteFieldPolynomial> for &FiniteFieldPolynomial {
    type Output = FiniteFieldPolynomial;

    fn sub(self, rhs: &FiniteFieldPolynomial) -> Self::Output {
        if self.field_mod != rhs.field_mod {
            panic!("Tried to add polynomials of different fields.")
        }
        let mut hm = HashMap::new();
        for (k, v) in self.coeffs.iter() {
            let e = hm.entry(*k).or_insert(FiniteField::new(0, self.field_mod));
            *e += v;
        }
        for (k, v) in rhs.coeffs.iter() {
            let e = hm.entry(*k).or_insert(FiniteField::new(0, self.field_mod));
            *e -= *v;
        }
        let mut p = FiniteFieldPolynomial {
            coeffs: hm,
            degree: self.degree.max(rhs.degree),
            field_mod: self.field_mod,
        };
        p.clean();
        p
    }
}

impl SubAssign for FiniteFieldPolynomial {
    fn sub_assign(&mut self, rhs: Self) {
        if self.field_mod != rhs.field_mod {
            panic!("Tried to add polynomials of different fields.")
        }
        for (k, v) in rhs.coeffs.into_iter() {
            let e = self
                .coeffs
                .entry(k)
                .or_insert((0_u32, self.field_mod).into());
            *e -= v;
        }
        self.clean();
    }
}

impl SubAssign<&FiniteFieldPolynomial> for FiniteFieldPolynomial {
    fn sub_assign(&mut self, rhs: &FiniteFieldPolynomial) {
        if self.field_mod != rhs.field_mod {
            panic!("Tried to add polynomials of different fields.")
        }
        for (k, v) in rhs.coeffs.iter() {
            let e = self
                .coeffs
                .entry(*k)
                .or_insert((0_u32, self.field_mod).into());
            *e -= v;
        }
        self.clean();
    }
}


impl From<&[(usize, FiniteField)]> for FiniteFieldPolynomial {
    fn from(value: &[(usize, FiniteField)]) -> Self {
        if value.len() == 0 {
            panic!("Cannot create polynomial without coefficients")
        }
        let mut max_degree = 0;
        let mut hm = HashMap::with_capacity(value.len());
        let p = value[0].1 .1;
        for (degree, coeff) in value {
            max_degree = max_degree.max(*degree);
            let e = hm.entry(*degree).or_insert(FiniteField::new(0, p));
            *e += coeff;
        }
        FiniteFieldPolynomial {
            coeffs: hm,
            degree: max_degree,
            field_mod: p,
        }
    }
}

// impl Rem for &FiniteFieldPolynomial {
//     type Output = FiniteFieldPolynomial;

//     fn rem(self, rhs: Self) -> Self::Output {
//         let (_, r) = self.division(rhs);
//         r
//     }
// }

fn clear_zeros(coeffs: HashMap<usize, FiniteField>) -> HashMap<usize, FiniteField> {
    if coeffs.len() == 0 {
        return coeffs;
    }
    let mut p = 0;
    let mut new_coeffs: HashMap<usize, FiniteField> = coeffs
        .into_iter()
        .filter(|(k, v)| {
            p = v.1;
            v.0 != 0
        })
        .collect();
    if new_coeffs.len() == 0 {
        new_coeffs.insert(0, FiniteField::new(0, p));
    }
    new_coeffs
}

fn remove_trailing_zeros(coeffs: &mut Vec<FiniteField>) {
    let n = coeffs.len();
    if n <= 1 {
        return;
    }
    let mut is_trailing_zeros = coeffs[n - 1].0 == 0;
    while is_trailing_zeros {
        if let Some(coeff) = coeffs.pop() {
            if coeff.0 != 0 || coeffs.is_empty() {
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

    use crate::math::{
        finite_field::FiniteField,
        polynomial::{get_divisors, get_smallest_divisor},
    };

    use super::{remove_trailing_zeros, FiniteFieldPolynomial};

    #[test]
    fn test_polynomial_arithmetic() {
        let coeffs_1: [(usize, FiniteField); 4] = [
            (0, (1, 199).into()),
            (1, (2, 199).into()),
            (2, (0, 199).into()),
            (3, (3, 199).into()),
        ];
        let coeffs_2: [(usize, FiniteField); 4] = [
            (0, (1, 199).into()),
            (1, (3, 199).into()),
            (2, (0, 199).into()),
            (3, (3, 199).into()),
        ];
        let p1 = FiniteFieldPolynomial::from(&coeffs_1[..]);
        let p2 = FiniteFieldPolynomial::from(&coeffs_2[..]);
        println!("p1 : {:}", p1);
        println!("p2 : {:}", p2);
        println!("p1 + p2: {:}", &p1 + &p2);
        println!("p1 - p2: {:}", &p1 - &p2);
        println!("p1 * p2: {:}", p1 * p2);
        // let z = FiniteField(0, 199);
        // let p2 = FiniteFieldPolynomial {
        //     coeffs: vec![z, z + 5, z + 7],
        //     field_mod: 199,
        // };
        // dbg!(p1 * p2);
    }

    #[test]
    fn test_remove_trailing_zeros() {
        let z = FiniteField(0, 199);
        let mut r = vec![z + 1, z + 2, z + 3, z + 4, z + 5, z, z, z + 2, z + 0, z + 0];
        remove_trailing_zeros(&mut r);
        println!("r capacity: {:}", r.capacity());
        dbg!(r);
    }

    #[test]
    fn test_constructors() {
        // let p = FiniteFieldPolynomial::monomial((0, 199).into(), 3, 199);
        // let poly = FiniteFieldPolynomial::new(vec![(0, 199).into()]);
        // println!("{:}", p);
        // println!(
        //     "p degree: {:}, leading_coeff: {:}",
        //     p.deg(),
        //     p.leading_coeff()
        // );
    }
    #[test]
    fn test_polynomial_division() {
        // let z = FiniteField(0, 199);
        // let b = FiniteFieldPolynomial {
        //     coeffs: vec![z + 1, z + 3, z + 21, z + 22, z + 23],
        //     field_mod: 199,
        // };
        // let q = FiniteFieldPolynomial {
        //     coeffs: vec![z + 17, z + 9, z + 13],
        //     field_mod: 199,
        // };
        // let r = FiniteFieldPolynomial {
        //     coeffs: vec![z + 5, z + 11],
        //     field_mod: 199,
        // };
        // let a = b.clone() * q.clone() + r.clone();
        // let (q_comp, r_comp) = a.division(&b);
        // assert_eq!(q, q_comp);
        // assert_eq!(r, r_comp);
        // println!("quotient - {:}", q_comp);
        // println!("remainder - {:}", r_comp);
    }

    #[test]
    fn test_prime_divisors() {
        let a = 7_u32.pow(3);
        let b = 3_u32.pow(4);
        let c = 2;
        let d = 199_u32.pow(2);
        let n = a * b * c * d;
        dbg!(n);
        dbg!(get_smallest_divisor(n));
        dbg!(get_divisors(n));
        dbg!(get_divisors(1));
        dbg!(get_divisors(199));
    }

    #[test]
    fn test_polynomial_irreducibility() {
        // let p = 199;
        // let coeffs = vec![
        //     (1, p).into(),
        //     (0, p).into(),
        //     (0, p).into(),
        //     (0, p).into(),
        //     (1, p).into(),
        // ];
        // dbg!(1);
        // let f = FiniteFieldPolynomial::new(coeffs);
        // dbg!(2);
        // let t1 = FiniteFieldPolynomial::monomial((1, 2).into(), 10, 2);
        // dbg!(3);
        // let t2 = FiniteFieldPolynomial::monomial((1, 2).into(), 3, 2);
        // dbg!(4);
        // let t3 = FiniteFieldPolynomial::monomial((1, 2).into(), 0, 2);
        // dbg!(5);
        // let g = t1 + t2 + t3;
        // println!("is {:} irreducible: {:}", f, f.is_irreducible());
        // println!("is {:} irreducible: {:}", g, g.is_irreducible());
    }

    #[test]
    fn test_is_primitive() {
        // let tester = FiniteFieldPolynomial::new(vec![(1, 3).into(), (0, 3).into(), (1, 3).into()]);
        // println!("tester: {:}", tester);
        // println!("irreducible: {:}", tester.is_irreducible());
        // println!("primitive: {:}", tester.is_primitive());
    }
}
