use fxhash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::{
    cmp::Ordering,
    collections::HashMap,
    fmt::Display,
    hash::Hash,
    ops::{Add, AddAssign, Div, Mul, MulAssign, Rem, Sub, SubAssign},
    str::FromStr,
};

use super::finite_field::{FFRep, FiniteField};

pub type PolyDegree = u32;

fn get_smallest_divisor(n: u32) -> Option<u32> {
    for x in 2..=(n / 2) {
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

/// Polynomial in single indeterminate. Uses a sparse representation.
/// Uses a small integer for the degree, so the maximum degree represented is 255
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct FFPolynomial {
    pub coeffs: FxHashMap<PolyDegree, FiniteField>,
    degree: PolyDegree,
    pub field_mod: u32,
}

impl FFPolynomial {
    /// Removes zeros and updates the degree
    pub fn clean(&mut self) {
        let mut new_degree = 0;
        let mut hm: FxHashMap<PolyDegree, FiniteField> = self
            .coeffs
            .clone()
            .into_iter()
            .filter(|(_, c)| c.0 != 0)
            .collect();
        if hm.len() == 0 {
            hm.insert(0, (0, self.field_mod).into());
        }
        for (d, _) in hm.iter() {
            new_degree = new_degree.max(*d);
        }
        self.degree = new_degree;
        self.coeffs = hm;
    }

    pub fn get_number(&self) -> FFRep {
        let mut num = 0;
        for (deg, coeff) in self.coeffs.iter() {
            num += coeff.0 * (self.field_mod.pow(*deg));
        }
        num
    }

    pub fn to_number(self) -> FFRep {
        let mut num = 0;
        for (deg, coeff) in self.coeffs.into_iter() {
            num += coeff.0 * (self.field_mod.pow(deg));
        }
        num
    }

    pub fn from_number(number: FFRep, field_mod: FFRep) -> Self {
        if number == 0 {
            return FFPolynomial::zero(field_mod);
        }
        let mut num = number;
        let mut deg_and_coeffs = Vec::new();
        let mut pow = 0;
        while num != 0 {
            let modder = field_mod.pow(pow + 1);
            let quotienter = field_mod.pow(pow);
            let remainder = num % modder;
            let coeff = remainder / quotienter;
            if remainder != 0 {
                deg_and_coeffs.push((pow, FiniteField::new(coeff, field_mod)));
                num -= remainder;
            }
            pow += 1;
        }
        FFPolynomial::from(&deg_and_coeffs[..])
    }

    pub fn evaluate(&self, x: &FiniteField) -> FiniteField {
        // This uses a very naive algorithm, it is unclear if Horner's
        // method would be better because this is using a sparse representation.
        let mut out = FiniteField::from((0, self.field_mod));
        for (d, c) in self.coeffs.iter() {
            out += *c * x.pow(*d as u32);
        }
        out
    }

    /// Returns the number of terms in the polynomial.
    /// `{x^1 + 2}.len() == 2`
    pub fn len(&self) -> usize {
        self.coeffs.len()
    }

    /// Creates the unique Lagrange interpolation polynomial through the provided (x_i, y_i) pairs.
    pub fn interpolation(points: Vec<(FiniteField, FiniteField)>) -> Self {
        if points.len() == 0 {
            panic!("Trying to interpolate on nothing?")
        }
        let p = points[0].0 .1;
        let linear_factors: Vec<FFPolynomial> = (0..points.len())
            .map(|k| {
                FFPolynomial::monomial((1, p).into(), 1)
                    - FFPolynomial::constant(points[k].0 .0.clone(), p)
            })
            .collect();
        let mut lagrange = FFPolynomial::zero(p);
        for k in 0..points.len() {
            let mut numerator = FFPolynomial::constant(1, p);
            let mut denominator = FiniteField::new(1, p);
            for m in 0..points.len() {
                if m != k {
                    numerator *= &linear_factors[m];
                    denominator *= points[k].0 - points[m].0;
                }
            }
            numerator.scale(&denominator.modular_inverse());
            numerator.scale(&points[k].1);
            lagrange += numerator;
        }
        lagrange
    }

    /// Multiplies all coefficients by `scalar`
    pub fn scale(&mut self, scalar: &FiniteField) {
        for (_, v) in self.coeffs.iter_mut() {
            *v *= scalar;
        }
    }

    /// Creates a new zero polynomial
    pub fn zero(field_mod: u32) -> Self {
        let z = FiniteField::new(0, field_mod);
        let mut hm = FxHashMap::default();
        hm.insert(0, z);
        FFPolynomial {
            coeffs: hm,
            degree: 0,
            field_mod,
        }
    }

    pub fn constant(c: u32, field_mod: u32) -> Self {
        let z = FiniteField::new(c, field_mod);
        let mut hm = FxHashMap::default();
        hm.insert(0, z);
        FFPolynomial {
            coeffs: hm,
            degree: 0,
            field_mod,
        }
    }

    pub fn leading_coeff(&self) -> FiniteField {
        if let Some(lc) = self.coeffs.get(&self.degree) {
            lc.clone()
        } else {
            println!("dirty degree polynomial: {:}", self);
            panic!("Dirty degree polynomial?")
        }
    }

    pub fn gcd(&self, rhs: &FFPolynomial) -> FFPolynomial {
        if rhs.is_zero() {
            self.clone()
        } else {
            let mut r = rhs.gcd(&(self % rhs));
            r.clean();
            r
        }
    }

    /// If `a = self` and `b = rhs`, returns `u, v, r` such that `u * a + v * b = r`, where we stop the process when the degree of `r` is less than `remainder degree`
    pub fn partial_gcd(
        &self,
        rhs: &FFPolynomial,
        remainder_degree: PolyDegree,
    ) -> (FFPolynomial, FFPolynomial, FFPolynomial) {
        let p = self.field_mod;
        let mut r = self.clone();
        let mut new_r = rhs.clone();
        let mut s = FFPolynomial::constant(1, self.field_mod);
        let mut new_s = FFPolynomial::zero(p);
        let mut t = new_s.clone();
        let mut new_t = s.clone();
        r.clean();
        new_r.clean();
        while r.degree >= remainder_degree {
            let (tmp, _) = &r / &new_r;
            let old_r = r.clone();
            r = new_r.clone();
            new_r = old_r - &tmp * &r;
            new_r.clean();

            let old_s = s.clone();
            s = new_s.clone();
            new_s = old_s - &tmp * &s;
            new_s.clean();

            let old_t = t.clone();
            t = new_t.clone();
            new_t = old_t - &tmp * &t;
            new_t.clean();
        }
        (s, t, r)
    }

    /// Implementation of the [Rabin Irreducibility Test ](https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields#Rabin's_test_of_irreducibility) to determine if `self` can be factored into two polynomials p(x) q(x) each of degree less than `self.deg()`.
    pub fn is_irreducible(&self) -> bool {
        let n = self.degree() as u32;
        let divisors = get_divisors(n);
        let powers: Vec<u32> = divisors
            .into_iter()
            .filter(|p| *p != n)
            .map(|p| n / p)
            .collect();

        for power in powers {
            let tmp_degree = self.field_mod.pow(power);
            let tmp_term_1 =
                FFPolynomial::monomial((1, self.field_mod).into(), tmp_degree.try_into().unwrap());
            let tmp_term_2 = FFPolynomial::monomial((-1, self.field_mod).into(), 1);
            let t = tmp_term_1 + tmp_term_2;
            let tmp = &t % &self;
            let g = self.gcd(&tmp);
            if g != FFPolynomial::monomial((1, self.field_mod).into(), 0) {
                return false;
            }
        }
        let coeff_1: FiniteField = (1, self.field_mod).into();
        let deg_1 = self.field_mod.pow(n).try_into().unwrap();
        let tmp_term_1 = FFPolynomial::monomial(coeff_1, deg_1);
        let tmp_term_2 = FFPolynomial::monomial((-1, self.field_mod).into(), 1);
        let t = tmp_term_1 + tmp_term_2;
        let g = &t % &self;
        g.is_zero()
    }

    /// This is based on a wikipedia entry [Primitive Polynomials in Finite Field](https://en.wikipedia.org/wiki/Primitive_polynomial_(field_theory))
    pub fn is_primitive(&self) -> bool {
        let d = self.degree;
        let upper_bound = self.field_mod.pow(d as u32) - 1;
        for n in 1..upper_bound {
            let buf = [
                (n as PolyDegree, (1, self.field_mod).into()),
                (0, (-1, self.field_mod).into()),
            ];
            let g = FFPolynomial::from(&buf[..]);
            let r = &g % self;
            if r.is_zero() {
                return false;
            }
        }
        let buf = [
            (upper_bound as PolyDegree, (1, self.field_mod).into()),
            (0, (-1, self.field_mod).into()),
        ];
        let g = FFPolynomial::from(&buf[..]);
        let r = &g % self;
        r.is_zero()
    }

    pub fn is_zero(&self) -> bool {
        if self.coeffs.len() == 0 {
            return true;
        }
        for (_, c) in self.coeffs.iter() {
            if c.0 != 0 {
                return false;
            }
        }
        true
    }

    pub fn is_one(&self) -> bool {
        let mut is_constant_one = false;
        let mut everything_else_zero = true;
        for (k, v) in self.coeffs.iter() {
            if *k == 0 {
                if v.0 == 1 {
                    is_constant_one = true;
                }
            } else {
                if v.0 != 0 {
                    everything_else_zero = false;
                }
            }
        }
        is_constant_one && everything_else_zero
    }

    pub fn degree(&self) -> PolyDegree {
        self.degree
    }

    pub fn monomial(coefficient: FiniteField, degree: PolyDegree) -> FFPolynomial {
        let buf = [(degree, coefficient)];
        FFPolynomial::from(&buf[..])
    }
}

impl PartialOrd for FFPolynomial {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.degree > other.degree {
            Some(Ordering::Greater)
        } else if self.degree < other.degree {
            Some(Ordering::Less)
        } else {
            let mut ret = Ordering::Equal;
            for ix in 0..=self.degree {
                let new_ix = self.degree - ix;
                let zero = FiniteField::new(0, self.field_mod);
                let left_coeff = self.coeffs.get(&new_ix).unwrap_or(&zero);
                let right_coeff = other.coeffs.get(&new_ix).unwrap_or(&zero);
                if new_ix > 0 {
                    if left_coeff > right_coeff {
                        ret = Ordering::Greater;
                        break;
                    }
                    if left_coeff < right_coeff {
                        ret = Ordering::Less;
                        break;
                    }
                } else {
                    ret = left_coeff.cmp(&right_coeff);
                    break;
                }
            }
            Some(ret)
        }
    }
}

impl Ord for FFPolynomial {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

impl Display for FFPolynomial {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.coeffs.len() == 0 {
            return f.write_str("zero polynomial encountered.");
        }
        let mut coeffs: Vec<(PolyDegree, FiniteField)> = self.coeffs.clone().into_iter().collect();
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
        // f.write_str(&format!(" mod {:}", self.field_mod))
        Ok(())
    }
}

impl Hash for FFPolynomial {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        let mut v: Vec<(PolyDegree, FiniteField)> = self.coeffs.clone().into_iter().collect();
        v.sort_by(|x, y| x.0.cmp(&y.0));
        for (d, c) in v {
            d.hash(state);
            c.hash(state);
        }
        self.degree.hash(state);
        self.field_mod.hash(state);
    }
}

impl Mul for FFPolynomial {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        &self * &rhs
    }
}

impl Mul<&FFPolynomial> for &FFPolynomial {
    type Output = FFPolynomial;

    fn mul(self, rhs: &FFPolynomial) -> Self::Output {
        if self.field_mod != rhs.field_mod {
            panic!("Tried to multiply polynomials of different fields.")
        }
        let mut hm = FxHashMap::default();
        for (d1, c1) in self.coeffs.iter() {
            for (d2, c2) in rhs.coeffs.iter() {
                let e = hm.entry(d1 + d2).or_insert((0_u32, self.field_mod).into());
                *e += *c1 * c2;
            }
        }
        let mut p = FFPolynomial {
            coeffs: hm,
            degree: self.degree + rhs.degree,
            field_mod: self.field_mod,
        };
        p.clean();
        p
    }
}

impl MulAssign for FFPolynomial {
    fn mul_assign(&mut self, rhs: Self) {
        *self *= &rhs;
    }
}

impl MulAssign<&FFPolynomial> for FFPolynomial {
    fn mul_assign(&mut self, rhs: &FFPolynomial) {
        if self.field_mod != rhs.field_mod {
            panic!("Tried to multiply polynomials of different fields.")
        }
        let mut hm = FxHashMap::default();
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

impl Div<&FFPolynomial> for &FFPolynomial {
    type Output = (FFPolynomial, FFPolynomial);

    fn div(self, rhs: &FFPolynomial) -> Self::Output {
        if rhs.degree > self.degree {
            return (FFPolynomial::zero(self.field_mod), self.clone());
        }
        if rhs.degree == 0 {
            let inv = rhs.leading_coeff().modular_inverse();
            let mut out = self.clone();
            out.scale(&inv);
            return (out, FFPolynomial::zero(self.field_mod));
        }
        let mut quotient = FFPolynomial::zero(self.field_mod);
        let mut remainder = self.clone();
        let d = rhs.degree;
        let lc = rhs.leading_coeff();
        while remainder.degree >= d && remainder.is_zero() == false {
            let coeff_s = remainder.leading_coeff() * lc.modular_inverse();
            let deg_s = remainder.degree - d;
            let s = FFPolynomial::monomial(coeff_s, deg_s);
            let remainder_sub = &s * &rhs;
            quotient = quotient + s;
            remainder -= remainder_sub;
            remainder.clean();
        }
        quotient.clean();
        remainder.clean();
        (quotient, remainder)
    }
}

impl Div for FFPolynomial {
    type Output = (FFPolynomial, FFPolynomial);

    fn div(self, rhs: Self) -> Self::Output {
        &self / &rhs
    }
}

impl Rem<&FFPolynomial> for &FFPolynomial {
    type Output = FFPolynomial;

    fn rem(self, rhs: &FFPolynomial) -> Self::Output {
        let (_, r) = self / rhs;
        r
    }
}

impl Add<&FFPolynomial> for &FFPolynomial {
    type Output = FFPolynomial;

    fn add(self, rhs: &FFPolynomial) -> Self::Output {
        if self.field_mod != rhs.field_mod {
            panic!("Tried to add polynomials of different fields.")
        }
        let mut hm = self.coeffs.clone();
        // for (k, v) in self.coeffs.iter() {
        //     let e = hm.entry(*k).or_insert(FiniteField::new(0, self.field_mod));
        //     *e += v;
        // }
        for (k, v) in rhs.coeffs.iter() {
            let mut cur_entry = hm
                .get(k)
                .cloned()
                .unwrap_or(FiniteField::new(0, self.field_mod));
            cur_entry += v;
            if cur_entry.0 == 0 {
                hm.remove(k);
            } else {
                hm.insert(*k, cur_entry);
            }
        }
        FFPolynomial {
            coeffs: hm,
            degree: self.degree.max(rhs.degree),
            field_mod: self.field_mod,
        }
    }
}

impl Add for FFPolynomial {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        &self + &rhs
    }
}

impl AddAssign for FFPolynomial {
    fn add_assign(&mut self, rhs: Self) {
        *self += &rhs;
    }
}

impl Add<u32> for &FFPolynomial {
    type Output = FFPolynomial;

    fn add(self, rhs: u32) -> Self::Output {
        let p = FFPolynomial::constant(rhs, self.field_mod);
        self + &p
    }
}

impl AddAssign<u32> for FFPolynomial {
    fn add_assign(&mut self, rhs: u32) {
        let c: FiniteField = (rhs, self.field_mod).into();
        let e = self.coeffs.entry(0).or_insert((0, self.field_mod).into());
        *e += c;
        if e.0 == 0 {
            self.coeffs.remove(&0);
        }
        self.clean();
    }
}

impl AddAssign<&FFPolynomial> for FFPolynomial {
    fn add_assign(&mut self, rhs: &FFPolynomial) {
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

impl Sub for FFPolynomial {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        &self - &rhs
    }
}

impl Sub<&FFPolynomial> for &FFPolynomial {
    type Output = FFPolynomial;

    fn sub(self, rhs: &FFPolynomial) -> Self::Output {
        if self.field_mod != rhs.field_mod {
            panic!("Tried to add polynomials of different fields.")
        }
        let mut hm = self.coeffs.clone();
        for (k, v) in rhs.coeffs.iter() {
            let e = hm.entry(*k).or_insert(FiniteField::new(0, self.field_mod));
            *e -= *v;
        }
        let mut p = FFPolynomial {
            coeffs: hm,
            degree: self.degree.max(rhs.degree),
            field_mod: self.field_mod,
        };
        p.clean();
        p
    }
}

impl SubAssign for FFPolynomial {
    fn sub_assign(&mut self, rhs: Self) {
        *self -= &rhs;
    }
}

impl SubAssign<&FFPolynomial> for FFPolynomial {
    fn sub_assign(&mut self, rhs: &FFPolynomial) {
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

impl From<&[(PolyDegree, FiniteField)]> for FFPolynomial {
    fn from(value: &[(PolyDegree, FiniteField)]) -> Self {
        if value.len() == 0 {
            panic!("Cannot create polynomial without coefficients")
        }
        let mut max_degree = 0;

        let mut hm = FxHashMap::default();
        let p = value[0].1 .1;
        for (degree, coeff) in value.iter().filter(|(_, c)| c.0 != 0) {
            max_degree = max_degree.max(*degree);
            let e = hm.entry(*degree).or_insert(FiniteField::new(0, p));
            *e += coeff;
        }
        FFPolynomial {
            coeffs: hm,
            degree: max_degree,
            field_mod: p,
        }
    }
}

impl From<(PolyDegree, FiniteField)> for FFPolynomial {
    fn from(value: (PolyDegree, FiniteField)) -> Self {
        let mut coeffs = FxHashMap::default();
        coeffs.insert(value.0, value.1);
        FFPolynomial {
            coeffs,
            degree: value.0,
            field_mod: value.1 .1,
        }
    }
}

#[derive(Debug)]
pub struct PolyParseError {}
impl FromStr for FFPolynomial {
    type Err = PolyParseError;

    /// String must be of the form: "a_d*x^d + a_{d-1}x^{d-1} + ... + a_0*x^0 % prime". terms are split at the plus and the percentage sign indicates the finite field being used.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let first_split: Vec<&str> = s.split("%").collect();
        if first_split.len() != 2 {
            println!("need to indicate the finite field being used.");
            return Err(PolyParseError {});
        }
        let p_str = first_split[1].trim();
        let p: u32 = p_str.parse().expect("could not parse prime field.");
        let poly_str = first_split[0].trim();
        let terms: Vec<&str> = poly_str.split("+").collect();
        let mut coefficients = Vec::new();
        for term in terms {
            let term_split: Vec<&str> = term.split("*").collect();
            let coeff: u32 = term_split[0]
                .trim()
                .parse()
                .expect("could not parse coefficient.");
            let deg_split: Vec<&str> = term_split[1].split("^").collect();
            let deg: PolyDegree = deg_split[1].trim().parse().expect("could not parse degree");
            coefficients.push((deg, FiniteField::new(coeff, p)));
        }
        let poly = FFPolynomial::from(&coefficients[..]);
        Ok(poly)
    }
}

#[cfg(test)]
mod tests {
    use std::{collections::HashSet, str::FromStr};

    use crate::math::{
        finite_field::FiniteField,
        polynomial::{get_divisors, PolyDegree},
    };

    use super::FFPolynomial;

    #[test]
    fn test_polynomial_arithmetic() {
        let coeffs_1: [(PolyDegree, FiniteField); 4] = [
            (0, (1, 199).into()),
            (1, (7, 199).into()),
            (2, (0, 199).into()),
            (3, (3, 199).into()),
        ];
        let coeffs_2: [(PolyDegree, FiniteField); 4] = [
            (0, (1, 199).into()),
            (1, (3, 199).into()),
            (2, (0, 199).into()),
            (3, (3, 199).into()),
        ];
        let coeffs_3: [(PolyDegree, FiniteField); 3] = [
            (0, (2, 199).into()),
            (1, (88, 199).into()),
            (2, (5, 199).into()),
        ];

        let mut p1 = FFPolynomial::from(&coeffs_1[..]);
        let mut p2 = FFPolynomial::from(&coeffs_2[..]);
        let mut p3 = FFPolynomial::from(&coeffs_3[..]);
        let mut p4 = &(&p1 * &p2) + &p3;

        p1.clean();
        p2.clean();
        p3.clean();
        p4.clean();
        let (q, r) = &p4 / &p2;

        let mut set = HashSet::new();
        set.insert(p1.clone());
        assert!(set.contains(&p1));
        assert_eq!(q, p1);
        assert_eq!(r, p3);
    }

    #[test]
    fn test_poly_cmp() {
        let p1 = FFPolynomial::monomial((1, 3).into(), 1) + FFPolynomial::constant(1, 3);
        let p2 = FFPolynomial::monomial((1, 3).into(), 1);
        println!("p1: {:}", p1);
        println!("p2: {:}", p2);
        println!("p1.cmp(&p2): {:?}", p1.cmp(&p2));
    }

    #[test]
    fn test_prime_divisors() {
        let a = 7_u32.pow(3);
        let b = 3_u32.pow(4);
        let c = 2;
        let d = 199_u32.pow(2);
        let n = a * b * c * d;
        let n_divisors = vec![2, 3, 7, 199];
        assert_eq!(n_divisors, get_divisors(n));
        assert_eq!(Vec::<u32>::new(), get_divisors(1));
        assert_eq!(vec![199], get_divisors(199));
        assert_eq!(vec![2], get_divisors(4));
    }

    #[test]
    fn test_polynomial_irreducibility() {
        let p = 53;
        let coeffs = vec![
            (0, (1, p).into()),
            (1, (0, p).into()),
            (2, (0, p).into()),
            (3, (0, p).into()),
            (4, (1, p).into()),
        ];

        let f = FFPolynomial::from(&coeffs[..]);
        let t1 = FFPolynomial::monomial((1, 2).into(), 10);

        let t2 = FFPolynomial::monomial((1, 2).into(), 3);

        let t3 = FFPolynomial::monomial((1, 2).into(), 0);

        let g = t1 + t2 + t3;
        assert!(f.is_irreducible() == false);
        assert!(g.is_irreducible());

        let poly = FFPolynomial::from_str("1*x^2 + 2*x^1 + 2*x^0 % 3").unwrap();
        println!("is poly irreducible? {:}", poly.is_irreducible());
    }

    #[test]
    fn test_is_primitive() {
        let buf = [(0, (1, 3).into()), (2, (1, 3).into())];
        let primitive_coeffs = [(2, (1, 3).into()), (1, (2, 3).into()), (0, (2, 3).into())];
        let tester = FFPolynomial::from(&buf[..]);
        let primitive_poly = FFPolynomial::from(&primitive_coeffs[..]);
        assert!(tester.is_irreducible());
        assert!(tester.is_primitive() == false);
        assert!(primitive_poly.is_irreducible());
        assert!(primitive_poly.is_primitive());
    }

    #[test]
    fn test_constant_addition() {
        let one = FFPolynomial::constant(1, 2);
        println!("one: {:}", one);
        println!("one + one = {:}", &one + 1);
    }

    #[test]
    fn test_lagrange_interpolation() {
        let p = 199_u32;
        let x_vals: Vec<FiniteField> =
            vec![(1, p).into(), (3, p).into(), (7, p).into(), (9, p).into()];
        let y_vals: Vec<FiniteField> = vec![
            (11, p).into(),
            (13, p).into(),
            (17, p).into(),
            (97, p).into(),
        ];
        let l = FFPolynomial::interpolation(
            x_vals
                .clone()
                .into_iter()
                .zip(y_vals.clone().into_iter())
                .collect(),
        );
        println!("l = {:}", l);
        for ix in 0..x_vals.len() {
            let y_computed = l.evaluate(&x_vals[ix]);
            assert_eq!(y_computed, y_vals[ix]);
        }
    }

    #[test]
    fn test_partial_gcd() {
        let p = 199_u32;
        let x_vals: Vec<FiniteField> =
            vec![(1, p).into(), (3, p).into(), (7, p).into(), (9, p).into()];
        let mut y_vals: Vec<FiniteField> = vec![
            (11, p).into(),
            (13, p).into(),
            (17, p).into(),
            (97, p).into(),
        ];
        let l = FFPolynomial::interpolation(
            x_vals
                .clone()
                .into_iter()
                .zip(y_vals.clone().into_iter())
                .collect(),
        );
        y_vals[2].0 = 31;
        let l2 = FFPolynomial::interpolation(x_vals.into_iter().zip(y_vals.into_iter()).collect());
        let out = l.partial_gcd(&l2, 3);
        println!("l = {:}", l);
        println!("l2 = {:}", l2);
        println!("out 1 = {:}", out.0);
        println!("out 2 = {:}", out.1);
        println!("out 3 = {:}", out.2);
        let checker = out.0 * l + out.1 * l2;
        println!("checker = {:}", checker);
    }

    #[test]
    fn test_serde() {
        let buf = [(0, (1, 3).into()), (2, (1, 3).into())];
        let primitive_coeffs = [(2, (1, 3).into()), (1, (2, 3).into()), (0, (2, 3).into())];
        let tester = FFPolynomial::from(&buf[..]);
        let _primitive_poly = FFPolynomial::from(&primitive_coeffs[..]);
        let s = serde_json::to_string(&tester).expect("serialized");
        println!("s: {:}", s);
        let _f: FFPolynomial = serde_json::from_str(&s).expect("could not deserialize.");
    }

    #[test]
    fn test_poly_from_str() {
        let s = "200 * x^3 + 7 * x ^ 2 + 13 * x^0 %199";
        let poly = FFPolynomial::from_str(s).unwrap();
        println!("parsed poly: {:} mod {:}", poly, poly.field_mod);
    }

    #[test]
    fn test_to_and_from_number() {
        // p = 7
        // x = 1 * 7^4 + 2 * 7^3 + 3 * 7^2 + 4 *7^1 + 5 * 7^0 == 3267
        let p = 7;
        let poly = FFPolynomial::from_number(3267, p);
        println!("poly: {:}", poly);
        let num = poly.get_number();
        println!("num: {:}", num);
    }
}
