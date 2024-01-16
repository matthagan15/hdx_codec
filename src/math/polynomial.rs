use core::time;
use deepsize::DeepSizeOf;
use gcd::binary_u32;
use serde::de;
use std::{
    hash::Hash,
    cmp::Ordering,
    collections::{HashMap, HashSet},
    fmt::Display,
    ops::{Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Rem, Sub, SubAssign},
    thread,
};

use super::{
    finite_field::FiniteField,
    group_ring_field::{self, Field, Ring},
};

#[derive(Debug, Clone, Hash, PartialEq, Eq, DeepSizeOf)]
pub struct QuotientPoly {
    pub poly: FiniteFieldPolynomial,
    /// the q in a = q * b + r
    pub quotient: FiniteFieldPolynomial,
    pub field_mod: u32,
}

impl QuotientPoly {
    /// creates a new zero polynomial with entries in F_{field_mod}. quotient is the polynomial you are modding by.
    pub fn zero(field_mod: u32, quotient: FiniteFieldPolynomial) -> Self {
        let poly = FiniteFieldPolynomial::zero(field_mod);
        // Currently do not run check to save time
        if poly.field_mod != quotient.field_mod {
            panic!("Polynomials defined over different fields.")
        }
        // TODO: need to make the divisor primitive or check if it is primitive.
        let p = poly.field_mod;
        let r = &poly % &quotient;
        QuotientPoly {
            poly: r,
            quotient,
            field_mod: p,
        }
    }

    pub fn monomial(coeff: u32, degree: usize, quotient: FiniteFieldPolynomial) -> Self {
        let p = FiniteFieldPolynomial::monomial((coeff, quotient.field_mod).into(), degree);
        let r = &p % &quotient;
        QuotientPoly { poly: r, quotient, field_mod: p.field_mod }
    }

    pub fn constant(coeff: u32, quotient: FiniteFieldPolynomial) -> QuotientPoly {
        let buf = [(0, (coeff, quotient.field_mod).into())];
        let n = quotient.field_mod;
        QuotientPoly {
            poly: FiniteFieldPolynomial::from(&buf[..]),
            quotient,
            field_mod: n,
        }
    }

    /// Returns the b such that a * b = 1 mod q. Panics if no inverse exists
    pub fn mul_inv(&self) -> Self {
        let g = self.poly.gcd(&self.quotient);
        if g.degree > 0 {
            panic!("Cannot find inverse as gcd is not constant.")
        }
        let g_inv = g.leading_coeff().modular_inverse();
        let p = self.field_mod;
        let mut r = self.poly.clone();
        let mut new_r = self.quotient.clone();
        let mut s = FiniteFieldPolynomial::constant(1, self.field_mod);
        let mut new_s = FiniteFieldPolynomial::zero(p);
        let mut t = new_s.clone();
        let mut new_t = s.clone();
        while new_r.is_zero() == false {
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
        s.scale(&g_inv);
        QuotientPoly::from((s, self.quotient.clone()))
    }
}

impl Display for QuotientPoly {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&format!("{:} % [{:}]", self.poly, self.quotient))
    }
}

impl Add<&FiniteFieldPolynomial> for &QuotientPoly {
    type Output = QuotientPoly;

    fn add(self, rhs: &FiniteFieldPolynomial) -> Self::Output {
        let r = &(&self.poly + rhs) % &self.quotient;
        QuotientPoly {
            poly: r,
            quotient: self.quotient.clone(),
            field_mod: self.field_mod,
        }
    }
}

impl Add<u32> for &QuotientPoly {
    type Output = QuotientPoly;

    fn add(self, rhs: u32) -> Self::Output {
        let new_rhs = QuotientPoly::constant(rhs, self.quotient.clone());
        self + &new_rhs
    }
}

impl AddAssign<u32> for QuotientPoly {
    fn add_assign(&mut self, rhs: u32) {
        self.poly += rhs;
    }
}

impl Add<&QuotientPoly> for &QuotientPoly {
    type Output = QuotientPoly;

    fn add(self, rhs: &QuotientPoly) -> Self::Output {
        let out_poly = &(&self.poly + &rhs.poly) % &self.quotient;
        QuotientPoly {
            poly: out_poly,
            quotient: self.quotient.clone(),
            field_mod: self.field_mod,
        }
    }
}

impl Add for QuotientPoly {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        &self + &rhs.poly
    }
}

impl AddAssign for QuotientPoly {
    fn add_assign(&mut self, rhs: Self) {
        let out = (self as &QuotientPoly) + &rhs.poly;
        self.poly = out.poly;
    }
}

impl Sub for QuotientPoly {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let out_poly = self.poly - rhs.poly;
        let out_quotiented = &out_poly % &self.quotient;
        QuotientPoly {
            poly: out_quotiented,
            quotient: self.quotient,
            field_mod: self.field_mod,
        }
    }
}

impl Mul for QuotientPoly {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let out_poly = self.poly * rhs.poly;
        let out_quotiented = &out_poly % &self.quotient;
        QuotientPoly {
            poly: out_quotiented,
            quotient: self.quotient,
            field_mod: self.field_mod,
        }
    }
}

impl Mul<&QuotientPoly> for &QuotientPoly {
    type Output = QuotientPoly;

    fn mul(self, rhs: &QuotientPoly) -> Self::Output {
        let out_poly = &(&self.poly * &rhs.poly) % &self.quotient;
        QuotientPoly {
            poly: out_poly,
            quotient: self.quotient.clone(),
            field_mod: self.field_mod,
        }
    }
}

impl MulAssign for QuotientPoly {
    fn mul_assign(&mut self, rhs: Self) {
        let out_poly = &(&self.poly * &rhs.poly) % &self.quotient;
        self.poly = out_poly;
    }
}

impl Mul<&FiniteFieldPolynomial> for &QuotientPoly {
    type Output = QuotientPoly;

    fn mul(self, rhs: &FiniteFieldPolynomial) -> Self::Output {
        let out_poly = &(&self.poly * rhs) % &self.quotient;
        QuotientPoly {
            poly: out_poly,
            quotient: self.quotient.clone(),
            field_mod: self.field_mod,
        }
    }
}

impl SubAssign for QuotientPoly {
    fn sub_assign(&mut self, rhs: Self) {
        let out_poly = &(&self.poly - &rhs.poly) % &self.quotient;
        self.poly = out_poly;
    }
}

impl Ring for QuotientPoly {
    fn zero(&self) -> Self {
        QuotientPoly {
            poly: FiniteFieldPolynomial::zero(self.field_mod),
            quotient: self.quotient.clone(),
            field_mod: self.field_mod,
        }
    }

    fn one(&self) -> Self {
        let buf = [(0, (1, self.field_mod).into())];
        QuotientPoly {
            poly: FiniteFieldPolynomial::from(&buf[..]),
            quotient: self.quotient.clone(),
            field_mod: self.field_mod,
        }
    }

    fn additive_inv(&self) -> Self {
        QuotientPoly::zero(self.field_mod, self.quotient.clone()) - self.clone()
    }
}

impl From<(&FiniteFieldPolynomial, &FiniteFieldPolynomial)> for QuotientPoly {
    fn from(value: (&FiniteFieldPolynomial, &FiniteFieldPolynomial)) -> Self {
        let mut p = value.0 % value.1;
        let mut q = value.1.clone();
        p.clean();
        q.clean();
        let n = p.field_mod;
        QuotientPoly {
            poly: p,
            quotient: q,
            field_mod: n,
        }
    }
}

impl From<(FiniteFieldPolynomial, FiniteFieldPolynomial)> for QuotientPoly {
    fn from(value: (FiniteFieldPolynomial, FiniteFieldPolynomial)) -> Self {
        let n = value.0.field_mod;
        let mut p = value.0;
        let mut q = value.1;
        p.clean();
        q.clean();
        QuotientPoly {
            poly: &p % &q,
            quotient: q,
            field_mod: n,
        }
    }
}

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

/// Polynomial in single indeterminate. Uses dense storage (vec)
#[derive(Debug, Clone, PartialEq, Eq, DeepSizeOf)]
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

    pub fn evaluate(&self, x: &FiniteField) -> FiniteField {
        let mut out = FiniteField::from((0, self.field_mod));
        for (d, c) in self.coeffs.iter() {
            out += *c * x.pow(*d as u32);
        }
        out
    }

    pub fn scale(&mut self, scalar: &FiniteField) {
        for (_, v) in self.coeffs.iter_mut() {
            *v *= scalar;
        }
    }

    /// Creates a new zero polynomial
    pub fn zero(field_mod: u32) -> Self {
        let z = FiniteField::new(0, field_mod);
        let hm = HashMap::from([(0, z)]);
        FiniteFieldPolynomial {
            coeffs: hm,
            degree: 0,
            field_mod,
        }
    }

    pub fn constant(c: u32, field_mod: u32) -> Self {
        FiniteFieldPolynomial {
            coeffs: HashMap::from([(0, (c, field_mod).into())]),
            degree: 0,
            field_mod,
        }
    }

    pub fn leading_coeff(&self) -> FiniteField {
        self.coeffs
            .get(&self.degree)
            .expect("Polynomial degree is dirty?")
            .clone()
    }

    pub fn gcd(&self, rhs: &FiniteFieldPolynomial) -> FiniteFieldPolynomial {
        if rhs.is_zero() {
            self.clone()
        } else {
            let mut r = rhs.gcd(&(self % rhs));
            r.clean();
            r
        }
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
                FiniteFieldPolynomial::monomial((1, self.field_mod).into(), tmp_degree as usize);
            let tmp_term_2 = FiniteFieldPolynomial::monomial((-1, self.field_mod).into(), 1);
            let t = tmp_term_1 + tmp_term_2;
            let tmp = &t % &self;
            let g = self.gcd(&tmp);
            if g != FiniteFieldPolynomial::monomial((1, self.field_mod).into(), 0) {
                return false;
            }
        }
        let coeff_1: FiniteField = (1, self.field_mod).into();
        let deg_1 = self.field_mod.pow(n) as usize;
        let tmp_term_1 = FiniteFieldPolynomial::monomial(coeff_1, deg_1);
        let tmp_term_2 = FiniteFieldPolynomial::monomial((-1, self.field_mod).into(), 1);
        let t = tmp_term_1 + tmp_term_2;
        let (q, r) = &t / &self;
        let x = self * &q;
        let g = &t % &self;
        g.is_zero()
    }

    /// This is based on a wikipedia entry [Primitive Polynomials in Finite Field](https://en.wikipedia.org/wiki/Primitive_polynomial_(field_theory))
    pub fn is_primitive(&self) -> bool {
        let d = self.degree;
        let upper_bound = self.field_mod.pow(d as u32) - 1;
        for n in 1..upper_bound {
            let buf = [
                (n as usize, (1, self.field_mod).into()),
                (0, (-1, self.field_mod).into()),
            ];
            let g = FiniteFieldPolynomial::from(&buf[..]);
            let r = &g % self;
            if r.is_zero() {
                return false;
            }
        }
        let buf = [
            (upper_bound as usize, (1, self.field_mod).into()),
            (0, (-1, self.field_mod).into()),
        ];
        let g = FiniteFieldPolynomial::from(&buf[..]);
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
        for (k,v) in self.coeffs.iter() {
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

    pub fn degree(&self) -> usize {
        self.degree
    }

    pub fn monomial(coefficient: FiniteField, degree: usize) -> FiniteFieldPolynomial {
        let buf = [(degree, coefficient)];
        FiniteFieldPolynomial::from(&buf[..])
    }
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
        // f.write_str(&format!(" mod {:}", self.field_mod))
        Ok(())
    }
}

impl Hash for FiniteFieldPolynomial {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        let mut v: Vec<(usize, FiniteField)> = self.coeffs.clone().into_iter().collect();
        v.sort_by(|x, y| {x.0.cmp(&y.0)});
        for (d, c) in v {
            d.hash(state);
            c.hash(state);
        }
        self.degree.hash(state);
        self.field_mod.hash(state);
    }
}

impl Mul for FiniteFieldPolynomial {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        &self * &rhs
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
        *self *= &rhs;
    }
}

impl MulAssign<&FiniteFieldPolynomial> for FiniteFieldPolynomial {
    fn mul_assign(&mut self, rhs: &FiniteFieldPolynomial) {
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

impl Div<&FiniteFieldPolynomial> for &FiniteFieldPolynomial {
    type Output = (FiniteFieldPolynomial, FiniteFieldPolynomial);

    fn div(self, rhs: &FiniteFieldPolynomial) -> Self::Output {
        if rhs.degree > self.degree {
            return (FiniteFieldPolynomial::zero(self.field_mod), self.clone());
        }
        let mut quotient = FiniteFieldPolynomial::zero(self.field_mod);
        let mut remainder = self.clone();
        let d = rhs.degree;
        let lc = rhs.leading_coeff();
        while remainder.degree >= d && remainder.is_zero() == false {
            let coeff_s = remainder.leading_coeff() * lc.modular_inverse();
            let deg_s = remainder.degree - d;
            let s = FiniteFieldPolynomial::monomial(coeff_s, deg_s);
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

impl Div for FiniteFieldPolynomial {
    type Output = (FiniteFieldPolynomial, FiniteFieldPolynomial);

    fn div(self, rhs: Self) -> Self::Output {
        &self / &rhs
    }
}

impl Rem<&FiniteFieldPolynomial> for &FiniteFieldPolynomial {
    type Output = FiniteFieldPolynomial;

    fn rem(self, rhs: &FiniteFieldPolynomial) -> Self::Output {
        let (_, r) = self / rhs;
        r
    }
}

impl Add<&FiniteFieldPolynomial> for &FiniteFieldPolynomial {
    type Output = FiniteFieldPolynomial;

    fn add(self, rhs: &FiniteFieldPolynomial) -> Self::Output {
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
        &self + &rhs
    }
}

impl AddAssign for FiniteFieldPolynomial {
    fn add_assign(&mut self, rhs: Self) {
        *self += &rhs;
    }
}

impl Add<u32> for &FiniteFieldPolynomial {
    type Output = FiniteFieldPolynomial;

    fn add(self, rhs: u32) -> Self::Output {
        let p = FiniteFieldPolynomial::constant(rhs, self.field_mod);
        self + &p
    }
}

impl AddAssign<u32> for FiniteFieldPolynomial {
    fn add_assign(&mut self, rhs: u32) {
        let c: FiniteField = (rhs, self.field_mod).into();
        let e = self.coeffs.entry(0).or_insert((0, self.field_mod).into());
        *e += c;
        if e.0 == 0 {
            self.coeffs.remove(&0);
        }
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
        &self - &rhs
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
        *self -= &rhs;
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
        let mut poly = FiniteFieldPolynomial {
            coeffs: hm,
            degree: max_degree,
            field_mod: p,
        };
        poly.clean();
        poly
    }
}

impl From<(usize, FiniteField)> for FiniteFieldPolynomial {
    fn from(value: (usize, FiniteField)) -> Self {
        FiniteFieldPolynomial {
            coeffs: HashMap::from([value]),
            degree: value.0,
            field_mod: value.1 .1,
        }
    }
}

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

mod tests {
    use std::collections::{HashMap, HashSet};

    use crate::math::{
        finite_field::FiniteField,
        group_ring_field::Ring,
        polynomial::{get_divisors, get_smallest_divisor},
    };

    use super::{remove_trailing_zeros, FiniteFieldPolynomial, QuotientPoly};

    #[test]
    fn test_polynomial_arithmetic() {
        let coeffs_1: [(usize, FiniteField); 4] = [
            (0, (1, 199).into()),
            (1, (7, 199).into()),
            (2, (0, 199).into()),
            (3, (3, 199).into()),
        ];
        let coeffs_2: [(usize, FiniteField); 4] = [
            (0, (1, 199).into()),
            (1, (3, 199).into()),
            (2, (0, 199).into()),
            (3, (3, 199).into()),
        ];
        let coeffs_3: [(usize, FiniteField); 3] = [
            (0, (2, 199).into()),
            (1, (88, 199).into()),
            (2, (5, 199).into()),
        ];

        let mut p1 = FiniteFieldPolynomial::from(&coeffs_1[..]);
        let mut p2 = FiniteFieldPolynomial::from(&coeffs_2[..]);
        let mut p3 = FiniteFieldPolynomial::from(&coeffs_3[..]);
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
    fn test_quotient_arithmetic() {
        let buff = [
            (0, (4, 199_u32).into()),
            (1, (7, 199_u32).into()),
            (2, (1, 199_u32).into()),
        ];
        let q = FiniteFieldPolynomial::from(&buff[..]);
        let p1 = FiniteFieldPolynomial::zero(199);
        let buff2 = [
            (0, (9, 199_u32).into()),
            (1, (17, 199_u32).into()),
            (4, (1, 199_u32).into()),
        ];
        let p2 = FiniteFieldPolynomial::from(&buff2[..]);
        let q1 = QuotientPoly::zero(199, q.clone());
        let added = &q1 + &p2;
        let multiplied = &added * &q;
        println!("added - {:}", added);
        println!("additive inverse - {:}", added.additive_inv());
        println!("multiplied - {:}", multiplied);
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

        let f = FiniteFieldPolynomial::from(&coeffs[..]);
        let t1 = FiniteFieldPolynomial::monomial((1, 2).into(), 10);

        let t2 = FiniteFieldPolynomial::monomial((1, 2).into(), 3);

        let t3 = FiniteFieldPolynomial::monomial((1, 2).into(), 0);

        let g = t1 + t2 + t3;
        assert!(f.is_irreducible() == false);
        assert!(g.is_irreducible());
    }

    #[test]
    fn test_is_primitive() {
        let buf = [(0, (1, 3).into()), (2, (1, 3).into())];
        let primitive_coeffs = [(2, (1, 3).into()), (1, (2, 3).into()), (0, (2, 3).into())];
        let tester = FiniteFieldPolynomial::from(&buf[..]);
        let primitive_poly = FiniteFieldPolynomial::from(&primitive_coeffs[..]);
        assert!(tester.is_irreducible());
        assert!(tester.is_primitive() == false);
        assert!(primitive_poly.is_irreducible());
        assert!(primitive_poly.is_primitive());
    }

    #[test]
    fn test_quotient_inverse() {
        let p = 3_u32;
        let buf = [(0, (1, p).into()), (2, (1, p).into())];
        let primitive_coeffs = [(2, (1, p).into()), (1, (2, p).into()), (0, (2, p).into())];
        let tester = FiniteFieldPolynomial::from(&buf[..]);
        let primitive_poly = FiniteFieldPolynomial::from(&primitive_coeffs[..]);
        let t = QuotientPoly::from((tester.clone(), primitive_poly.clone()));
        let g = tester.gcd(&primitive_poly);
        println!("tester: {:}", tester);
        println!("primitive: {:}", primitive_poly);
        println!("quotiented: {:}", t);
        println!("gcd: {:}", g);
        let t_inv = t.mul_inv();
        println!("t_inv: {:}", t_inv);
        println!("product: {:}", t_inv * t);
    }

    #[test]
    fn test_constant_addition() {
        let one = FiniteFieldPolynomial::constant(1, 2);
        println!("one: {:}", one);
        println!("one + one = {:}", &one + 1);
    }


}
