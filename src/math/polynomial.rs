use regex::Regex;
use serde::{Deserialize, Serialize};
use std::{
    cmp::Ordering,
    fmt::Display,
    hash::Hash,
    ops::{Add, AddAssign, Div, Mul, MulAssign, Rem, Sub, SubAssign},
    str::FromStr,
};

use super::finite_field::{FFRep, FiniteField};

pub type PolyDegree = usize;

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
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct FFPolynomial {
    pub coeffs: Vec<(PolyDegree, FFRep)>,
    pub field_mod: u32,
}

impl FFPolynomial {
    pub fn new(
        degrees_and_coefficients: impl AsRef<[(PolyDegree, FFRep)]>,
        field_mod: FFRep,
    ) -> Self {
        if degrees_and_coefficients.as_ref().len() == 0 {
            panic!("Cannot create polynomial without a field_mod");
        }

        let mut coeffs_with_possible_duplicates: Vec<(usize, FFRep)> = degrees_and_coefficients
            .as_ref()
            .iter()
            .filter(|(_d, c)| *c > 0)
            .cloned()
            .collect();
        coeffs_with_possible_duplicates.sort_by_key(|x| x.0);
        let mut coeffs = Vec::new();
        for (d, c) in coeffs_with_possible_duplicates.into_iter() {
            if coeffs.is_empty() {
                coeffs.push((d, c));
            } else {
                let last_coeff = coeffs.last().unwrap();
                if last_coeff.0 == d {
                    let last_coeff = coeffs.pop().unwrap();
                    let new_coeff = (last_coeff.1 + c) % field_mod;
                    if new_coeff > 0 {
                        coeffs.push((d, new_coeff));
                    }
                } else {
                    coeffs.push((d, c));
                }
            }
        }
        FFPolynomial { coeffs, field_mod }
    }

    pub fn monomial(coefficient: FiniteField, degree: PolyDegree) -> FFPolynomial {
        if coefficient.0 == 0 {
            FFPolynomial {
                coeffs: Vec::new(),
                field_mod: coefficient.1,
            }
        } else {
            let buf = [(degree, coefficient)];
            FFPolynomial::from(&buf[..])
        }
    }

    pub fn get_number(&self) -> FFRep {
        let mut num = 0;
        for (deg, coeff) in self.coeffs.iter() {
            num += coeff
                * (self.field_mod.pow(
                    (*deg)
                        .try_into()
                        .expect("Too large a polynomial degree encountered."),
                ));
        }
        num
    }

    pub fn to_number(self) -> FFRep {
        let mut num = 0;
        for (deg, coeff) in self.coeffs.into_iter() {
            num += coeff
                * (self.field_mod.pow(
                    deg.try_into()
                        .expect("Too large a polynomial degree encountered."),
                ));
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
        let mut max_deg = 0;
        while num != 0 {
            let modder = field_mod.pow(pow + 1);
            let quotienter = field_mod.pow(pow);
            let remainder = num % modder;
            let coeff = remainder / quotienter;
            if remainder != 0 {
                max_deg = max_deg.max(pow);
                deg_and_coeffs.push((pow as usize, coeff));
                num -= remainder;
            }
            pow += 1;
        }
        FFPolynomial {
            coeffs: deg_and_coeffs,
            field_mod,
        }
    }

    pub fn get_coeff_of_degree(&self, degree: PolyDegree) -> FiniteField {
        if let Ok(ix) = self.coeffs.binary_search_by_key(&degree, |x| x.0) {
            FiniteField::new(self.coeffs[ix].1, self.field_mod)
        } else {
            FiniteField::new(0, self.field_mod)
        }
    }

    pub fn evaluate(&self, x: &FiniteField) -> FiniteField {
        // This uses a very naive algorithm, it is unclear if Horner's
        // method would be better because this is using a sparse representation.
        let mut out = 0;
        for (d, c) in self.coeffs.iter() {
            out += *c * x.pow(*d as u32).0;
        }
        FiniteField::new(out, self.field_mod)
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
        if scalar.0 == 0 {
            self.coeffs = Vec::new();
        } else {
            for (_, v) in self.coeffs.iter_mut() {
                *v *= scalar.0;
                *v %= self.field_mod;
            }
        }
    }

    /// Creates a new zero polynomial
    pub fn zero(field_mod: u32) -> Self {
        FFPolynomial {
            coeffs: Vec::new(),
            field_mod,
        }
    }

    pub fn constant(c: u32, field_mod: u32) -> Self {
        FFPolynomial {
            coeffs: vec![(0, c)],
            field_mod,
        }
    }

    pub fn leading_coeff(&self) -> FiniteField {
        self.coeffs
            .last()
            .map_or(FiniteField(0, self.field_mod), |coeff| {
                FiniteField::new(coeff.1, self.field_mod)
            })
    }

    pub fn gcd(&self, rhs: &FFPolynomial) -> FFPolynomial {
        if rhs.is_zero() {
            self.clone()
        } else {
            rhs.gcd(&(self % rhs))
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
        while r.degree() >= remainder_degree {
            let (tmp, _) = &r / &new_r;
            let old_r = r.clone();
            r = new_r.clone();
            new_r = old_r - &tmp * &r;

            let old_s = s.clone();
            s = new_s.clone();
            new_s = old_s - &tmp * &s;

            let old_t = t.clone();
            t = new_t.clone();
            new_t = old_t - &tmp * &t;
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

        for power in powers.iter() {
            let tmp_degree = self.field_mod.pow(*power);
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
        let d = self.degree();
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
        self.coeffs.is_empty()
    }

    pub fn is_one(&self) -> bool {
        let mut is_constant_one = false;
        let mut everything_else_zero = true;
        for (k, v) in self.coeffs.iter() {
            if *k == 0 {
                if *v == 1 {
                    is_constant_one = true;
                }
            } else {
                if *v != 0 {
                    everything_else_zero = false;
                }
            }
        }
        is_constant_one && everything_else_zero
    }

    pub fn degree(&self) -> PolyDegree {
        self.coeffs.last().map_or(0, |(deg, _coeff)| *deg)
    }
}

impl PartialOrd for FFPolynomial {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.degree() > other.degree() {
            Some(Ordering::Greater)
        } else if self.degree() < other.degree() {
            Some(Ordering::Less)
        } else {
            // These conditions can still occur even given the above checks,
            // ex: 1*x^0 compared to 0 polynomial.
            if self.coeffs.len() == 0 && other.coeffs.len() > 0 {
                return Some(Ordering::Less);
            } else if self.coeffs.len() > 0 && other.coeffs.len() == 0 {
                return Some(Ordering::Greater);
            }
            let mut lhs_ix = self.coeffs.len() - 1;
            let mut rhs_ix = other.coeffs.len() - 1;
            while lhs_ix > 0 && rhs_ix > 0 {
                if self.coeffs[lhs_ix].0 < other.coeffs[rhs_ix].0 {
                    return Some(Ordering::Less);
                } else if self.coeffs[lhs_ix].0 > other.coeffs[rhs_ix].0 {
                    return Some(Ordering::Greater);
                } else {
                    if self.coeffs[lhs_ix].1 < other.coeffs[rhs_ix].1 {
                        return Some(Ordering::Less);
                    } else if self.coeffs[lhs_ix].1 > other.coeffs[rhs_ix].1 {
                        return Some(Ordering::Greater);
                    } else {
                        lhs_ix -= 1;
                        rhs_ix -= 1;
                    }
                }
            }
            match (lhs_ix == 0, rhs_ix == 0) {
                (true, true) => return self.coeffs[lhs_ix].1.partial_cmp(&other.coeffs[rhs_ix].1),
                (true, false) => {
                    return Some(Ordering::Less);
                }
                (false, true) => return Some(Ordering::Greater),
                (false, false) => {
                    panic!("Loop invariant should guarantee one of lhs_ix or rhs_ix are zero.")
                }
            }
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
            return f.write_str("0");
        }
        let mut display = self
            .coeffs
            .iter()
            .rev()
            .map(|(deg, coeff)| {
                if *deg > 1 {
                    let mut s = if *coeff == 1 {
                        String::new()
                    } else {
                        coeff.to_string()
                    };
                    s.push_str("x^");
                    s.push_str(&deg.to_string()[..]);
                    s
                } else if *deg == 1 {
                    let mut s = coeff.to_string();
                    s.push('x');
                    s
                } else {
                    coeff.to_string()
                }
            })
            .collect::<Vec<String>>()
            .join(" + ");
        display.push_str(" mod ");
        display.push_str(&self.field_mod.to_string());
        f.write_str(&format!("{display}"))
    }
}

impl Hash for FFPolynomial {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        for (d, c) in self.coeffs.iter() {
            (*d).hash(state);
            (*c).hash(state);
        }
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
        if self.coeffs.is_empty() || rhs.coeffs.is_empty() {
            return FFPolynomial::zero(self.field_mod);
        }
        let mut new_polys = Vec::with_capacity(self.coeffs.len());
        for (deg, coeff) in self.coeffs.iter() {
            let new_poly: Vec<(usize, u32)> = rhs
                .coeffs
                .iter()
                .map(|(d, c)| ((d + *deg), (*c * *coeff) % self.field_mod))
                .collect();
            new_polys.push(new_poly);
        }
        let mut out = FFPolynomial {
            coeffs: new_polys.pop().unwrap(),
            field_mod: self.field_mod,
        };
        for coeffs in new_polys.into_iter() {
            let poly = FFPolynomial {
                coeffs,
                field_mod: self.field_mod,
            };
            out += poly;
        }
        out
    }
}

impl MulAssign for FFPolynomial {
    fn mul_assign(&mut self, rhs: Self) {
        *self *= &rhs;
    }
}

impl MulAssign<&FFPolynomial> for FFPolynomial {
    fn mul_assign(&mut self, rhs: &FFPolynomial) {
        *self = &*self * rhs;
    }
}

impl Div<&FFPolynomial> for &FFPolynomial {
    type Output = (FFPolynomial, FFPolynomial);

    fn div(self, rhs: &FFPolynomial) -> Self::Output {
        if rhs.is_zero() {
            panic!("Cannot divide by zero.");
        }
        if rhs.degree() > self.degree() {
            return (FFPolynomial::zero(self.field_mod), self.clone());
        }
        if rhs.degree() == 0 {
            let inv = rhs.leading_coeff().modular_inverse();
            let mut out = self.clone();
            out.scale(&inv);
            return (out, FFPolynomial::zero(self.field_mod));
        }
        let mut quotient = FFPolynomial::zero(self.field_mod);
        let mut remainder = self.clone();
        let d = rhs.degree();
        let lc = rhs.leading_coeff();
        while remainder.degree() >= d && remainder.is_zero() == false {
            let coeff_s = remainder.leading_coeff() * lc.modular_inverse();
            let deg_s = remainder.degree() - d;
            let s = FFPolynomial::monomial(coeff_s, deg_s);
            let remainder_sub = &s * &rhs;
            quotient = quotient + s;
            remainder -= remainder_sub;
        }
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
        let mut out_coeffs = Vec::with_capacity(self.len().max(rhs.len()));
        let mut lhs_ix = 0;
        let mut rhs_ix = 0;
        for _ in 0..(self.coeffs.len() + rhs.coeffs.len()) {
            let l = self.coeffs.get(lhs_ix);
            let r = rhs.coeffs.get(rhs_ix);
            match (l, r) {
                (None, None) => break,
                (None, Some((rhs_deg, rhs_coeff))) => {
                    out_coeffs.push((*rhs_deg, *rhs_coeff));
                    rhs_ix += 1;
                }
                (Some((lhs_deg, lhs_coeff)), None) => {
                    out_coeffs.push((*lhs_deg, *lhs_coeff));
                    lhs_ix += 1;
                }
                (Some((lhs_deg, lhs_coeff)), Some((rhs_deg, rhs_coeff))) => {
                    if lhs_deg < rhs_deg {
                        out_coeffs.push((*lhs_deg, *lhs_coeff));
                        lhs_ix += 1;
                    } else if lhs_deg > rhs_deg {
                        out_coeffs.push((*rhs_deg, *rhs_coeff));
                        rhs_ix += 1;
                    } else {
                        let coeff = (*rhs_coeff + *lhs_coeff) % self.field_mod;
                        if coeff > 0 {
                            out_coeffs.push((*rhs_deg, coeff));
                        }
                        rhs_ix += 1;
                        lhs_ix += 1;
                    }
                }
            }
        }
        FFPolynomial {
            coeffs: out_coeffs,
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

impl AddAssign<&FFPolynomial> for FFPolynomial {
    fn add_assign(&mut self, rhs: &FFPolynomial) {
        let out = &*self + rhs;
        *self = out;
    }
}

impl AddAssign for FFPolynomial {
    fn add_assign(&mut self, rhs: Self) {
        *self += &rhs;
    }
}

impl Add<FFRep> for &FFPolynomial {
    type Output = FFPolynomial;

    fn add(self, rhs: FFRep) -> Self::Output {
        let p = FFPolynomial::constant(rhs, self.field_mod);
        self + &p
    }
}

impl AddAssign<u32> for FFPolynomial {
    fn add_assign(&mut self, rhs: u32) {
        let coeff = rhs % self.field_mod;
        if self.coeffs.is_empty() {
            self.coeffs.push((0, coeff));
        } else {
            if self.coeffs[0].0 == 0 {
                self.coeffs[0].1 += coeff;
                self.coeffs[0].1 %= self.field_mod;
                if self.coeffs[0].1 == 0 {
                    self.coeffs.remove(0);
                }
            } else {
                self.coeffs.insert(0, (0, coeff));
            }
        }
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
        if self.coeffs.is_empty() && rhs.coeffs.is_empty() == false {
            let mut ret = rhs.clone();
            for ix in 0..ret.coeffs.len() {
                let mut new_coeff = ret.coeffs[ix].1;
                new_coeff = (self.field_mod - new_coeff) % self.field_mod;
                ret.coeffs[ix].1 = new_coeff;
            }
            return ret;
        } else if self.coeffs.is_empty() == false && rhs.coeffs.is_empty() {
            return self.clone();
        } else if self.coeffs.is_empty() && rhs.coeffs.is_empty() {
            return self.clone();
        }
        let mut out_coeffs = Vec::with_capacity(self.len().max(rhs.len()));
        let mut lhs_ix = 0;
        let mut rhs_ix = 0;
        for _ in 0..(self.coeffs.len() + rhs.coeffs.len()) {
            let l = self.coeffs.get(lhs_ix);
            let r = rhs.coeffs.get(rhs_ix);
            match (l, r) {
                (None, None) => break,
                (None, Some((rhs_deg, rhs_coeff))) => {
                    let new_coeff = (self.field_mod - *rhs_coeff) % self.field_mod;
                    if new_coeff > 0 {
                        out_coeffs.push((*rhs_deg, new_coeff));
                    }
                    rhs_ix += 1;
                }
                (Some((lhs_deg, lhs_coeff)), None) => {
                    out_coeffs.push((*lhs_deg, *lhs_coeff));
                    lhs_ix += 1;
                }
                (Some((lhs_deg, lhs_coeff)), Some((rhs_deg, rhs_coeff))) => {
                    if lhs_deg < rhs_deg {
                        out_coeffs.push((*lhs_deg, *lhs_coeff));
                        lhs_ix += 1;
                    } else if lhs_deg > rhs_deg {
                        let new_coeff = (self.field_mod - *rhs_coeff) % self.field_mod;
                        if new_coeff > 0 {
                            out_coeffs.push((*rhs_deg, new_coeff));
                        }
                        rhs_ix += 1;
                    } else {
                        let out_coeff = if *rhs_coeff > *lhs_coeff {
                            let mut tmp = *lhs_coeff + self.field_mod;
                            tmp -= *rhs_coeff;
                            tmp % self.field_mod
                        } else {
                            (*lhs_coeff - *rhs_coeff) % self.field_mod
                        };
                        if out_coeff > 0 {
                            out_coeffs.push((*rhs_deg, out_coeff))
                        };
                        rhs_ix += 1;
                        lhs_ix += 1;
                    }
                }
            }
        }
        FFPolynomial {
            coeffs: out_coeffs,
            field_mod: self.field_mod,
        }
    }
}

impl SubAssign for FFPolynomial {
    fn sub_assign(&mut self, rhs: Self) {
        *self -= &rhs;
    }
}

impl SubAssign<&FFPolynomial> for FFPolynomial {
    fn sub_assign(&mut self, rhs: &FFPolynomial) {
        *self = &*self - rhs;
    }
}

impl From<&[(PolyDegree, FiniteField)]> for FFPolynomial {
    fn from(value: &[(PolyDegree, FiniteField)]) -> Self {
        if value.len() == 0 {
            panic!("Cannot create polynomial without a field_mod");
        }

        let field_mod = value[0].1 .1;
        FFPolynomial::new(
            value
                .iter()
                .map(|(d, c)| (*d, c.0))
                .collect::<Vec<(PolyDegree, FFRep)>>(),
            field_mod,
        )
    }
}

impl From<(PolyDegree, FiniteField)> for FFPolynomial {
    fn from(value: (PolyDegree, FiniteField)) -> Self {
        FFPolynomial {
            coeffs: vec![(value.0, value.1 .0)],
            field_mod: value.1 .1,
        }
    }
}

impl FromStr for FFPolynomial {
    type Err = ();

    /// String must be of the form: "a_d*x^d + a_{d-1}x^{d-1} + ... + a_0*x^0 % prime". terms are split at the plus and the percentage sign indicates the finite field being used.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let field_mod_split = match (s.contains('%'), s.contains("mod")) {
            (true, true) => return Err(()),
            (true, false) => s.split('%').collect::<Vec<&str>>(),
            (false, true) => s.split("mod").collect::<Vec<&str>>(),
            (false, false) => return Err(()),
        };
        if field_mod_split.len() != 2 {
            return Err(());
        }
        let field_mod = field_mod_split[1].trim().parse::<u32>();
        if field_mod.is_err() {
            return Err(());
        }
        let field_mod = field_mod.unwrap();
        let mut poly = field_mod_split[0].split_whitespace().collect::<String>();
        if poly.is_empty() {
            return Ok(FFPolynomial::zero(field_mod));
        }
        if poly.chars().nth(0) != Some('+') && poly.chars().nth(0) != Some('-') {
            poly.insert(0, '+');
        }
        let mut coeffs = Vec::new();
        let re =
            Regex::new(r"(?<sign>[\+,\-])(?<coeff>\d*)(x\^(?<exp>\d+)|(?<lone_x>x?))").unwrap();
        for poly_term in re.captures_iter(&poly[..]) {
            let negate: bool = poly_term.name("sign").unwrap().as_str() == "-";
            let coeff: u32 = match poly_term.name("coeff").unwrap().as_str() {
                "" => {
                    if negate {
                        field_mod - 1
                    } else {
                        1
                    }
                }
                c => {
                    if let Ok(coeff) = c.parse::<u32>() {
                        if negate {
                            (field_mod - coeff) % field_mod
                        } else {
                            coeff % field_mod
                        }
                    } else {
                        return Err(());
                    }
                }
            };
            let degree: usize = match (
                poly_term.name("exp").is_some(),
                poly_term.name("lone_x").is_some(),
            ) {
                (true, true) => return Err(()),
                (true, false) => {
                    if let Ok(degree) = poly_term.name("exp").unwrap().as_str().parse::<usize>() {
                        degree
                    } else {
                        return Err(());
                    }
                }
                (false, true) => {
                    if poly_term.name("lone_x").unwrap().as_str() == "x" {
                        1
                    } else {
                        if poly_term.name("coeff").unwrap().as_str().is_empty() {
                            return Err(());
                        }
                        0
                    }
                }
                (false, false) => return Err(()),
            };
            coeffs.push((degree, coeff));
        }
        Ok(FFPolynomial::new(coeffs, field_mod))
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use regex::Regex;

    use crate::math::{
        finite_field::FiniteField,
        polynomial::{get_divisors, PolyDegree},
    };

    use super::FFPolynomial;

    #[test]
    fn test_polynomial_arithmetic() {
        let p = 199;
        let coeffs_1: [(PolyDegree, FiniteField); 3] = [
            (0, (1, 199).into()),
            (1, (7, 199).into()),
            (3, (3, 199).into()),
        ];
        let coeffs_2: [(PolyDegree, FiniteField); 3] = [
            (0, (1, 199).into()),
            (1, (3, 199).into()),
            (3, (3, 199).into()),
        ];
        let coeffs_3: [(PolyDegree, FiniteField); 3] = [
            (0, (2, 199).into()),
            (1, (88, 199).into()),
            (2, (5, 199).into()),
        ];
        let coeffs_5: [(PolyDegree, FiniteField); 3] = [
            (0, (197, 199).into()),
            (1, (111, 199).into()),
            (2, (194, 199).into()),
        ];
        let coeffs_6: Vec<(usize, FiniteField)> = vec![
            (0, (2, 199).into()),
            (1, (88 + 2 * 7, 199).into()),
            (2, (1 * 5 + 7 * 88, 199).into()),
            (3, (2 * 3 + 7 * 5, 199).into()),
            (4, (3 * 88, 199).into()),
            (5, (3 * 5, 199).into()),
        ];

        let coeffs_7: Vec<(usize, FiniteField)> = vec![(0, (1, p).into())];
        let coeffs_8: Vec<(usize, FiniteField)> = vec![(70, (2, p).into())];
        let coeffs_9: Vec<(usize, FiniteField)> = vec![(140, (4, p).into())];

        let p1 = FFPolynomial::from(&coeffs_1[..]);
        let p2 = FFPolynomial::from(&coeffs_2[..]);
        let p3 = FFPolynomial::from(&coeffs_3[..]);
        let p4 = &(&p1 * &p2) + &p3;
        let p5 = FFPolynomial::from(&coeffs_5[..]);
        let p6 = &p1 * &p3;
        let p7 = FFPolynomial::from(&coeffs_7[..]);
        let p8 = FFPolynomial::from(&coeffs_8[..]);
        let p9 = FFPolynomial::from(&coeffs_9[..]);
        assert_eq!(p8, &p7 * &p8);
        assert_eq!(p9, &p8 * &p8);

        let zero = &p3 + &p5;
        assert!(zero.coeffs.is_empty());
        assert_eq!(p6, FFPolynomial::from(&coeffs_6[..]));
        let (q, r) = &p4 / &p2;
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
        let p = 5;
        let coeffs = vec![
            (0, (1, p).into()),
            // (1, (0, p).into()),
            // (2, (0, p).into()),
            // (3, (0, p).into()),
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
        assert!(poly.is_irreducible());
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
        let known_interpolation = FFPolynomial::from(
            &vec![
                (0, (80, p).into()),
                (1, (163, p).into()),
                (2, (103, p).into()),
                (3, (63, p).into()),
            ][..],
        );
        assert_eq!(known_interpolation, l);
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
        let tester = FFPolynomial::from(&buf[..]);
        let s = serde_json::to_string(&tester).expect("serialized");
        let f: FFPolynomial = serde_json::from_str(&s).expect("could not deserialize.");
        assert_eq!(f, tester);
    }

    #[test]
    fn test_poly_from_str() {
        let s1 = "1 + 2x +3x^2 + 4x^3 mod 3";
        let f1 = FFPolynomial::from_str(s1);
        let v1 = vec![(0, 1), (1, 2), (3, 1)];
        let g1 = FFPolynomial::new(v1, 3);
        assert_eq!(Ok(g1), f1);

        let s2 = "+1 + 2x     \n  +3x^2 + 4x^3 %       3";
        let f2 = FFPolynomial::from_str(s2);
        let v2 = vec![(0, 1), (1, 2), (3, 1)];
        let g2 = FFPolynomial::new(v2, 3);
        assert_eq!(Ok(g2), f2);

        let s3 = "1 + 2 + 3 + 4 mod 5";
        let f3 = FFPolynomial::from_str(s3);
        assert_eq!(Ok(FFPolynomial::zero(5)), f3);

        let s4 = "1 + x";
        let f4 = FFPolynomial::from_str(s4);
        assert!(f4.is_err());

        let s5 = "++ x^2 mod 3";
        let f5 = FFPolynomial::from_str(s5);
        assert!(f5.is_err());

        let s6 = "2 + x + x^2 - x + 2 x^2  - 2 % 3";
        let f6 = FFPolynomial::from_str(s6);
        assert_eq!(Ok(FFPolynomial::zero(3)), f6);

        let s7 = "3 ++ x % 5";
        let f7 = FFPolynomial::from_str(s7);
        assert!(f7.is_err());

        let s8 = "3 -- x % 5";
        let f8 = FFPolynomial::from_str(s8);
        assert!(f8.is_err());
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
