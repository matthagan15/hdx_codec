use std::{
    fmt::Display,
    ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign},
};

use super::{
    group_ring_field::Ring,
    polynomial::{FFPolynomial, PolyDegree},
};

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct QuotientPoly {
    pub poly: FFPolynomial,
    /// the q in a = q * b + r
    pub quotient: FFPolynomial,
    pub field_mod: u32,
}

impl QuotientPoly {
    /// creates a new zero polynomial with entries in F_{field_mod}. quotient is the polynomial you are modding by.
    pub fn zero(field_mod: u32, quotient: FFPolynomial) -> Self {
        let poly = FFPolynomial::zero(field_mod);
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

    pub fn monomial(coeff: u32, degree: PolyDegree, quotient: FFPolynomial) -> Self {
        let p = FFPolynomial::monomial((coeff, quotient.field_mod).into(), degree);
        let r = &p % &quotient;
        QuotientPoly {
            poly: r,
            quotient,
            field_mod: p.field_mod,
        }
    }

    pub fn constant(coeff: u32, quotient: FFPolynomial) -> QuotientPoly {
        let buf = [(0, (coeff, quotient.field_mod).into())];
        let n = quotient.field_mod;
        QuotientPoly {
            poly: FFPolynomial::from(&buf[..]),
            quotient,
            field_mod: n,
        }
    }

    /// Returns the b such that a * b = 1 mod q. Panics if no inverse exists
    pub fn mul_inv(&self) -> Self {
        let g = self.poly.gcd(&self.quotient);
        if g.degree() > 0 {
            panic!("Cannot find inverse as gcd is not constant.")
        }
        let g_inv = g.leading_coeff().modular_inverse();
        let p = self.field_mod;
        let mut r = self.poly.clone();
        let mut new_r = self.quotient.clone();
        let mut s = FFPolynomial::constant(1, self.field_mod);
        let mut new_s = FFPolynomial::zero(p);
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

impl Add<&FFPolynomial> for &QuotientPoly {
    type Output = QuotientPoly;

    fn add(self, rhs: &FFPolynomial) -> Self::Output {
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

impl Mul<&FFPolynomial> for &QuotientPoly {
    type Output = QuotientPoly;

    fn mul(self, rhs: &FFPolynomial) -> Self::Output {
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
            poly: FFPolynomial::zero(self.field_mod),
            quotient: self.quotient.clone(),
            field_mod: self.field_mod,
        }
    }

    fn one(&self) -> Self {
        let buf = [(0, (1, self.field_mod).into())];
        QuotientPoly {
            poly: FFPolynomial::from(&buf[..]),
            quotient: self.quotient.clone(),
            field_mod: self.field_mod,
        }
    }

    fn additive_inv(&self) -> Self {
        QuotientPoly::zero(self.field_mod, self.quotient.clone()) - self.clone()
    }
}

impl From<(&FFPolynomial, &FFPolynomial)> for QuotientPoly {
    fn from(value: (&FFPolynomial, &FFPolynomial)) -> Self {
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

impl From<(FFPolynomial, FFPolynomial)> for QuotientPoly {
    fn from(value: (FFPolynomial, FFPolynomial)) -> Self {
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

mod tests {
    use crate::math::{
        group_ring_field::Ring, polynomial::FFPolynomial, quotient_polynomial::QuotientPoly,
    };

    #[test]
    fn test_quotient_arithmetic() {
        let buff = [
            (0, (4, 199_u32).into()),
            (1, (7, 199_u32).into()),
            (2, (1, 199_u32).into()),
        ];
        let q = FFPolynomial::from(&buff[..]);
        let p1 = FFPolynomial::zero(199);
        let buff2 = [
            (0, (9, 199_u32).into()),
            (1, (17, 199_u32).into()),
            (4, (1, 199_u32).into()),
        ];
        let p2 = FFPolynomial::from(&buff2[..]);
        let q1 = QuotientPoly::zero(199, q.clone());
        let added = &q1 + &p2;
        let multiplied = &added * &q;
        println!("added - {:}", added);
        println!("additive inverse - {:}", added.additive_inv());
        println!("multiplied - {:}", multiplied);
    }

    #[test]
    fn test_quotient_inverse() {
        let p = 3_u32;
        let buf = [(0, (1, p).into()), (2, (1, p).into())];
        let primitive_coeffs = [(2, (1, p).into()), (1, (2, p).into()), (0, (2, p).into())];
        let tester = FFPolynomial::from(&buf[..]);
        let primitive_poly = FFPolynomial::from(&primitive_coeffs[..]);
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
}
