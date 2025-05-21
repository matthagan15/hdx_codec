use fxhash::FxHashMap;
use serde::{Deserialize, Serialize};

use super::{finite_field::FFRep, polynomial::FFPolynomial};

/// A lookup table for computing additions and multiplications of elements in a finite
/// field. Represents the elements of the finite field
/// as integers mod p^d and converts to/from polynomials over $F_p$ if an element
/// is not contained in the multiplication lookup. a hashmap with the `Fx` hash function
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GaloisField {
    pub lookup: FxHashMap<(FFRep, FFRep), FFRep>,
    pub field_mod: FFRep,
    pub quotient: FFPolynomial,
}

impl GaloisField {
    pub fn new(quotient: FFPolynomial) -> Self {
        let lookup: FxHashMap<(FFRep, FFRep), FFRep> = FxHashMap::default();
        GaloisField {
            lookup,
            field_mod: quotient.field_mod,
            quotient,
        }
    }

    /// Computes the multiplication of two integer encoded finite fields and returns the
    /// result as an integer.
    ///
    /// `&mut self` is needed due to the fact that if the answer is not in the lookup table
    /// then it will be added.
    pub fn mul(&mut self, lhs: FFRep, rhs: FFRep) -> FFRep {
        let query = (lhs.min(rhs), lhs.max(rhs));
        *self.lookup.entry(query).or_insert({
            let p1 = FFPolynomial::from_number(lhs.min(rhs), self.field_mod);
            let p2 = FFPolynomial::from_number(lhs.max(rhs), self.field_mod);
            let out = &(&p1 * &p2) % &self.quotient;
            let res = out.get_number();
            res
        })
    }

    pub fn add(&self, lhs: &FFRep, rhs: &FFRep) -> FFRep {
        let lhs_poly = FFPolynomial::from_number(*lhs, self.field_mod);
        let rhs_poly = FFPolynomial::from_number(*rhs, self.field_mod);
        let sum = lhs_poly + rhs_poly;
        sum.to_number()
    }

    pub fn mul_poly(&mut self, lhs: &FFPolynomial, rhs: &FFPolynomial) -> FFPolynomial {
        let n1 = lhs.get_number();
        let n2 = rhs.get_number();
        let out = self.mul(n1, n2);
        FFPolynomial::from_number(out, lhs.field_mod)
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use crate::math::polynomial::FFPolynomial;

    use super::GaloisField;

    #[test]
    fn test_multiplication_table() {
        let p = 3;
        let quotient = FFPolynomial::from_str("1*x^2 + 2*x^1 + 2*x^0 % 3").expect("parse?");
        let mut galois = GaloisField::new(quotient.clone());
        let p1 = FFPolynomial::monomial((1, p).into(), 1);
        let p2 = FFPolynomial::constant(2, p);
        let p3 = FFPolynomial::from_str("1*x^1 + 2 * x^0 % 3").expect("parse?");
        let n1 = p1.get_number();
        let n2 = p2.get_number();
        println!("{:}", galois.mul_poly(&p1, &p2));
        dbg!(galois.mul(n1, n2));
        println!("{:}", galois.mul_poly(&p1, &p1));
        dbg!(galois.mul(n1, n1));
        println!("{:}", galois.mul_poly(&p1, &p3));
    }
}
