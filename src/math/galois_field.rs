use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{Read, Write},
    path::Path,
};

use fxhash::FxHashMap;
use serde::{Deserialize, Serialize};

use super::{
    finite_field::{FFRep, FiniteField},
    polynomial::FFPolynomial,
};

fn generate_all_polys(field_mod: u32, max_degree: usize) -> HashMap<usize, HashSet<FFPolynomial>> {
    if max_degree == 0 {
        let mut ret = HashSet::new();
        for a in 0..field_mod {
            ret.insert(FFPolynomial::constant(a, field_mod));
        }
        HashMap::from([(0, ret)])
    } else {
        let smalls = generate_all_polys(field_mod, max_degree - 1);
        let mut ret = HashMap::new();
        let monomials: Vec<FFPolynomial> = (1..field_mod)
            .into_iter()
            .map(|x| FFPolynomial::monomial(FiniteField::new(x, field_mod), max_degree as u32))
            .collect();
        for (deg, polys) in smalls.into_iter() {
            for poly in polys.iter() {
                for mono in monomials.iter() {
                    let with_mono: &mut HashSet<FFPolynomial> = ret.entry(max_degree).or_default();
                    with_mono.insert(mono + poly);

                    let without_mono = ret.entry(poly.degree() as usize).or_default();
                    without_mono.insert(poly.clone());
                }
            }
        }
        ret
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GaloisField {
    lookup: FxHashMap<(FFRep, FFRep), FFRep>,
    pub field_mod: FFRep,
    pub quotient: FFPolynomial,
}

impl GaloisField {
    pub fn new(quotient: FFPolynomial) -> Self {
        // Generate all polys
        // Compute all pairs
        let all_polys = generate_all_polys(quotient.field_mod, (quotient.degree() - 1) as usize);
        let mut all_polys_vec = Vec::new();
        for (_, set) in all_polys.into_iter() {
            for poly in set.into_iter() {
                all_polys_vec.push(poly);
            }
        }
        let mut lookup: FxHashMap<(FFRep, FFRep), FFRep> = FxHashMap::default();
        for p1 in all_polys_vec.iter() {
            for p2 in all_polys_vec.iter() {
                if p1 <= p2 {
                    let out = &(p1 * p2) % &quotient;
                    let n1 = p1.clone().get_number();
                    let n2 = p2.clone().get_number();
                    let res = out.get_number();
                    lookup.insert((n1, n2), res);
                }
            }
        }
        GaloisField {
            lookup,
            field_mod: quotient.field_mod,
            quotient,
        }
    }

    pub fn mul(&self, lhs: FFRep, rhs: FFRep) -> FFRep {
        let query = (lhs.min(rhs), lhs.max(rhs));
        *self.lookup.get(&query).expect("Not precomputed.")
    }

    pub fn add(&self, lhs: &FFRep, rhs: &FFRep) -> FFRep {
        // TODO: implement this in FFRep arithmetic instead of converting, adding, then converting back. This would be a minor speedup I think but not asymptotically advantageous.
        let lhs_poly = FFPolynomial::from_number(*lhs, self.field_mod);
        let rhs_poly = FFPolynomial::from_number(*rhs, self.field_mod);
        let sum = lhs_poly + rhs_poly;
        sum.to_number()
    }

    pub fn mul_poly(&self, lhs: &FFPolynomial, rhs: &FFPolynomial) -> FFPolynomial {
        let n1 = lhs.get_number();
        let n2 = rhs.get_number();
        let out = self.mul(n1, n2);
        FFPolynomial::from_number(out, lhs.field_mod)
    }

    // pub fn to_disk(&self, filename: &Path) {
    //     let s = serde_json::to_string(&self).expect("Could not serialize GaloisField");
    //     let mut file = File::create(filename).expect("Coudl not open file for GaloisField");
    //     file.write_all(s.as_bytes()).expect("Could not write to file for GaloisField.");
    // }

    // pub fn from_disk(filename: &Path) -> Self {
    //     let mut file = File::open(filename).expect("Could not open file for read for GaloisField");
    //     let mut buf = String::new();
    //     file.read_to_string(&mut buf).expect("could not read file for GaloisField.");
    //     serde_json::from_str(&buf).expect("Could not deserialize GaloisField")
    // }
}

mod tests {
    use std::str::FromStr;

    use crate::math::{finite_field::FFRep, polynomial::FFPolynomial};

    use super::{generate_all_polys, GaloisField};

    #[test]
    fn test_polynomial_generation() {
        let p: FFRep = 7;
        let out = generate_all_polys(p, 3);
        for (d, set) in out {
            println!("{:}", "*".repeat(50));
            println!("deg: {:}", d);
            for poly in set {
                println!("{:}", poly);
                let num = poly.get_number();
                let computed = FFPolynomial::from_number(num, p);
                assert_eq!(poly, computed)
            }
        }
    }

    #[test]
    fn test_multiplication_table() {
        let p = 3;
        let quotient = FFPolynomial::from_str("1*x^2 + 2*x^1 + 2*x^0 % 3").expect("parse?");
        let all_polys = generate_all_polys(p, 1);
        let mut all_polys_flat = Vec::new();
        for (deg, set) in all_polys.into_iter() {
            for p in set.into_iter() {
                all_polys_flat.push(p);
            }
        }
        let galois = GaloisField::new(quotient.clone());
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
