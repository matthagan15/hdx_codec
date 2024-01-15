use std::collections::{HashSet, HashMap};

use mhgl::HGraph;

use super::{polynomial::{FiniteFieldPolynomial, QuotientPoly}, matrix::PolyMatrix, finite_field::FiniteField};

fn compute_generators(dim: usize, quotient: FiniteFieldPolynomial) -> CosetGenerators {
    let e = PolyMatrix::id(dim, quotient.clone());
    let p = quotient.field_mod;
    let mut generators: HashMap<usize, Vec<PolyMatrix>> = HashMap::new();
    for j in 0..dim {
        let mut v = Vec::new();
        for a in 0..p {
            for b in 0..p {
                if (a, b) == (0, 0) {
                    continue;
                }
                let coeffs: Vec<(usize, FiniteField)> = vec![(0, (b, p).into()), (1, (a, p).into())];
                let poly = FiniteFieldPolynomial::from(&coeffs[..]);
                let mut mat = e.clone();
                let entry = mat.get_mut(j, j + 1);
                *entry = QuotientPoly::from((poly, quotient.clone()));
                v.push(mat);
            }
        }
        generators.insert(j, v);
    }
    generators
}

type CosetGenerators = HashMap<usize, Vec<PolyMatrix>>;

fn compute_group(generators: CosetGenerators) -> HashSet<PolyMatrix> {
    let mut ret = HashSet::new();
    if generators.is_empty() {
        return ret;
    }

    ret
}

/// dim is the dimension of the resulting coset complex. So `dim = 2`` would give a graph and `dim = 3` would give a triangle complex and so on.
fn generate_complex(quotient: FiniteFieldPolynomial, dim: usize) -> HGraph {
    HGraph::new()
}