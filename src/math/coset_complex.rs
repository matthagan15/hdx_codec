use std::collections::HashSet;

use mhgl::HGraph;

use super::{polynomial::FiniteFieldPolynomial, matrix::PolyMatrix};

fn generate_group(dim: usize, quotient: FiniteFieldPolynomial) -> HashSet<PolyMatrix> {
    let e = PolyMatrix::id(dim, quotient);
    HashSet::new()
}

/// dim is the dimension of the resulting coset complex. So `dim = 2`` would give a graph and `dim = 3` would give a triangle complex and so on.
fn generate_complex(quotient: FiniteFieldPolynomial, dim: usize) -> HGraph {
    HGraph::new()
}