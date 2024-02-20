use std::{collections::{HashMap, HashSet}, path::PathBuf};

use mhgl::HGraph;

use crate::{hdx_code::HDXCodeConfig, lps::compute_generators};
use crate::math::coset_complex::*;

use super::{finite_field::FiniteField, polymatrix::PolyMatrix, polynomial::FiniteFieldPolynomial};

struct Coset {
    type_ix: usize,
    set: Vec<PolyMatrix>,
}

struct CosetGenerators {
    pub type_to_generators: HashMap<usize, HashSet<PolyMatrix>>,
    dim: usize,
    quotient: FiniteFieldPolynomial,
}

pub struct CosetComplex {
    file_base: String,
    dim: usize,
    quotient: FiniteFieldPolynomial,
    group: Option<HashSet<PolyMatrix>>,
    subgroups: Option<CosetGenerators>,
    hgraph: HGraph,
    node_to_coset: Option<HashMap<u32, Coset>>,
}

struct DFNode {
    matrix: PolyMatrix,
    visited: bool,
    distance: usize,
}
struct DFSurfer {
    base_dir: PathBuf,
    id_to_matrix: HashMap<bool, bool>
}

struct Triangle {
    type_zero_coset: Vec<PolyMatrix>,
    type_one_coset: Vec<PolyMatrix>,
    type_two_coset: Vec<PolyMatrix>,
    distance_from_origin: usize,
}

fn h_type_subgroup(type_ix: usize, quotient: FiniteFieldPolynomial) -> Vec<PolyMatrix> {
    let mut ret = Vec::new();
    let dim = 3;
    let p = quotient.field_mod;
    let id = PolyMatrix::id(dim, quotient.clone());
    let mut row_ix = type_ix as i32 - 1;
    while row_ix <= 0 {
        row_ix += dim as i32;
    }
    row_ix %= dim as i32;
    let mut col_ix = type_ix;
    for a in 0..p {
        let mut tmp = id.clone();
        let e = tmp.get_mut(row_ix as usize, col_ix);
        *e = FiniteFieldPolynomial::monomial(FiniteField::new(a, p), 1);
        ret.push(tmp);
    }
    ret
}

impl Triangle {
    fn compute_neighbors(&self, subgroups: CosetGenerators) -> Vec<Triangle> {
        // first type zero neighbors:
        let v1 = &self.type_one_coset;
        let v2 = &self.type_two_coset;

        Vec::new()
    }
}

pub fn compute_hgraph(hdx_conf: HDXCodeConfig) -> HGraph {
    let subgroups = compute_subgroups(hdx_conf.dim, hdx_conf.quotient_poly.clone());
    let mut hg = HGraph::new();

    hg
}