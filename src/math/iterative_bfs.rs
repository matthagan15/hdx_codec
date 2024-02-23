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

impl Triangle {
    pub fn new(g_1: PolyMatrix, g_2: PolyMatrix, g_3: PolyMatrix) -> Self {
        if g_1.is_square() | g_2.is_square() | g_3.is_square() == false {
            panic!("Need all square matrices")
        }
        if (g_1.n_rows, g_1.n_cols) != (g_3.n_rows, g_3.n_cols) {
            panic!("Shapes do not match.")
        }
        if (g_1.n_rows, g_1.n_cols) != (g_2.n_rows, g_2.n_cols) {
            panic!("Shapes do not match.")
        }
        if (g_2.n_rows, g_2.n_cols) != (g_3.n_rows, g_3.n_cols) {
            panic!("Shapes do not match.")
        }
        let e = PolyMatrix::id(g_1.n_rows, g_1.quotient.clone());
        let cg = CosetGenerators::new();
        let c0 = compute_coset(&e, 0, coset_gens: &);
    }
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

/// Currently comptes the entire group using Breadth-First-Search
/// starting at the identity matrix over the generators provided.
fn compute_group(generators: &CosetGenerators, verbose: bool) -> HashSet<PolyMatrix> {
    if verbose {
        println!("Computing group with the following parameters:");
        println!("quotient: {:}", generators.quotient);
        println!("dim: {:}", generators.dim);
        for (t, set) in generators.type_to_generators.iter() {
            println!("{:}", "*".repeat(50));
            println!("Type: {:}", t);
            for s in set.iter() {
                println!("{:}", s);
            }
        }
    }
    let e = PolyMatrix::id(generators.dim, generators.quotient.clone());
    let starting_triangle = Triangle::new(e, e, e);

    // TODO: Currently a matrix is being stored twice while it is in
    // the frontier as we also put it in visited. Instead just keep track
    // of an extra bit if the matrix is visited or not.
    // Also, do not store completed in RAM. come up with some way of storing
    // them on disk in the meanwhile.
    let mut completed = HashSet::new();
    let mut frontier = VecDeque::from([e.clone()]);
    let mut visited = HashSet::from([e.clone()]);
    let gens = generators.type_to_generators.clone();
    
    if gens.is_empty() {
        completed.insert(e);
        return completed;
    }
   

    let mut counter = 0;
    while frontier.len() > 0 {
        counter += 1;
        if (counter % 100) == 0 {
            println!("{:}", ".".repeat(50));
            println!(
                "frontier length: {:} , {:.4}% of visitied.",
                frontier.len(),
                frontier.len() as f64 / visited.len() as f64
            );
            println!("completed length: {:}", completed.len());
            println!("visited length: {:}", visited.len());
        }
        let x = frontier.pop_front().expect("no frontier?");
        for (j, gen_list) in gens.iter() {
            for g in gen_list {
                let new = g * &x;
                if visited.contains(&new) == false {
                    visited.insert(new.clone());
                    frontier.push_back(new);
                }
            }
        }
        completed.insert(x);
    }
    completed
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