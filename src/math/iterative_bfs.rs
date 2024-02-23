use std::{collections::{HashMap, HashSet, VecDeque}, path::PathBuf};

use mhgl::HGraph;

use crate::{hdx_code::HDXCodeConfig, lps::compute_generators};
use crate::math::coset_complex::*;

use super::{finite_field::FiniteField, polymatrix::PolyMatrix, polynomial::FiniteFieldPolynomial};

struct Coset {
    type_ix: usize,
    set: Vec<PolyMatrix>,
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

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
struct Triangle {
    type_zero_coset: Vec<PolyMatrix>,
    type_one_coset: Vec<PolyMatrix>,
    type_two_coset: Vec<PolyMatrix>,
    distance_from_origin: usize,
}

impl Triangle {
    pub fn new(g: &PolyMatrix, subgroups: &CosetGenerators) -> Self {
        let type_zero_coset = subgroups.type_to_generators.get(&0).expect("no zero").clone();
        let type_one_coset = subgroups.type_to_generators.get(&1).expect("no zero").clone();
        let type_two_coset = subgroups.type_to_generators.get(&2).expect("no zero").clone();
        Self {
            type_zero_coset: type_zero_coset.into_iter().map(|k_0| g * &k_0).collect(),
            type_one_coset: type_one_coset.into_iter().map(|k_1| g * &k_1).collect(),
            type_two_coset: type_two_coset.into_iter().map(|k_2| g * &k_2).collect(),
            distance_from_origin: 0
        }
    }

    pub fn neighbors(&self, cg: &CosetGenerators) -> Vec<Self> {
        let mut ret = Vec::new();
        for (type_ix, subgroup) in cg.type_to_generators.iter() {
            for new_generator in subgroup.iter() {
                let mut new_type_one_coset = self.type_one_coset
                    .clone()
                    .into_iter()
                    .map(|g| {
                        new_generator * &g
                    })
                    .collect();
                let mut new_type_two_coset = self.type_two_coset
                    .clone()
                    .into_iter()
                    .map(|g| {
                        new_generator * &g
                    })
                    .collect();
                let mut new_type_zero_coset = self.type_zero_coset
                    .clone()
                    .into_iter()
                    .map(|g| {
                        new_generator * &g
                    })
                    .collect();
                ret.push(Triangle { type_zero_coset: new_type_zero_coset, type_one_coset: new_type_one_coset, type_two_coset: new_type_two_coset, distance_from_origin: self.distance_from_origin + 1 });
            }
        } 
        ret
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
fn triangle_based_bfs(generators: &CosetGenerators, verbose: bool) {
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
    let starting_triangle = Triangle::new(&e, generators);

    // TODO: Currently a matrix is being stored twice while it is in
    // the frontier as we also put it in visited. Instead just keep track
    // of an extra bit if the matrix is visited or not.
    // Also, do not store completed in RAM. come up with some way of storing
    // them on disk in the meanwhile.
    let mut completed = HashSet::new();
    let mut frontier = VecDeque::from([starting_triangle.clone()]);
    let mut visited = HashSet::from([starting_triangle.clone()]);
    let gens = generators.type_to_generators.clone();
   
    let mut counter = 0;
    let mut last_flushed_distance = 0;
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
                // let new = g * &x;
                // if visited.contains(&new) == false {
                //     visited.insert(new.clone());
                //     frontier.push_back(new);
                // }
            }
        }
        completed.insert(x);
    }
}

pub fn compute_hgraph(hdx_conf: HDXCodeConfig) -> HGraph {
    let subgroups = CosetGenerators::new(hdx_conf.dim, &hdx_conf.quotient_poly);
    let mut hg = HGraph::new();

    hg
}