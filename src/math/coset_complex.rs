
use std::{collections::{HashSet, HashMap, VecDeque}, time};

use mhgl::HGraph;

use super::{polynomial::{FiniteFieldPolynomial, QuotientPoly}, matrix::PolyMatrix, finite_field::FiniteField};

fn compute_generators(dim: usize, quotient: FiniteFieldPolynomial) -> CosetGenerators {
    let q = &quotient;

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
                *entry = &poly % q;
                v.push(mat);
            }
        }
        generators.insert(j, v);
    }
    CosetGenerators {
        type_to_generators: generators,
        dim,
        quotient,
    }
}

struct CosetGenerators {
    type_to_generators: HashMap<usize, Vec<PolyMatrix>>,
    dim: usize,
    quotient: FiniteFieldPolynomial,
}
struct HDXCode {
    
}

/// Currently comptes the entire group using Breadth-First-Search 
/// starting at the identity matrix over the generators provided.
fn compute_group(generators: &CosetGenerators) -> HashSet<PolyMatrix> {

    let num_matrices_upper_bound = (generators.quotient.field_mod.pow(generators.quotient.degree() as u32)).pow((generators.dim * generators.dim) as u32);
    let start_time = std::time::Instant::now();

    let mut completed = HashSet::new();
    let gens = generators.type_to_generators.clone();
    let e = PolyMatrix::id(generators.dim, generators.quotient.clone());
    if gens.is_empty() {
        completed.insert(e);
        return completed;
    }
    let mut frontier = VecDeque::from([e.clone()]);
    let mut visited = HashSet::from([e.clone()]);

    while frontier.len() > 0 {
        let x = frontier.pop_front().expect("no frontier?");
        // println!("Visiting: {:}", x);
        for (j, gen_list) in gens.iter() {
            for g in gen_list {
                let new = g * &x;
                // println!("visiting: {:}", new);
                if visited.contains(&new) == false {
                    visited.insert(new.clone());
                    frontier.push_back(new);
                }
            }
        }
        // if (completed.len() + 1) % 1000 == 0 {
        //     let percent = (visited.len() as f64) / num_matrices_upper_bound as f64;
        //     let duration = start_time.elapsed().as_secs_f64();
        //     let tot_time = duration / percent;
        //     let remaining_time = tot_time - duration;
        //     println!("upper bound on time remaining: {:} min", remaining_time / 60.);
        //     println!("upper bound on time remaining: {:} hours", remaining_time / (60. * 60.));
        //     println!("completed this many: {:}", completed.len() + 1);
        // }
        completed.insert(x);
    }
    completed
}

fn compute_coset(start: &PolyMatrix, gens: &Vec<PolyMatrix>) -> HashSet<PolyMatrix> {
    let mut completed = HashSet::new();
    let mut frontier = VecDeque::from([start.clone()]);
    let mut visited = HashSet::from([start.clone()]);
    while frontier.len() > 0 {
        let x = frontier.pop_front().expect("no frontier?");
        for g in gens {
            let new = &x * g;
            // println!("visiting: {:}", new);
            if visited.contains(&new) == false {
                visited.insert(new.clone());
                frontier.push_back(new);
            }
        }
        completed.insert(x);
    }
    completed
}

fn compute_vertices(gens: &CosetGenerators, group: HashSet<PolyMatrix>) {
    // let mut v = HashSet::new();
    // for (t, gen_vec) in gens.type_to_generators.iter() {
    //     let mut coset = HashSet::new();
    //     for gen_mat in gen_vec.iter() {

    //     }
    // }
}

/// dim is the dimension of the resulting coset complex. So `dim = 2`` would give a graph and `dim = 3` would give a triangle complex and so on.
fn generate_complex(quotient: FiniteFieldPolynomial, dim: usize) -> HGraph {
    HGraph::new()
}

mod tests {
    use std::collections::HashSet;

    use crate::math::{polynomial::{QuotientPoly, FiniteFieldPolynomial}, matrix::PolyMatrix, coset_complex::compute_coset};

    use super::{compute_generators, compute_group, CosetGenerators};

    use deepsize::DeepSizeOf;

    fn simplest_group() -> (CosetGenerators, HashSet<PolyMatrix>) {
        let p = 2_u32;
        let primitive_coeffs = [
            (2, (1, p).into()),
            (1, (2, p).into()),
            (0, (2, p).into())];
        let primitive_poly = FiniteFieldPolynomial::from(&primitive_coeffs[..]);
        let dim = 3;
        let gens = compute_generators(dim, primitive_poly);
        let g = compute_group(&gens);
        (gens, g)
    }

    #[test]
    fn test_compute_group() {
        let (gens, g) = simplest_group();
        let size = g.deep_size_of();
        println!("size of g: {:}", size);
        // for m in g.into_iter() {
        //     println!("{:}", m);
        // }
    }

    #[test]
    fn test_coset() {
        let (gens, g) = simplest_group();
        println!("len of g: {:}", g.len());
        for x in g.iter() {
            println!("x: {:}", x);
            let coset = compute_coset(&x, gens.type_to_generators.get(&0).unwrap());
            println!("coset.");
            for h in coset.iter() {
                println!("{:}", h);
            }
            break;
        }
    }
}