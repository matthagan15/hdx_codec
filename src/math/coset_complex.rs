
use std::{collections::{HashSet, HashMap, VecDeque}, time};

use mhgl::HGraph;

use super::{polynomial::{FiniteFieldPolynomial, QuotientPoly}, matrix::PolyMatrix, finite_field::FiniteField};

fn compute_generators(dim: usize, quotient: FiniteFieldPolynomial) -> CosetGenerators {
    let mut q = quotient.clone();
    println!("irreducible? {:}", q.is_irreducible());
    println!("primitive? {:}", q.is_primitive());

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
                let alpha = QuotientPoly::from((poly, quotient.clone()));
                *entry = alpha;
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
    // println!("generators:");
    // for (k, v) in gens.iter() {
    //     println!("k = {:}", k);
    //     for m in v.iter() {
    //         println!("{:}", m);
    //     }
    //     println!("{:}", "*".repeat(70));
    // }
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
        // println!("frontier:\n");
        // for v in frontier.iter() {
        //     println!("{:}", v);
        // }
        // println!("completed:");
        // for c in completed.iter() {
        //     println!("{:}", c);
        // }
        // println!("completed vertex: {:}", x);
        if (completed.len() + 1) % 1000 == 0 {
            let percent = (visited.len() as f64) / num_matrices_upper_bound as f64;
            let duration = start_time.elapsed().as_secs_f64();
            let tot_time = duration / percent;
            let remaining_time = tot_time - duration;
            println!("upper bound on time remaining: {:} min", remaining_time / 60.);
            println!("upper bound on time remaining: {:} hours", remaining_time / (60. * 60.));
            println!("completed this many: {:}", completed.len() + 1);
        }
        completed.insert(x);
        // std::thread::sleep(time::Duration::from_millis(1000));
        // break;
    }
    completed
}

/// dim is the dimension of the resulting coset complex. So `dim = 2`` would give a graph and `dim = 3` would give a triangle complex and so on.
fn generate_complex(quotient: FiniteFieldPolynomial, dim: usize) -> HGraph {
    HGraph::new()
}

mod tests {
    use crate::math::polynomial::{QuotientPoly, FiniteFieldPolynomial};

    use super::{compute_generators, compute_group};

    use deepsize::DeepSizeOf;

    #[test]
    fn test_compute_group() {
        let p = 2_u32;
        let primitive_coeffs = [
            (2, (1, p).into()),
            (1, (2, p).into()),
            (0, (2, p).into())];
        let primitive_poly = FiniteFieldPolynomial::from(&primitive_coeffs[..]);
        let dim = 3;
        let gens = compute_generators(dim, primitive_poly);
        let g = compute_group(&gens);
        let size = g.deep_size_of();
        println!("size of g: {:}", size);
        // for m in g.into_iter() {
        //     println!("{:}", m);
        // }
    }
}