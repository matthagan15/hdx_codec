use std::collections::{HashMap, HashSet};

use mhgl::Graph;
use ndarray::prelude::*;
use rand::prelude::*;

use crate::math::pauli::PauliString;

struct ReedSolomon {}

#[derive(Debug)]
struct LDPC {
    factor_graph: Graph<u32>,
    // This is supposed to be a precomputed encoder
    // encoder: HashMap<u32, PauliString>,
    bit_nodes: HashSet<u32>,
    check_nodes: HashSet<u32>,
    weight: usize,
    k: usize,
    n: usize,
}

impl LDPC {
    pub fn random(n: usize, k: usize, weight: usize) -> Self {
        let m = n - k;
        let mut g = Graph::new();
        let left_nodes: Vec<u32> = (0..n as u32).collect();
        let right_nodes: Vec<u32> = (n as u32..n as u32 + m as u32).collect();
        g.add_nodes(left_nodes.clone());
        g.add_nodes(right_nodes.clone());
        // Sample m times of w different selections
        let mut seen_connections = HashSet::new();
        for ix in 0..right_nodes.len() {
            let mut connections = left_nodes.clone();
            shuffle(&mut connections);
            connections.resize(weight, 0);
            connections.sort();
            if seen_connections.contains(&connections) {
                continue;
            }
            for node in connections.iter() {
                g.add_edge(right_nodes[ix], *node);
            }
            seen_connections.insert(connections);
        }
        LDPC {
            factor_graph: g,
            bit_nodes: left_nodes.into_iter().collect(),
            check_nodes: right_nodes.into_iter().collect(),
            weight,
            k,
            n,
        }
    }
}

struct ExpanderCode {}

/// Shuffles a vector uniformly. Probably not at all efficient,
/// samples a permutation uniformly at random.
fn shuffle<T: Copy>(v: &mut Vec<T>) {
    let n = v.len();
    let mut rng = thread_rng();
    for ix in 0..n - 1 {
        for jx in ix..n {
            if rng.gen_bool(0.5) {
                let tmp = v[ix];
                v[ix] = v[jx];
                v[jx] = tmp;
            }
        }
    }
}

mod tests {
    use super::LDPC;

    #[test]
    fn randomized_construction() {
        let ldpc = LDPC::random(10, 3, 4);
        dbg!(ldpc);
    }
}
