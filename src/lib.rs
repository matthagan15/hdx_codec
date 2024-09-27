//! A library for codes over High Dimensional Expanders (HDX)
//!
//! Intended to be a relatively self-contained implementation of the
//! error correcting codes from (New Codes on High Dimensional Expanders)[https://arxiv.org/abs/2308.15563] and the various implementations of c^3
//! locally testable codes, which can be found scattered across the following papers
//! - (Asymptotically Good Quantum and Locally Testable Classical LDPC Codes)[https://arxiv.org/abs/2111.03654]
//! - (Good Quantum LDPC Codes with Linear Time Decoders)[https://arxiv.org/abs/2206.07750]
//! - (Quantum Tanner codes)[https://arxiv.org/abs/2202.13641]
pub mod code;
pub mod hdx_code;
pub type HDXCode = hdx_code::NewHDXCode;
pub mod math;
pub mod matrices;
pub mod quantum;
pub mod rank_estimator_sparse;
pub mod rank_estimator_sparse_sparse;
pub mod reed_solomon;
pub mod tanner_code;

use math::polynomial::FFPolynomial;
use matrices::mat_trait::RankMatrix;

use std::{io::Write, path::PathBuf};

pub use math::lps;
use mhgl::{HGraph, HyperGraph};

pub struct Conf {
    quotient_poly: FFPolynomial,
    dim: usize,
    reed_solomon_degree: usize,
    cache_file: Option<PathBuf>,
    hgraph_file: PathBuf,
    num_threads: usize,
}

pub fn factorial(n: usize) -> usize {
    (2..=n).fold(1, |a, i| a * i)
}

pub fn coset_complex_to_disk(hgraph: &HGraph<u16, ()>, filename: PathBuf) {
    if let Ok(mut file) = std::fs::File::create(filename) {
        write!(&mut file, "nodes\n").expect("Can't write?");
        for node in hgraph.nodes() {
            let type_ix = *hgraph.get_node(&node).unwrap();
            write!(&mut file, "{:},{:}\n", node, type_ix).expect("Can't write");
        }
        write!(&mut file, "edges\n").expect("Cannot write");
        for edge in hgraph.edges() {
            let nodes = hgraph.query_edge(&edge).unwrap();
            let mut s = String::new();
            for ix in 0..nodes.len() - 1 {
                s.push_str(&format!("{:}", nodes[ix])[..]);
                s.push(',');
            }
            s.push_str(&format!("{:}", nodes[nodes.len() - 1])[..]);
            write!(&mut file, "{:}", s).expect("Cannot write");
        }
    }
}

pub fn binomial(n: usize, k: usize) -> usize {
    let k = k.min(n - k);
    let top = ((n - k + 1)..=n).fold(1, |a, i| a * i);
    let bot = factorial(k);
    top / bot
}

mod tests {
    use super::{binomial, factorial};
    #[test]
    fn factorial_and_binomial() {
        assert_eq!(factorial(4), 24);
        assert_eq!(factorial(10), 3_628_800);
        assert_eq!(binomial(10, 4), 210);
    }
}
