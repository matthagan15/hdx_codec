//! A library for error correcting codes over High Dimensional Expanders (HDX).
//! The code is defined over a simplicial complex $X$  
//! Constructs parity check matrices extracted from boundary operators.
//!
//! Intended to be a relatively self-contained implementation of the
// ! error correcting codes from (New Codes on High Dimensional Expanders)[https://arxiv.org/abs/2308.15563]
pub mod code;
pub mod first_node;
pub mod math;
pub mod matrices;
pub mod quantum;
pub mod rank_estimator_sparse;
pub mod rate_and_dist_estimator;
pub mod reed_solomon;
pub mod tanner_code;

use std::{collections::HashSet, path::PathBuf, time::Instant};

use math::finite_field::FFRep;
pub use math::lps;
use matrices::{ffmatrix::FFMatrix, sparse_ffmatrix::MemoryLayout};
use serde::Deserialize;

#[inline]
pub fn factorial(n: usize) -> usize {
    (2..=n).fold(1, |a, i| a * i)
}

pub fn binomial(n: usize, k: usize) -> usize {
    let k = k.min(n - k);
    let top = ((n - k + 1)..=n).fold(1, |a, i| a * i);
    let bot = factorial(k);
    top / bot
}

pub fn powerset(input: Vec<u32>) -> Vec<Vec<u32>> {
    if input.len() == 0 {
        vec![vec![]]
    } else {
        let mut ret = HashSet::new();
        for ix in 0..input.len() {
            let mut new_input = input.clone();
            new_input.remove(ix);
            new_input.sort();
            for new_vec in powerset(new_input.clone()) {
                ret.insert(new_vec);
            }
        }
        let mut final_input = input.clone();
        final_input.sort();
        ret.insert(final_input);
        ret.into_iter().collect()
    }
}

pub fn row_reduction_benchmark(file: PathBuf) {
    #[derive(Deserialize)]
    struct SampleData {
        input: Vec<FFRep>,
        output: Vec<FFRep>,
        rank: usize,
    }
    #[derive(Deserialize)]
    struct BenchData {
        num_samples: usize,
        n_rows: usize,
        n_cols: usize,
        finite_field: FFRep,
        #[allow(dead_code)]
        memory_layout: MemoryLayout,
        time_per_matrix: f64,
        data: Vec<SampleData>,
    }
    if let Ok(s) = std::fs::read_to_string(file) {
        if let Ok(bench_data) = serde_json::from_str::<BenchData>(&s[..]) {
            let mut tot_time = 0.0;
            let mut num_correct = 0;
            let mut correct_ranks = 0;
            for sample in bench_data.data {
                let mut input_matrix = FFMatrix::new(
                    sample.input,
                    bench_data.n_rows,
                    bench_data.n_cols,
                    bench_data.finite_field,
                );
                let start = Instant::now();
                input_matrix.rref();
                tot_time += start.elapsed().as_secs_f64();
                let computed_rank = input_matrix.rank();

                let python_output = FFMatrix::new(
                    sample.output,
                    bench_data.n_rows,
                    bench_data.n_cols,
                    bench_data.finite_field,
                );
                if computed_rank == sample.rank {
                    correct_ranks += 1;
                }

                if python_output == input_matrix {
                    num_correct += 1;
                }
            }
            let time_per_matrix = tot_time / bench_data.num_samples as f64;
            let ratio = bench_data.time_per_matrix / time_per_matrix;
            println!("python rref / rust: {:}", ratio);
            println!(
                "num_correct / num_samples: {:}",
                num_correct as f64 / bench_data.num_samples as f64
            );
            println!("correct ranks: {:}", correct_ranks as f64);
        } else {
            println!("couldnt deserialize.");
        }
    } else {
        println!("couldn't open file.");
    }
}

#[cfg(test)]
mod tests {
    use crate::{binomial, factorial, powerset};

    #[test]
    fn factorial_and_binomial() {
        assert_eq!(factorial(4), 24);
        assert_eq!(factorial(10), 3_628_800);
        assert_eq!(binomial(10, 4), 210);
    }

    #[test]
    fn powerset_small() {
        let v = vec![0, 1, 2, 3, 4];
        let mut count = 0;
        for p in powerset(v.clone()) {
            count += 1;
            let mut s = String::new();

            for n in p {
                s.push_str(&n.to_string()[..]);
                s.push_str(", ");
            }
            s.pop();
            let x = s.pop();
            if x.is_none() {
                println!("empty set");
            } else {
                println!("{s}");
            }
        }
        assert_eq!(count, 2_usize.pow(v.len() as u32));
    }
}
