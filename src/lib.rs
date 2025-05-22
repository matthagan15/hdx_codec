//! A library for error correcting codes over High Dimensional Expanders (HDX).
//! The code is defined over a simplicial complex $X$  
//! Constructs parity check matrices extracted from boundary operators.
//!
//! Intended to be a relatively self-contained implementation of the
// ! error correcting codes from (New Codes on High Dimensional Expanders)[https://arxiv.org/abs/2308.15563]
pub mod code;
pub mod math;
pub mod matrices;
pub mod quantum;
pub mod rank_estimator_sparse;
pub mod reed_solomon;
pub mod tanner_code;

use std::{path::PathBuf, time::Instant};

use math::finite_field::{FFRep, FiniteField};
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
            let mut printed = false;
            for sample in bench_data.data {
                let mut input_matrix = FFMatrix::new(
                    sample
                        .input
                        .into_iter()
                        .map(|x| FiniteField::new(x, bench_data.finite_field))
                        .collect(),
                    bench_data.n_rows,
                    bench_data.n_cols,
                );
                let start = Instant::now();
                input_matrix.rref();
                tot_time += start.elapsed().as_secs_f64();
                let computed_rank = input_matrix.rank();

                let python_output = FFMatrix::new(
                    sample
                        .output
                        .into_iter()
                        .map(|x| FiniteField::new(x, bench_data.finite_field))
                        .collect(),
                    bench_data.n_rows,
                    bench_data.n_cols,
                );
                if computed_rank == sample.rank {
                    correct_ranks += 1;
                }
                if printed == false {
                    printed = true;
                    println!("Python output:\n{:}", python_output);
                    println!("Rust output:\n{:}", input_matrix);
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
    use crate::{binomial, factorial};

    #[test]
    fn factorial_and_binomial() {
        assert_eq!(factorial(4), 24);
        assert_eq!(factorial(10), 3_628_800);
        assert_eq!(binomial(10, 4), 210);
    }
}
