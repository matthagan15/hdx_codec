//! A library for error correcting codes over High Dimensional Expanders (HDX).
//! The code is defined over a simplicial complex $X$  
//! Constructs parity check matrices extracted from boundary operators.
//!
//! Intended to be a relatively self-contained implementation of the
// ! error correcting codes from (New Codes on High Dimensional Expanders)[https://arxiv.org/abs/2308.15563]
pub mod code;
pub mod hdx_code;
pub type HDXCode = hdx_code::NewHDXCode;
pub mod math;
pub mod matrices;
pub mod quantum;
pub mod rank_estimator_sparse;
pub mod reed_solomon;
pub mod tanner_code;

pub use math::lps;

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
