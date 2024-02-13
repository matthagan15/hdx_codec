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
pub type HDXCode = hdx_code::HDXCode;
pub mod math;
pub mod reed_solomon;
mod tanner_code;
pub mod quantum;

pub use math::lps;