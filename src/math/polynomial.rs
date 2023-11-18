use std::collections::HashMap;

use super::group_ring_field::Ring;

/// Polynomial in single indeterminate
pub struct Polynomial<R> where R: Ring {
    coeffs: HashMap<i32, R>,
}

impl<R: Ring> Polynomial<R> {
    
}
