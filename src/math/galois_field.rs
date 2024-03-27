use std::collections::HashMap;

use super::{finite_field::FFRep, polynomial::FiniteFieldPolynomial};

pub struct GaloisField {
    lookup: HashMap<(FFRep, FFRep), FFRep>,
}

impl GaloisField {
    pub fn new(quotient: FiniteFieldPolynomial) -> Self {
        todo!()
    }
}