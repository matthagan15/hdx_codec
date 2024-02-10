use serde::{Deserialize, Serialize};

use crate::math::finite_field::FiniteField;

#[derive(Debug, Clone, Hash, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct FFMatrix {
    pub entries: Vec<FiniteField>,
    pub n_rows: usize,
    pub n_cols: usize,
    pub field_mod: u32,
}

impl FFMatrix {
    pub fn new(entries: Vec<FiniteField>,
        n_rows: usize,
        n_cols: usize,
        field_mod: u32) -> Self {
            Self {
                entries,
                n_rows,
                n_cols,
                field_mod,
            }
        }
}