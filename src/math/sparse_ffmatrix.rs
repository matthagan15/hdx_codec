use core::num;
use std::{
    collections::HashMap, fmt::{Display, Write}, ops::{Index, Mul}
};

use serde::{Deserialize, Serialize};

use crate::math::finite_field::FiniteField as FF;
use crate::math::finite_field::FiniteFieldExt as FFX;

use super::finite_field::FFRep;

// fn make_orthogonal_to_single(start_vec: Vec<FFX>, ortho_basis: &Vec<FFX>) -> Vec<FFX> {
//     let dot = dot_prod(&start_vec, ortho_basis);
//     Vec::new()
// }

// fn dot_prod(left: &Vec<FFX>, right: &Vec<FFX>) -> FFX {
//     if left.len() != right.len() {
//         panic!("Cannot take dot product of two unequal length vectors.")
//     }
//     // let mut ret = FFX::new(0, left[0].1, left[0].2);
//     // for ix in 0..left.len() {
//     //     ret = ret + &(left[ix] * &right[ix]);
//     // }
//     // ret
//     todo!()
// }

// /// Constructs the Vandermonde matrix V such that V_{i,j} = elements[j]^i.
// pub fn vandermonde(elements: &Vec<FF>, n_rows: usize) -> SparseFFMatrix {
//     let mut entries = Vec::new();
//     for k in 0..n_rows {
//         let mut tmp = elements
//             .clone()
//             .into_iter()
//             .map(|x| x.pow(k as u32))
//             .collect();
//         entries.append(&mut tmp);
//     }
//     SparseFFMatrix::new(entries, n_rows, elements.len())
// }

#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
struct SparseSection(Vec<(usize, FFRep)>);
impl SparseSection {
    
}

/// Used to determine if a memory layout for a matrix (posibly extend to tensor?)
/// is RowMajor or ColMajor, used to determine if the sparse sections
/// retrieved are rows or columns.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum MemoryLayout {
    RowMajor,
    ColMajor,
}

// TODO: Would be nice to convert this and PolyMatrix to be instantiations
// of the same generic type, but I'm not sure that would work due to the
// quotient behavior of PolyMatrix. Either way there should be a trait Matrix
// which captures behavior of both of these and allows for some generic behavior, like RREF and stuff.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SparseFFMatrix {
    ix_to_section: HashMap<u32, SparseSection>,
    pub n_rows: usize,
    pub n_cols: usize,
    pub field_mod: FFRep,
    memory_layout: MemoryLayout,
}

impl SparseFFMatrix {
    /// Creates an all zeros sparse matrix with the specified memory layout.
    /// If your matrix is row-sparse then use `MemoryLayout::RowMajor`, if it is
    /// column sparse use `MemoryLayout::ColMajor`.
    pub fn new(n_rows: usize, n_cols: usize, field_mod: FFRep, layout: MemoryLayout) -> Self {
        Self {
            ix_to_section: HashMap::new(),
            n_rows,
            n_cols,
            field_mod,
            memory_layout: layout,
        }
    }

    pub fn is_square(&self) -> bool {
        self.n_rows == self.n_cols
    }
}