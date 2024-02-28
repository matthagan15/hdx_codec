use std::{clone, collections::HashMap};

use mhgl::HGraph;
use uuid::Uuid;

use crate::{
    code::Code,
    math::{
        ffmatrix::FFMatrix, finite_field::FiniteField as FF, polynomial::FiniteFieldPolynomial,
    },
};

struct TannerCode<C: Code> {
    data_bits: HashMap<u32, FF>,
    parity_bits: HashMap<u32, FF>,
    local_code: C,
    graph: HGraph,
}

impl<C: Code> TannerCode<C> {}

pub fn get_parity_check_matrix_from(generator_matrix: &FFMatrix) -> FFMatrix {
    // The idea behind here is to clone the generator matrix and put it into
    // reduced row echelon form. We can then use the grassmanian/leftover matrix
    // on the right to compute the parity check matrix
    let n = generator_matrix.n_rows;
    let k = generator_matrix.n_cols;
    let mut cloned = generator_matrix.clone();
    if n > k {
        cloned.transpose();
    }
    cloned.rref();

    let corner1 = (0, k);
    let corner2 = (cloned.n_rows - 1, cloned.n_cols - 1);
    let mut parity_sub_matrix = cloned.clone_block(corner1, corner2);
    parity_sub_matrix.transpose();
    parity_sub_matrix.scale(FF::from((-1, cloned.field_mod)));
    let parity_identity = FFMatrix::id(n - k, generator_matrix.field_mod);
    let mut entries = Vec::new();
    for row_ix in 0..parity_sub_matrix.n_rows {
        let mut row_left = parity_sub_matrix.get_row(row_ix);
        let mut row_right = parity_identity.get_row(row_ix);
        entries.append(&mut row_left);
        entries.append(&mut row_right);
    }
    FFMatrix::new(entries, n - k, n)
}
