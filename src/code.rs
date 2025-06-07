use crate::math::finite_field::FiniteField as FF;
use crate::matrices::ffmatrix::FFMatrix;

pub trait Code {
    /// Takes ownership of a logical message and outputs the physical bits.
    fn encode(&self, message: &Vec<FF>) -> Vec<FF>;

    /// Takes in a logical string and returns the closest codeword,
    /// ties broken arbitrarily.
    fn decode(&self, encrypted: &Vec<FF>) -> Vec<FF>;

    /// Returns `true` if `word` is in the codespace, `false` otherwise.
    fn is_codeword(&self, encrypted: &Vec<FF>) -> bool;

    /// Returns the parity check of the provided message
    fn parity_check(&self, encrypted: &Vec<FF>) -> Vec<FF>;

    fn parity_check_matrix(&self) -> FFMatrix;

    fn parity_check_len(&self) -> usize;

    fn message_len(&self) -> usize;
}

pub fn get_parity_check_matrix_from(generator: &FFMatrix) -> FFMatrix {
    let left_inverse = generator.rref_left_inverse();
    let rref = &left_inverse * generator;
    let id = FFMatrix::id(rref.num_zero_rows(), generator.field_mod);
    let mut parity_check_entries = Vec::new();
    let num_zero_cols = generator.n_rows - id.n_cols;
    for row_ix in 0..id.n_rows {
        for _ in 0..num_zero_cols {
            parity_check_entries.push(0_u32);
        }
        parity_check_entries.append(&mut id.get_row(row_ix));
    }
    &FFMatrix::new(
        parity_check_entries,
        id.n_rows,
        generator.n_rows,
        generator.field_mod,
    ) * &left_inverse
}

pub fn get_generator_from_parity_check(parity_check: &FFMatrix) -> FFMatrix {
    let n = parity_check.n_cols;
    let _k = n - parity_check.n_rows;
    let mut cloned = parity_check.clone();
    if cloned.n_rows > cloned.n_cols {
        cloned.transpose();
    }
    cloned.rref();

    let corner1 = (0, cloned.n_rows);
    let corner2 = (cloned.n_rows - 1, cloned.n_cols - 1);
    let mut p_sub_matrix = cloned.clone_block(corner1, corner2);
    p_sub_matrix.scale(FF::from((-1, cloned.field_mod)));

    let needed_identity = FFMatrix::id(p_sub_matrix.n_cols, cloned.field_mod);
    let mut entries = Vec::new();
    for row_ix in 0..p_sub_matrix.n_rows {
        let mut new_row = p_sub_matrix.get_row(row_ix);
        entries.append(&mut new_row);
    }
    for row_ix in 0..needed_identity.n_rows {
        let mut new_row = needed_identity.get_row(row_ix);
        entries.append(&mut new_row);
    }
    FFMatrix::new(
        entries,
        p_sub_matrix.n_rows + needed_identity.n_rows,
        p_sub_matrix.n_cols,
        parity_check.field_mod,
    )
}
