use crate::math::{ffmatrix::FFMatrix, finite_field::FiniteField as FF};

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

pub fn get_parity_check_matrix_from(generator_matrix: &FFMatrix) -> FFMatrix {
    // The idea behind here is to clone the generator matrix and put it into
    // reduced row echelon form. We can then use the grassmanian/leftover matrix
    // on the right to compute the parity check matrix
    let mut cloned = generator_matrix.clone();
    if cloned.n_rows > cloned.n_cols {
        cloned.transpose();
    }
    cloned.rref();
    if cloned.n_rows < cloned.n_cols {
        cloned.transpose();
    }
    if cloned.n_rows == cloned.n_cols {
        panic!("Why do I have a square generator matrix?.");
    }
    let p_mat_rows = cloned.n_rows - cloned.n_cols;
    let p_mat_cols = cloned.n_cols;
    let corner1 = (cloned.n_cols, 0);
    let corner2 = (cloned.n_cols + p_mat_rows - 1, cloned.n_cols - 1);
    let mut parity_sub_matrix = cloned.clone_block(corner1, corner2);
    parity_sub_matrix.transpose();
    parity_sub_matrix.scale(FF::from((-1, cloned.field_mod)));

    let parity_identity = FFMatrix::id(cloned.n_cols, generator_matrix.field_mod);
    let mut entries = Vec::new();
    for row_ix in 0..parity_sub_matrix.n_rows {
        let mut row_left = parity_sub_matrix.get_row(row_ix);
        let mut row_right = parity_identity.get_row(row_ix);
        entries.append(&mut row_left);
        entries.append(&mut row_right);
    }
    FFMatrix::new(
        entries,
        parity_sub_matrix.n_rows,
        parity_sub_matrix.n_cols + parity_identity.n_cols,
    )
}

pub fn get_generator_from_parity_check(parity_check: &FFMatrix) -> FFMatrix {
    let n = parity_check.n_cols;
    let k = n - parity_check.n_rows;
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
    )
}
