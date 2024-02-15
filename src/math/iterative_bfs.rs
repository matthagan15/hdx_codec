use std::{collections::HashMap, path::PathBuf};

use super::polymatrix::PolyMatrix;

struct BFNode {
    matrix: PolyMatrix,
    visited: bool,
    distance: usize,
}


struct BFSurfer {
    base_dir: PathBuf,
    id_to_matrix: HashMap<bool, bool>
}