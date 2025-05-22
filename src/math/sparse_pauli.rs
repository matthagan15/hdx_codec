use crate::matrices::sparse_vec::SparseVector;

#[allow(dead_code)]
enum Clifford {
    Hadamard(usize),
    Phase(usize),
    CNOT(usize, usize),
}

#[allow(dead_code)]
struct SparsePauli {
    x_ixs: Vec<usize>,
    z_ixs: Vec<usize>,
}

#[allow(dead_code)]
struct CliffordTableau {
    phases: SparseVector,
}

impl SparsePauli {
    #[allow(dead_code)]
    pub fn new() -> Self {
        SparsePauli {
            x_ixs: Vec::new(),
            z_ixs: Vec::new(),
        }
    }

    #[allow(dead_code)]
    /// flips the bit associated with `ix`
    pub fn add_x(&mut self, ix: usize) {
        match self.x_ixs.binary_search(&ix) {
            Ok(ix_storage) => {
                self.x_ixs.remove(ix_storage);
            }
            Err(ix_storage) => {
                self.x_ixs.insert(ix_storage, ix);
            }
        }
    }
    #[allow(dead_code)]
    pub fn add_xs(&mut self, ixs: impl AsRef<[usize]>) {
        for ix in ixs.as_ref() {
            self.add_x(*ix);
        }
    }

    #[allow(dead_code)]
    pub fn add_z(&mut self, ix: usize) {
        match self.z_ixs.binary_search(&ix) {
            Ok(ix_storage) => {
                self.z_ixs.remove(ix_storage);
            }
            Err(ix_storage) => {
                self.z_ixs.insert(ix_storage, ix);
            }
        }
    }

    #[allow(dead_code)]
    pub fn add_zs(&mut self, ixs: impl AsRef<[usize]>) {
        for ix in ixs.as_ref() {
            self.add_z(*ix);
        }
    }

    #[allow(dead_code)]
    pub fn multiply(&mut self, other: &SparsePauli) {
        let mut ret = Vec::with_capacity(self.x_ixs.len().min(other.x_ixs.len()));
        let mut self_ix = 0;
        let mut other_ix = 0;
        for _ix in 0..(self.x_ixs.len() + other.x_ixs.len()) {
            if self_ix == self.x_ixs.len() || other_ix == other.x_ixs.len() {
                break;
            }
            if self.x_ixs[self_ix] < other.x_ixs[other_ix] {
                self_ix += 1;
            } else if self.x_ixs[self_ix] > other.x_ixs[other_ix] {
                other_ix += 1;
            } else {
                ret.push(self.x_ixs[self_ix]);
                self_ix += 1;
                other_ix += 1;
            }
        }
        self.x_ixs = ret;
    }
}
