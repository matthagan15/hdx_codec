enum Clifford {
    Hadamard(usize),
    Phase(usize),
    CNOT(usize, usize),
}

struct SparsePauli {
    ixs: Vec<usize>,
}

impl SparsePauli {
    pub fn new() -> Self {
        SparsePauli { ixs: Vec::new() }
    }

    /// flips the bit associated with `ix`
    pub fn add(&mut self, ix: usize) {
        match self.ixs.binary_search(&ix) {
            Ok(ix_storage) => {
                self.ixs.remove(ix_storage);
            }
            Err(ix_storage) => {
                self.ixs.insert(ix_storage, ix);
            }
        }
    }
    pub fn add_a_lot(&mut self, ixs: impl AsRef<[usize]>) {
        for ix in ixs.as_ref() {
            self.add(*ix);
        }
    }

    pub fn multiply(&mut self, other: &SparsePauli) {
        let mut ret = Vec::with_capacity(self.ixs.len().min(other.ixs.len()));
        let mut self_ix = 0;
        let mut other_ix = 0;
        for _ix in 0..(self.ixs.len() + other.ixs.len()) {
            if self_ix == self.ixs.len() || other_ix == other.ixs.len() {
                break;
            }
            if self.ixs[self_ix] < other.ixs[other_ix] {
                self_ix += 1;
            } else if self.ixs[self_ix] > other.ixs[other_ix] {
                other_ix += 1;
            } else {
                ret.push(self.ixs[self_ix]);
                self_ix += 1;
                other_ix += 1;
            }
        }
        self.ixs = ret;
    }
}
