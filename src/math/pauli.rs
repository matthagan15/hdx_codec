use std::fmt::Display;

use bitvec::vec::BitVec;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum Phase {
    pr,
    pi,
    nr,
    ni,
}

impl Phase {
    pub fn mul(&self, rhs: &Phase) -> Phase {
        match (self, rhs) {
            (Phase::pr, Phase::pr) => Phase::pr,
            (Phase::pr, Phase::pi) => Phase::pi,
            (Phase::pr, Phase::nr) => Phase::nr,
            (Phase::pr, Phase::ni) => Phase::ni,
            (Phase::pi, Phase::pr) => Phase::pi,
            (Phase::pi, Phase::pi) => Phase::nr,
            (Phase::pi, Phase::nr) => Phase::ni,
            (Phase::pi, Phase::ni) => Phase::pr,
            (Phase::nr, Phase::pr) => Phase::nr,
            (Phase::nr, Phase::pi) => Phase::ni,
            (Phase::nr, Phase::nr) => Phase::pr,
            (Phase::nr, Phase::ni) => Phase::pi,
            (Phase::ni, Phase::pr) => Phase::ni,
            (Phase::ni, Phase::pi) => Phase::pr,
            (Phase::ni, Phase::nr) => Phase::pi,
            (Phase::ni, Phase::ni) => Phase::nr,
        }
    }
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum Pauli {
    I,
    X,
    Y,
    Z,
}

impl Pauli {
    pub fn mul(&self, rhs: &Pauli) -> (Phase, Pauli) {
        match (self, rhs) {
            (Pauli::I, Pauli::I) => (Phase::pr, Pauli::I),
            (Pauli::I, Pauli::X) => (Phase::pr, Pauli::X),
            (Pauli::I, Pauli::Y) => (Phase::pr, Pauli::Y),
            (Pauli::I, Pauli::Z) => (Phase::pr, Pauli::Z),
            (Pauli::X, Pauli::I) => (Phase::pr, Pauli::X),
            (Pauli::X, Pauli::X) => (Phase::pr, Pauli::I),
            (Pauli::X, Pauli::Y) => (Phase::pi, Pauli::Z),
            (Pauli::X, Pauli::Z) => (Phase::ni, Pauli::Y),
            (Pauli::Y, Pauli::I) => (Phase::pr, Pauli::Y),
            (Pauli::Y, Pauli::X) => (Phase::ni, Pauli::Z),
            (Pauli::Y, Pauli::Y) => (Phase::pr, Pauli::I),
            (Pauli::Y, Pauli::Z) => (Phase::pi, Pauli::X),
            (Pauli::Z, Pauli::I) => (Phase::pr, Pauli::Z),
            (Pauli::Z, Pauli::X) => (Phase::pi, Pauli::Y),
            (Pauli::Z, Pauli::Y) => (Phase::ni, Pauli::X),
            (Pauli::Z, Pauli::Z) => (Phase::pr, Pauli::I),
        }
    }
}

#[derive(Debug, Clone)]
pub struct ExperimentalPauliString {
    phase: Phase,
    pub x_string: BitVec,
    pub z_string: BitVec,
}

#[derive(Debug, Clone)]
pub struct PauliString {
    pub phase: Phase,
    pub string: Vec<Pauli>, // TODO: how to store pauli string? currently using a string, should use a bitvec
}

impl PauliString {
    // TODO: Should really use regex for this?
    // TODO: Add phase?
    pub fn from(input: String) -> Self {
        let paulis = input
            .into_bytes()
            .into_iter()
            .map(|x| match x {
                b'I' => Pauli::I,
                b'X' => Pauli::X,
                b'Y' => Pauli::Y,
                b'Z' => Pauli::Z,
                _ => panic!("Encountered non-pauli character in input string"),
            })
            .collect();
        PauliString {
            phase: Phase::pr,
            string: paulis,
        }
    }

    pub fn index(&self, index: usize) -> Pauli {
        self.string
            .get(index)
            .expect("Bad PauliString index")
            .clone()
    }

    pub fn len(&self) -> usize {
        self.string.len()
    }

    pub fn new() -> Self {
        PauliString {
            phase: Phase::pr,
            string: Vec::new(),
        }
    }

    pub fn push(&mut self, p: Pauli) {
        self.string.push(p);
    }

    pub fn set_phase(&mut self, phase: Phase) {
        self.phase = phase;
    }

    pub fn mul(&self, rhs: &PauliString) -> PauliString {
        let mut tot_phase = Phase::pr;
        if self.len() != rhs.len() {
            PauliString::new()
        } else {
            let mut ps = PauliString::new();
            for ix in 0..self.len() {
                let (tmp_phase, pauli) = self.index(ix).mul(&rhs.index(ix));
                tot_phase = tot_phase.mul(&tmp_phase);
                ps.push(pauli);
            }
            ps.set_phase(tot_phase);
            ps
        }
    }

    pub fn commutes(&self, rhs: &PauliString) -> bool {
        let self_on_left = self.mul(rhs);
        let self_on_right = rhs.mul(self);
        for ix in 0..self.len() {
            if self_on_left.index(ix) != self_on_right.index(ix) {
                return false;
            }
        }
        if self_on_left.phase != self_on_right.phase {
            return false;
        }
        true
    }

    pub fn anti_commutes(&self, rhs: &PauliString) -> bool {
        let self_on_left = self.mul(rhs);
        let self_on_right = rhs.mul(self);
        for ix in 0..self.len() {
            if self_on_left.index(ix) != self_on_right.index(ix) {
                return false;
            }
        }
        if self_on_left.phase != self_on_right.phase.mul(&Phase::nr) {
            return false;
        }
        true
    }
}

impl Display for PauliString {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let _ = f.write_str(match self.phase {
            Phase::pr => "+",
            Phase::pi => "+i",
            Phase::nr => "-",
            Phase::ni => "-i",
        });
        for pauli in self.string.iter() {
            let _ = f.write_str(match pauli {
                Pauli::I => "I",
                Pauli::X => "X",
                Pauli::Y => "Y",
                Pauli::Z => "Z",
            });
        }
        Ok(())
    }
}
