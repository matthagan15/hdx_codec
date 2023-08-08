use std::collections::HashMap;

use mhgl::PGraph;
use rand::prelude::*;

struct CssCode {}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
enum Phase {
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
            (Phase::pi, Phase::pr) => Phase::pr,
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
enum Pauli {
    I,
    X,
    Y,
    Z,
}

#[derive(Debug, Clone)]
struct PauliString {
    phase: Phase,
    string: Vec<Pauli>,
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
        self.string.get(index).expect("Bad PauliString index")
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
impl Pauli {
    fn mul(&self, rhs: &Pauli) -> (Phase, Pauli) {
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
            (Pauli::Z, Pauli::I) => (Phase::pr, Pauli::I),
            (Pauli::Z, Pauli::X) => (Phase::pi, Pauli::Y),
            (Pauli::Z, Pauli::Y) => (Phase::ni, Pauli::X),
            (Pauli::Z, Pauli::Z) => (Phase::pr, Pauli::I),
        }
    }
}

// fn sample_errors(x_prob: f64, y_prob: f64, z_prob: f64, length: usize) -> PauliString {
//     let mut rng = thread_rng();
//     let paulis = Vec::with_capacity(length);
//     for ix in 0..length {
//         let mut p = Pauli::I;
//         if rng.gen_bool(x_prob) {
//             p = p.mul(&Pauli::X);
//         }
//     }
// }

struct FiveOneThree {
    stabilizers: Vec<PauliString>,
    syndrome: Option<Vec<bool>>,
}

impl FiveOneThree {
    pub fn new() -> Self {
        let s = ["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"];
        let mut stabs = Vec::new();
        for x in s {
            stabs.push(PauliString::from(x.to_string()));
        }
        FiveOneThree {
            stabilizers: stabs,
            syndrome: None,
        }
    }

    fn generate_syndrome(&self) {
        // TODO: For now I am just checking X errors
        let p = 0.05_f64;
        let mut rng = thread_rng();
        let mut e = PauliString::new();
        for _ in 0..5 {
            e.push(rng.gen_bool(p));
        }
        let mut syndrome = Vec::new();
        for stab in self.stabilizers.iter() {
            if e.anti_commutes(stab) {
                syndrome.push(true);
            } else {
                syndrome.push(false);
            }
        }
        self.syndrome = Some(syndrome);
    }
}

mod tests {
    use crate::quantum::PauliString;

    #[test]
    fn test_pauli_string() {
        let s = "IIXXIIZZY".to_string();
        let p = PauliString::from(s);
        dbg!(p);
    }
}
