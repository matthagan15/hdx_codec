use std::collections::HashMap;

use mhgl::PGraph;
use rand::prelude::*;

struct CssCode {}

#[derive(Debug, Clone, Copy)]
enum Phase {
    pr,
    pi,
    nr,
    ni,
}

#[derive(Debug, Clone, Copy)]
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

struct SurfaceCode {
    // assume lattice for now
    length: usize,
    height: usize,
    distance: f64, // TODO: might only need usize
    graph: PGraph<u32>,
    lattice_to_node: HashMap<(usize, usize), u32>,
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
