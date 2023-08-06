use std::collections::HashMap;

use rand::prelude::*;
use mhgl::PGraph;

struct CssCode {
    
}

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
    Z
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
        let paulis = input.into_bytes().into_iter().map(|x| {
            match x {
                b'I' => Pauli::I,
                b'X' => Pauli::X,
                b'Y' => Pauli::Y,
                b'Z' => Pauli::Z,
                _ => panic!("Encountered non-pauli character in input string"),
            }
        }).collect();
        PauliString { phase: Phase::pr, string: paulis }
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

struct ErrorTracker {
    pub x: bool,
    pub z: bool,
    error_probability: f64,
}

impl ErrorTracker {
    fn new(error_probability: f64) -> Self {
        ErrorTracker { x: false, z: false, error_probability }
    }
    fn flip_x(&mut self) {
        self.x = true ^ self.x;
    }

    fn flip_z(&mut self) {
        self.z = true ^ self.z;
    }

    fn meas(&self) -> (bool, bool) {
        (self.x, self.z)
    }

    fn sample_errors(&mut self) {
        let mut rng = thread_rng();
        self.x = rng.gen_bool(self.error_probability);
        self.z = rng.gen_bool(self.error_probability);
    }
}

struct SurfaceCode {
    // assume lattice for now
    length: usize,
    height: usize,
    distance: f64, // TODO: might only need usize
    graph: PGraph<u32>,
    lattice_to_node: HashMap<(usize, usize), u32>,
    lattice_index_to_errors: HashMap<(usize, usize), ErrorTracker>,
}

mod tests {
    use super::PauliString;

    #[test]
    fn test_pauli_string() {
        let s = "IIXXIIZZY".to_string();
        let p = PauliString::from(s);
        dbg!(p);
    }
}