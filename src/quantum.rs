use std::{
    collections::{HashMap, HashSet},
    fmt::Display,
};

use mhgl::{HGraph, HyperGraph, PGraph, SparseBasis};
use rand::prelude::*;

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
enum Pauli {
    I,
    X,
    Y,
    Z,
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
            (Pauli::Z, Pauli::I) => (Phase::pr, Pauli::Z),
            (Pauli::Z, Pauli::X) => (Phase::pi, Pauli::Y),
            (Pauli::Z, Pauli::Y) => (Phase::ni, Pauli::X),
            (Pauli::Z, Pauli::Z) => (Phase::pr, Pauli::I),
        }
    }
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

#[derive(Debug)]
pub struct SurfaceCode {
    // TODO: this is currently based on indexes, not good!
    x_checks: Vec<PauliString>,
    z_checks: Vec<PauliString>,
    qubits: Vec<u128>,
    hgraph: HGraph,
}

impl SurfaceCode {
    pub fn from_hgraph(hgraph: HGraph) -> Self {
        let mut x_pauli_strings = Vec::new();
        let mut z_pauli_strings = Vec::new();
        let nodes = hgraph.nodes();
        let qubits = hgraph.edges_of_size(2);
        let mut qubit_to_index = HashMap::new();
        for ix in 0..qubits.len() as usize {
            qubit_to_index.insert(qubits[ix], ix);
        }
        let z_checks = hgraph.edges_of_size(4);
        for x_check in nodes.iter() {
            let edges = hgraph.get_containing_edges(&[*x_check]);
            let mut pauli_indices = HashSet::new();
            for edge in edges {
                if qubit_to_index.contains_key(&edge) {
                    pauli_indices.insert(qubit_to_index.get(&edge).unwrap());
                }
            }
            let mut pauli_string = PauliString::new();
            for ix in 0..qubits.len() as usize {
                if pauli_indices.contains(&ix) {
                    pauli_string.push(Pauli::X);
                } else {
                    pauli_string.push(Pauli::I);
                }
            }
            x_pauli_strings.push(pauli_string);
        }
        for z_check in z_checks {
            let nodes = hgraph.query_nodes_in_edge(&z_check);
            let mut pauli_indices = HashSet::new();
            for ix in 0..nodes.len() {
                for jx in 0..nodes.len() {
                    if ix == jx {
                        continue;
                    }
                    let e = hgraph.query_edge_id(&[nodes[ix], nodes[jx]]);
                    if e.is_none() {
                        continue;
                    }
                    let id = e.unwrap();
                    if qubit_to_index.contains_key(&id) {
                        pauli_indices.insert(*qubit_to_index.get(&id).unwrap());
                    }
                }
            }
            let mut ps = PauliString::new();
            for ix in 0..qubits.len() {
                if pauli_indices.contains(&ix) {
                    ps.push(Pauli::Z);
                } else {
                    ps.push(Pauli::I);
                }
            }
            z_pauli_strings.push(ps);
        }
        SurfaceCode {
            x_checks: x_pauli_strings,
            z_checks: z_pauli_strings,
            qubits,
            hgraph,
        }
    }

    pub fn print_stabilizers(&self) {
        println!("X-Stabilizers\n");
        for ix in 0..self.x_checks.len() {
            println!("{:}: {:}", ix, self.x_checks[ix]);
        }
        println!("Z-Stabilizers\n");
        for ix in 0..self.z_checks.len() {
            println!("{:}: {:}", ix, self.z_checks[ix]);
        }
    }

    fn sample_error(&self) {
        let prob_x_error = 0.01_f64;
        let prob_z_error = 0.01_f64;
        let mut rng = thread_rng();
        let mut error_pauli_string = PauliString::new();
        for ix in 0..self.qubits.len() {
            let x_err = rng.gen_bool(prob_x_error);
            let z_err = rng.gen_bool(prob_z_error);
            error_pauli_string.push(match (x_err, z_err) {
                (true, true) => Pauli::Y,
                (true, false) => Pauli::X,
                (false, true) => Pauli::Z,
                (false, false) => Pauli::I,
            });
        }
        println!("sampled error: {:}", error_pauli_string);
    }
}

#[derive(Debug, Clone)]
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

    fn sample_error(&mut self) {
        let p_x = 0.2_f64;
        let p_z = 0.2_f64;
        let mut rng = thread_rng();
        let mut e = PauliString::new();
        for _ in 0..5 {
            let x_sample = rng.gen_bool(p_x);
            let z_sample = rng.gen_bool(p_z);
            match (x_sample, z_sample) {
                (false, false) => e.push(Pauli::I),
                (false, true) => e.push(Pauli::Z),
                (true, false) => e.push(Pauli::X),
                (true, true) => e.push(Pauli::Y),
            }
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

    /// For now this should just return the pauli string
    /// corresponding to the best guess of the error that
    /// occurred.
    fn decode(&self) {}
}

mod tests {
    use super::{FiveOneThree, SurfaceCode};
    use crate::{
        left_right_cayley::surface_code_hgraph,
        quantum::{Pauli, PauliString, Phase},
    };

    #[test]
    fn test_pauli_string() {
        let s = "IIXXIIZZY".to_string();
        let t = "XIXZIIYZX".to_string();
        let p = PauliString::from(s);
        let q = PauliString::from(t);
        dbg!(p.mul(&q));
        dbg!(q.mul(&p));
        dbg!(p);
        let left = PauliString {
            phase: Phase::pr,
            string: [Pauli::X, Pauli::I, Pauli::I, Pauli::X, Pauli::X].to_vec(),
        };
        let rhs = PauliString {
            phase: Phase::pr,
            string: [Pauli::I, Pauli::X, Pauli::Z, Pauli::Z, Pauli::X].to_vec(),
        };
        let l_then_r = left.mul(&rhs);
        let r_then_l = rhs.mul(&left);
        dbg!(l_then_r);
        dbg!(r_then_l);
    }

    #[test]
    fn test_5_1_3_syndrome() {
        let mut code = FiveOneThree::new();
        code.sample_error();

        dbg!(code);
    }

    #[test]
    fn test_surface_from_lr_cayley() {
        let hg = surface_code_hgraph();
        let sc = SurfaceCode::from_hgraph(hg);
        sc.print_stabilizers();
    }

    #[test]
    fn test_surface_code_error_model() {
        let sc = SurfaceCode::from_hgraph(surface_code_hgraph());
        sc.sample_error();
    }

    #[test]
    fn test_pauli_string_display() {
        let mut ps = PauliString::new();
        ps.set_phase(Phase::ni);
        ps.push(Pauli::I);
        ps.push(Pauli::X);
        ps.push(Pauli::Y);
        ps.push(Pauli::Z);
        println!("pauli displayed: {:}", ps);
    }
}
