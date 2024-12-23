use std::{
    collections::{HashMap, HashSet},
    fmt::Display,
};

use crate::math::pauli::*;
use mhgl::{ConGraph, HyperGraph};
use rand::prelude::*;

struct CSSCode {}

#[derive(Debug)]
pub struct SurfaceCode {
    // TODO: this is currently based on indexes, not good!
    x_checks: Vec<PauliString>,
    z_checks: Vec<PauliString>,
    qubits: Vec<u64>,
    hgraph: ConGraph,
}

struct PauliTableau {}

trait CSS {
    fn x_checks(&self) -> PauliTableau;
    fn z_checks(&self) -> PauliTableau;
    fn stabilizers(&self) -> PauliTableau;
    fn logical_dim(&self) -> usize;
    fn physical_dim(&self) -> usize;
}

impl SurfaceCode {
    pub fn from_hgraph(hgraph: ConGraph) -> Self {
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
            let edges = hgraph.containing_edges_of_nodes(&[*x_check]);
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
            let nodes = hgraph.query_edge(&z_check).expect("z_check has no nodes?");
            let mut pauli_indices = HashSet::new();
            for ix in 0..nodes.len() {
                for jx in 0..nodes.len() {
                    if ix == jx {
                        continue;
                    }
                    let e = hgraph.find_id(&[nodes[ix], nodes[jx]]);
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

struct LDPC {
    nodes: Vec<u32>,
    lr_complex: ConGraph,
}

#[derive(Debug, Clone)]
struct FiveOneThree {
    stabilizers: Vec<PauliString>,
    syndrome: Option<Vec<bool>>,
}

impl Display for FiveOneThree {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("FiveOneThree code with stabilizers:\n");
        for s in self.stabilizers.iter() {
            f.write_fmt(format_args!("{}\n", s));
        }
        Ok(())
    }
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
        math::left_right_cayley::surface_code_hgraph,
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

        println!("{}", code);
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
