use std::{
    collections::{HashMap, HashSet},
    fmt::Display,
    hash::Hash,
    ops::{Add, AddAssign, Mul, MulAssign, Sub},
};

use crate::lps::*;
use ff::*;
use mhgl::{HGraph, HyperGraph};
use uuid::Uuid;

use crate::math::pauli::*;

use crate::math::finite_field::CyclicGroup as CyclicGroup;

/// Returns a hypergraph representing a left-right Cayley complex
/// provided by the given set of generators.
fn left_right_cayley<G>(
    nodes: HashSet<G>,
    left_generators: HashSet<G>,
    right_generators: HashSet<G>,
) -> HGraph
where
    G: Clone + Eq + Hash + Mul<G, Output = G>,
{
    let mut hg = HGraph::new();
    let node_ids = hg.add_nodes(nodes.len());
    let node_ids: HashMap<G, u32> =
        HashMap::from_iter(nodes.clone().into_iter().zip(node_ids.into_iter()));
    for g in nodes.into_iter() {
        for a in left_generators.clone() {
            for b in right_generators.clone() {
                // TODO: replace this with reference multiplication
                let u = a.clone() * g.clone();
                let v = g.clone() * b.clone();
                let w = a.clone() * (g.clone() * b);

                let g_u_slice = [node_ids[&g], node_ids[&u]];
                let g_v_slice = [node_ids[&g], node_ids[&v]];
                let g_u_v_w_slice = [node_ids[&g], node_ids[&v], node_ids[&u], node_ids[&w]];
                hg.create_edge(&g_u_slice);
                hg.create_edge(&g_v_slice);
                hg.create_edge(&g_u_v_w_slice);
            }
        }
    }
    hg
}

struct QLDPC {
    hgraph: HGraph,
    x_stabilizers: Vec<PauliString>,
    z_stabilizers: Vec<PauliString>,
    qubits: Vec<Uuid>,
}

impl QLDPC {
    pub fn from(p: u32, q: u32) -> Option<Self> {
        let mut hg = HGraph::new();
        match legendre_symbol(p as i32, q as i32) {
            1 => {
                if p + 1 > ((q.pow(3) - q) - 1) {
                    println!("Degree cannot exceed the group order.");
                    return None;
                }
                let generators: HashSet<PGL2> = compute_generators(p, q)
                    .into_iter()
                    .filter_map(|x| PGL2::from_coeffs(x))
                    .collect();
                let matrices: HashSet<PGL2> = generate_all_pgl2(q).into_iter().collect();
                let hg = left_right_cayley(matrices, generators.clone(), generators.clone());
                None
            }
            -1 => None,
            _ => None,
        }
    }
}

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
struct Lattice {
    x_coord: CyclicGroup,
    y_coord: CyclicGroup,
}

impl From<(CyclicGroup, CyclicGroup)> for Lattice {
    fn from(value: (CyclicGroup, CyclicGroup)) -> Self {
        Lattice {
            x_coord: value.0,
            y_coord: value.1,
        }
    }
}
impl Mul for Lattice {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        (self.x_coord + rhs.x_coord, self.y_coord + rhs.y_coord).into()
    }
}

/// Goal here is to generate the stabilizers for the surface code from a left-right cayley complex. Should return a StabilizerCode object?
pub fn surface_code_hgraph() -> HGraph {
    let lattice_length = 3_u32;
    let mut points: Vec<Lattice> = Vec::new();
    for ix in 0..lattice_length {
        for jx in 0..lattice_length {
            let x = CyclicGroup(ix as u32, lattice_length);
            let y = CyclicGroup(jx as u32, lattice_length);
            points.push((x, y).into());
        }
    }

    let a_gens: HashSet<Lattice> = HashSet::from([
        (
            CyclicGroup(1, lattice_length),
            CyclicGroup(0, lattice_length),
        )
            .into(),
        (
            CyclicGroup(lattice_length - 1, lattice_length),
            CyclicGroup(0, lattice_length),
        )
            .into(),
    ]);
    let b_gens: HashSet<Lattice> = HashSet::from([
        (
            CyclicGroup(0, lattice_length),
            CyclicGroup(1, lattice_length),
        )
            .into(),
        (
            CyclicGroup(0, lattice_length),
            CyclicGroup(lattice_length - 1, lattice_length),
        )
            .into(),
    ]);
    left_right_cayley(points.into_iter().collect(), a_gens, b_gens)
}

mod tests {
    use ff::Field;

    use super::{surface_code_hgraph, CyclicGroup};

    #[test]
    #[should_panic]
    fn test_finite_field_add_nonequal() {
        let a = CyclicGroup(1, 7);
        let b = CyclicGroup(3, 8);
        let _ = a + &b;
    }

    #[test]
    #[should_panic]
    fn test_finite_field_mul_nonequal() {
        let a = CyclicGroup(1, 7);
        let b = CyclicGroup(3, 8);
        let _ = a * &b;
    }

    #[test]
    fn test_finite_field_add() {
        let a = CyclicGroup(5, 7);
        let b = CyclicGroup(3, 7);
        assert_eq!(CyclicGroup(1, 7), a + &b);
    }

    #[test]
    fn test_finite_field_mul() {
        let a = CyclicGroup(5, 7);
        let b = CyclicGroup(3, 7);
        assert_eq!(CyclicGroup(1, 7), a * &b);
        println!("a = {:}", a);
    }

    #[test]
    fn test_surface_code_graph() {
        dbg!(surface_code_hgraph());
    }
}
