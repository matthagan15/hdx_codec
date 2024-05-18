use std::{
    collections::{HashMap, HashSet},
    fmt::Display,
    hash::Hash,
    ops::{Add, AddAssign, Mul, MulAssign, Sub},
};

use crate::lps::*;

use mhgl::{ConGraph, HGraph, HyperGraph};
use uuid::Uuid;

use crate::math::pauli::*;

use crate::math::finite_field::FiniteField;

struct LeftRightCayley<N>
where
    N: Hash,
{
    nodes_to_type: HashMap<N, u8>,
}

/// Returns a hypergraph representing a left-right Cayley complex
/// provided by the given set of generators.
fn left_right_cayley<G>(
    nodes: HashSet<G>,
    left_generators: HashSet<G>,
    right_generators: HashSet<G>,
) -> ConGraph
where
    G: Clone + Eq + Hash + Mul<G, Output = G>,
{
    // TODO: also need to add information about the vertex type and edge type.
    let mut hg = ConGraph::new();
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
                hg.add_edge(&g_u_slice);
                hg.add_edge(&g_v_slice);
                hg.add_edge(&g_u_v_w_slice);
            }
        }
    }
    hg
}

struct QLDPC {
    hgraph: ConGraph,
    x_stabilizers: Vec<PauliString>,
    z_stabilizers: Vec<PauliString>,
    qubits: Vec<Uuid>,
}

impl QLDPC {
    pub fn from(p: u32, q: u32) -> Option<Self> {
        let mut hg = ConGraph::new();
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
    x_coord: FiniteField,
    y_coord: FiniteField,
}

impl From<(FiniteField, FiniteField)> for Lattice {
    fn from(value: (FiniteField, FiniteField)) -> Self {
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

// TODO: This is garbage, but it gives the right stabilizers. The whole point
// of the A, B generators for the left-right cayley complex is they both
// generate the entire group, whereas here they don't.
/// Goal here is to generate the stabilizers for the surface code from a left-right cayley complex. Should return a StabilizerCode object?
pub fn surface_code_hgraph() -> ConGraph {
    let lattice_length = 3_u32;
    let mut points: Vec<Lattice> = Vec::new();
    for ix in 0..lattice_length {
        for jx in 0..lattice_length {
            let x = FiniteField(ix as u32, lattice_length);
            let y = FiniteField(jx as u32, lattice_length);
            points.push((x, y).into());
        }
    }

    let a_gens: HashSet<Lattice> = HashSet::from([
        (
            FiniteField(1, lattice_length),
            FiniteField(0, lattice_length),
        )
            .into(),
        (
            FiniteField(lattice_length - 1, lattice_length),
            FiniteField(0, lattice_length),
        )
            .into(),
    ]);
    let b_gens: HashSet<Lattice> = HashSet::from([
        (
            FiniteField(0, lattice_length),
            FiniteField(1, lattice_length),
        )
            .into(),
        (
            FiniteField(0, lattice_length),
            FiniteField(lattice_length - 1, lattice_length),
        )
            .into(),
    ]);
    left_right_cayley(points.into_iter().collect(), a_gens, b_gens)
}

mod tests {

    use super::{surface_code_hgraph, FiniteField};

    #[test]
    #[should_panic]
    fn test_finite_field_add_nonequal() {
        let a = FiniteField(1, 7);
        let b = FiniteField(3, 8);
        let _ = a + &b;
    }

    #[test]
    #[should_panic]
    fn test_finite_field_mul_nonequal() {
        let a = FiniteField(1, 7);
        let b = FiniteField(3, 8);
        let _ = a * &b;
    }

    #[test]
    fn test_finite_field_add() {
        let a = FiniteField(5, 7);
        let b = FiniteField(3, 7);
        assert_eq!(FiniteField(1, 7), a + &b);
    }

    #[test]
    fn test_finite_field_mul() {
        let a = FiniteField(5, 7);
        let b = FiniteField(3, 7);
        assert_eq!(FiniteField(1, 7), a * &b);
        println!("a = {:}", a);
    }

    #[test]
    fn test_surface_code_graph() {
        dbg!(surface_code_hgraph());
    }
}
