use std::{
    collections::{HashMap, HashSet},
    fmt::Display,
    hash::Hash,
    ops::{Add, AddAssign, Mul, MulAssign, Sub},
};

use crate::lps::*;

use mhgl::{HGraph, HyperGraph};
use uuid::Uuid;

use crate::math::pauli::*;

use crate::math::finite_field::FiniteField;

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
    // TODO: also need to add information about the vertex type and edge type.
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
pub fn surface_code_hgraph() -> HGraph {
    let lattice_length = 3;
    let mut points: Vec<Lattice> = Vec::new();
    for ix in 0..lattice_length {
        for jx in 0..lattice_length {
            let x = FiniteField(ix, lattice_length);
            let y = FiniteField(jx, lattice_length);
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
