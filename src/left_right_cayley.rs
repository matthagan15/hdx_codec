use std::{
    collections::{HashMap, HashSet},
    hash::Hash,
    ops::{Add, AddAssign, Mul, MulAssign},
};

use ff::*;
use mhgl::{HGraph, HyperGraph};
use ndarray::Array2;


#[derive(PrimeField)]
#[PrimeFieldModulus = "199"]
#[PrimeFieldGenerator = "1"]
#[PrimeFieldReprEndianness = "little"]
struct Fp([u64; 1]);

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Eq, Ord, Hash)]
struct FiniteField(u32, u32);

// TODO: Eliminating these checks could introduce bugs but might be a lot faster.
impl Add<&FiniteField> for FiniteField {
    type Output = FiniteField;

    fn add(self, rhs: &FiniteField) -> Self::Output {
        if self.1 != rhs.1 {
            panic!("[FiniteField] addition not defined for different fields.");
        }
        FiniteField((self.0 + rhs.0) % self.1, self.1)
    }
}

impl Add<&FiniteField> for &FiniteField {
    type Output = FiniteField;

    fn add(self, rhs: &FiniteField) -> Self::Output {
        if self.1 != rhs.1 {
            panic!("[FiniteField] addition not defined for different fields.");
        }
        FiniteField((self.0 + rhs.0) % self.1, self.1)
    }
}

impl AddAssign<&FiniteField> for FiniteField {
    fn add_assign(&mut self, rhs: &FiniteField) {
        if self.1 != rhs.1 {
            panic!("[FiniteField] Addition not defined for different fields.")
        }
        self.0 = (self.0 + rhs.0) % self.1;
    }
}

impl Mul<&FiniteField> for FiniteField {
    type Output = FiniteField;

    fn mul(self, rhs: &FiniteField) -> Self::Output {
        if self.1 != rhs.1 {
            panic!("[FiniteField] Multiplication not defined for different fields.")
        }
        FiniteField((self.0 * rhs.0) % self.1, self.1)
    }
}

impl MulAssign<&FiniteField> for FiniteField {
    fn mul_assign(&mut self, rhs: &FiniteField) {
        if self.1 != rhs.1 {
            panic!("[FiniteField] Multiplication not defined for different fields.")
        }
        self.0 = (self.0 * rhs.0) % self.1;
    }
}

//
struct PGL {
    mat: Array2<FiniteField>,
}

fn left_right_cayley<G: Clone + Eq + Hash, F>(
    nodes: HashSet<G>,
    left_generators: HashSet<G>,
    right_generators: HashSet<G>,
    f: F,
) -> HGraph
where
    F: Fn(&G, &G) -> G,
{
    let mut hg = HGraph::new();
    let node_ids = hg.create_nodes(nodes.len());
    let node_ids: HashMap<G, u128> =
        HashMap::from_iter(nodes.clone().into_iter().zip(node_ids.into_iter()));
    for g in nodes.clone().into_iter() {
        for a in left_generators.clone() {
            for b in right_generators.clone() {
                let u = f(&a, &g);
                let v = f(&g, &b);
                let w = f(&f(&a, &g), &b);
                // TODO: Check if the edge is already present in the graph
                hg.create_edge(&[node_ids[&g], node_ids[&u]], 1.);
                hg.create_edge(&[node_ids[&g], node_ids[&v]], 1.0);
                hg.create_edge(
                    &[node_ids[&g], node_ids[&u], node_ids[&v], node_ids[&w]],
                    1.,
                );
            }
        }
    }
    hg
}

/// Goal here is to generate the stabilizers for the surface code from a left-right cayley complex. Should return a StabilizerCode object?
pub fn surface_code_hgraph() -> HGraph {
    let lattice_length = 3_u32;
    let mut points = Vec::new();
    for ix in 0..lattice_length {
        for jx in 0..lattice_length {
            let x = FiniteField(ix as u32, lattice_length);
            let y = FiniteField(jx as u32, lattice_length);
            points.push((x, y));
        }
    }
    let x = |a: &(FiniteField, FiniteField), b: &(FiniteField, FiniteField)| {
        let x_pos = &a.0 + &b.0;
        let y_pos = &a.1 + &b.1;
        (x_pos, y_pos)
    };
    let a_gens = HashSet::from([
        (
            FiniteField(1, lattice_length),
            FiniteField(0, lattice_length),
        ),
        (
            FiniteField(lattice_length - 1, lattice_length),
            FiniteField(0, lattice_length),
        ),
    ]);
    let b_gens = HashSet::from([
        (
            FiniteField(0, lattice_length),
            FiniteField(1, lattice_length),
        ),
        (
            FiniteField(0, lattice_length),
            FiniteField(lattice_length - 1, lattice_length),
        ),
    ]);
    left_right_cayley(points.into_iter().collect(), a_gens, b_gens, x)
}

mod tests {
    use ff::Field;

    use super::{ FiniteField, Fp, surface_code_hgraph};

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
    }

    #[test]
    fn test_surface_code_graph() {
        dbg!(surface_code_hgraph());
    }
}
