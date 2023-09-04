use std::{
    collections::{HashMap, HashSet},
    hash::Hash,
    ops::{Add, AddAssign, Mul, MulAssign, Sub}, fmt::Display,
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
pub struct CyclicGroup(pub u32, pub u32);

// TODO: Eliminating these checks could introduce bugs but might be a lot faster.
impl Add<&CyclicGroup> for CyclicGroup {
    type Output = CyclicGroup;

    fn add(self, rhs: &CyclicGroup) -> Self::Output {
        if self.1 != rhs.1 {
            panic!("[FiniteField] addition not defined for different fields.");
        }
        CyclicGroup((self.0 + rhs.0) % self.1, self.1)
    }
}

impl Mul<CyclicGroup> for CyclicGroup {
    type Output = CyclicGroup;

    fn mul(self, rhs: CyclicGroup) -> Self::Output {
        self * &rhs
    }
}

impl Mul<CyclicGroup> for i32 {
    type Output = CyclicGroup;

    fn mul(self, rhs: CyclicGroup) -> Self::Output {
        let mut a = (rhs.0 as i32) * self;
        while a < 0 {
            a += rhs.1 as i32;
        }
        CyclicGroup::from((a as u32, rhs.1))
    }
}

impl Mul<i32> for CyclicGroup {
    type Output = CyclicGroup;

    fn mul(self, rhs: i32) -> Self::Output {
        let mut a = rhs * (self.0 as i32);
        while a < 0 {
            a += self.1 as i32;
        }
        CyclicGroup::from((a as u32 % self.1, self.1))
    }
}
impl Add<i32> for CyclicGroup {
    type Output = CyclicGroup;
    fn add(self, rhs: i32) -> Self::Output {
        let mut rhs_mod_n = rhs;
        while rhs_mod_n < 0 {
            rhs_mod_n += self.1 as i32;
        }
        let a = (rhs_mod_n as u32) + self.0;
        CyclicGroup::from((a % self.1, self.1))
    }
}

impl Sub<CyclicGroup> for CyclicGroup {
    type Output = CyclicGroup;

    fn sub(self, rhs: CyclicGroup) -> Self::Output {
        if self.1 != rhs.1 {
            panic!("Subtraction among non-equal CyclicGroups.");
        }
        let mut a = (self.0 as i32) - (rhs.0 as i32);
        while a < 0 {
            a += self.1 as i32;
        }
        CyclicGroup::from(((a as u32) % self.1, self.1))
    }
}
impl Add<CyclicGroup> for CyclicGroup {
    type Output = Self;

    fn add(self, rhs: CyclicGroup) -> Self::Output {
        if self.1 != rhs.1 {
            panic!("Addition among non-equal CyclicGroups")
        }
        CyclicGroup((self.0 + rhs.0) % self.1, self.1)
    }
}

impl From<(u32, u32)> for CyclicGroup {
    fn from(value: (u32, u32)) -> Self {
        CyclicGroup(value.0, value.1)
    }
}

impl From<(i32, u32)> for CyclicGroup {
    fn from(value: (i32, u32)) -> Self {
        let mut a = value.0;
        while a < 0 {
            a += value.1 as i32;
        }
        CyclicGroup((a as u32) % value.1, value.1)
    }
}

impl From<(i32, i32)> for CyclicGroup {
    fn from(value: (i32, i32)) -> Self {
        (value.0 as u32, value.1 as u32).into()
    }
}

impl Add<&CyclicGroup> for &CyclicGroup {
    type Output = CyclicGroup;

    fn add(self, rhs: &CyclicGroup) -> Self::Output {
        if self.1 != rhs.1 {
            panic!("[FiniteField] addition not defined for different fields.");
        }
        CyclicGroup((self.0 + rhs.0) % self.1, self.1)
    }
}

impl AddAssign<&CyclicGroup> for CyclicGroup {
    fn add_assign(&mut self, rhs: &CyclicGroup) {
        if self.1 != rhs.1 {
            panic!("[FiniteField] Addition not defined for different fields.")
        }
        self.0 = (self.0 + rhs.0) % self.1;
    }
}

impl Mul<&CyclicGroup> for CyclicGroup {
    type Output = CyclicGroup;

    fn mul(self, rhs: &CyclicGroup) -> Self::Output {
        if self.1 != rhs.1 {
            panic!("[FiniteField] Multiplication not defined for different fields.")
        }
        CyclicGroup((self.0 * rhs.0) % self.1, self.1)
    }
}

impl MulAssign<&CyclicGroup> for CyclicGroup {
    fn mul_assign(&mut self, rhs: &CyclicGroup) {
        if self.1 != rhs.1 {
            panic!("[FiniteField] Multiplication not defined for different fields.")
        }
        self.0 = (self.0 * rhs.0) % self.1;
    }
}

impl CyclicGroup {
    pub const ZERO: CyclicGroup = CyclicGroup(0, 0);
}

impl Display for CyclicGroup {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&format!("{:} mod {:}", self.0, self.1))
    }
}

struct FFPolynomial {
    coeffs: Vec<CyclicGroup>,
}

// TODO: What representation do I use for the equivalency classes?
struct PGL {
    mat: Array2<CyclicGroup>,
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
    let node_ids = hg.add_nodes(nodes.len());
    let node_ids: HashMap<G, u32> =
        HashMap::from_iter(nodes.clone().into_iter().zip(node_ids.into_iter()));
    for g in nodes.clone().into_iter() {
        for a in left_generators.clone() {
            for b in right_generators.clone() {
                let u = f(&a, &g);
                let v = f(&g, &b);
                let w = f(&f(&a, &g), &b);
                // TODO: Check if the edge is already present in the graph
                let first_slice = [node_ids[&g], node_ids[&u]];
                // hg.create_edge(node_ids[&g, &u]);
                // hg.create_edge(&[node_ids[&g], node_ids[&v]], 1.0);
                // hg.create_edge(
                //     &[node_ids[&g], node_ids[&u], node_ids[&v], node_ids[&w]],
                //     1.,
                // );
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
            let x = CyclicGroup(ix as u32, lattice_length);
            let y = CyclicGroup(jx as u32, lattice_length);
            points.push((x, y));
        }
    }
    let x = |a: &(CyclicGroup, CyclicGroup), b: &(CyclicGroup, CyclicGroup)| {
        let x_pos = &a.0 + &b.0;
        let y_pos = &a.1 + &b.1;
        (x_pos, y_pos)
    };
    let a_gens = HashSet::from([
        (
            CyclicGroup(1, lattice_length),
            CyclicGroup(0, lattice_length),
        ),
        (
            CyclicGroup(lattice_length - 1, lattice_length),
            CyclicGroup(0, lattice_length),
        ),
    ]);
    let b_gens = HashSet::from([
        (
            CyclicGroup(0, lattice_length),
            CyclicGroup(1, lattice_length),
        ),
        (
            CyclicGroup(0, lattice_length),
            CyclicGroup(lattice_length - 1, lattice_length),
        ),
    ]);
    left_right_cayley(points.into_iter().collect(), a_gens, b_gens, x)
}

mod tests {
    use ff::Field;

    use super::{surface_code_hgraph, CyclicGroup, Fp};

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
