use std::{
    collections::{HashMap, HashSet},
    fmt::Display,
    hash::Hash,
    ops::{Add, AddAssign, Mul, MulAssign, Sub},
};

use ff::*;
use mhgl::{HGraph, HyperGraph};

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

/// Returns a hypergraph representing a left-right Cayley complex
/// provided by the given set of generators.
fn left_right_cayley<G: Clone + Eq + Hash>(
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



#[derive(Debug, Clone, Hash, PartialEq, Eq)]
struct Lattice {
x_coord: CyclicGroup,
y_coord: CyclicGroup,
}

impl From<(CyclicGroup, CyclicGroup)> for Lattice {
fn from(value: (CyclicGroup, CyclicGroup)) -> Self {
    Lattice { x_coord: value.0, y_coord: value.1 }
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
        ).into(),
        (
            CyclicGroup(lattice_length - 1, lattice_length),
            CyclicGroup(0, lattice_length),
        ).into(),
    ]);
    let b_gens:HashSet<Lattice> = HashSet::from([
        (
            CyclicGroup(0, lattice_length),
            CyclicGroup(1, lattice_length),
        ).into(),
        (
            CyclicGroup(0, lattice_length),
            CyclicGroup(lattice_length - 1, lattice_length),
        ).into(),
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
