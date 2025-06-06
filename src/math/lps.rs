use std::{
    collections::{HashMap, HashSet},
    fmt::Display,
    ops::{Index, Mul},
};

use mhgl::ConGraph;
use serde::{Deserialize, Serialize};

use crate::math::finite_field::FiniteField;

use crate::matrices::ffmatrix::FFMatrix;

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Hash, PartialEq, Eq)]
pub struct GeneralSquaresSolution(pub i32, pub i32, pub i32, pub i32);

impl From<(i32, i32, i32, i32)> for GeneralSquaresSolution {
    fn from(value: (i32, i32, i32, i32)) -> Self {
        GeneralSquaresSolution(value.0, value.1, value.2, value.3)
    }
}

impl From<Vec<i32>> for GeneralSquaresSolution {
    fn from(value: Vec<i32>) -> Self {
        if value.len() == 0 {
            (0, 0, 0, 0).into()
        } else if value.len() == 1 {
            (value[0], 0, 0, 0).into()
        } else if value.len() == 2 {
            (value[0], value[1], 0, 0).into()
        } else if value.len() == 3 {
            (value[0], value[1], value[2], 0).into()
        } else {
            (value[0], value[1], value[2], value[3]).into()
        }
    }
}

impl GeneralSquaresSolution {
    fn equivalent_forms(&self) -> Vec<Self> {
        let mut ret = Vec::new();
        let og = vec![self.0, self.1, self.2, self.3];
        let shuffled = shuffle(og);
        for v in shuffled {
            let mut signed = generate_signs(v);
            ret.append(&mut signed);
        }
        ret.into_iter().map(|v| v.into()).collect()
    }
}

// TODO: Change the PGL2 code over to use FFMatrix. This will reduce codespace down by at least 100-200 lines and just be nicer to read.
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct FFMatrixPGL2 {
    matrix: FFMatrix,
}

#[derive(Debug, Clone, Copy, Hash, PartialEq, Eq)]
pub struct PGL2 {
    pub coeffs: [FiniteField; 4],
}

impl PGL2 {
    pub fn det(&self) -> FiniteField {
        self.coeffs[0] * self.coeffs[3] - self.coeffs[1] * self.coeffs[2]
    }
    /// Returns the identity matrix [[1, 0], [0, 1]];
    pub fn identity(mod_n: u32) -> Self {
        PGL2 {
            coeffs: [
                FiniteField(1, mod_n),
                FiniteField(0, mod_n),
                FiniteField(0, mod_n),
                FiniteField(1, mod_n),
            ],
        }
    }

    /// Returns the canonical form of the member of PGL2 from the given
    /// coefficients. If the determinant is zero then this returns `None`.
    pub fn from_coeffs(coeffs: [FiniteField; 4]) -> Option<Self> {
        let det = coeffs[0] * &coeffs[3] - coeffs[1] * &coeffs[2];
        if det.0 == 0 {
            None
        } else {
            let mod_inv = if coeffs[0].0 != 0 {
                modular_inverse(coeffs[0].0 as i32, coeffs[0].1 as i32)
            } else {
                modular_inverse(coeffs[2].0 as i32, coeffs[2].1 as i32)
            };
            if mod_inv.is_none() {
                None
            } else {
                let mod_inv = mod_inv.unwrap();
                let mat = PGL2 {
                    coeffs: [
                        coeffs[0] * mod_inv,
                        coeffs[1] * mod_inv,
                        coeffs[2] * mod_inv,
                        coeffs[3] * mod_inv,
                    ],
                };
                Some(mat)
            }
        }
    }
}

impl Mul<i32> for PGL2 {
    type Output = PGL2;

    fn mul(self, rhs: i32) -> Self::Output {
        let new_coeffs = [
            self.coeffs[0] * rhs,
            self.coeffs[1] * rhs,
            self.coeffs[2] * rhs,
            self.coeffs[3] * rhs,
        ];
        PGL2 { coeffs: new_coeffs }
    }
}

impl Mul<&Self> for PGL2 {
    type Output = PGL2;

    fn mul(self, rhs: &Self) -> Self::Output {
        // TODO: fill in with the matrix mul alg.
        let new: [FiniteField; 4] = [
            self[[0, 0]] * rhs[[0, 0]] + self[[0, 1]] * rhs[[1, 0]],
            self[[0, 0]] * rhs[[0, 1]] + self[[0, 1]] * rhs[[1, 1]],
            self[[1, 0]] * rhs[[0, 0]] + self[[1, 1]] * rhs[[1, 0]],
            self[[1, 0]] * rhs[[0, 1]] + self[[1, 1]] * rhs[[1, 1]],
        ];
        PGL2::from_coeffs(new).expect("Could not multiply matrices.")
    }
}

impl Mul<Self> for PGL2 {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        self * &rhs
    }
}

impl Index<[usize; 2]> for PGL2 {
    type Output = FiniteField;

    /// Matrices are zero-indexed
    fn index(&self, index: [usize; 2]) -> &Self::Output {
        let ix = index[0] * 2 + index[1];
        &self.coeffs[ix]
    }
}

impl From<[FiniteField; 4]> for PGL2 {
    /// Note: if an invalid matrix (i.e. a non-invertible matrix) is given
    /// this will panic.
    fn from(value: [FiniteField; 4]) -> Self {
        PGL2::from_coeffs(value).expect("[PGL2] Tried creating PGL2 with zero determinant matrix.")
    }
}

impl Display for PGL2 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = format!(
            "[[{:}, {:}]\n[{:}, {:}]] mod {:}",
            self.coeffs[0].0,
            self.coeffs[1].0,
            self.coeffs[2].0,
            self.coeffs[3].0,
            self.coeffs[0].1
        );
        f.write_str(&s)
    }
}

pub fn generate_all_pgl2(mod_p: u32) -> Vec<PGL2> {
    let mut ret = Vec::with_capacity(mod_p.pow(4) as usize);
    // x_{1,1} is nonzero
    for b in 0..mod_p {
        for c in 0..mod_p {
            for d in 0..mod_p {
                if let Some(mat) = PGL2::from_coeffs([
                    FiniteField(1, mod_p),
                    FiniteField(b, mod_p),
                    FiniteField(c, mod_p),
                    FiniteField(d, mod_p),
                ]) {
                    ret.push(mat);
                }
            }
        }
    }
    // x_{1,1} is zero, x_{2,1} is non-zero
    for b in 0..mod_p {
        for d in 0..mod_p {
            if let Some(mat) = PGL2::from_coeffs([
                FiniteField(0, mod_p),
                FiniteField(b, mod_p),
                FiniteField(1, mod_p),
                FiniteField(d, mod_p),
            ]) {
                ret.push(mat);
            }
        }
    }
    ret
}

/// 2 x 2 matrix with unit determinant in canonical form
/// where the first non-zero entry in the first column is
/// in the range {1, ..., (p - 1) / 2}
#[derive(Debug, Clone, Copy, Hash, Eq, PartialEq)]
pub struct PSL2 {
    matrix: PGL2,
}

impl PSL2 {
    pub fn from(value: PGL2) -> Option<Self> {
        let det = value.det();
        let prime_sqrt_sols = prime_mod_sqrt(det.0 as i32, det.1 as i32);
        if prime_sqrt_sols.len() == 0 {
            return None;
        }
        let sqrt_det = prime_sqrt_sols[0];
        let sqrt_det_inv =
            modular_inverse(sqrt_det, det.1 as i32).expect("No inverse for PSL determinant");
        let det_normalized_matrix = value * sqrt_det_inv;
        let standard_range: HashSet<u32> = (1..=(det.1 - 1) / 2).collect();
        let first_nonzero = if det_normalized_matrix.coeffs[0].0 != 0 {
            det_normalized_matrix.coeffs[0].0
        } else {
            det_normalized_matrix.coeffs[2].0
        };
        if standard_range.contains(&first_nonzero) {
            Some(PSL2 {
                matrix: det_normalized_matrix,
            })
        } else {
            Some(PSL2 {
                matrix: det_normalized_matrix * -1,
            })
        }
    }

    pub fn from_coeffs(coeffs: [FiniteField; 4]) -> Option<Self> {
        PGL2::from_coeffs(coeffs)
            .map(|pgl2| PSL2::from(pgl2))
            .flatten()
    }
}

impl Mul for PSL2 {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let pgl2 = self.matrix * rhs.matrix;
        PSL2::from(pgl2).expect("PSL2 should be closed under multiplication.")
    }
}

impl Display for PSL2 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = format!(
            "[[{:}, {:}],\n[{:}, {:}]] mod {:}",
            self.matrix.coeffs[0].0,
            self.matrix.coeffs[1].0,
            self.matrix.coeffs[2].0,
            self.matrix.coeffs[3].0,
            self.matrix.coeffs[0].1
        );
        f.write_str(&s)
    }
}

pub fn generate_all_psl2(mod_p: u32) -> Vec<PSL2> {
    let pgl2_matrices = generate_all_pgl2(mod_p);
    pgl2_matrices
        .into_iter()
        .filter_map(|m| PSL2::from(m))
        .collect()
}

/// Returns all shuffled versions of the vec
fn shuffle(input: Vec<i32>) -> Vec<Vec<i32>> {
    if input.len() == 0 {
        Vec::new()
    } else if input.len() == 1 {
        vec![input]
    } else {
        let mut ret = Vec::new();
        for ix in 0..input.len() {
            let mut removed = input.clone();
            let removed_number = removed.remove(ix);
            let mut remainder = shuffle(removed);
            for v in remainder.iter_mut() {
                v.insert(0, removed_number);
            }
            ret.append(&mut remainder);
        }
        ret
    }
}

// TODO: use the implementation in math::finite_field::FiniteField.
/// Computes the modular inverse, returns `a^{-1}` such taht
/// `a * a^{-1} mod n == 1 mod n`
fn modular_inverse(a: i32, n: i32) -> Option<i32> {
    let mut t = 0;
    let mut r = n;
    let mut new_t = 1;
    let mut new_r = a;
    while new_r != 0 {
        if let Some(q) = r.checked_div(new_r) {
            (t, new_t) = (new_t, t - q * new_t);
            (r, new_r) = (new_r, r - q * new_r);
        } else {
            return None;
        }
    }
    if r > 1 {
        None
    } else {
        while t < 0 {
            t += n;
        }
        Some(t)
    }
}

// TODO: Either use an implementation in crate::math::finite_field::FiniteField or move this one there.
/// Taken from `[Rust Programming](https://rustp.org/number-theory/modular-exponentiation/)`
/// computes `n^x mod p` using repeated squares
fn modular_exponent(mut n: i32, mut x: i32, p: i32) -> i32 {
    let mut ans = 1;
    if x <= 0 {
        return 1;
    }
    if n == 0 {
        return 0;
    }
    if n < 0 {
        while n < 0 {
            n += p;
        }
    }
    loop {
        if x == 1 {
            return (ans * n) % p;
        }
        if x & 1 == 0 {
            n = (n * n) % p;
            x >>= 1;
            continue;
        } else {
            ans = (ans * n) % p;
            x -= 1;
        }
    }
}

// TODO: Also move this to the FF
pub fn legendre_symbol(a: i32, p: i32) -> i32 {
    let ls = modular_exponent(a, (p - 1) / 2, p);
    if ls == p - 1 {
        -1
    } else {
        ls
    }
}

// TODO: Also move this to the FF
/// Returns solutions to `x^2 = a mod p`.
fn prime_mod_sqrt(a: i32, p: i32) -> Vec<i32> {
    let reduced_a = a % p;
    if reduced_a == 0 {
        return vec![0];
    }
    if p == 2 {
        return vec![reduced_a];
    }
    if legendre_symbol(a, p) != 1 {
        return Vec::new();
    }

    if p % 4 == 3 {
        let x = modular_exponent(reduced_a, (p + 1) / 4, p);
        return vec![x, p - x];
    }

    let mut q = p - 1;
    let mut s = 0;
    while q % 2 == 0 {
        s += 1;
        q /= 2;
    }

    let mut z = 1;
    while legendre_symbol(z, p) != -1 {
        z += 1;
    }
    let mut c = modular_exponent(z, q, p);

    let mut x = modular_exponent(reduced_a, (q + 1) / 2, p);
    let mut t = modular_exponent(reduced_a, q, p);
    while t != 1 {
        let mut i = 0_i32;
        for ix in 1..s as u32 {
            if modular_exponent(t, 2_i32.pow(ix), p) == 1 {
                i = ix as i32;
                break;
            }
        }
        let b = if s - i - 1 < 0 {
            1
        } else {
            modular_exponent(c, 2_i32.pow((s - i - 1) as u32), p)
        };
        x = (x * b) % p;
        t = (t * b * b) % p;
        c = (b * b) % p;
        s = i;
    }
    vec![x, p - x]
}

pub fn generate_signs(mut input: Vec<i32>) -> Vec<Vec<i32>> {
    if input.len() == 1 {
        return vec![vec![input[0]], vec![-input[0]]];
    }
    let first_elem = input.remove(0);
    let mut remaining_elems = generate_signs(input);
    let mut copy = remaining_elems.clone();
    for v in remaining_elems.iter_mut() {
        v.insert(0, first_elem);
    }
    for v in copy.iter_mut() {
        v.insert(0, -1 * first_elem);
    }
    remaining_elems.append(&mut copy);
    remaining_elems
}

/// Solves a_1^2 + a_2^2 + a_3^2 + a_4^2 - q == 0
/// If limit is `None` then returns all possible solutions, else limit contains the maximum number of solutions to return
pub fn diophantine_squares_solver(q: i32, limit: Option<usize>) -> Vec<GeneralSquaresSolution> {
    let mut ret = HashSet::new();
    // Store the canonical form of visited sols.
    // let mut visited = HashSet::new();
    'outer_loop: for a_1 in 0..=q {
        for a_2 in a_1..=q - a_1 {
            for a_3 in a_2..=q - a_1 - a_2 {
                for a_4 in a_3..=q - a_1 - a_2 - a_3 {
                    let sum = a_1.pow(2) + a_2.pow(2) + a_3.pow(2) + a_4.pow(2);
                    if sum - q == 0 {
                        let sol: GeneralSquaresSolution = (a_1, a_2, a_3, a_4).into();
                        let all_sols = sol.equivalent_forms();
                        for sol in all_sols {
                            ret.insert(sol);
                        }
                        if ret.len() == limit.map_or(usize::MAX, |v| v) {
                            break 'outer_loop;
                        }
                    }
                }
            }
        }
    }
    ret.into_iter().collect()
}

fn reduce_diophantine_solutions(
    sols: Vec<GeneralSquaresSolution>,
    mod_p: u32,
) -> Vec<GeneralSquaresSolution> {
    sols.into_iter()
        .filter(|x| {
            if mod_p % 4 == 1 {
                x.0 > 0 && x.0 % 2 == 1
            } else if mod_p % 4 == 3 {
                let cond_1 = x.0 > 0 && x.0 % 2 == 0;
                let cond_2 = x.0 == 0 && x.1 > 0;
                cond_1 || cond_2
            } else {
                false
            }
        })
        .collect()
}

pub fn compute_generators(p: u32, q: u32) -> Vec<[FiniteField; 4]> {
    let solutions = diophantine_squares_solver(p as i32, None);
    let reduced_sols = reduce_diophantine_solutions(solutions, p);
    let mut keepers = HashSet::new();
    for sol in reduced_sols {
        let negated_sol = GeneralSquaresSolution(sol.0, -1 * sol.1, -1 * sol.2, -1 * sol.3);
        if keepers.contains(&negated_sol) == false {
            keepers.insert(sol);
        }
    }
    if let Some((x, y)) = solve_mod(q) {
        let cyclic_x = FiniteField::from((x, q));
        let cyclic_y = FiniteField::from((y, q));
        keepers
            .into_iter()
            .map(|sol| {
                let (a, b, c, d) = (sol.0, sol.1, sol.2, sol.3);
                let coeffs = [
                    (cyclic_x * b + a) + cyclic_y * d,
                    (cyclic_y * -b) + c + cyclic_x * d,
                    (cyclic_y * -b) + (-1 * c) + cyclic_x * d,
                    (FiniteField::from((a, q)) - cyclic_x * b - cyclic_y * d),
                ];
                coeffs
            })
            .collect()
    } else {
        Vec::new()
    }
}

fn solve_mod(q: u32) -> Option<(u32, u32)> {
    for x in 0..q {
        for y in 0..q {
            if (x * x + y * y + 1) % q == 0 {
                return Some((x, y));
            }
        }
    }
    None
}

/// Creates a new LPS graph given `p` and `q`. `q` should be chosen as a prime number as it is used as the finite field and therefore dictates the number of nodes in the graph and `p` determines the degree of the graph.
/// Computation fails if these numbers are not chosen properly.
pub fn compute_lps_graph(p: u32, q: u32) -> Option<ConGraph> {
    let mut hg = ConGraph::new();
    match legendre_symbol(p as i32, q as i32) {
        // PGL
        -1 => {
            if p + 1 > ((q.pow(3) - q) - 1) {
                println!("Degree cannot exceed the group order.");
                return None;
            }
            let generators = compute_generators(p, q);
            let matrices = generate_all_pgl2(q);
            let nodes = hg.add_nodes(matrices.len());
            let mat_iter = matrices.into_iter();
            let nodes_iter = nodes.into_iter();
            let mat_to_node: HashMap<PGL2, u32> =
                HashMap::from_iter(Iterator::zip(mat_iter, nodes_iter));
            for generator in generators
                .iter()
                .map(|coeffs| PGL2::from_coeffs(coeffs.clone()).unwrap())
            {
                for (mat, node) in mat_to_node.iter() {
                    let out = *mat * generator;
                    let out_node = mat_to_node
                        .get(&out)
                        .expect("Multiplication is supposed to be closed");
                    hg.add_edge(&[*node, *out_node]);
                }
            }
            Some(hg)
        }
        // PSL
        1 => {
            if p + 1 > ((q.pow(3) - q) / 2) - 1 {
                println!("Degree cannot exceed the group order.");
                return None;
            }
            let generators = compute_generators(p, q);
            let matrices = generate_all_psl2(q);
            let nodes = hg.add_nodes(matrices.len());
            let mat_iter = matrices.into_iter();
            let node_iter = nodes.into_iter();
            let mat_to_node: HashMap<PSL2, u32> =
                HashMap::from_iter(Iterator::zip(mat_iter, node_iter));
            for generator in generators.iter().map(|coeffs| {
                PSL2::from_coeffs(*coeffs).expect("cannot convert coefficients to matrix.")
            }) {
                for (mat, node) in mat_to_node.iter() {
                    let out = *mat * generator;
                    let out_node = mat_to_node
                        .get(&out)
                        .expect("Multiplication is supposed to be closed");
                    hg.add_edge(&[*node, *out_node]);
                }
            }
            Some(hg)
        }
        _ => None,
    }
}

#[cfg(test)]
mod tests {

    use crate::{
        lps::{generate_all_pgl2, modular_inverse, prime_mod_sqrt, reduce_diophantine_solutions},
        math::finite_field::FiniteField,
    };

    use super::{
        compute_generators, compute_lps_graph, diophantine_squares_solver, generate_all_psl2,
        modular_exponent, PGL2, PSL2,
    };

    #[test]
    fn test_mod_inverse() {
        let a = 7;
        let n = 18;
        println!("mod inverse: {:?}", modular_inverse(a, n));
    }

    #[test]
    fn test_pgl2_multiplication() {
        let a = PGL2::from_coeffs([
            (1, 11).into(),
            (2, 11).into(),
            (2, 11).into(),
            (1, 11).into(),
        ])
        .unwrap();
        let b = PGL2::from_coeffs([
            (3, 11).into(),
            (2, 11).into(),
            (2, 11).into(),
            (8, 11).into(),
        ])
        .unwrap();
        println!("a = {:?}", a);
        println!("b = {:?}", b);
        println!("a * b = {:?}", a * b)
    }

    #[test]
    fn test_graph_creation_small() {
        let hg = compute_lps_graph(5, 3).expect("no graph?");
        println!("graph:\n{:}", hg);
    }

    #[test]
    fn test_diophantine_squares_solver() {
        let sols = diophantine_squares_solver(3, None);
        let num_total_sols = sols.len();
        let reduced_sols = reduce_diophantine_solutions(sols, 3);
        println!("number solutions: {:}", num_total_sols);
        println!("reduced solutions: {:?}", reduced_sols)
    }

    #[test]
    fn test_mod_exponent() {
        let n = -3;
        let p = 57;
        let x = 1;
        println!("n^x mod p = {:}", modular_exponent(n, x, p));
    }

    #[test]
    fn test_prime_mod_sqrt() {
        let p = 13;
        let a = 4;
        println!("x^2 = a mod p : {:?}", prime_mod_sqrt(a, p));
    }

    #[test]
    fn test_pgl2_from_coeffs() {
        let n = 5;
        let coeffs = [
            FiniteField(1, n),
            FiniteField(2, n),
            FiniteField(2, n),
            FiniteField(0, n),
        ];
        let mat = PGL2::from_coeffs(coeffs);
        dbg!(mat);
        println!("mat: {:}", mat.unwrap());
    }

    #[test]
    fn test_generate_all_pgl2() {
        let p = 5;
        let vertices = generate_all_pgl2(p);
        for vertex in vertices {
            println!("{:}\n", vertex);
        }
    }

    #[test]
    fn test_psl2_construction() {
        let p = 7;
        let m = PGL2::from_coeffs([
            FiniteField(1, p),
            FiniteField(2, p),
            FiniteField(3, p),
            FiniteField(4, p),
        ])
        .expect("aint no zero det.");
        let psl = PSL2::from(m);
        dbg!(psl);
    }

    #[test]
    fn test_generate_all_psl2() {
        let p = 3;
        let mats = generate_all_psl2(p);
        for mat in mats {
            println!("{:}", mat);
        }
    }

    #[test]
    fn test_get_generators() {
        let gens = compute_generators(7, 3);
        for gen in gens {
            println!("{:?}", gen);
        }
    }
}
