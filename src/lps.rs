use std::{
    collections::{HashMap, HashSet},
    hash::Hash,
    ops::{Add, AddAssign, Mul}, fmt::Display,
};

use ff::Field;
use ndarray::{Array2, ShapeBuilder};
use serde::{Deserialize, Serialize};

use crate::left_right_cayley::CyclicGroup;

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Hash)]
struct GeneralSquaresSolution(i32, i32, i32, i32);

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
    fn canonical_form(&self) -> Self {
        let mut v = vec![self.0.abs(), self.1.abs(), self.2.abs(), self.3.abs()];
        v.sort();
        GeneralSquaresSolution(v[0], v[1], v[2], v[3])
    }

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

#[derive(Debug, Clone, Copy, Hash, Eq)]
struct PGL2 {
    pub coeffs: [CyclicGroup; 4],
}

impl PartialEq for PGL2 {
    fn eq(&self, other: &Self) -> bool {
        let a = self.generate_equivalence_class();
        a.contains(other)
    }
}

impl PGL2 {
    /// Returns the identity matrix [[1, 0], [0, 1]];
    pub fn identity(mod_n: u32) -> Self {
        PGL2 {
            coeffs: [
                CyclicGroup(1, mod_n),
                CyclicGroup(0, mod_n),
                CyclicGroup(0, mod_n),
                CyclicGroup(1, mod_n),
            ],
        }
    }

    /// Returns the canonical form of the member of PGL2 from the given
    /// coefficients. If the determinant is zero then this returns `None`.
    pub fn from_coeffs(coeffs: [CyclicGroup; 4]) -> Option<Self> {
        let det = coeffs[0] * &coeffs[3] - coeffs[1] * &coeffs[2];
        if det == CyclicGroup::from((0, det.1)) {
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

    fn generate_equivalence_class(&self) -> HashSet<Self> {
        let minus_identity = PGL2 {
            coeffs: [
                -1 * self.coeffs[0],
                -1 * self.coeffs[1],
                -1 * self.coeffs[2],
                -1 * self.coeffs[3],
            ],
        };
        todo!()
    }
}

impl Mul<i32> for PGL2 {
    type Output = PGL2;

    fn mul(self, rhs: i32) -> Self::Output {
        let new_coeffs: Vec<CyclicGroup> =
            self.coeffs.clone().into_iter().map(|x| x * rhs).collect();
        PGL2 {
            coeffs: [new_coeffs[0], new_coeffs[1], new_coeffs[2], new_coeffs[3]],
        }
    }
}

impl Mul<&Self> for PGL2 {
    type Output = PGL2;

    fn mul(self, rhs: &Self) -> Self::Output {
        // TODO: fill in with the matrix mul alg.
        let new: [CyclicGroup; 4] = [
            CyclicGroup::ZERO,
            CyclicGroup::ZERO,
            CyclicGroup::ZERO,
            CyclicGroup::ZERO,
        ];
        PGL2 { coeffs: new }
    }
}

impl From<[CyclicGroup; 4]> for PGL2 {
    /// Note: if an invalid matrix (i.e. a non-invertible matrix) is given 
    /// this will return an identity matrix.
    fn from(value: [CyclicGroup; 4]) -> Self {
        PGL2 { coeffs: value }
    }
}

impl Display for PGL2 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = format!("[[{:}, {:}]\n[{:}, {:}]]", self.coeffs[0], self.coeffs[1], self.coeffs[2], self.coeffs[3]);
        f.write_str(&s)
    }
}

fn generate_all_pgl2(mod_p: u32) -> Vec<PGL2> {
    let mut ret = Vec::with_capacity(mod_p.pow(4) as usize);
    // x_{1,1} is nonzero
    for b in 0..mod_p {
        for c in 0..mod_p {
            for d in 0..mod_p {
                if let Some(mat) = PGL2::from_coeffs([
                    CyclicGroup(1, mod_p),
                    CyclicGroup(b, mod_p),
                    CyclicGroup(c, mod_p),
                    CyclicGroup(d, mod_p),
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
                CyclicGroup(0, mod_p),
                CyclicGroup(b, mod_p),
                CyclicGroup(1, mod_p),
                CyclicGroup(d, mod_p),
            ]) {
                ret.push(mat);
            }
        }
    }
    ret
}

/// Returns all shuffled versions of the vec
fn shuffle(input: Vec<i32>) -> Vec<Vec<i32>> {
    if input.len() == 1 {
        return vec![input];
    }
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

/// Computes the modular inverse
/// `a * a^{-1} mod n == 1 mod n
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

fn legendre_symbol(a: i32, p: i32) -> i32 {
    let ls = modular_exponent(a, (p - 1) / 2, p);
    if ls == p - 1 {
        -1
    } else {
        ls
    }
}

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

pub trait Group: for<'a> Mul<&'a Self> + Hash {
    const ID: Self;
    fn inv(&self) -> Self;
}

pub trait AbelianGroup: for<'a> Add<&'a Self> + for<'a> AddAssign<&'a Self> + Hash {
    const ZERO: Self;
}

pub trait Ring: Group {
    const ZERO: Self;
}

impl Group for crate::left_right_cayley::CyclicGroup {
    const ID: Self = CyclicGroup::ZERO;

    fn inv(&self) -> Self {
        self.inv()
    }
}

struct CoxeterComplex {
    letters: HashSet<u32>,
    m_matrix: HashMap<(u32, u32), usize>,
}

fn generate_signs(mut input: Vec<i32>) -> Vec<Vec<i32>> {
    if input.len() == 1 {
        return vec![vec![input[0]], vec![-input[0]]];
    }
    let first_elem = input.remove(0);
    let mut remaining_elems = generate_signs(input);
    let mut copy = remaining_elems.clone();
    for mut v in remaining_elems.iter_mut() {
        v.insert(0, first_elem);
    }
    for mut v in copy.iter_mut() {
        v.insert(0, -1 * first_elem);
    }
    remaining_elems.append(&mut copy);
    remaining_elems
}

fn general_squares_solver(q: i32, limit: Option<usize>) -> Vec<GeneralSquaresSolution> {
    let range: Vec<i32> = (0..=q).collect();
    let mut ret = Vec::new();
    // Store the canonical form of visited sols.
    // let mut visited = HashSet::new();
    'outer_loop: for a_1 in range.clone().into_iter() {
        for a_2 in a_1..=q {
            for a_3 in a_2..=q {
                for a_4 in a_3..=q {
                    let sum = a_1.pow(2) + a_2.pow(2) + a_3.pow(2) + a_4.pow(2);
                    if sum - q == 0 {
                        let sol: GeneralSquaresSolution = (a_1, a_2, a_3, a_4).into();
                        let mut all_sols = sol.equivalent_forms();
                        ret.append(&mut all_sols);
                        if ret.len() == limit.map_or(1, |v| v) {
                            break 'outer_loop;
                        }
                    }
                }
            }
        }
    }
    ret
}

mod tests {
    use crate::{lps::{modular_inverse, prime_mod_sqrt, generate_all_pgl2}, left_right_cayley::CyclicGroup};

    use super::{modular_exponent, PGL2};

    #[test]
    fn test_mod_inverse() {
        let a = 6;
        let n = 18;
        println!("mod inverse: {:?}", modular_inverse(a, n));
    }

    #[test]
    fn test_mod_exponent() {
        let mut n = -3;
        let p = 57;
        let mut x = 1;
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
        let n = 53;
        let coeffs = [
            CyclicGroup(52, n),
            CyclicGroup(6, n),
            CyclicGroup(7, n),
            CyclicGroup(8, n),
        ];
        let mat = PGL2::from_coeffs(coeffs);
        dbg!(mat);
    }

    #[test]
    fn test_generate_all_pgl2() {
        let p = 3;
        let vertices = generate_all_pgl2(p);
        for vertex in vertices {
            println!("{:}\n", vertex);
        }
    }
}
