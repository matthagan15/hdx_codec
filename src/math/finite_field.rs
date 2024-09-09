use core::panic;
use std::{
    fmt::Display,
    hash::Hash,
    ops::{Add, AddAssign, Div, Mul, MulAssign, Sub, SubAssign},
};

use rand::prelude::*;
use serde::{Deserialize, Serialize};

use super::group_ring_field::{Group, Ring};

pub fn random_message(message_len: usize, field_mod: u32) -> Vec<FiniteField> {
    let mut ret = Vec::with_capacity(message_len);
    for _ in 0..message_len {
        ret.push(0_u32);
    }
    let mut rng = thread_rng();
    rng.fill(&mut ret[..]);
    ret.into_iter()
        .map(|x| FiniteField::new(x, field_mod))
        .collect()
}

pub fn zero_vec(length: usize, field_mod: u32) -> Vec<FiniteField> {
    (0..length)
        .map(|_| FiniteField::new(0, field_mod))
        .collect()
}

pub type FFRep = u32;
/// A member of a finite field extension, where the elements are:
/// (value, prime_base, prime_power). So 4 mod 3^2 would be (4, 3, 9), if you want to get the power idk I dont want to store too many numbers but I also want to be able to do multiplications faster (don't want to compute pow every time)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct FiniteFieldExt(pub FFRep, FFRep, FFRep);

impl FiniteFieldExt {
    pub fn new(value: FFRep, prime_base: FFRep, prime_power: FFRep) -> Self {
        // todo: technically should do a `.checked_pow(prime_power).unwrap()`
        // but if these objects are getting made a lot we don't want that
        // overhead.
        FiniteFieldExt(
            value % prime_power.pow(prime_power),
            prime_power,
            prime_power.pow(prime_power),
        )
    }
    pub fn get_basis(degree: u32, prime_base: FFRep, prime_power: FFRep) -> Self {
        let a = prime_base.pow(degree);
        FiniteFieldExt::new(a, prime_base, prime_power)
    }
}

/// THIS IS WRONG! Multiplication in F_{p^k} is not the
/// same as multiplication in Z / p^k Z.
impl Mul<&FiniteFieldExt> for FiniteFieldExt {
    type Output = Self;

    fn mul(self, rhs: &FiniteFieldExt) -> Self::Output {
        FiniteFieldExt((self.0 * rhs.0) % self.2, self.1, self.2)
    }
}

impl Add<FFRep> for FiniteFieldExt {
    type Output = FiniteFieldExt;

    fn add(self, rhs: FFRep) -> Self::Output {
        FiniteFieldExt::new(self.0 + rhs, self.1, self.2)
    }
}

impl AddAssign<FFRep> for FiniteFieldExt {
    fn add_assign(&mut self, rhs: FFRep) {
        let a = (self.0 + rhs) % self.2;
        self.0 = a;
    }
}

impl Mul<i32> for &FiniteFieldExt {
    type Output = FiniteFieldExt;

    fn mul(self, rhs: i32) -> Self::Output {
        let mut a = (self.0 as i32) * rhs;
        while a < 0 {
            a += self.2 as i32;
        }
        FiniteFieldExt(a as FFRep, self.1, self.2)
    }
}

impl Ord for FiniteFieldExt {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.cmp(&other.0)
    }
}

impl PartialOrd for FiniteFieldExt {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

impl Hash for FiniteFieldExt {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.0.hash(state);
    }
}

impl Add<&FiniteFieldExt> for FiniteFieldExt {
    type Output = Self;

    fn add(self, rhs: &FiniteFieldExt) -> Self::Output {
        let a = (self.0 + rhs.0) % self.2;
        FiniteFieldExt(a, self.1, self.2)
    }
}

impl Mul<&FiniteFieldExt> for i32 {
    type Output = FiniteFieldExt;

    fn mul(self, rhs: &FiniteFieldExt) -> Self::Output {
        let mut a = self * (rhs.0 as i32);
        while a < 0 {
            a += rhs.2 as i32;
        }
        FiniteFieldExt::new(a as u32, rhs.1, rhs.2)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct FiniteField(pub u32, pub u32);

impl Ord for FiniteField {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.cmp(&other.0)
    }
}

impl PartialOrd for FiniteField {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

// TODO: Eliminating these checks could introduce bugs but might be a lot faster.
impl Add<&FiniteField> for FiniteField {
    type Output = FiniteField;

    fn add(self, rhs: &FiniteField) -> Self::Output {
        // if self.1 != rhs.1 {
        //     panic!("[FiniteField] addition not defined for different fields.");
        // }
        FiniteField((self.0 + rhs.0) % self.1, self.1)
    }
}

impl Mul<FiniteField> for FiniteField {
    type Output = FiniteField;

    #[inline(always)]
    fn mul(self, rhs: FiniteField) -> Self::Output {
        self * &rhs
    }
}

impl Mul<FiniteField> for i32 {
    type Output = FiniteField;

    fn mul(self, rhs: FiniteField) -> Self::Output {
        let mut a = (rhs.0 as i32) * self;
        while a < 0 {
            a += rhs.1 as i32;
        }
        FiniteField::from((a as u32, rhs.1))
    }
}

impl Mul<i32> for FiniteField {
    type Output = FiniteField;

    fn mul(self, rhs: i32) -> Self::Output {
        let mut a = rhs * (self.0 as i32);
        while a < 0 {
            a += self.1 as i32;
        }
        FiniteField::from((a as u32 % self.1, self.1))
    }
}
impl Mul<u32> for &FiniteField {
    type Output = FiniteField;

    fn mul(self, rhs: u32) -> Self::Output {
        let a = rhs * self.0;
        FiniteField(a % self.1, self.1)
    }
}

impl crate::math::group_ring_field::Group for FiniteField {
    fn id() -> Self {
        FiniteField(1, 0)
    }
    fn inv(&self) -> Self {
        FiniteField::from((self.0 as i32 * -1, self.1))
    }
}

impl Add<i32> for FiniteField {
    type Output = FiniteField;
    fn add(self, rhs: i32) -> Self::Output {
        let mut rhs_mod_n = rhs;
        while rhs_mod_n < 0 {
            rhs_mod_n += self.1 as i32;
        }
        let a = (rhs_mod_n as u32) + self.0;
        FiniteField::from((a % self.1, self.1))
    }
}

impl Add<u32> for &FiniteField {
    type Output = FiniteField;
    fn add(self, rhs: u32) -> Self::Output {
        FiniteField((rhs + self.0) % self.1, self.1)
    }
}

impl Sub<FiniteField> for FiniteField {
    type Output = FiniteField;

    fn sub(self, rhs: FiniteField) -> Self::Output {
        // if self.1 != rhs.1 {
        //     panic!("Subtraction among non-equal CyclicGroups.");
        // }
        let mut a = (self.0 as i32) - (rhs.0 as i32);
        while a < 0 {
            a += self.1 as i32;
        }
        FiniteField::from(((a as u32) % self.1, self.1))
    }
}
impl Add<FiniteField> for FiniteField {
    type Output = Self;

    fn add(self, rhs: FiniteField) -> Self::Output {
        // if self.1 != rhs.1 {
        //     panic!("Addition among non-equal CyclicGroups")
        // }
        FiniteField((self.0 + rhs.0) % self.1, self.1)
    }
}

impl From<(u32, u32)> for FiniteField {
    fn from(value: (u32, u32)) -> Self {
        FiniteField(value.0 % value.1, value.1)
    }
}

impl From<(i32, u32)> for FiniteField {
    fn from(value: (i32, u32)) -> Self {
        let mut a = value.0;
        while a < 0 {
            a += value.1 as i32;
        }
        FiniteField((a as u32) % value.1, value.1)
    }
}

impl From<(i32, i32)> for FiniteField {
    fn from(value: (i32, i32)) -> Self {
        let p = value.1.abs() as u32;
        let mut x = value.0;
        while x < 0 {
            x += p as i32;
        }
        FiniteField(x as u32, p)
    }
}

impl Add<&FiniteField> for &FiniteField {
    type Output = FiniteField;

    fn add(self, rhs: &FiniteField) -> Self::Output {
        // if self.1 != rhs.1 {
        //     panic!("[FiniteField] addition not defined for different fields.");
        // }
        FiniteField((self.0 + rhs.0) % self.1, self.1)
    }
}

impl AddAssign<&FiniteField> for FiniteField {
    fn add_assign(&mut self, rhs: &FiniteField) {
        // if self.1 != rhs.1 {
        //     panic!("[FiniteField] Addition not defined for different fields.")
        // }
        self.0 = (self.0 + rhs.0) % self.1;
    }
}

impl Mul<&FiniteField> for FiniteField {
    type Output = FiniteField;

    fn mul(self, rhs: &FiniteField) -> Self::Output {
        // if self.1 != rhs.1 {
        //     panic!("[FiniteField] Multiplication not defined for different fields.")
        // }
        FiniteField((self.0 * rhs.0) % self.1, self.1)
    }
}

impl MulAssign<&FiniteField> for FiniteField {
    fn mul_assign(&mut self, rhs: &FiniteField) {
        // if self.1 != rhs.1 {
        //     panic!("[FiniteField] Multiplication not defined for different fields.")
        // }
        self.0 = (self.0 * rhs.0) % self.1;
    }
}

impl FiniteField {
    pub const ZERO: FiniteField = FiniteField(0, 0);

    pub fn pow(&self, exponent: u32) -> Self {
        if exponent == 0 {
            (1, self.1).into()
        } else if exponent == 1 {
            self.clone()
        } else {
            let mut acc = 2;
            let mut out = self.clone() * self.clone();
            while acc * acc <= exponent {
                out = out * out;
                acc *= acc;
            }
            let leftover = self.pow(exponent - acc);
            leftover * out
        }
    }

    /// gives you element mod characterstic
    pub fn new(element: u32, characteristic: u32) -> Self {
        FiniteField(element % characteristic, characteristic)
    }

    /// Performs Euclidean remainder to get the multiplicative inverse. Panics if one cannot be found.
    pub fn modular_inverse(&self) -> FiniteField {
        let mut t = 0 as i64;
        let mut r = self.1 as i64;
        let mut new_t = 1 as i64;
        let mut new_r = self.0 as i64;
        while new_r != 0 {
            if let Some(q) = r.checked_div(new_r) {
                (t, new_t) = (new_t, t - q * new_t);
                (r, new_r) = (new_r, r - q * new_r);
            } else {
                panic!("Could not find modular inverse.")
            }
        }
        if r > 1 {
            println!("self = {:}", self);
            println!("r = {:}", r);
            println!("t = {:}", t);
            println!("new_r = {:}", new_r);
            println!("new_t = {:}", new_t);
            panic!("Could not find modular inverse.")
        } else {
            while t < 0 {
                t += self.1 as i64;
            }
            (t as u32, self.1).into()
        }
    }
}

impl Display for FiniteField {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&format!("{:} mod {:}", self.0, self.1))
    }
}

impl SubAssign for FiniteField {
    fn sub_assign(&mut self, rhs: Self) {
        // if self.1 != rhs.1 {
        //     panic!("Subtraction among non-equal CyclicGroups.");
        // }
        let mut a = (self.0 as i32) - (rhs.0 as i32);
        while a < 0 {
            a += self.1 as i32;
        }
        self.0 = (a as u32) % self.1;
    }
}

impl SubAssign<&FiniteField> for FiniteField {
    fn sub_assign(&mut self, rhs: &FiniteField) {
        // if self.1 != rhs.1 {
        //     panic!("Subtraction among non-equal CyclicGroups.");
        // }
        let mut a = (self.0 as i32) - (rhs.0 as i32);
        while a < 0 {
            a += self.1 as i32;
        }
        self.0 = (a as u32) % self.1;
    }
}

impl AddAssign for FiniteField {
    fn add_assign(&mut self, rhs: Self) {
        // if self.1 != rhs.1 {
        //     panic!("Addition among non-equal CyclicGroups")
        // }
        self.0 = (self.0 + rhs.0) % self.1;
    }
}

impl MulAssign for FiniteField {
    fn mul_assign(&mut self, rhs: Self) {
        // if self.1 != rhs.1 {
        //     panic!("non-equal modulus")
        // }
        self.0 = (self.0 * rhs.0) % self.1;
    }
}

impl Ring for FiniteField {
    fn zero(&self) -> Self {
        (0, self.1).into()
    }

    fn one(&self) -> Self {
        (1, self.1).into()
    }

    fn additive_inv(&self) -> Self {
        self.inv()
    }
}

mod tests {
    use super::FiniteField;

    #[test]
    fn test_modular_inverse() {
        let a = FiniteField(23, 199);
        let a_inv = a.modular_inverse();
        let b = a * a_inv;
        dbg!(&a_inv);
        assert_eq!(b.0, 1);
    }

    #[test]
    fn test_pow() {
        let a = FiniteField::from((2, 199));
        for k in 0..10 {
            println!("2 ^ {:} % 199 = {:}", k, a.pow(k));
        }
    }

    #[test]
    fn test_multiplication() {
        let f = FiniteField::new(7, 9);
        println!("{:}", -1 * f);
    }
}
