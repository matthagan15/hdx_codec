use core::panic;
use std::{ops::{Add, Mul, Sub, AddAssign, MulAssign, SubAssign}, fmt::Display};

use super::group_ring_field::{Ring, Group};


#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Eq, Ord, Hash)]
pub struct CyclicGroup(pub u64, pub u64);

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
        CyclicGroup::from((a as u64, rhs.1))
    }
}

impl Mul<i32> for CyclicGroup {
    type Output = CyclicGroup;

    fn mul(self, rhs: i32) -> Self::Output {
        let mut a = rhs * (self.0 as i32);
        while a < 0 {
            a += self.1 as i32;
        }
        CyclicGroup::from((a as u64 % self.1, self.1))
    }
}

impl crate::math::group_ring_field::Group for CyclicGroup {
    fn id() -> Self {
        CyclicGroup(1, 0)
    }
    fn inv(&self) -> Self {
        CyclicGroup::from((self.0 as i32 * -1, self.1))
    }
}

impl Add<i32> for CyclicGroup {
    type Output = CyclicGroup;
    fn add(self, rhs: i32) -> Self::Output {
        let mut rhs_mod_n = rhs;
        while rhs_mod_n < 0 {
            rhs_mod_n += self.1 as i32;
        }
        let a = (rhs_mod_n as u64) + self.0;
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
        CyclicGroup::from(((a as u64) % self.1, self.1))
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

impl From<(u64, u64)> for CyclicGroup {
    fn from(value: (u64, u64)) -> Self {
        CyclicGroup(value.0 % value.1, value.1)
    }
}

impl From<(i32, u64)> for CyclicGroup {
    fn from(value: (i32, u64)) -> Self {
        let mut a = value.0;
        while a < 0 {
            a += value.1 as i32;
        }
        CyclicGroup((a as u64) % value.1, value.1)
    }
}

impl From<(i32, i32)> for CyclicGroup {
    fn from(value: (i32, i32)) -> Self {
        (value.0 as u64, value.1 as u64).into()
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

    /// Performs Euclidean remainder to get the multiplicative inverse. Panics if one cannot be found
    pub fn mul_inv(&self) -> CyclicGroup {
        let mut t = 0_i32;
        let mut r = self.1 as i32;
        let mut new_t = 1_i32;
        let mut new_r = self.0 as i32;
        while new_r != 0 {
            if let Some(q) = r.checked_div(new_r) {
                (t, new_t) = (new_t, t - q * new_t);
                (r, new_r) = (new_r, r - q * new_r);
            } else {
                panic!("Could not find modular inverse.")
            }
        }
        if r > 1 {
            panic!("Could not find modular inverse.")
        } else {
            while t < 0 {
                t += self.1 as i32;
            }
            (t, self.1).into()
        }
    }
}

impl Display for CyclicGroup {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&format!("{:} mod {:}", self.0, self.1))
    }
}

impl SubAssign for CyclicGroup {
    fn sub_assign(&mut self, rhs: Self) {
        if self.1 != rhs.1 {
            panic!("Subtraction among non-equal CyclicGroups.");
        }
        let mut a = (self.0 as i32) - (rhs.0 as i32);
        while a < 0 {
            a += self.1 as i32;
        }
        self.0 = (a as u64) % self.1;
    }
}

impl AddAssign for CyclicGroup {
    fn add_assign(&mut self, rhs: Self) {
        if self.1 != rhs.1 {
            panic!("Addition among non-equal CyclicGroups")
        }
        self.0 = (self.0 + rhs.0) % self.1;
    }
}

impl MulAssign for CyclicGroup {
    fn mul_assign(&mut self, rhs: Self) {
        if self.1 != rhs.1 {
            panic!("non-equal modulus")
        }
        self.0 = (self.0 * rhs.0) % self.1;
    }
}

impl Ring for CyclicGroup {
    fn zero() -> Self {
        (0, u64::MAX).into()
    }

    fn one() -> Self {
        (1, u64::MAX).into()
    }

    fn additive_inv(&self) -> Self {
        self.inv()
    }
}

mod tests {
    use super::CyclicGroup;

    #[test]
    fn test_modular_inverse() {
        let a = CyclicGroup(23, 199);
        let a_inv = a.mul_inv();
        let b = a * a_inv;
        dbg!(&a_inv);
        assert_eq!(b.0, 1);
    }
}