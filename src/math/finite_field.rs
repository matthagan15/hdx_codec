use core::panic;
use std::{ops::{Add, Mul, Sub, AddAssign, MulAssign, SubAssign}, fmt::Display};

use super::group_ring_field::{Ring, Group};


#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Eq, Ord, Hash)]
pub struct FiniteField(pub u32, pub u32);

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

impl Mul<FiniteField> for FiniteField {
    type Output = FiniteField;

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

impl Sub<FiniteField> for FiniteField {
    type Output = FiniteField;

    fn sub(self, rhs: FiniteField) -> Self::Output {
        if self.1 != rhs.1 {
            panic!("Subtraction among non-equal CyclicGroups.");
        }
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
        if self.1 != rhs.1 {
            panic!("Addition among non-equal CyclicGroups")
        }
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
        (value.0 as u32, value.1 as u32).into()
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

impl FiniteField {
    pub const ZERO: FiniteField = FiniteField(0, 0);

    /// Performs Euclidean remainder to get the multiplicative inverse. Panics if one cannot be found
    pub fn mul_inv(&self) -> FiniteField {
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

impl Display for FiniteField {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&format!("{:} mod {:}", self.0, self.1))
    }
}

impl SubAssign for FiniteField {
    fn sub_assign(&mut self, rhs: Self) {
        if self.1 != rhs.1 {
            panic!("Subtraction among non-equal CyclicGroups.");
        }
        let mut a = (self.0 as i32) - (rhs.0 as i32);
        while a < 0 {
            a += self.1 as i32;
        }
        self.0 = (a as u32) % self.1;
    }
}

impl AddAssign for FiniteField {
    fn add_assign(&mut self, rhs: Self) {
        if self.1 != rhs.1 {
            panic!("Addition among non-equal CyclicGroups")
        }
        self.0 = (self.0 + rhs.0) % self.1;
    }
}

impl MulAssign for FiniteField {
    fn mul_assign(&mut self, rhs: Self) {
        if self.1 != rhs.1 {
            panic!("non-equal modulus")
        }
        self.0 = (self.0 * rhs.0) % self.1;
    }
}

impl Ring for FiniteField {
    fn zero() -> Self {
        (0, u32::MAX).into()
    }

    fn one() -> Self {
        (1, u32::MAX).into()
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
        let a_inv = a.mul_inv();
        let b = a * a_inv;
        dbg!(&a_inv);
        assert_eq!(b.0, 1);
    }
}