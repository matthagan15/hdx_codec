use std::{ops::{Add, Mul, Sub, AddAssign, MulAssign}, fmt::Display};


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
