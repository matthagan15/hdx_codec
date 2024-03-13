use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};

pub trait Ring:
    Mul<Self, Output = Self> + MulAssign + Add + AddAssign + Sub + SubAssign + Clone + Sized
{
    fn zero(&self) -> Self;
    fn one(&self) -> Self;
    fn additive_inv(&self) -> Self;
}

pub trait Group: Mul + Sized {
    fn id() -> Self;
    fn inv(&self) -> Self;
}

pub trait AbelianGroup: Add + Sized {
    fn zero() -> Self;
}

pub trait Module: AbelianGroup + Sized {
    type Scalar: Ring;
}

pub trait Field: Ring + std::fmt::Debug + std::cmp::PartialEq + Sized {
    fn mul_inv(&self) -> Self;
    // fn division();
    fn zero() -> Self;
    fn one() -> Self;
}
