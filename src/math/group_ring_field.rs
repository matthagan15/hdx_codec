use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};

pub trait Ring:
    Mul<Self, Output = Self> + MulAssign + Add + AddAssign + Sub + SubAssign + Copy + Sized
{
    fn zero() -> Self;
    fn one() -> Self;
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

pub struct Matrix<R: Ring> {
    data: Vec<R>,
    dimensions: [u16; 2],
}
