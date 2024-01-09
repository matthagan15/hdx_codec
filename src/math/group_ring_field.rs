use std::ops::{Mul, Add, AddAssign};

pub trait Ring:
Mul<Self, Output = Self> + Add + AddAssign + Copy + Sized 
{
    fn zero() -> Self;
    fn one() -> Self;
    fn additive_inv(&self) -> Self;
}

pub trait Group:
Mul + Sized
{
    fn one() -> Self;
    fn inv(&self) -> Self;
}

pub trait AbelianGroup:
Add + Sized
{
    fn zero() -> Self;
}

pub trait Module:
AbelianGroup + Sized
{
    type Scalar: Ring;
}

pub trait EuclideanDomain {
    fn eul_division(&self, rhs: &Self) -> Self;
}

pub trait Field:
Ring + std::fmt::Debug + std::cmp::PartialEq + Sized
{
    fn mul_inv(&self) -> Self;
    // fn division();
    fn zero() -> Self;
    fn one() -> Self;
}

pub struct Matrix<R: Ring> {
    data: Vec<R>,
    dimensions: [u16; 2],
}

