use std::ops::{Mul, Add};

pub trait Ring:
Mul + Add + Sized
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

pub trait Field:
Ring + Sized
{

}

pub struct Matrix<R: Ring> {
    data: Vec<R>,
    dimensions: [u16; 2],
}

