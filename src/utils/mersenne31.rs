use std::fmt::Display;
use std::ops::{
    Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Rem, RemAssign, Sub, SubAssign,
};

use super::{ComplexConjugate, FieldExpOps};
use crate::impl_field;

pub const MODULUS_BITS: u32 = 31;
pub const N_BYTES_FELT: usize = 4;
pub const P: u32 = 2147483647; // 2 ** 31 -1

#[repr(transparent)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct M31(pub u32);
pub type BaseField = M31;

impl_field!(M31, P);

impl M31 {
    pub fn sqrt(&self) -> Option<Self> {
        let result = self.pow(1 << 29);
        (result.square() == *self).then_some(result)
    }

    /// Assumes that `val` is in the range [0, 2 * `P`) and returns `val` % `P`.
    pub fn partial_reduce(val: u32) -> Self {
        Self(val.checked_sub(P).unwrap_or(val))
    }

    /// Assumes that `val` is in the range [0, `P`.pow(2)) and returns `val` % `P`.
    pub fn reduce(val: u64) -> Self {
        Self((((((val >> MODULUS_BITS) + val + 1) >> MODULUS_BITS) + val) & (P as u64)) as u32)
    }

    pub const fn from_u32_unchecked(arg: u32) -> Self {
        Self(arg)
    }
}

impl Display for M31 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl Add for M31 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::partial_reduce(self.0 + rhs.0)
    }
}

impl Neg for M31 {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::partial_reduce(P - self.0)
    }
}

impl Sub for M31 {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::partial_reduce(self.0 + P - rhs.0)
    }
}

impl Mul for M31 {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self::reduce((self.0 as u64) * (rhs.0 as u64))
    }
}

impl ComplexConjugate for M31 {}

impl One for M31 {
    fn one() -> Self {
        Self(1)
    }
}

impl Zero for M31 {
    fn zero() -> Self {
        Self(0)
    }

    fn is_zero(&self) -> bool {
        *self == Self::zero()
    }
}

impl From<usize> for M31 {
    fn from(value: usize) -> Self {
        M31::reduce(value.try_into().unwrap())
    }
}

impl From<u32> for M31 {
    fn from(value: u32) -> Self {
        M31::reduce(value.into())
    }
}

impl From<i32> for M31 {
    fn from(value: i32) -> Self {
        M31::reduce(value.try_into().unwrap())
    }
}
