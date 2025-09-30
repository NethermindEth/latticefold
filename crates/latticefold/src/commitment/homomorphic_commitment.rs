use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::ops::{AddAssign, MulAssign, SubAssign};
use stark_rings::Ring;

use crate::{ark_base::*, impl_additive_ops, impl_multiplicative_ops, impl_subtractive_ops};

/// The Ajtai commitment type. Meant to contain the output of the
/// matrix-vector multiplication `A \cdot x`.
/// Since Ajtai commitment is bounded-additively homomorphic
/// one can add commitments and multiply them by a scalar.
#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize)]
pub struct Commitment<R1: Ring> {
    val: Vec<R1>,
}

impl<R: Ring> Commitment<R> {
    pub(crate) fn from_vec_raw(vec: Vec<R>) -> Self {
        Self { val: vec }
    }

    pub(crate) fn zeroed(kappa: usize) -> Self {
        Self {
            val: vec![R::zero(); kappa],
        }
    }

    pub fn len(&self) -> usize {
        self.val.len()
    }

    pub fn is_empty(&self) -> bool {
        self.val.is_empty()
    }
}

impl<'a, R: Ring> From<&'a [R]> for Commitment<R> {
    fn from(slice: &'a [R]) -> Self {
        Self { val: slice.into() }
    }
}

impl<R: Ring> From<Vec<R>> for Commitment<R> {
    fn from(vec: Vec<R>) -> Self {
        Self::from_vec_raw(vec)
    }
}

impl<R: Ring> AsRef<[R]> for Commitment<R> {
    fn as_ref(&self) -> &[R] {
        &self.val
    }
}

impl<'a, R: Ring> AddAssign<&'a Commitment<R>> for Commitment<R> {
    fn add_assign(&mut self, rhs: &'a Commitment<R>) {
        self.val
            .iter_mut()
            .zip(rhs.val.iter())
            .for_each(|(a, b)| *a += b);
    }
}

impl<'a, R: Ring> SubAssign<&'a Commitment<R>> for Commitment<R> {
    fn sub_assign(&mut self, rhs: &'a Commitment<R>) {
        self.val
            .iter_mut()
            .zip(rhs.val.iter())
            .for_each(|(a, b)| *a -= b)
    }
}

impl<'a, R: Ring> MulAssign<&'a R> for Commitment<R> {
    fn mul_assign(&mut self, rhs: &'a R) {
        self.val.iter_mut().for_each(|a| *a *= rhs)
    }
}

impl_additive_ops!(Commitment, Ring, usize);
impl_subtractive_ops!(Commitment, Ring, usize);
impl_multiplicative_ops!(Commitment, Ring, usize);
