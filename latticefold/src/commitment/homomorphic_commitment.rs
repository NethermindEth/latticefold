use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{
    ops::{Add, Mul, Sub},
    Zero,
};
use lattirust_ring::Ring;

use crate::{
    commitment::CommitmentError, impl_additive_ops_from_ref, impl_multiplicative_ops_from_ref,
    impl_subtractive_ops_from_ref,
};

/// The Ajtai commitment type. Meant to contain the output of the
/// matrix-vector multiplication `A \cdot x`.
/// Enforced to have the length `C`.
/// Since Ajtai commitment is bounded-additively homomorphic
/// one can add commitments and multiply them by a scalar.
#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize)]
pub struct Commitment<const C: usize, R1: Ring> {
    val: Vec<R1>,
}

impl<const C: usize, R: Ring> Commitment<C, R> {
    pub(super) fn from_vec_raw(vec: Vec<R>) -> Self {
        Self { val: vec }
    }
}

impl<const C: usize, R: Ring> From<[R; C]> for Commitment<C, R> {
    fn from(val: [R; C]) -> Self {
        Self { val: val.into() }
    }
}

impl<'a, const C: usize, R: Ring> From<&'a [R]> for Commitment<C, R> {
    fn from(slice: &'a [R]) -> Self {
        Self { val: slice.into() }
    }
}

impl<const C: usize, R: Ring> TryFrom<Vec<R>> for Commitment<C, R> {
    type Error = CommitmentError;

    fn try_from(vec: Vec<R>) -> Result<Self, Self::Error> {
        if vec.len() != C {
            return Err(CommitmentError::WrongCommitmentLength(vec.len(), C));
        }

        Ok(Self { val: vec })
    }
}

impl<const C: usize, R: Ring> AsRef<[R]> for Commitment<C, R> {
    fn as_ref(&self) -> &[R] {
        &self.val
    }
}

impl<'a, 'b, const C: usize, R: Ring> Add<&'a Commitment<C, R>> for &'b Commitment<C, R> {
    type Output = Commitment<C, R>;

    fn add(self, rhs: &'a Commitment<C, R>) -> Self::Output {
        let mut res_vec = vec![R::zero(); C];

        res_vec
            .iter_mut()
            .zip(self.val.iter())
            .zip(rhs.val.iter())
            .for_each(|((res, &a), &b)| *res = a + b);

        Commitment::from_vec_raw(res_vec)
    }
}

impl<'a, 'b, const C: usize, R: Ring> Sub<&'a Commitment<C, R>> for &'b Commitment<C, R> {
    type Output = Commitment<C, R>;

    fn sub(self, rhs: &'a Commitment<C, R>) -> Self::Output {
        let mut res_vec = vec![R::zero(); C];

        res_vec
            .iter_mut()
            .zip(self.val.iter())
            .zip(rhs.val.iter())
            .for_each(|((res, &a), &b)| *res = a - b);

        Commitment::from_vec_raw(res_vec)
    }
}

impl<'a, 'b, const C: usize, R: Ring> Mul<&'a R> for &'b Commitment<C, R> {
    type Output = Commitment<C, R>;

    fn mul(self, rhs: &'a R) -> Self::Output {
        let mut res_vec = vec![R::zero(); C];

        res_vec
            .iter_mut()
            .zip(self.val.iter())
            .for_each(|(res, &a)| *res = a * rhs);

        Commitment::from_vec_raw(res_vec)
    }
}

impl_additive_ops_from_ref!(Commitment, Ring, usize);
impl_subtractive_ops_from_ref!(Commitment, Ring, usize);
impl_multiplicative_ops_from_ref!(Commitment, Ring, usize);

impl<const C: usize, R: Ring> Zero for Commitment<C, R> {
    fn zero() -> Self {
        Self::from([R::zero(); C])
    }

    fn is_zero(&self) -> bool {
        self.val == [R::zero(); C]
    }
}
