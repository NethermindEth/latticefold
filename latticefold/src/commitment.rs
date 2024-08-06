use std::{
    marker::PhantomData,
    ops::{Add, Mul, Sub},
};

use lattirust_arithmetic::balanced_decomposition::decompose_balanced_slice_polyring;
use lattirust_arithmetic::ring::{PolyRing, Ring};
use lattirust_arithmetic::{
    challenge_set::latticefold_challenge_set::OverField,
    linear_algebra::{Matrix, Vector},
};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum CommitmentError {
    #[error("Wrong length of the witness: {0}, expected: {1}")]
    WrongWitnessLength(usize, usize),
}

/// Ajtai commitment parameters.
/// Convenient to enforce them compile-time.
pub trait AjtaiParams: Clone {
    /// The MSIS bound.
    const B: u128;
    /// The ring modulus should be < B^L.
    const L: usize;
    /// The number of rows of the Ajtai matrix.
    const WITNESS_SIZE: usize;
    /// The number of columns of the Ajtai matrix.
    const OUTPUT_SIZE: usize;
}

/// The Ajtai commitment type. Meant to contain the output of the
/// matrix-vector multiplication `A \cdot x`.
/// Enforced to have the length `P::OUTPUT_SIZE`.
/// Since Ajtai commitment is bounded-additively homomorphic
/// one can add commitments and multiply them by a scalar.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Commitment<R: Ring, P: AjtaiParams> {
    _phantom: PhantomData<P>,
    val: Vec<R>,
}

impl<R: Ring, P: AjtaiParams> Commitment<R, P> {
    fn from_vec_raw(vec: Vec<R>) -> Self {
        Self {
            _phantom: PhantomData,
            val: vec,
        }
    }
}

impl<'a, R: Ring, P: AjtaiParams> TryFrom<&'a [R]> for Commitment<R, P> {
    type Error = CommitmentError;

    fn try_from(slice: &'a [R]) -> Result<Self, Self::Error> {
        if slice.len() != P::WITNESS_SIZE {
            return Err(CommitmentError::WrongWitnessLength(
                slice.len(),
                P::WITNESS_SIZE,
            ));
        }

        Ok(Self {
            _phantom: PhantomData,
            val: Vec::from(slice),
        })
    }
}

impl<R: Ring, P: AjtaiParams> TryFrom<Vec<R>> for Commitment<R, P> {
    type Error = CommitmentError;

    fn try_from(vec: Vec<R>) -> Result<Self, Self::Error> {
        if vec.len() != P::WITNESS_SIZE {
            return Err(CommitmentError::WrongWitnessLength(
                vec.len(),
                P::WITNESS_SIZE,
            ));
        }

        Ok(Self {
            _phantom: PhantomData,
            val: vec,
        })
    }
}

impl<'a, R: Ring, P: AjtaiParams> Add<&'a Commitment<R, P>> for &'a Commitment<R, P> {
    type Output = Commitment<R, P>;

    fn add(self, rhs: &'a Commitment<R, P>) -> Self::Output {
        let mut res_vec = Vec::<R>::with_capacity(P::OUTPUT_SIZE);

        res_vec
            .iter_mut()
            .zip(self.val.iter())
            .zip(rhs.val.iter())
            .for_each(|((res, &a), &b)| *res = a + b);

        Commitment::from_vec_raw(res_vec)
    }
}

impl<'a, R: Ring, P: AjtaiParams> Add<Commitment<R, P>> for &'a Commitment<R, P> {
    type Output = Commitment<R, P>;

    fn add(self, rhs: Commitment<R, P>) -> Self::Output {
        let mut res_vec = Vec::<R>::with_capacity(P::OUTPUT_SIZE);

        res_vec
            .iter_mut()
            .zip(self.val.iter())
            .zip(rhs.val.iter())
            .for_each(|((res, &a), &b)| *res = a + b);

        Commitment::from_vec_raw(res_vec)
    }
}

impl<'a, R: Ring, P: AjtaiParams> Add<&'a Commitment<R, P>> for Commitment<R, P> {
    type Output = Commitment<R, P>;

    fn add(self, rhs: &'a Commitment<R, P>) -> Self::Output {
        let mut res_vec = Vec::<R>::with_capacity(P::OUTPUT_SIZE);

        res_vec
            .iter_mut()
            .zip(self.val.iter())
            .zip(rhs.val.iter())
            .for_each(|((res, &a), &b)| *res = a + b);

        Commitment::from_vec_raw(res_vec)
    }
}

impl<R: Ring, P: AjtaiParams> Add<Commitment<R, P>> for Commitment<R, P> {
    type Output = Commitment<R, P>;

    fn add(self, rhs: Commitment<R, P>) -> Self::Output {
        let mut res_vec = Vec::<R>::with_capacity(P::OUTPUT_SIZE);

        res_vec
            .iter_mut()
            .zip(self.val.iter())
            .zip(rhs.val.iter())
            .for_each(|((res, &a), &b)| *res = a + b);

        Commitment::from_vec_raw(res_vec)
    }
}

impl<'a, R: Ring, P: AjtaiParams> Sub<&'a Commitment<R, P>> for &'a Commitment<R, P> {
    type Output = Commitment<R, P>;

    fn sub(self, rhs: &'a Commitment<R, P>) -> Self::Output {
        let mut res_vec = Vec::<R>::with_capacity(P::OUTPUT_SIZE);

        res_vec
            .iter_mut()
            .zip(self.val.iter())
            .zip(rhs.val.iter())
            .for_each(|((res, &a), &b)| *res = a - b);

        Commitment::from_vec_raw(res_vec)
    }
}

impl<'a, R: Ring, P: AjtaiParams> Sub<Commitment<R, P>> for &'a Commitment<R, P> {
    type Output = Commitment<R, P>;

    fn sub(self, rhs: Commitment<R, P>) -> Self::Output {
        let mut res_vec = Vec::<R>::with_capacity(P::OUTPUT_SIZE);

        res_vec
            .iter_mut()
            .zip(self.val.iter())
            .zip(rhs.val.iter())
            .for_each(|((res, &a), &b)| *res = a - b);

        Commitment::from_vec_raw(res_vec)
    }
}

impl<'a, R: Ring, P: AjtaiParams> Sub<&'a Commitment<R, P>> for Commitment<R, P> {
    type Output = Commitment<R, P>;

    fn sub(self, rhs: &'a Commitment<R, P>) -> Self::Output {
        let mut res_vec = Vec::<R>::with_capacity(P::OUTPUT_SIZE);

        res_vec
            .iter_mut()
            .zip(self.val.iter())
            .zip(rhs.val.iter())
            .for_each(|((res, &a), &b)| *res = a - b);

        Commitment::from_vec_raw(res_vec)
    }
}

impl<R: Ring, P: AjtaiParams> Sub<Commitment<R, P>> for Commitment<R, P> {
    type Output = Commitment<R, P>;

    fn sub(self, rhs: Commitment<R, P>) -> Self::Output {
        let mut res_vec = Vec::<R>::with_capacity(P::OUTPUT_SIZE);

        res_vec
            .iter_mut()
            .zip(self.val.iter())
            .zip(rhs.val.iter())
            .for_each(|((res, &a), &b)| *res = a - b);

        Commitment::from_vec_raw(res_vec)
    }
}

impl<'a, 'b, R: Ring, P: AjtaiParams> Mul<&'b R> for &'a Commitment<R, P> {
    type Output = Commitment<R, P>;

    fn mul(self, rhs: &'b R) -> Self::Output {
        let mut res_vec = Vec::<R>::with_capacity(P::OUTPUT_SIZE);

        res_vec
            .iter_mut()
            .zip(self.val.iter())
            .for_each(|(res, &a)| *res = a * rhs);

        Commitment::from_vec_raw(res_vec)
    }
}

impl<'b, R: Ring, P: AjtaiParams> Mul<&'b R> for Commitment<R, P> {
    type Output = Commitment<R, P>;

    fn mul(self, rhs: &'b R) -> Self::Output {
        let mut res_vec = Vec::<R>::with_capacity(P::OUTPUT_SIZE);

        res_vec
            .iter_mut()
            .zip(self.val.iter())
            .for_each(|(res, &a)| *res = a * rhs);

        Commitment::from_vec_raw(res_vec)
    }
}

impl<R: Ring, P: AjtaiParams> Mul<R> for Commitment<R, P> {
    type Output = Commitment<R, P>;

    fn mul(self, rhs: R) -> Self::Output {
        let mut res_vec = Vec::<R>::with_capacity(P::OUTPUT_SIZE);

        res_vec
            .iter_mut()
            .zip(self.val.iter())
            .for_each(|(res, &a)| *res = a * rhs);

        Commitment::from_vec_raw(res_vec)
    }
}

impl<'a, R: Ring, P: AjtaiParams> Mul<R> for &'a Commitment<R, P> {
    type Output = Commitment<R, P>;

    fn mul(self, rhs: R) -> Self::Output {
        let mut res_vec = Vec::<R>::with_capacity(P::OUTPUT_SIZE);

        res_vec
            .iter_mut()
            .zip(self.val.iter())
            .for_each(|(res, &a)| *res = a * rhs);

        Commitment::from_vec_raw(res_vec)
    }
}

// TODO: use macros to implement the other operations

/// A concrete instantiation of the Ajtai commitment scheme.
/// Contains a random Ajtai matrix for the coresssponding Ajtai parameters
/// `CR` is the type parameter for the coefficient representation of the ring
/// `NTT` is the NTT representation of the same ring.
#[derive(Clone, Debug)]
pub struct AjtaiCommitmentScheme<CR: PolyRing, NTT: OverField, P: AjtaiParams>
where
    NTT: Into<CR> + From<CR>,
{
    _cr: PhantomData<CR>,
    _p: PhantomData<P>,
    matrix: Matrix<NTT>,
}

impl<CR: PolyRing, NTT: OverField, P: AjtaiParams> AjtaiCommitmentScheme<CR, NTT, P>
where
    NTT: Into<CR> + From<CR>,
{
    pub fn rand<Rng: rand::Rng + ?Sized>(rng: &mut Rng) -> Self {
        Self {
            _cr: PhantomData,
            _p: PhantomData,
            matrix: Matrix::rand(P::WITNESS_SIZE, P::OUTPUT_SIZE, rng),
        }
    }
}

impl<CR: PolyRing, NTT: OverField, P: AjtaiParams> AjtaiCommitmentScheme<CR, NTT, P>
where
    NTT: Into<CR> + From<CR>,
{
    pub fn commit_ntt(&self, f: &[NTT]) -> Result<Commitment<NTT, P>, CommitmentError> {
        // TODO: a lot of clones and copies. Can we optimise this somehow?
        if f.len() != P::WITNESS_SIZE {
            return Err(CommitmentError::WrongWitnessLength(
                f.len(),
                P::WITNESS_SIZE,
            ));
        }

        let commitment_vec = self.matrix.clone() * Vector::from(Vec::from(f));

        Commitment::try_from(commitment_vec.iter().copied().collect::<Vec<_>>())
    }

    pub fn commit_coeff(&self, f: &[CR]) -> Result<Commitment<NTT, P>, CommitmentError> {
        if f.len() != P::WITNESS_SIZE {
            return Err(CommitmentError::WrongWitnessLength(
                f.len(),
                P::WITNESS_SIZE,
            ));
        }

        self.commit_ntt(&f.iter().map(|&x| x.into()).collect::<Vec<NTT>>())
    }

    pub fn decompose_and_commit_coeff(
        &self,
        f: &[CR],
    ) -> Result<Commitment<NTT, P>, CommitmentError> {
        let f = decompose_balanced_slice_polyring(f, P::B, Some(P::L))
            .into_iter()
            .flatten()
            .collect::<Vec<_>>();

        self.commit_coeff(&f)
    }

    /// Commits the gadgeted witness, i.e. w = G_B f, for some f.
    pub fn decompose_and_commit_ntt(
        &self,
        w: &[NTT],
    ) -> Result<Commitment<NTT, P>, CommitmentError> {
        let f: Vec<NTT> = decompose_balanced_slice_polyring(
            &w.iter().map(|&x| x.into()).collect::<Vec<CR>>(),
            P::B,
            Some(P::L),
        )
        .iter()
        .flatten()
        .map(|&x| x.into())
        .collect();

        self.commit_ntt(&f)
    }
}
