use std::{
    marker::PhantomData,
    ops::{Add, Mul, Sub},
};

use lattirust_arithmetic::ring::Ring;
use lattirust_arithmetic::{
    balanced_decomposition::decompose_balanced_vec,
    challenge_set::latticefold_challenge_set::OverField,
    linear_algebra::{Matrix, Vector},
    ring::ConvertibleRing,
};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum CommitmentError {
    #[error("Wrong length of the witness {0}")]
    WrongWitnessLength(usize),
}

pub trait AjtaiParams<R: Ring>: Clone {
    // The MSIS bound.
    const B: u128;
    // The ring modulus should be < B^L.
    const L: usize;
    const WITNESS_SIZE: usize;
    const OUTPUT_SIZE: usize;
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Commitment<R: Ring, P: AjtaiParams<R>> {
    _phantom: PhantomData<P>,
    val: Vec<R>,
}

impl<R: Ring, P: AjtaiParams<R>> Commitment<R, P> {
    fn from_vec_raw(vec: Vec<R>) -> Self {
        Self {
            _phantom: PhantomData,
            val: vec,
        }
    }
}

impl<'a, R: Ring, P: AjtaiParams<R>> TryFrom<&'a [R]> for Commitment<R, P> {
    type Error = CommitmentError;

    fn try_from(slice: &'a [R]) -> Result<Self, Self::Error> {
        if slice.len() != P::WITNESS_SIZE {
            return Err(CommitmentError::WrongWitnessLength(slice.len()));
        }

        Ok(Self {
            _phantom: PhantomData,
            val: Vec::from(slice),
        })
    }
}

impl<R: Ring, P: AjtaiParams<R>> TryFrom<Vec<R>> for Commitment<R, P> {
    type Error = CommitmentError;

    fn try_from(vec: Vec<R>) -> Result<Self, Self::Error> {
        if vec.len() != P::WITNESS_SIZE {
            return Err(CommitmentError::WrongWitnessLength(vec.len()));
        }

        Ok(Self {
            _phantom: PhantomData,
            val: vec,
        })
    }
}

impl<'a, R: Ring, P: AjtaiParams<R>> Add<&'a Commitment<R, P>> for &'a Commitment<R, P> {
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

impl<'a, R: Ring, P: AjtaiParams<R>> Add<Commitment<R, P>> for &'a Commitment<R, P> {
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

impl<'a, R: Ring, P: AjtaiParams<R>> Add<&'a Commitment<R, P>> for Commitment<R, P> {
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

impl<R: Ring, P: AjtaiParams<R>> Add<Commitment<R, P>> for Commitment<R, P> {
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

impl<'a, R: Ring, P: AjtaiParams<R>> Sub<&'a Commitment<R, P>> for &'a Commitment<R, P> {
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

impl<'a, R: Ring, P: AjtaiParams<R>> Sub<Commitment<R, P>> for &'a Commitment<R, P> {
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

impl<'a, R: Ring, P: AjtaiParams<R>> Sub<&'a Commitment<R, P>> for Commitment<R, P> {
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

impl<R: Ring, P: AjtaiParams<R>> Sub<Commitment<R, P>> for Commitment<R, P> {
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

impl<'a, 'b, R: Ring, P: AjtaiParams<R>> Mul<&'b R> for &'a Commitment<R, P> {
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

impl<'b, R: Ring, P: AjtaiParams<R>> Mul<&'b R> for Commitment<R, P> {
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

impl<R: Ring, P: AjtaiParams<R>> Mul<R> for Commitment<R, P> {
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

impl<'a, R: Ring, P: AjtaiParams<R>> Mul<R> for &'a Commitment<R, P> {
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

#[derive(Clone, Debug)]
pub struct AjtaiCommitmentScheme<
    CR: ConvertibleRing,
    R: OverField + Into<CR> + From<CR>,
    P: AjtaiParams<R>,
> {
    _cr: PhantomData<CR>,
    _p: PhantomData<P>,
    matrix: Matrix<R>,
}

impl<CR: ConvertibleRing, R: OverField + Into<CR> + From<CR>, P: AjtaiParams<R>>
    AjtaiCommitmentScheme<CR, R, P>
{
    pub fn rand<Rng: rand::Rng + ?Sized>(rng: &mut Rng) -> Self {
        Self {
            _cr: PhantomData,
            _p: PhantomData,
            matrix: Matrix::rand(P::WITNESS_SIZE, P::OUTPUT_SIZE, rng),
        }
    }
}

impl<CR: ConvertibleRing, R: OverField + Into<CR> + From<CR>, P: AjtaiParams<R>>
    AjtaiCommitmentScheme<CR, R, P>
{
    pub fn commit_ntt(&self, f: &[R]) -> Result<Commitment<R, P>, CommitmentError> {
        // TODO: a lot of clones and copies. Can we optimise this somehow?
        if f.len() != P::WITNESS_SIZE {
            return Err(CommitmentError::WrongWitnessLength(f.len()));
        }

        let commitment_vec = self.matrix.clone() * Vector::from(Vec::from(f));

        Commitment::try_from(commitment_vec.iter().copied().collect::<Vec<_>>())
    }

    pub fn commit_coeff(&self, f: &[CR]) -> Result<Commitment<R, P>, CommitmentError> {
        if f.len() != P::WITNESS_SIZE {
            return Err(CommitmentError::WrongWitnessLength(f.len()));
        }

        self.commit_ntt(&f.iter().map(|&x| x.into()).collect::<Vec<R>>())
    }

    pub fn decompose_and_commit_coeff(
        &self,
        f: &[CR],
    ) -> Result<Commitment<R, P>, CommitmentError> {
        let f = decompose_balanced_vec(f, P::B, Some(P::L))
            .into_iter()
            .flatten()
            .collect::<Vec<_>>();

        self.commit_coeff(&f)
    }

    /// Commits the gadgeted witness, i.e. w = G_B f, for some f.
    pub fn decompose_and_commit_ntt(&self, w: &[R]) -> Result<Commitment<R, P>, CommitmentError> {
        let f: Vec<R> = decompose_balanced_vec(
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