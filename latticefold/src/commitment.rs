use ark_std::UniformRand;
use std::{
    marker::PhantomData,
    ops::{Add, Mul, Sub},
};

use lattirust_arithmetic::{
    balanced_decomposition::decompose_balanced_vec,
    linear_algebra::{Matrix, Vector},
    ring::ConvertibleRing,
};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum CommitmentError {
    #[error("Wrong length of the witness {0}")]
    WrongWitnessLength(usize),
}

pub trait AjtaiParams<R: ConvertibleRing> {
    // The MSIS bound.
    const B: u128;
    // The ring modulus should be < B^L.
    const L: usize;
    const WITNESS_SIZE: usize;
    const OUTPUT_SIZE: usize;
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Commitment<R: ConvertibleRing, P: AjtaiParams<R>> {
    _phantom: PhantomData<P>,
    val: Vec<R>,
}

impl<R: ConvertibleRing, P: AjtaiParams<R>> Commitment<R, P> {
    fn from_vec_raw(vec: Vec<R>) -> Self {
        Self {
            _phantom: PhantomData,
            val: vec,
        }
    }
}

impl<'a, R: ConvertibleRing, P: AjtaiParams<R>> TryFrom<&'a [R]> for Commitment<R, P> {
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

impl<R: ConvertibleRing, P: AjtaiParams<R>> TryFrom<Vec<R>> for Commitment<R, P> {
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

impl<'a, R: ConvertibleRing, P: AjtaiParams<R>> Add<&'a Commitment<R, P>> for &'a Commitment<R, P> {
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

impl<'a, R: ConvertibleRing, P: AjtaiParams<R>> Add<Commitment<R, P>> for &'a Commitment<R, P> {
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

impl<'a, R: ConvertibleRing, P: AjtaiParams<R>> Add<&'a Commitment<R, P>> for Commitment<R, P> {
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

impl<R: ConvertibleRing, P: AjtaiParams<R>> Add<Commitment<R, P>> for Commitment<R, P> {
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

impl<'a, R: ConvertibleRing, P: AjtaiParams<R>> Sub<&'a Commitment<R, P>> for &'a Commitment<R, P> {
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

impl<'a, R: ConvertibleRing, P: AjtaiParams<R>> Sub<Commitment<R, P>> for &'a Commitment<R, P> {
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

impl<'a, R: ConvertibleRing, P: AjtaiParams<R>> Sub<&'a Commitment<R, P>> for Commitment<R, P> {
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

impl<R: ConvertibleRing, P: AjtaiParams<R>> Sub<Commitment<R, P>> for Commitment<R, P> {
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

impl<'a, 'b, R: ConvertibleRing, P: AjtaiParams<R>> Mul<&'b R> for &'a Commitment<R, P> {
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

impl<'b, R: ConvertibleRing, P: AjtaiParams<R>> Mul<&'b R> for Commitment<R, P> {
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

impl<R: ConvertibleRing, P: AjtaiParams<R>> Mul<R> for Commitment<R, P> {
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

impl<'a, R: ConvertibleRing, P: AjtaiParams<R>> Mul<R> for &'a Commitment<R, P> {
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
pub struct AjtaiCommitmentScheme<R: ConvertibleRing, P: AjtaiParams<R>> {
    _phantom: PhantomData<P>,
    matrix: Matrix<R>,
}

impl<R: ConvertibleRing + UniformRand, P: AjtaiParams<R>> AjtaiCommitmentScheme<R, P> {
    pub fn rand<Rng: rand::Rng + ?Sized>(rng: &mut Rng) -> Self {
        Self {
            _phantom: PhantomData,
            matrix: Matrix::rand(P::WITNESS_SIZE, P::OUTPUT_SIZE, rng),
        }
    }
}

impl<R: ConvertibleRing, P: AjtaiParams<R>> AjtaiCommitmentScheme<R, P> {
    pub fn commit_pre_gadget(&self, f: &[R]) -> Result<Commitment<R, P>, CommitmentError> {
        // TODO: a lot of clones and copies. Can we optimise this somehow?
        if f.len() != P::WITNESS_SIZE {
            return Err(CommitmentError::WrongWitnessLength(f.len()));
        }

        let commitment_vec = self.matrix.clone() * Vector::from(Vec::from(f));

        Commitment::try_from(commitment_vec.iter().copied().collect::<Vec<_>>())
    }

    /// Commits the gadgeted witness, i.e. w = G_B f, for some f.
    pub fn commit(&self, w: &[R]) -> Result<Commitment<R, P>, CommitmentError> {
        let f: Vec<R> = decompose_balanced_vec(w, P::B, Some(P::L))
            .iter()
            .flatten()
            .copied()
            .collect();

        self.commit_pre_gadget(&f)
    }
}
