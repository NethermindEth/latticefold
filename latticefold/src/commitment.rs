use lattirust_arithmetic::{
    balanced_decomposition::{decompose_balanced_slice_polyring, pad_and_transpose},
    challenge_set::latticefold_challenge_set::OverField,
    ring::{PolyRing, Pow2CyclotomicPolyRing, Pow2CyclotomicPolyRingNTT, Ring, Zq},
};
use std::{
    fmt::Display,
    marker::PhantomData,
    ops::{Add, Mul, Sub},
};
use thiserror::Error;

/// A concrete instantiation of the Ajtai commitment scheme.
/// Contains a random Ajtai matrix for the coresssponding Ajtai parameters
/// `CR` is the type parameter for the coefficient representation of the ring
/// `NTT` is the NTT representation of the same ring.
#[derive(Clone, Debug)]
pub struct AjtaiCommitmentScheme<CR: PolyRing, NTT: OverField, P: AjtaiParams>
where
    CR: Into<NTT> + From<NTT>,
{
    _cr: PhantomData<CR>,
    _p: PhantomData<P>,
    matrix: Vec<Vec<NTT>>,
}

impl<CR: PolyRing, NTT: OverField, P: AjtaiParams> TryFrom<Vec<Vec<NTT>>>
    for AjtaiCommitmentScheme<CR, NTT, P>
where
    CR: Into<NTT> + From<NTT>,
{
    type Error = CommitmentError;

    fn try_from(matrix: Vec<Vec<NTT>>) -> Result<Self, Self::Error> {
        if matrix.len() != P::OUTPUT_SIZE || matrix[0].len() != P::WITNESS_SIZE {
            return Err(CommitmentError::WrongAjtaiMatrixDimensions(
                matrix.len(),
                matrix[0].len(),
                P::OUTPUT_SIZE,
                P::WITNESS_SIZE,
            ));
        }

        Ok(Self {
            _cr: PhantomData,
            _p: PhantomData,
            matrix,
        })
    }
}

impl<CR: PolyRing, NTT: OverField, P: AjtaiParams> AjtaiCommitmentScheme<CR, NTT, P>
where
    CR: Into<NTT> + From<NTT>,
{
    pub fn rand<Rng: rand::Rng + ?Sized>(rng: &mut Rng) -> Self {
        let mut matrix = Vec::<Vec<NTT>>::with_capacity(P::OUTPUT_SIZE);

        for _i in 0..P::OUTPUT_SIZE {
            let mut row = Vec::<NTT>::with_capacity(P::WITNESS_SIZE);
            for _j in 0..P::WITNESS_SIZE {
                row.push(NTT::rand(rng));
            }
            matrix.push(row);
        }

        Self {
            _cr: PhantomData,
            _p: PhantomData,
            matrix,
        }
    }

    /// Commit to a witness in the NTT form.
    /// The most basic one just multiplies by the matrix.
    pub fn commit_ntt(&self, f: &[NTT]) -> Result<Commitment<NTT, P>, CommitmentError> {
        if f.len() != P::WITNESS_SIZE {
            return Err(CommitmentError::WrongWitnessLength(
                f.len(),
                P::WITNESS_SIZE,
            ));
        }

        Ok(Commitment::from_vec_raw(
            self.matrix
                .iter()
                .map(|row| row.iter().zip(f).map(|(&m, &x)| m * x).sum())
                .collect(),
        ))
    }

    /// Commit to a witness in the coefficient form.
    /// Performs NTT on each component of the witness and then does Ajtai commitment.
    pub fn commit_coeff(&self, f: &[CR]) -> Result<Commitment<NTT, P>, CommitmentError> {
        if f.len() != P::WITNESS_SIZE {
            return Err(CommitmentError::WrongWitnessLength(
                f.len(),
                P::WITNESS_SIZE,
            ));
        }

        self.commit_ntt(&f.iter().map(|&x| x.into()).collect::<Vec<NTT>>())
    }

    /// Takes a coefficient form witness, decomposes it vertically in radix-B
    /// and Ajtai commits to the result.
    pub fn decompose_and_commit_coeff(
        &self,
        f: &[CR],
    ) -> Result<Commitment<NTT, P>, CommitmentError> {
        let f = pad_and_transpose(decompose_balanced_slice_polyring(f, P::B, Some(P::L)))
            .into_iter()
            .flatten()
            .collect::<Vec<_>>();

        self.commit_coeff(&f)
    }

    /// Takes an NTT form witness, transforms it into the coefficient form,
    /// decomposes it vertically in radix-B and Ajtai commits to the result.
    pub fn decompose_and_commit_ntt(
        &self,
        w: &[NTT],
    ) -> Result<Commitment<NTT, P>, CommitmentError> {
        let f: Vec<NTT> = pad_and_transpose(decompose_balanced_slice_polyring(
            &w.iter().map(|&x| x.into()).collect::<Vec<CR>>(),
            P::B,
            Some(P::L),
        ))
        .iter()
        .flatten()
        .map(|&x| x.into())
        .collect();

        self.commit_ntt(&f)
    }
}

#[derive(Debug, Error)]
pub enum CommitmentError {
    #[error("Wrong length of the witness: {0}, expected: {1}")]
    WrongWitnessLength(usize, usize),
    #[error("Ajtai matrix has dimensions: {0}x{1}, expected: {2}x{3}")]
    WrongAjtaiMatrixDimensions(usize, usize, usize, usize),
}

/// Ajtai commitment parameters.
/// Convenient to enforce them compile-time.
pub trait AjtaiParams: Clone + Display + Eq + PartialEq {
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

impl<R: Ring, P: AjtaiParams> AsRef<[R]> for Commitment<R, P> {
    fn as_ref(&self) -> &[R] {
        &self.val
    }
}

impl<'a, 'b, R: Ring, P: AjtaiParams> Add<&'a Commitment<R, P>> for &'b Commitment<R, P> {
    type Output = Commitment<R, P>;

    fn add(self, rhs: &'a Commitment<R, P>) -> Self::Output {
        let mut res_vec = vec![R::zero(); P::OUTPUT_SIZE];

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
        self + &rhs
    }
}

impl<'a, R: Ring, P: AjtaiParams> Add<&'a Commitment<R, P>> for Commitment<R, P> {
    type Output = Commitment<R, P>;

    fn add(self, rhs: &'a Commitment<R, P>) -> Self::Output {
        &self + rhs
    }
}

impl<R: Ring, P: AjtaiParams> Add<Commitment<R, P>> for Commitment<R, P> {
    type Output = Commitment<R, P>;

    fn add(self, rhs: Commitment<R, P>) -> Self::Output {
        &self + &rhs
    }
}

impl<'a, 'b, R: Ring, P: AjtaiParams> Sub<&'a Commitment<R, P>> for &'b Commitment<R, P> {
    type Output = Commitment<R, P>;

    fn sub(self, rhs: &'a Commitment<R, P>) -> Self::Output {
        let mut res_vec = vec![R::zero(); P::OUTPUT_SIZE];

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
        self - &rhs
    }
}

impl<'a, R: Ring, P: AjtaiParams> Sub<&'a Commitment<R, P>> for Commitment<R, P> {
    type Output = Commitment<R, P>;

    fn sub(self, rhs: &'a Commitment<R, P>) -> Self::Output {
        &self - rhs
    }
}

impl<R: Ring, P: AjtaiParams> Sub<Commitment<R, P>> for Commitment<R, P> {
    type Output = Commitment<R, P>;

    fn sub(self, rhs: Commitment<R, P>) -> Self::Output {
        &self - &rhs
    }
}

impl<'a, 'b, R: Ring, P: AjtaiParams> Mul<&'a R> for &'b Commitment<R, P> {
    type Output = Commitment<R, P>;

    fn mul(self, rhs: &'a R) -> Self::Output {
        let mut res_vec = vec![R::zero(); P::OUTPUT_SIZE];
        res_vec
            .iter_mut()
            .zip(self.val.iter())
            .for_each(|(res, &a)| *res = a * rhs);

        Commitment::from_vec_raw(res_vec)
    }
}

impl<'a, R: Ring, P: AjtaiParams> Mul<&'a R> for Commitment<R, P> {
    type Output = Commitment<R, P>;

    fn mul(self, rhs: &'a R) -> Self::Output {
        &self * rhs
    }
}

impl<R: Ring, P: AjtaiParams> Mul<R> for Commitment<R, P> {
    type Output = Commitment<R, P>;

    #[allow(clippy::op_ref)]
    fn mul(self, rhs: R) -> Self::Output {
        &self * &rhs
    }
}

impl<'a, R: Ring, P: AjtaiParams> Mul<R> for &'a Commitment<R, P> {
    type Output = Commitment<R, P>;

    #[allow(clippy::op_ref)]
    fn mul(self, rhs: R) -> Self::Output {
        self * &rhs
    }
}

// TODO: use macros to implement the other operations

// Some classic lattice parameter sets.

pub const DILITHIUM_PRIME: u64 = 0x00000000_007FE001;

pub type DilithiumCR = Pow2CyclotomicPolyRing<Zq<DILITHIUM_PRIME>, 256>;
pub type DilithiumNTT = Pow2CyclotomicPolyRingNTT<DILITHIUM_PRIME, 256>;

#[derive(Clone, Copy, Eq, PartialEq)]
pub struct DilithiumTestParams;

// TODO: Revise this later
impl AjtaiParams for DilithiumTestParams {
    const B: u128 = 1 << 13;
    const L: usize = 2;
    const WITNESS_SIZE: usize = 1 << 15;
    const OUTPUT_SIZE: usize = 9;
}

impl Display for DilithiumTestParams {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "B={}, l={}, m={}, n={}",
            Self::B,
            Self::L,
            Self::OUTPUT_SIZE,
            Self::WITNESS_SIZE
        )
    }
}

#[cfg(test)]
mod tests {
    use lattirust_arithmetic::{
        challenge_set::latticefold_challenge_set::OverField, ring::PolyRing,
    };

    use super::{
        AjtaiCommitmentScheme, AjtaiParams, CommitmentError, DilithiumCR, DilithiumNTT,
        DilithiumTestParams,
    };

    pub(crate) fn generate_ajtai<
        CR: PolyRing + From<NTT> + Into<NTT>,
        NTT: OverField,
        P: AjtaiParams,
    >(
        m: usize,
        n: usize,
    ) -> Result<AjtaiCommitmentScheme<CR, NTT, P>, CommitmentError> {
        let mut matrix = Vec::<Vec<NTT>>::new();

        for i in 0..m {
            let mut row = Vec::<NTT>::new();
            for j in 0..n {
                row.push(NTT::from((i * n + j) as u128));
            }
            matrix.push(row)
        }

        AjtaiCommitmentScheme::<CR, NTT, P>::try_from(matrix)
    }

    #[test]
    fn test_commit_ntt() -> Result<(), CommitmentError> {
        let ajtai_data: AjtaiCommitmentScheme<DilithiumCR, DilithiumNTT, DilithiumTestParams> =
            generate_ajtai(
                DilithiumTestParams::OUTPUT_SIZE,
                DilithiumTestParams::WITNESS_SIZE,
            )?;
        let input: Vec<_> = (0..(1 << 15)).map(|_| 2_u128.into()).collect();

        let committed = ajtai_data.commit_ntt(&input)?;

        for (i, &x) in committed.as_ref().iter().enumerate() {
            let expected: u128 = ((DilithiumTestParams::WITNESS_SIZE)
                * (2 * i * DilithiumTestParams::WITNESS_SIZE
                    + (DilithiumTestParams::WITNESS_SIZE - 1)))
                as u128;
            assert_eq!(x, expected.into());
        }

        Ok(())
    }
}
