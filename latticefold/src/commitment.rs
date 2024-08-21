use lattirust_arithmetic::{
    balanced_decomposition::decompose_balanced_slice_polyring,
    challenge_set::latticefold_challenge_set::OverField,
    ring::{PolyRing, Ring},
};
use std::ops::{Add, Mul, Sub};
use thiserror::Error;

use crate::parameters::DecompositionParams;

/// A concrete instantiation of the Ajtai commitment scheme.
/// Contains a random Ajtai matrix for the corresponding Ajtai parameters
/// `C` is the length of commitment vectors or, equivalently, the number of rows of the Ajtai matrix.
/// `W` is the length of witness vectors or, equivalently, the number of columns of the Ajtai matrix.
/// `NTT` is a cyclotomic ring, for better results it should be in the NTT form.
#[derive(Clone, Debug)]
pub struct AjtaiCommitmentScheme<const C: usize, const W: usize, NTT: Ring> {
    matrix: Vec<Vec<NTT>>,
}

impl<const C: usize, const W: usize, NTT: OverField> TryFrom<Vec<Vec<NTT>>>
    for AjtaiCommitmentScheme<C, W, NTT>
{
    type Error = CommitmentError;

    fn try_from(matrix: Vec<Vec<NTT>>) -> Result<Self, Self::Error> {
        if matrix.len() != C || matrix[0].len() != W {
            return Err(CommitmentError::WrongAjtaiMatrixDimensions(
                matrix.len(),
                matrix[0].len(),
                C,
                W,
            ));
        }

        let mut ajtai_matrix: Vec<Vec<NTT>> = Vec::with_capacity(C);

        for row in matrix.into_iter() {
            let len = row.len();

            if len != W {
                return Err(CommitmentError::WrongAjtaiMatrixDimensions(C, len, C, W));
            }
            ajtai_matrix.push(row)
        }

        Ok(Self {
            matrix: ajtai_matrix,
        })
    }
}

impl<const C: usize, const W: usize, NTT: OverField> AjtaiCommitmentScheme<C, W, NTT> {
    pub fn rand<Rng: rand::Rng + ?Sized>(rng: &mut Rng) -> Self {
        Self {
            matrix: vec![vec![NTT::rand(rng); W]; C],
        }
    }
}

impl<const C: usize, const W: usize, NTT: OverField> AjtaiCommitmentScheme<C, W, NTT> {
    /// Commit to a witness in the NTT form.
    /// The most basic one just multiplies by the matrix.
    pub fn commit_ntt(&self, f: &[NTT]) -> Result<Commitment<C, NTT>, CommitmentError> {
        if f.len() != W {
            return Err(CommitmentError::WrongWitnessLength(f.len(), W));
        }

        let mut commitment: Vec<NTT> = vec![NTT::zero(); C];

        commitment
            .iter_mut()
            .zip(&self.matrix)
            .for_each(|(x, row)| *x = row.iter().zip(f).map(|(&m, &x)| m * x).sum());

        Ok(Commitment::from_vec_raw(commitment))
    }

    /// Commit to a witness in the coefficient form.
    /// Performs NTT on each component of the witness and then does Ajtai commitment.
    pub fn commit_coeff<CR: PolyRing + From<NTT> + Into<NTT>, P: DecompositionParams>(
        &self,
        f: &[CR],
    ) -> Result<Commitment<C, NTT>, CommitmentError> {
        if f.len() != W {
            return Err(CommitmentError::WrongWitnessLength(f.len(), W));
        }

        self.commit_ntt(&f.iter().map(|&x| x.into()).collect::<Vec<NTT>>())
    }

    /// Takes a coefficient form witness, decomposes it vertically in radix-B
    /// and Ajtai commits to the result.
    pub fn decompose_and_commit_coeff<
        CR: PolyRing + From<NTT> + Into<NTT>,
        P: DecompositionParams,
    >(
        &self,
        f: &[CR],
    ) -> Result<Commitment<C, NTT>, CommitmentError> {
        let f = decompose_balanced_slice_polyring(f, P::B, Some(P::L))
            .into_iter()
            .flatten()
            .collect::<Vec<_>>();

        self.commit_coeff::<_, P>(&f)
    }

    /// Takes an NTT form witness, transforms it into the coefficient form,
    /// decomposes it vertically in radix-B and Ajtai commits to the result.
    pub fn decompose_and_commit_ntt<
        CR: PolyRing + From<NTT> + Into<NTT>,
        P: DecompositionParams,
    >(
        &self,
        w: &[NTT],
    ) -> Result<Commitment<C, NTT>, CommitmentError> {
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

#[derive(Debug, Error)]
pub enum CommitmentError {
    #[error("Wrong length of the witness: {0}, expected: {1}")]
    WrongWitnessLength(usize, usize),
    #[error("Wrong length of the commitment: {0}, expected: {1}")]
    WrongCommitmentLength(usize, usize),
    #[error("Ajtai matrix has dimensions: {0}x{1}, expected: {2}x{3}")]
    WrongAjtaiMatrixDimensions(usize, usize, usize, usize),
}

/// The Ajtai commitment type. Meant to contain the output of the
/// matrix-vector multiplication `A \cdot x`.
/// Enforced to have the length `C`.
/// Since Ajtai commitment is bounded-additively homomorphic
/// one can add commitments and multiply them by a scalar.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Commitment<const C: usize, R: Ring> {
    val: Vec<R>,
}

impl<const C: usize, R: Ring> Commitment<C, R> {
    fn from_vec_raw(vec: Vec<R>) -> Self {
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
        let mut res_vec = Vec::<R>::with_capacity(C);

        res_vec
            .iter_mut()
            .zip(self.val.iter())
            .zip(rhs.val.iter())
            .for_each(|((res, &a), &b)| *res = a + b);

        Commitment::from_vec_raw(res_vec)
    }
}

impl<'a, const C: usize, R: Ring> Add<Commitment<C, R>> for &'a Commitment<C, R> {
    type Output = Commitment<C, R>;

    fn add(self, rhs: Commitment<C, R>) -> Self::Output {
        self + &rhs
    }
}

impl<'a, const C: usize, R: Ring> Add<&'a Commitment<C, R>> for Commitment<C, R> {
    type Output = Commitment<C, R>;

    fn add(self, rhs: &'a Commitment<C, R>) -> Self::Output {
        &self + rhs
    }
}

impl<const C: usize, R: Ring> Add<Commitment<C, R>> for Commitment<C, R> {
    type Output = Commitment<C, R>;

    fn add(self, rhs: Commitment<C, R>) -> Self::Output {
        &self + &rhs
    }
}

impl<'a, 'b, const C: usize, R: Ring> Sub<&'a Commitment<C, R>> for &'b Commitment<C, R> {
    type Output = Commitment<C, R>;

    fn sub(self, rhs: &'a Commitment<C, R>) -> Self::Output {
        let mut res_vec = Vec::<R>::with_capacity(C);

        res_vec
            .iter_mut()
            .zip(self.val.iter())
            .zip(rhs.val.iter())
            .for_each(|((res, &a), &b)| *res = a - b);

        Commitment::from_vec_raw(res_vec)
    }
}

impl<'a, const C: usize, R: Ring> Sub<Commitment<C, R>> for &'a Commitment<C, R> {
    type Output = Commitment<C, R>;

    fn sub(self, rhs: Commitment<C, R>) -> Self::Output {
        self - &rhs
    }
}

impl<'a, const C: usize, R: Ring> Sub<&'a Commitment<C, R>> for Commitment<C, R> {
    type Output = Commitment<C, R>;

    fn sub(self, rhs: &'a Commitment<C, R>) -> Self::Output {
        &self - rhs
    }
}

impl<const C: usize, R: Ring> Sub<Commitment<C, R>> for Commitment<C, R> {
    type Output = Commitment<C, R>;

    fn sub(self, rhs: Commitment<C, R>) -> Self::Output {
        &self - &rhs
    }
}

impl<'a, 'b, const C: usize, R: Ring> Mul<&'a R> for &'b Commitment<C, R> {
    type Output = Commitment<C, R>;

    fn mul(self, rhs: &'a R) -> Self::Output {
        let mut res_vec = Vec::<R>::with_capacity(C);

        res_vec
            .iter_mut()
            .zip(self.val.iter())
            .for_each(|(res, &a)| *res = a * rhs);

        Commitment::from_vec_raw(res_vec)
    }
}

impl<'a, const C: usize, R: Ring> Mul<&'a R> for Commitment<C, R> {
    type Output = Commitment<C, R>;

    fn mul(self, rhs: &'a R) -> Self::Output {
        &self * rhs
    }
}

impl<const C: usize, R: Ring> Mul<R> for Commitment<C, R> {
    type Output = Commitment<C, R>;

    #[allow(clippy::op_ref)]
    fn mul(self, rhs: R) -> Self::Output {
        &self * &rhs
    }
}

impl<'a, const C: usize, R: Ring> Mul<R> for &'a Commitment<C, R> {
    type Output = Commitment<C, R>;

    #[allow(clippy::op_ref)]
    fn mul(self, rhs: R) -> Self::Output {
        self * &rhs
    }
}

// TODO: use macros to implement the other operations

#[cfg(test)]
mod tests {
    use lattirust_arithmetic::challenge_set::latticefold_challenge_set::OverField;

    use super::{AjtaiCommitmentScheme, CommitmentError};
    use crate::parameters::DilithiumNTT;

    pub(crate) fn generate_ajtai<const C: usize, const W: usize, NTT: OverField>(
    ) -> Result<AjtaiCommitmentScheme<C, W, NTT>, CommitmentError> {
        let mut matrix = Vec::<Vec<NTT>>::new();

        for i in 0..C {
            let mut row = Vec::<NTT>::new();
            for j in 0..W {
                row.push(NTT::from((i * W + j) as u128));
            }
            matrix.push(row)
        }

        AjtaiCommitmentScheme::try_from(matrix)
    }

    #[test]
    fn test_commit_ntt() -> Result<(), CommitmentError> {
        const WITNESS_SIZE: usize = 1 << 15;
        const OUTPUT_SIZE: usize = 9;

        let ajtai_data: AjtaiCommitmentScheme<OUTPUT_SIZE, WITNESS_SIZE, DilithiumNTT> =
            generate_ajtai()?;
        let input: Vec<_> = (0..(1 << 15)).map(|_| 2_u128.into()).collect();

        let committed = ajtai_data.commit_ntt(&input)?;

        for (i, &x) in committed.as_ref().iter().enumerate() {
            let expected: u128 =
                ((WITNESS_SIZE) * (2 * i * WITNESS_SIZE + (WITNESS_SIZE - 1))) as u128;
            assert_eq!(x, expected.into());
        }

        Ok(())
    }
}
