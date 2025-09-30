use cyclotomic_rings::rings::SuitableRing;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use stark_rings::{
    balanced_decomposition::DecomposeToVec,
    cyclotomic_ring::{CRT, ICRT},
    Ring,
};
use stark_rings_linalg::Matrix;

use super::homomorphic_commitment::Commitment;
use crate::{
    ark_base::*, commitment::CommitmentError, decomposition_parameters::DecompositionParams,
};

/// A concrete instantiation of the Ajtai commitment scheme.
/// Contains a random Ajtai matrix (kappa x n).
#[derive(Clone, Debug)]
pub struct AjtaiCommitmentScheme<R> {
    matrix: Matrix<R>,
}

impl<R> AjtaiCommitmentScheme<R> {
    /// Create a new scheme using the provided Ajtai matrix
    pub fn new(matrix: Matrix<R>) -> Self {
        Self { matrix }
    }
}

impl<R: Ring> AjtaiCommitmentScheme<R> {
    /// Returns a random Ajtai commitment matrix
    pub fn rand<Rng: rand::Rng + ?Sized>(kappa: usize, n: usize, rng: &mut Rng) -> Self {
        Self::new(vec![vec![R::rand(rng); n]; kappa].into())
    }
}

impl<R: Ring> AjtaiCommitmentScheme<R> {
    /// Commit to a witness
    pub fn commit(&self, f: &[R]) -> Result<Commitment<R>, CommitmentError> {
        if f.len() != self.matrix.ncols {
            return Err(CommitmentError::WrongWitnessLength(f.len(), self.matrix.ncols));
        }

        let commitment = self.matrix.checked_mul_vec(f).ok_or(CommitmentError::WrongWitnessLength(f.len(), self.matrix.ncols))?;

        Ok(Commitment::from_vec_raw(commitment))
    }

    /// Ajtai matrix number of rows
    ///
    /// This value affects the security of the scheme.
    pub fn kappa(&self) -> usize {
        self.matrix.nrows
    }

    /// Ajtai matrix number of columns
    ///
    /// The size of the witness must be equal to this value.
    pub fn width(&self) -> usize {
        self.matrix.ncols
    }
}

// SuitableRing helpers
impl<NTT: SuitableRing> AjtaiCommitmentScheme<NTT> {
    /// Commit to a witness in the NTT form.
   /// The most basic one just multiplies by the matrix.
    pub fn commit_ntt(&self, f: &[NTT]) -> Result<Commitment<NTT>, CommitmentError> {
        self.commit(f)
     }

    /// Commit to a witness in the coefficient form.
    /// Performs NTT on each component of the witness and then does Ajtai commitment.
    pub fn commit_coeff<P: DecompositionParams>(
        &self,
        f: Vec<NTT::CoefficientRepresentation>,
    ) -> Result<Commitment<NTT>, CommitmentError> {
        self.commit_ntt(&CRT::elementwise_crt(f))
    }

    /// Takes a coefficient form witness, decomposes it vertically in radix-B,
    /// i.e. computes a preimage G_B^{-1}(w), and Ajtai commits to the result.
    pub fn decompose_and_commit_coeff<P: DecompositionParams>(
        &self,
        f: &[NTT::CoefficientRepresentation],
    ) -> Result<Commitment<NTT>, CommitmentError> {
        let f = f
            .decompose_to_vec(P::B, P::L)
            .into_iter()
            .flatten()
            .collect::<Vec<_>>();

        self.commit_coeff::<P>(f)
    }

    /// Takes an NTT form witness, transforms it into the coefficient form,
    /// decomposes it vertically in radix-B, i.e.
    /// computes a preimage G_B^{-1}(w), and Ajtai commits to the result.
    pub fn decompose_and_commit_ntt<P: DecompositionParams>(
        &self,
        w: Vec<NTT>,
    ) -> Result<Commitment<NTT>, CommitmentError> {
        let coeff: Vec<NTT::CoefficientRepresentation> = ICRT::elementwise_icrt(w);

        self.decompose_and_commit_coeff::<P>(&coeff)
    }
}

#[cfg(test)]
mod tests {
    use cyclotomic_rings::rings::GoldilocksRingNTT;
    use stark_rings::OverField;

    use super::{AjtaiCommitmentScheme, CommitmentError};
    use crate::ark_base::*;

    pub(crate) fn generate_ajtai<NTT: OverField>(
        kappa: usize,
        n: usize,
    ) -> AjtaiCommitmentScheme<NTT> {
        let mut matrix = Vec::<Vec<NTT>>::new();

        for i in 0..kappa {
            let mut row = Vec::<NTT>::new();
            for j in 0..n {
                row.push(NTT::from((i * n + j) as u128));
            }
            matrix.push(row)
        }

        AjtaiCommitmentScheme::new(matrix.into())
    }

    #[test]
    fn test_commit_ntt() -> Result<(), CommitmentError> {
        const WITNESS_SIZE: usize = 1 << 15;
        const OUTPUT_SIZE: usize = 9;

        let ajtai_data: AjtaiCommitmentScheme<GoldilocksRingNTT> =
            generate_ajtai(OUTPUT_SIZE, WITNESS_SIZE);
        let witness: Vec<_> = (0..(1 << 15)).map(|_| 2_u128.into()).collect();

        let committed = ajtai_data.commit_ntt(&witness)?;

        for (i, &x) in committed.as_ref().iter().enumerate() {
            let expected: u128 =
                ((WITNESS_SIZE) * (2 * i * WITNESS_SIZE + (WITNESS_SIZE - 1))) as u128;
            assert_eq!(x, expected.into());
        }

        Ok(())
    }
}
