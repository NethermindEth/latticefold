use crate::arith::{Witness, CCS, LCCCS};
use crate::decomposition_parameters::DecompositionParams;
use crate::nifs::error::FoldingError;
use crate::utils::sumcheck;
use crate::{ark_base::Vec, transcript::TranscriptWithShortChallenges};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::marker::PhantomData;
use cyclotomic_rings::rings::SuitableRing;
use lattirust_poly::mle::DenseMultilinearExtension;
use lattirust_ring::OverField;

#[derive(Clone, Debug, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct FoldingProof<NTT: OverField> {
    // Step 2.
    pub pointshift_sumcheck_proof: sumcheck::Proof<NTT>,
    // Step 3
    pub theta_s: Vec<Vec<NTT>>,
    pub eta_s: Vec<Vec<NTT>>,
}

pub trait FoldingProver<NTT: SuitableRing, T: TranscriptWithShortChallenges<NTT>> {
    fn prove<const C: usize, P: DecompositionParams>(
        cm_i_s: &[LCCCS<C, NTT>],
        w_s: Vec<Witness<NTT>>,
        transcript: &mut impl TranscriptWithShortChallenges<NTT>,
        ccs: &CCS<NTT>,
        mz_mles: &[Vec<DenseMultilinearExtension<NTT>>],
    ) -> Result<(LCCCS<C, NTT>, Witness<NTT>, FoldingProof<NTT>), FoldingError<NTT>>;
}

pub trait FoldingVerifier<NTT: SuitableRing, T: TranscriptWithShortChallenges<NTT>> {
    fn verify<const C: usize, P: DecompositionParams>(
        cm_i_s: &[LCCCS<C, NTT>],
        proof: &FoldingProof<NTT>,
        transcript: &mut impl TranscriptWithShortChallenges<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<LCCCS<C, NTT>, FoldingError<NTT>>;
}

pub struct LFFoldingProver<NTT, T> {
    _ntt: PhantomData<NTT>,
    _t: PhantomData<T>,
}

pub struct LFFoldingVerifier<NTT, T> {
    _ntt: PhantomData<NTT>,
    _t: PhantomData<T>,
}
