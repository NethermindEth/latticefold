#![allow(non_snake_case)]
use ark_std::marker::PhantomData;
use cyclotomic_rings::SuitableRing;
use lattirust_ring::OverField;

use crate::nifs::error::FoldingError;
use crate::transcript::TranscriptWithSmallChallenges;
use crate::{
    arith::{Witness, CCS, LCCCS},
    decomposition_parameters::DecompositionParams,
    utils::sumcheck,
};

#[derive(Clone)]
pub struct FoldingProof<NTT: OverField> {
    // Step 2.
    pub pointshift_sumcheck_proof: sumcheck::Proof<NTT>,
    // Step 3
    pub theta_s: Vec<NTT>,
    pub eta_s: Vec<Vec<NTT>>,
}

pub trait FoldingProver<NTT: SuitableRing, T: TranscriptWithSmallChallenges<NTT>> {
    fn prove<const C: usize, P: DecompositionParams>(
        cm_i_s: &[LCCCS<C, NTT>],
        w_s: &[Witness<NTT>],
        transcript: &mut impl TranscriptWithSmallChallenges<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(LCCCS<C, NTT>, Witness<NTT>, FoldingProof<NTT>), FoldingError<NTT>>;
}

pub trait FoldingVerifier<NTT: SuitableRing, T: TranscriptWithSmallChallenges<NTT>> {
    fn verify<const C: usize, P: DecompositionParams>(
        cm_i_s: &[LCCCS<C, NTT>],
        proof: &FoldingProof<NTT>,
        transcript: &mut impl TranscriptWithSmallChallenges<NTT>,
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
