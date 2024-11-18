use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::marker::PhantomData;
use cyclotomic_rings::rings::SuitableRing;

use crate::{
    arith::{Witness, CCCS, CCS, LCCCS},
    nifs::error::LinearizationError,
    transcript::Transcript,
    utils::sumcheck,
};

use lattirust_poly::mle::DenseMultilinearExtension;
use lattirust_ring::OverField;

#[derive(Clone, Debug, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct LinearizationProof<NTT: OverField> {
    // Sent in the step 2. of the linearization subprotocol
    pub linearization_sumcheck: sumcheck::Proof<NTT>,
    // Sent in the step 3.
    pub v: Vec<NTT>,
    pub u: Vec<NTT>,
}

pub(crate) trait ChallengeGenerator<NTT: OverField> {
    fn generate_challenges(transcript: &mut impl Transcript<NTT>, log_m: usize) -> Vec<NTT>;
}

pub(crate) struct BetaChallengeGenerator<NTT> {
    _ntt: PhantomData<NTT>,
}

pub(crate) struct ProverState<NTT: OverField> {
    pub beta_s: Vec<NTT>,
    pub z_ccs: Vec<NTT>,
    pub Mz_mles: Vec<DenseMultilinearExtension<NTT>>,
}

pub(crate) struct VerifierState<NTT: OverField> {
    pub beta_s: Vec<NTT>,
    pub point_r: Vec<NTT>,
    pub s: NTT,
}

pub trait LinearizationProver<NTT: SuitableRing, T: Transcript<NTT>> {
    fn prove<const C: usize>(
        cm_i: &CCCS<C, NTT>,
        wit: &Witness<NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(LCCCS<C, NTT>, LinearizationProof<NTT>), LinearizationError<NTT>>;
}

pub trait LinearizationVerifier<NTT: OverField, T: Transcript<NTT>> {
    fn verify<const C: usize>(
        cm_i: &CCCS<C, NTT>,
        proof: &LinearizationProof<NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<LCCCS<C, NTT>, LinearizationError<NTT>>;
}

pub struct LFLinearizationProver<NTT, T> {
    _ntt: PhantomData<NTT>,
    _t: PhantomData<T>,
}

pub struct LFLinearizationVerifier<NTT, T> {
    _ntt: PhantomData<NTT>,
    _t: PhantomData<T>,
}
