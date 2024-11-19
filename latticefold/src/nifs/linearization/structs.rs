use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::marker::PhantomData;
use cyclotomic_rings::rings::SuitableRing;

use crate::{
    arith::{Witness, CCCS, CCS, LCCCS},
    nifs::error::LinearizationError,
    transcript::Transcript,
    utils::sumcheck,
};

use crate::utils::sumcheck::Proof;
use lattirust_poly::polynomials::VirtualPolynomial;
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

pub trait LinearizationProver<NTT: SuitableRing, T: Transcript<NTT>> {
    fn prove<const C: usize>(
        cm_i: &CCCS<C, NTT>,
        wit: &Witness<NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(LCCCS<C, NTT>, LinearizationProof<NTT>), LinearizationError<NTT>>;

    fn compute_z_ccs<const C: usize>(
        wit: &Witness<NTT>,
        x_ccs: &[NTT],
    ) -> Result<Vec<NTT>, LinearizationError<NTT>>;

    fn construct_polynomial_g(
        z_state: &[NTT],
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(VirtualPolynomial<NTT>, Vec<NTT>), LinearizationError<NTT>>;

    fn generate_sumcheck_proof(
        g: &VirtualPolynomial<NTT>,
        beta_s: &[NTT],
        transcript: &mut impl Transcript<NTT>,
    ) -> Result<(Proof<NTT>, Vec<NTT>), LinearizationError<NTT>>;

    fn compute_evaluation_vectors(
        wit: &Witness<NTT>,
        point_r: &[NTT],
        ccs: &CCS<NTT>,
        z_state: &[NTT],
    ) -> Result<(Vec<NTT>, Vec<NTT>, Vec<NTT>), LinearizationError<NTT>>;
}

pub trait LinearizationVerifier<NTT: OverField, T: Transcript<NTT>> {
    fn verify<const C: usize>(
        cm_i: &CCCS<C, NTT>,
        proof: &LinearizationProof<NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<LCCCS<C, NTT>, LinearizationError<NTT>>;

    fn verify_sumcheck_proof(
        proof: &LinearizationProof<NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(Vec<NTT>, NTT), LinearizationError<NTT>>;

    fn verify_evaluation_claim(
        beta_s: &[NTT],
        point_r: &[NTT],
        s: NTT,
        proof: &LinearizationProof<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(), LinearizationError<NTT>>;

    fn prepare_final_state<const C: usize>(
        cm_i: &CCCS<C, NTT>,
        point_r: Vec<NTT>,
        proof: &LinearizationProof<NTT>,
    ) -> LCCCS<C, NTT>;
}

pub struct LFLinearizationProver<NTT, T> {
    _ntt: PhantomData<NTT>,
    _t: PhantomData<T>,
}

pub struct LFLinearizationVerifier<NTT, T> {
    _ntt: PhantomData<NTT>,
    _t: PhantomData<T>,
}
