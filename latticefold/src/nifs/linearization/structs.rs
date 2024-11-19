use ark_ff::PrimeField;
use ark_ff::Field;
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
use lattirust_poly::mle::DenseMultilinearExtension;
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

pub(crate) struct VerifierState<NTT: OverField> {
    pub beta_s: Vec<NTT>,
    pub point_r: Vec<NTT>,
    pub s: NTT,
}

pub(crate) struct EvaluationState<NTT: OverField> {
    pub point_r: Vec<NTT>,
    pub v: Vec<NTT>,
    pub u: Vec<NTT>,
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
    ) -> Result<EvaluationState<NTT>, LinearizationError<NTT>>;
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
