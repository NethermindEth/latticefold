use ark_std::marker::PhantomData;

use crate::{
    arith::{Witness, CCCS, CCS, LCCCS},
    nifs::error::LinearizationError,
    transcript::Transcript,
    utils::sumcheck,
};

use lattirust_ring::OverField;
#[derive(Clone)]
pub struct LinearizationProof<NTT: OverField> {
    // Sent in the step 2. of the linearization subprotocol
    pub linearization_sumcheck: sumcheck::Proof<NTT>,
    // Sent in the step 3.
    pub v: Vec<NTT>,
    pub u: Vec<NTT>,
}

pub trait LinearizationProver<NTT: OverField, T: Transcript<NTT>> {
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
