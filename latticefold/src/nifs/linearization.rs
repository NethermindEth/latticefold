use ark_crypto_primitives::sponge::Absorb;
use ark_ff::Field;
use lattirust_arithmetic::challenge_set::latticefold_challenge_set::OverField;

use crate::{
    arith::{Witness, CCCS, CCS, LCCCS},
    transcript::Transcript,
    utils::sumcheck::SumCheckProof,
};

use super::{error::LinearizationError, NIFSProver, NIFSVerifier};

#[derive(Clone)]
pub struct LinearizationProof<F: Field, R: OverField<F>>
where
    F: Absorb,
{
    // Sent in the step 2. of the linearization subprotocol
    pub linearization_sumcheck: SumCheckProof<F, R>,
    // Sent in the step 3.
    pub v: R,
    pub u: Vec<R>,
}

pub trait LinearizationProver<F: Field, R: OverField<F>, T: Transcript<F, R>> {
    type Proof: Clone;
    type Error: std::error::Error;

    fn prove(
        cm_i: &CCCS<R>,
        wit: &Witness<R>,
        transcript: &mut impl Transcript<F, R>,
        ccs: &CCS<R>,
    ) -> Result<(LCCCS<R>, Self::Proof), Self::Error>;
}

pub trait LinearizationVerifier<F: Field, R: OverField<F>, T: Transcript<F, R>> {
    type Prover: LinearizationProver<F, R, T>;
    type Error = <Self::Prover as LinearizationProver<F, R, T>>::Error;

    fn verify(
        cm_i: &CCCS<R>,
        proof: &<Self::Prover as LinearizationProver<F, R, T>>::Proof,
        transcript: &mut impl Transcript<F, R>,
        ccs: &CCS<R>,
    ) -> Result<LCCCS<R>, Self::Error>;
}

impl<F: Field, R: OverField<F>, T: Transcript<F, R>> LinearizationProver<F, R, T>
    for NIFSProver<F, R, T>
where
    F: Absorb,
{
    type Proof = LinearizationProof<F, R>;
    type Error = LinearizationError<R>;

    fn prove(
        _cm_i: &CCCS<R>,
        _wit: &Witness<R>,
        _transcript: &mut impl Transcript<F, R>,
        _ccs: &CCS<R>,
    ) -> Result<(LCCCS<R>, LinearizationProof<F, R>), LinearizationError<R>> {
        todo!()
    }
}

impl<F: Field, R: OverField<F>, T: Transcript<F, R>> LinearizationVerifier<F, R, T>
    for NIFSVerifier<F, R, T>
where
    F: Absorb,
{
    type Prover = NIFSProver<F, R, T>;

    fn verify(
        _cm_i: &CCCS<R>,
        _proof: &<Self::Prover as LinearizationProver<F, R, T>>::Proof,
        _transcript: &mut impl Transcript<F, R>,
        _ccs: &CCS<R>,
    ) -> Result<LCCCS<R>, LinearizationError<R>> {
        todo!()
    }
}
