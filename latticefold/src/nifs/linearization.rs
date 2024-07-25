use ark_crypto_primitives::sponge::Absorb;
use ark_ff::Field;
use lattirust_arithmetic::{challenge_set::latticefold_challenge_set::OverField, ring::Ring};
use thiserror_no_std::Error;

use crate::{
    arith::{Witness, CCCS, LCCCS},
    transcript::Transcript,
    utils::sumcheck::{SumCheckError, SumCheckProof},
};

use super::{NIFSProver, NIFSVerifier};

#[derive(Debug, Error)]
pub enum LinearizationError<R: Ring> {
    SumCheckError(#[from] SumCheckError<R>),
}

pub trait LinearizationProver<F: Field, R: OverField<F>, T: Transcript<F, R>> {
    type Proof: Clone;
    fn prove(
        cm_i: &CCCS<R>,
        wit: &Witness<R>,
        transcript: &mut impl Transcript<F, R>,
    ) -> Result<(LCCCS<R>, Self::Proof), LinearizationError<R>>;
}

pub trait LinearizationVerifier<F: Field, R: OverField<F>, T: Transcript<F, R>> {
    type Prover: LinearizationProver<F, R, T>;

    fn verify(
        cm_i: &CCCS<R>,
        proof: &<Self::Prover as LinearizationProver<F, R, T>>::Proof,
        transcript: &mut impl Transcript<F, R>,
    ) -> Result<LCCCS<R>, LinearizationError<R>>;
}

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

impl<F: Field, R: OverField<F>, T: Transcript<F, R>> LinearizationProver<F, R, T>
    for NIFSProver<F, R, T>
where
    F: Absorb,
{
    type Proof = LinearizationProof<F, R>;

    fn prove(
        cm_i: &CCCS<R>,
        wit: &Witness<R>,
        transcript: &mut impl Transcript<F, R>,
    ) -> Result<(LCCCS<R>, Self::Proof), LinearizationError<R>> {
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
        cm_i: &CCCS<R>,
        proof: &<Self::Prover as LinearizationProver<F, R, T>>::Proof,
        transcript: &mut impl Transcript<F, R>,
    ) -> Result<LCCCS<R>, LinearizationError<R>> {
        todo!()
    }
}
