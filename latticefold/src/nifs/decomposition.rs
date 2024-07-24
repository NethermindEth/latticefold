use std::marker::PhantomData;

use ark_crypto_primitives::sponge::Absorb;
use ark_ff::Field;
use lattirust_arithmetic::{challenge_set::latticefold_challenge_set::OverField, ring::Ring};
use thiserror_no_std::Error;

use crate::{
    arith::{Witness, CCCS, LCCCS},
    transcript::Transcript,
};

use super::{NIFSProver, NIFSVerifier};

#[derive(Debug, Error)]
pub enum DecompositionError<R: Ring> {
    PhantomRRemoveThisLater(R),
    IncorrectLength,
}

pub trait DecompositionProver<F: Field, R: OverField<F>, T: Transcript<F, R>> {
    type Proof: Clone;
    fn prove(
        cm_i: &LCCCS,
        wit: &Witness,
        transcript: &mut impl Transcript<F, R>,
    ) -> Result<(Vec<(LCCCS, Witness)>, Self::Proof), DecompositionError<R>>;
}

pub trait DecompositionVerifier<F: Field, R: OverField<F>, T: Transcript<F, R>> {
    type Prover: DecompositionProver<F, R, T>;

    fn verify(
        cm_i: &LCCCS,
        proof: &<Self::Prover as DecompositionProver<F, R, T>>::Proof,
        transcript: &mut impl Transcript<F, R>,
    ) -> Result<Vec<LCCCS>, DecompositionError<R>>;
}

#[derive(Clone)]
pub struct DecompositionProof<F: Field, R: OverField<F>>
where
    F: Absorb,
{
    _f: PhantomData<F>,
    pub u_s: Vec<Vec<R>>,
    pub v_s: Vec<R>,
    pub x_s: Vec<Vec<R>>,
    pub y_s: Vec<Vec<R>>,
}

impl<F: Field, R: OverField<F>, T: Transcript<F, R>> DecompositionProver<F, R, T>
    for NIFSProver<F, R, T>
where
    F: Absorb,
{
    type Proof = DecompositionProof<F, R>;

    fn prove(
        cm_i: &LCCCS,
        wit: &Witness,
        transcript: &mut impl Transcript<F, R>,
    ) -> Result<(Vec<(LCCCS, Witness)>, Self::Proof), DecompositionError<R>> {
        todo!()
    }
}

impl<F: Field, R: OverField<F>, T: Transcript<F, R>> DecompositionVerifier<F, R, T>
    for NIFSVerifier<F, R, T>
where
    F: Absorb,
{
    type Prover = NIFSProver<F, R, T>;

    fn verify(
        cm_i: &LCCCS,
        proof: &<Self::Prover as DecompositionProver<F, R, T>>::Proof,
        transcript: &mut impl Transcript<F, R>,
    ) -> Result<Vec<LCCCS>, DecompositionError<R>> {
        todo!()
    }
}
