use std::marker::PhantomData;

use ark_crypto_primitives::sponge::Absorb;
use ark_ff::Field;
use lattirust_arithmetic::{challenge_set::latticefold_challenge_set::OverField, ring::Ring};
use thiserror_no_std::Error;

use crate::{
    arith::{Witness, CCCS, LCCCS},
    transcript::Transcript,
    utils::sumcheck::SumCheckProof,
};

use super::{NIFSProver, NIFSVerifier};

#[derive(Debug, Error)]
pub enum FoldingError<R: Ring> {
    PhantomRRemoveThisLater(R),
}

pub trait FoldingProver<F: Field, R: OverField<F>, T: Transcript<F, R>> {
    type Proof: Clone;
    fn prove(
        cm_i_s: &[LCCCS],
        w_s: &[Witness],
        transcript: &mut impl Transcript<F, R>,
    ) -> Result<(LCCCS, Witness, Self::Proof), FoldingError<R>>;
}

pub trait FoldingVerifier<F: Field, R: OverField<F>, T: Transcript<F, R>> {
    type Prover: FoldingProver<F, R, T>;

    fn verify(
        cm_i_s: &[LCCCS],
        proof: &<Self::Prover as FoldingProver<F, R, T>>::Proof,
        transcript: &mut impl Transcript<F, R>,
    ) -> Result<LCCCS, FoldingError<R>>;
}

#[derive(Clone)]
pub struct FoldingProof<F: Field, R: OverField<F>>
where
    F: Absorb,
{
    // Step 2.
    pub pointshift_sumcheck_proof: SumCheckProof<F, R>,
    // Step 3
    pub theta_s: Vec<R>,
    pub eta_s: Vec<R>,
}

impl<F: Field, R: OverField<F>, T: Transcript<F, R>> FoldingProver<F, R, T> for NIFSProver<F, R, T>
where
    F: Absorb,
{
    type Proof = FoldingProof<F, R>;

    fn prove(
        cm_i_s: &[LCCCS],
        w_s: &[Witness],
        transcript: &mut impl Transcript<F, R>,
    ) -> Result<(LCCCS, Witness, Self::Proof), FoldingError<R>> {
        todo!()
    }
}

impl<F: Field, R: OverField<F>, T: Transcript<F, R>> FoldingVerifier<F, R, T>
    for NIFSVerifier<F, R, T>
where
    F: Absorb,
{
    type Prover = NIFSProver<F, R, T>;

    fn verify(
        cm_i_s: &[LCCCS],
        proof: &<Self::Prover as FoldingProver<F, R, T>>::Proof,
        transcript: &mut impl Transcript<F, R>,
    ) -> Result<LCCCS, FoldingError<R>> {
        todo!()
    }
}
