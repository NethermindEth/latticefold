#![allow(non_snake_case, clippy::upper_case_acronyms)]
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{marker::PhantomData, ops::MulAssign};
use lattirust_ring::{OverField, Ring};

use crate::{
    arith::{Witness, CCS, LCCCS},
    commitment::AjtaiCommitmentScheme,
    commitment::Commitment,
    decomposition_parameters::DecompositionParams,
    nifs::error::DecompositionError,
    transcript::Transcript,
};
use cyclotomic_rings::rings::SuitableRing;

#[derive(Clone, Debug, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct DecompositionProof<const C: usize, NTT: Ring> {
    pub u_s: Vec<Vec<NTT>>,
    pub v_s: Vec<Vec<NTT>>,
    pub x_s: Vec<Vec<NTT>>,
    pub y_s: Vec<Commitment<C, NTT>>,
}

pub trait DecompositionProver<NTT: SuitableRing, T: Transcript<NTT>> {
    fn prove<const W: usize, const C: usize, P: DecompositionParams>(
        cm_i: &LCCCS<C, NTT>,
        wit: &Witness<NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
        scheme: &AjtaiCommitmentScheme<C, W, NTT>,
    ) -> Result<
        (
            Vec<LCCCS<C, NTT>>,
            Vec<Witness<NTT>>,
            DecompositionProof<C, NTT>,
        ),
        DecompositionError,
    >
    where
        for<'a> <NTT as SuitableRing>::CoefficientRepresentation: MulAssign<&'a u128>;
}

pub trait DecompositionVerifier<NTT: OverField, T: Transcript<NTT>> {
    fn verify<const C: usize, P: DecompositionParams>(
        cm_i: &LCCCS<C, NTT>,
        proof: &DecompositionProof<C, NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<Vec<LCCCS<C, NTT>>, DecompositionError>;
}

pub struct LFDecompositionProver<NTT, T> {
    _ntt: PhantomData<NTT>,
    _t: PhantomData<T>,
}

pub struct LFDecompositionVerifier<NTT, T> {
    _ntt: PhantomData<NTT>,
    _t: PhantomData<T>,
}
