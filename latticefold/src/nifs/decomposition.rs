use lattirust_arithmetic::{challenge_set::latticefold_challenge_set::OverField, ring::PolyRing};
use std::marker::PhantomData;

use super::error::DecompositionError;
use crate::{
    arith::{Witness, CCS, LCCCS},
    commitment::Commitment,
    parameters::DecompositionParams,
    transcript::Transcript,
};

#[derive(Clone)]
pub struct DecompositionProof<const C: usize, NTT: OverField> {
    pub u_s: Vec<Vec<NTT>>,
    pub v_s: Vec<NTT>,
    pub x_s: Vec<Vec<NTT>>,
    pub y_s: Vec<Commitment<C, NTT>>,
}

pub trait DecompositionProver<NTT: OverField, T: Transcript<NTT>> {
    fn prove<const C: usize, CR: PolyRing, P: DecompositionParams>(
        cm_i: &LCCCS<C, NTT>,
        wit: &Witness<NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<
        (
            Vec<LCCCS<C, NTT>>,
            Vec<Witness<NTT>>,
            DecompositionProof<C, NTT>,
        ),
        DecompositionError<NTT>,
    >;
}

pub trait DecompositionVerifier<NTT: OverField, T: Transcript<NTT>> {
    fn verify<const C: usize, P: DecompositionParams>(
        cm_i: &LCCCS<C, NTT>,
        proof: &DecompositionProof<C, NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<Vec<LCCCS<C, NTT>>, DecompositionError<NTT>>;
}

pub struct LFDecompositionProver<NTT, T> {
    _ntt: PhantomData<NTT>,
    _t: PhantomData<T>,
}

pub struct LFDecompositionVerifier<NTT, T> {
    _ntt: PhantomData<NTT>,
    _t: PhantomData<T>,
}

impl<NTT: OverField, T: Transcript<NTT>> DecompositionProver<NTT, T>
    for LFDecompositionProver<NTT, T>
{
    fn prove<const C: usize, CR: PolyRing, P: DecompositionParams>(
        _cm_i: &LCCCS<C, NTT>,
        _wit: &Witness<NTT>,
        _transcript: &mut impl Transcript<NTT>,
        _ccs: &CCS<NTT>,
    ) -> Result<
        (
            Vec<LCCCS<C, NTT>>,
            Vec<Witness<NTT>>,
            DecompositionProof<C, NTT>,
        ),
        DecompositionError<NTT>,
    > {
        todo!()
    }
}

impl<NTT: OverField, T: Transcript<NTT>> DecompositionVerifier<NTT, T>
    for LFDecompositionVerifier<NTT, T>
{
    fn verify<const C: usize, P: DecompositionParams>(
        _cm_i: &LCCCS<C, NTT>,
        _proof: &DecompositionProof<C, NTT>,
        _transcript: &mut impl Transcript<NTT>,
        _ccs: &CCS<NTT>,
    ) -> Result<Vec<LCCCS<C, NTT>>, DecompositionError<NTT>> {
        todo!()
    }
}
