use ark_std::marker::PhantomData;

use lattirust_arithmetic::challenge_set::latticefold_challenge_set::OverField;
use lattirust_arithmetic::ring::PolyRing;

use super::error::FoldingError;
use crate::{
    arith::{Witness, CCS, LCCCS},
    parameters::DecompositionParams,
    transcript::Transcript,
    utils::sumcheck,
};

#[derive(Clone)]
pub struct FoldingProof<NTT: OverField> {
    // Step 2.
    pub pointshift_sumcheck_proof: sumcheck::Proof<NTT>,
    // Step 3
    pub theta_s: Vec<NTT>,
    pub eta_s: Vec<NTT>,
}

pub trait FoldingProver<NTT: OverField, T: Transcript<NTT>> {
    fn prove<const C: usize, CR: PolyRing, P: DecompositionParams>(
        cm_i_s: &[LCCCS<C, NTT>],
        w_s: &[Witness<NTT>],
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(LCCCS<C, NTT>, Witness<NTT>, FoldingProof<NTT>), FoldingError<NTT>>;
}

pub trait FoldingVerifier<NTT: OverField, T: Transcript<NTT>> {
    fn verify<const C: usize, P: DecompositionParams>(
        cm_i_s: &[LCCCS<C, NTT>],
        proof: &FoldingProof<NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<LCCCS<C, NTT>, FoldingError<NTT>>;
}

pub struct LFFoldingProver<NTT, T> {
    _ntt: PhantomData<NTT>,
    _t: PhantomData<T>,
}

pub struct LFFoldingVerifier<NTT, T> {
    _ntt: PhantomData<NTT>,
    _t: PhantomData<T>,
}

impl<NTT: OverField, T: Transcript<NTT>> FoldingProver<NTT, T> for LFFoldingProver<NTT, T> {
    fn prove<const C: usize, CR: PolyRing, P: DecompositionParams>(
        _cm_i_s: &[LCCCS<C, NTT>],
        _w_s: &[Witness<NTT>],
        _transcript: &mut impl Transcript<NTT>,
        _ccs: &CCS<NTT>,
    ) -> Result<(LCCCS<C, NTT>, Witness<NTT>, FoldingProof<NTT>), FoldingError<NTT>> {
        todo!()
    }
}

impl<NTT: OverField, T: Transcript<NTT>> FoldingVerifier<NTT, T> for LFFoldingVerifier<NTT, T> {
    fn verify<const C: usize, P: DecompositionParams>(
        _cm_i_s: &[LCCCS<C, NTT>],
        _proof: &FoldingProof<NTT>,
        _transcript: &mut impl Transcript<NTT>,
        _ccs: &CCS<NTT>,
    ) -> Result<LCCCS<C, NTT>, FoldingError<NTT>> {
        todo!()
    }
}
