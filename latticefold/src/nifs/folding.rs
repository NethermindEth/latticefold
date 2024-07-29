use lattirust_arithmetic::challenge_set::latticefold_challenge_set::OverField;

use crate::{
    arith::{Witness, CCS, LCCCS},
    transcript::Transcript,
    utils::sumcheck::SumCheckProof,
};

use super::{error::FoldingError, NIFSProver, NIFSVerifier};

#[derive(Clone)]
pub struct FoldingProof<R: OverField> {
    // Step 2.
    pub pointshift_sumcheck_proof: SumCheckProof<R>,
    // Step 3
    pub theta_s: Vec<R>,
    pub eta_s: Vec<R>,
}

pub trait FoldingProver<R: OverField, T: Transcript<R>> {
    type Proof: Clone;
    type Error: std::error::Error;

    fn prove(
        cm_i_s: &[LCCCS<R>],
        w_s: &[Witness<R>],
        transcript: &mut impl Transcript<R>,
        ccs: &CCS<R>,
    ) -> Result<(LCCCS<R>, Witness<R>, Self::Proof), Self::Error>;
}

pub trait FoldingVerifier<R: OverField, T: Transcript<R>> {
    type Prover: FoldingProver<R, T>;
    type Error = <Self::Prover as FoldingProver<R, T>>::Error;

    fn verify(
        cm_i_s: &[LCCCS<R>],
        proof: &<Self::Prover as FoldingProver<R, T>>::Proof,
        transcript: &mut impl Transcript<R>,
        ccs: &CCS<R>,
    ) -> Result<LCCCS<R>, Self::Error>;
}

impl<R: OverField, T: Transcript<R>> FoldingProver<R, T> for NIFSProver<R, T> {
    type Proof = FoldingProof<R>;
    type Error = FoldingError<R>;

    fn prove(
        _cm_i_s: &[LCCCS<R>],
        _w_s: &[Witness<R>],
        _transcript: &mut impl Transcript<R>,
        _ccs: &CCS<R>,
    ) -> Result<(LCCCS<R>, Witness<R>, FoldingProof<R>), FoldingError<R>> {
        todo!()
    }
}

impl<R: OverField, T: Transcript<R>> FoldingVerifier<R, T> for NIFSVerifier<R, T> {
    type Prover = NIFSProver<R, T>;

    fn verify(
        _cm_i_s: &[LCCCS<R>],
        _proof: &<Self::Prover as FoldingProver<R, T>>::Proof,
        _transcript: &mut impl Transcript<R>,
        _ccs: &CCS<R>,
    ) -> Result<LCCCS<R>, FoldingError<R>> {
        todo!()
    }
}
