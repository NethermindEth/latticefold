use lattirust_arithmetic::challenge_set::latticefold_challenge_set::OverField;

use crate::{
    arith::{Witness, CCS, LCCCS},
    transcript::Transcript,
};

use super::{error::DecompositionError, NIFSProver, NIFSVerifier};

#[derive(Clone)]
pub struct DecompositionProof<R: OverField> {
    pub u_s: Vec<Vec<R>>,
    pub v_s: Vec<R>,
    pub x_s: Vec<Vec<R>>,
    pub y_s: Vec<Vec<R>>,
}

pub trait DecompositionProver<R: OverField, T: Transcript<R>> {
    type Proof: Clone;
    type Error: std::error::Error;

    fn prove(
        cm_i: &LCCCS<R>,
        wit: &Witness<R>,
        transcript: &mut impl Transcript<R>,
        ccs: &CCS<R>,
    ) -> Result<(Vec<LCCCS<R>>, Vec<Witness<R>>, Self::Proof), Self::Error>;
}

pub trait DecompositionVerifier<R: OverField, T: Transcript<R>> {
    type Prover: DecompositionProver<R, T>;
    type Error = <Self::Prover as DecompositionProver<R, T>>::Error;

    fn verify(
        cm_i: &LCCCS<R>,
        proof: &<Self::Prover as DecompositionProver<R, T>>::Proof,
        transcript: &mut impl Transcript<R>,
        ccs: &CCS<R>,
    ) -> Result<Vec<LCCCS<R>>, Self::Error>;
}

impl<R: OverField, T: Transcript<R>> DecompositionProver<R, T> for NIFSProver<R, T> {
    type Proof = DecompositionProof<R>;
    type Error = DecompositionError<R>;

    fn prove(
        _cm_i: &LCCCS<R>,
        _wit: &Witness<R>,
        _transcript: &mut impl Transcript<R>,
        _ccs: &CCS<R>,
    ) -> Result<(Vec<LCCCS<R>>, Vec<Witness<R>>, DecompositionProof<R>), DecompositionError<R>>
    {
        todo!()
    }
}

impl<R: OverField, T: Transcript<R>> DecompositionVerifier<R, T> for NIFSVerifier<R, T> {
    type Prover = NIFSProver<R, T>;

    fn verify(
        _cm_i: &LCCCS<R>,
        _proof: &<Self::Prover as DecompositionProver<R, T>>::Proof,
        _transcript: &mut impl Transcript<R>,
        _ccs: &CCS<R>,
    ) -> Result<Vec<LCCCS<R>>, DecompositionError<R>> {
        todo!()
    }
}
