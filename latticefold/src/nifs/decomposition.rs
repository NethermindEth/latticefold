use lattirust_arithmetic::challenge_set::latticefold_challenge_set::OverField;
use lattirust_arithmetic::ring::PolyRing;

use crate::commitment::{AjtaiParams, Commitment};

use crate::{
    arith::{Witness, CCS, LCCCS},
    transcript::Transcript,
};

use super::{error::DecompositionError, NIFSProver, NIFSVerifier};

#[derive(Clone)]
pub struct DecompositionProof<NTT: OverField, P: AjtaiParams> {
    pub u_s: Vec<Vec<NTT>>,
    pub v_s: Vec<NTT>,
    pub x_s: Vec<Vec<NTT>>,
    pub y_s: Vec<Commitment<NTT, P>>,
}

pub trait DecompositionProver<CR: PolyRing, NTT: OverField, P: AjtaiParams, T: Transcript<NTT>> {
    type Proof: Clone;
    type Error: std::error::Error;

    fn prove(
        cm_i: &LCCCS<NTT, P>,
        wit: &Witness<NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(Vec<LCCCS<NTT, P>>, Vec<Witness<NTT>>, Self::Proof), Self::Error>;
}

pub trait DecompositionVerifier<CR: PolyRing, NTT: OverField, P: AjtaiParams, T: Transcript<NTT>> {
    type Prover: DecompositionProver<CR, NTT, P, T>;
    type Error = <Self::Prover as DecompositionProver<CR, NTT, P, T>>::Error;

    fn verify(
        cm_i: &LCCCS<NTT, P>,
        proof: &<Self::Prover as DecompositionProver<CR, NTT, P, T>>::Proof,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<Vec<LCCCS<NTT, P>>, Self::Error>;
}

impl<CR: PolyRing, NTT: OverField, P: AjtaiParams, T: Transcript<NTT>>
    DecompositionProver<CR, NTT, P, T> for NIFSProver<CR, NTT, P, T>
{
    type Proof = DecompositionProof<NTT, P>;
    type Error = DecompositionError<NTT>;

    fn prove(
        _cm_i: &LCCCS<NTT, P>,
        _wit: &Witness<NTT>,
        _transcript: &mut impl Transcript<NTT>,
        _ccs: &CCS<NTT>,
    ) -> Result<
        (
            Vec<LCCCS<NTT, P>>,
            Vec<Witness<NTT>>,
            DecompositionProof<NTT, P>,
        ),
        DecompositionError<NTT>,
    > {
        todo!()
    }
}

impl<CR: PolyRing, NTT: OverField, P: AjtaiParams, T: Transcript<NTT>>
    DecompositionVerifier<CR, NTT, P, T> for NIFSVerifier<CR, NTT, P, T>
{
    type Prover = NIFSProver<CR, NTT, P, T>;

    fn verify(
        _cm_i: &LCCCS<NTT, P>,
        _proof: &<Self::Prover as DecompositionProver<CR, NTT, P, T>>::Proof,
        _transcript: &mut impl Transcript<NTT>,
        _ccs: &CCS<NTT>,
    ) -> Result<Vec<LCCCS<NTT, P>>, DecompositionError<NTT>> {
        todo!()
    }
}
