use lattirust_arithmetic::challenge_set::latticefold_challenge_set::OverField;
use lattirust_arithmetic::ring::ConvertibleRing;

use crate::commitment::AjtaiParams;

use crate::{
    arith::{Witness, CCS, LCCCS},
    transcript::Transcript,
};

use super::{error::DecompositionError, NIFSProver, NIFSVerifier};

#[derive(Clone)]
pub struct DecompositionProof<NTT: OverField> {
    pub u_s: Vec<Vec<NTT>>,
    pub v_s: Vec<NTT>,
    pub x_s: Vec<Vec<NTT>>,
    pub y_s: Vec<Vec<NTT>>,
}

pub trait DecompositionProver<
    CR: ConvertibleRing,
    NTT: OverField,
    P: AjtaiParams,
    T: Transcript<NTT>,
>
{
    type Proof: Clone;
    type Error: std::error::Error;

    fn prove(
        cm_i: &LCCCS<CR, NTT, P>,
        wit: &Witness<CR, NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(Vec<LCCCS<CR, NTT, P>>, Vec<Witness<CR, NTT>>, Self::Proof), Self::Error>;
}

pub trait DecompositionVerifier<
    CR: ConvertibleRing,
    NTT: OverField,
    P: AjtaiParams,
    T: Transcript<NTT>,
>
{
    type Prover: DecompositionProver<CR, NTT, P, T>;
    type Error = <Self::Prover as DecompositionProver<CR, NTT, P, T>>::Error;

    fn verify(
        cm_i: &LCCCS<CR, NTT, P>,
        proof: &<Self::Prover as DecompositionProver<CR, NTT, P, T>>::Proof,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<Vec<LCCCS<CR, NTT, P>>, Self::Error>;
}

impl<CR: ConvertibleRing, NTT: OverField, P: AjtaiParams, T: Transcript<NTT>>
    DecompositionProver<CR, NTT, P, T> for NIFSProver<CR, NTT, P, T>
{
    type Proof = DecompositionProof<NTT>;
    type Error = DecompositionError<NTT>;

    fn prove(
        _cm_i: &LCCCS<CR, NTT, P>,
        _wit: &Witness<CR, NTT>,
        _transcript: &mut impl Transcript<NTT>,
        _ccs: &CCS<NTT>,
    ) -> Result<
        (
            Vec<LCCCS<CR, NTT, P>>,
            Vec<Witness<CR, NTT>>,
            DecompositionProof<NTT>,
        ),
        DecompositionError<NTT>,
    > {
        todo!()
    }
}

impl<CR: ConvertibleRing, NTT: OverField, P: AjtaiParams, T: Transcript<NTT>>
    DecompositionVerifier<CR, NTT, P, T> for NIFSVerifier<CR, NTT, P, T>
{
    type Prover = NIFSProver<CR, NTT, P, T>;

    fn verify(
        _cm_i: &LCCCS<CR, NTT, P>,
        _proof: &<Self::Prover as DecompositionProver<CR, NTT, P, T>>::Proof,
        _transcript: &mut impl Transcript<NTT>,
        _ccs: &CCS<NTT>,
    ) -> Result<Vec<LCCCS<CR, NTT, P>>, DecompositionError<NTT>> {
        todo!()
    }
}
