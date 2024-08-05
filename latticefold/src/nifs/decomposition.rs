use lattirust_arithmetic::challenge_set::latticefold_challenge_set::OverField;
use lattirust_arithmetic::ring::ConvertibleRing;

use crate::commitment::AjtaiParams;

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

pub trait DecompositionProver<
    CR: ConvertibleRing,
    R: OverField,
    P: AjtaiParams<CR>,
    T: Transcript<R>,
>
{
    type Proof: Clone;
    type Error: std::error::Error;

    fn prove(
        cm_i: &LCCCS<CR, R, P>,
        wit: &Witness<CR, R>,
        transcript: &mut impl Transcript<R>,
        ccs: &CCS<R>,
    ) -> Result<(Vec<LCCCS<CR, R, P>>, Vec<Witness<CR, R>>, Self::Proof), Self::Error>;
}

pub trait DecompositionVerifier<
    CR: ConvertibleRing,
    R: OverField,
    P: AjtaiParams<CR>,
    T: Transcript<R>,
>
{
    type Prover: DecompositionProver<CR, R, P, T>;
    type Error = <Self::Prover as DecompositionProver<CR, R, P, T>>::Error;

    fn verify(
        cm_i: &LCCCS<CR, R, P>,
        proof: &<Self::Prover as DecompositionProver<CR, R, P, T>>::Proof,
        transcript: &mut impl Transcript<R>,
        ccs: &CCS<R>,
    ) -> Result<Vec<LCCCS<CR, R, P>>, Self::Error>;
}

impl<CR: ConvertibleRing, R: OverField, P: AjtaiParams<CR>, T: Transcript<R>>
    DecompositionProver<CR, R, P, T> for NIFSProver<CR, R, P, T>
{
    type Proof = DecompositionProof<R>;
    type Error = DecompositionError<R>;

    fn prove(
        _cm_i: &LCCCS<CR, R, P>,
        _wit: &Witness<CR, R>,
        _transcript: &mut impl Transcript<R>,
        _ccs: &CCS<R>,
    ) -> Result<
        (
            Vec<LCCCS<CR, R, P>>,
            Vec<Witness<CR, R>>,
            DecompositionProof<R>,
        ),
        DecompositionError<R>,
    > {
        todo!()
    }
}

impl<CR: ConvertibleRing, R: OverField, P: AjtaiParams<CR>, T: Transcript<R>>
    DecompositionVerifier<CR, R, P, T> for NIFSVerifier<CR, R, P, T>
{
    type Prover = NIFSProver<CR, R, P, T>;

    fn verify(
        _cm_i: &LCCCS<CR, R, P>,
        _proof: &<Self::Prover as DecompositionProver<CR, R, P, T>>::Proof,
        _transcript: &mut impl Transcript<R>,
        _ccs: &CCS<R>,
    ) -> Result<Vec<LCCCS<CR, R, P>>, DecompositionError<R>> {
        todo!()
    }
}
