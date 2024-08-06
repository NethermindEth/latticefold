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

pub trait DecompositionParams: Clone {
    type AP: AjtaiParams;
    // the small b such that we decompose a witness into vectors of the norm b.
    const SMALL_B: u128;
    // k such that b^k = B
    const K: usize;
}

pub trait DecompositionProver<
    CR: PolyRing,
    NTT: OverField,
    P: DecompositionParams,
    T: Transcript<NTT>,
>
{
    type Proof: Clone;
    type Error: std::error::Error;

    fn prove(
        cm_i: &LCCCS<NTT, P::AP>,
        wit: &Witness<NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(Vec<LCCCS<NTT, P::AP>>, Vec<Witness<NTT>>, Self::Proof), Self::Error>;
}

pub trait DecompositionVerifier<
    CR: PolyRing,
    NTT: OverField,
    P: DecompositionParams,
    T: Transcript<NTT>,
>
{
    type Prover: DecompositionProver<CR, NTT, P, T>;
    type Error = <Self::Prover as DecompositionProver<CR, NTT, P, T>>::Error;

    fn verify(
        cm_i: &LCCCS<NTT, P::AP>,
        proof: &<Self::Prover as DecompositionProver<CR, NTT, P, T>>::Proof,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<Vec<LCCCS<NTT, P::AP>>, Self::Error>;
}

impl<
        CR: PolyRing,
        NTT: OverField,
        P: AjtaiParams,
        DP: DecompositionParams<AP = P>,
        T: Transcript<NTT>,
    > DecompositionProver<CR, NTT, DP, T> for NIFSProver<CR, NTT, P, DP, T>
{
    type Proof = DecompositionProof<NTT, P>;
    type Error = DecompositionError<NTT>;

    fn prove(
        cm_i: &LCCCS<NTT, P>,
        wit: &Witness<NTT>,
        transcript: &mut impl Transcript<NTT>,
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

impl<
        CR: PolyRing,
        NTT: OverField,
        P: AjtaiParams,
        DP: DecompositionParams<AP = P>,
        T: Transcript<NTT>,
    > DecompositionVerifier<CR, NTT, DP, T> for NIFSVerifier<CR, NTT, P, DP, T>
{
    type Prover = NIFSProver<CR, NTT, P, DP, T>;

    fn verify(
        _cm_i: &LCCCS<NTT, P>,
        _proof: &<Self::Prover as DecompositionProver<CR, NTT, DP, T>>::Proof,
        _transcript: &mut impl Transcript<NTT>,
        _ccs: &CCS<NTT>,
    ) -> Result<Vec<LCCCS<NTT, P>>, DecompositionError<NTT>> {
        todo!()
    }
}
