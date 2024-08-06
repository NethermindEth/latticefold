use crate::commitment::AjtaiParams;
use lattirust_arithmetic::challenge_set::latticefold_challenge_set::OverField;
use lattirust_arithmetic::ring::ConvertibleRing;

use crate::{
    arith::{Witness, CCS, LCCCS},
    transcript::Transcript,
    utils::sumcheck::SumCheckProof,
};

use super::{error::FoldingError, NIFSProver, NIFSVerifier};

#[derive(Clone)]
pub struct FoldingProof<NTT: OverField> {
    // Step 2.
    pub pointshift_sumcheck_proof: SumCheckProof<NTT>,
    // Step 3
    pub theta_s: Vec<NTT>,
    pub eta_s: Vec<NTT>,
}

pub trait FoldingProver<CR: ConvertibleRing, NTT: OverField, P: AjtaiParams<CR>, T: Transcript<NTT>>
{
    type Proof: Clone;
    type Error: std::error::Error;

    fn prove(
        cm_i_s: &[LCCCS<CR, NTT, P>],
        w_s: &[Witness<CR, NTT>],
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(LCCCS<CR, NTT, P>, Witness<CR, NTT>, Self::Proof), Self::Error>;
}

pub trait FoldingVerifier<
    CR: ConvertibleRing,
    NTT: OverField,
    P: AjtaiParams<CR>,
    T: Transcript<NTT>,
>
{
    type Prover: FoldingProver<CR, NTT, P, T>;
    type Error = <Self::Prover as FoldingProver<CR, NTT, P, T>>::Error;

    fn verify(
        cm_i_s: &[LCCCS<CR, NTT, P>],
        proof: &<Self::Prover as FoldingProver<CR, NTT, P, T>>::Proof,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<LCCCS<CR, NTT, P>, Self::Error>;
}

impl<CR: ConvertibleRing, NTT: OverField, P: AjtaiParams<CR>, T: Transcript<NTT>>
    FoldingProver<CR, NTT, P, T> for NIFSProver<CR, NTT, P, T>
{
    type Proof = FoldingProof<NTT>;
    type Error = FoldingError<NTT>;

    fn prove(
        _cm_i_s: &[LCCCS<CR, NTT, P>],
        _w_s: &[Witness<CR, NTT>],
        _transcript: &mut impl Transcript<NTT>,
        _ccs: &CCS<NTT>,
    ) -> Result<(LCCCS<CR, NTT, P>, Witness<CR, NTT>, FoldingProof<NTT>), FoldingError<NTT>> {
        todo!()
    }
}

impl<CR: ConvertibleRing, NTT: OverField, P: AjtaiParams<CR>, T: Transcript<NTT>>
    FoldingVerifier<CR, NTT, P, T> for NIFSVerifier<CR, NTT, P, T>
{
    type Prover = NIFSProver<CR, NTT, P, T>;

    fn verify(
        _cm_i_s: &[LCCCS<CR, NTT, P>],
        _proof: &<Self::Prover as FoldingProver<CR, NTT, P, T>>::Proof,
        _transcript: &mut impl Transcript<NTT>,
        _ccs: &CCS<NTT>,
    ) -> Result<LCCCS<CR, NTT, P>, FoldingError<NTT>> {
        todo!()
    }
}
