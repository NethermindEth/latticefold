use crate::commitment::AjtaiParams;
use lattirust_arithmetic::challenge_set::latticefold_challenge_set::OverField;
use lattirust_arithmetic::ring::PolyRing;

use crate::{
    arith::{Witness, CCS, LCCCS},
    transcript::Transcript,
    utils::sumcheck::SumCheckProof,
};

use super::{decomposition::DecompositionParams, error::FoldingError, NIFSProver, NIFSVerifier};

#[derive(Clone)]
pub struct FoldingProof<NTT: OverField> {
    // Step 2.
    pub pointshift_sumcheck_proof: SumCheckProof<NTT>,
    // Step 3
    pub theta_s: Vec<NTT>,
    pub eta_s: Vec<NTT>,
}

pub trait FoldingProver<CR: PolyRing, NTT: OverField, P: AjtaiParams, T: Transcript<NTT>> {
    type Proof: Clone;
    type Error: std::error::Error;

    fn prove(
        cm_i_s: &[LCCCS<NTT, P>],
        w_s: &[Witness<NTT>],
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(LCCCS<NTT, P>, Witness<NTT>, Self::Proof), Self::Error>;
}

pub trait FoldingVerifier<CR: PolyRing, NTT: OverField, P: AjtaiParams, T: Transcript<NTT>> {
    type Prover: FoldingProver<CR, NTT, P, T>;
    type Error = <Self::Prover as FoldingProver<CR, NTT, P, T>>::Error;

    fn verify(
        cm_i_s: &[LCCCS<NTT, P>],
        proof: &<Self::Prover as FoldingProver<CR, NTT, P, T>>::Proof,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<LCCCS<NTT, P>, Self::Error>;
}

impl<CR: PolyRing, NTT: OverField, P: AjtaiParams, DP: DecompositionParams, T: Transcript<NTT>>
    FoldingProver<CR, NTT, P, T> for NIFSProver<CR, NTT, P, DP, T>
{
    type Proof = FoldingProof<NTT>;
    type Error = FoldingError<NTT>;

    fn prove(
        _cm_i_s: &[LCCCS<NTT, P>],
        _w_s: &[Witness<NTT>],
        _transcript: &mut impl Transcript<NTT>,
        _ccs: &CCS<NTT>,
    ) -> Result<(LCCCS<NTT, P>, Witness<NTT>, FoldingProof<NTT>), FoldingError<NTT>> {
        todo!()
    }
}

impl<CR: PolyRing, NTT: OverField, P: AjtaiParams, DP: DecompositionParams, T: Transcript<NTT>>
    FoldingVerifier<CR, NTT, P, T> for NIFSVerifier<CR, NTT, P, DP, T>
{
    type Prover = NIFSProver<CR, NTT, P, DP, T>;

    fn verify(
        _cm_i_s: &[LCCCS<NTT, P>],
        _proof: &<Self::Prover as FoldingProver<CR, NTT, P, T>>::Proof,
        _transcript: &mut impl Transcript<NTT>,
        _ccs: &CCS<NTT>,
    ) -> Result<LCCCS<NTT, P>, FoldingError<NTT>> {
        todo!()
    }
}
