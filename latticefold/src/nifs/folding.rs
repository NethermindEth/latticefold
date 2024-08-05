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
pub struct FoldingProof<R: OverField> {
    // Step 2.
    pub pointshift_sumcheck_proof: SumCheckProof<R>,
    // Step 3
    pub theta_s: Vec<R>,
    pub eta_s: Vec<R>,
}

pub trait FoldingProver<CR: ConvertibleRing, R: OverField, P: AjtaiParams<CR>, T: Transcript<R>> {
    type Proof: Clone;
    type Error: std::error::Error;

    fn prove(
        cm_i_s: &[LCCCS<CR, R, P>],
        w_s: &[Witness<CR, R>],
        transcript: &mut impl Transcript<R>,
        ccs: &CCS<R>,
    ) -> Result<(LCCCS<CR, R, P>, Witness<CR, R>, Self::Proof), Self::Error>;
}

pub trait FoldingVerifier<CR: ConvertibleRing, R: OverField, P: AjtaiParams<CR>, T: Transcript<R>> {
    type Prover: FoldingProver<CR, R, P, T>;
    type Error = <Self::Prover as FoldingProver<CR, R, P, T>>::Error;

    fn verify(
        cm_i_s: &[LCCCS<CR, R, P>],
        proof: &<Self::Prover as FoldingProver<CR, R, P, T>>::Proof,
        transcript: &mut impl Transcript<R>,
        ccs: &CCS<R>,
    ) -> Result<LCCCS<CR, R, P>, Self::Error>;
}

impl<CR: ConvertibleRing, R: OverField, P: AjtaiParams<CR>, T: Transcript<R>>
    FoldingProver<CR, R, P, T> for NIFSProver<CR, R, P, T>
{
    type Proof = FoldingProof<R>;
    type Error = FoldingError<R>;

    fn prove(
        _cm_i_s: &[LCCCS<CR, R, P>],
        _w_s: &[Witness<CR, R>],
        _transcript: &mut impl Transcript<R>,
        _ccs: &CCS<R>,
    ) -> Result<(LCCCS<CR, R, P>, Witness<CR, R>, FoldingProof<R>), FoldingError<R>> {
        todo!()
    }
}

impl<CR: ConvertibleRing, R: OverField, P: AjtaiParams<CR>, T: Transcript<R>>
    FoldingVerifier<CR, R, P, T> for NIFSVerifier<CR, R, P, T>
{
    type Prover = NIFSProver<CR, R, P, T>;

    fn verify(
        _cm_i_s: &[LCCCS<CR, R, P>],
        _proof: &<Self::Prover as FoldingProver<CR, R, P, T>>::Proof,
        _transcript: &mut impl Transcript<R>,
        _ccs: &CCS<R>,
    ) -> Result<LCCCS<CR, R, P>, FoldingError<R>> {
        todo!()
    }
}
