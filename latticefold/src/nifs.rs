use ark_std::marker::PhantomData;

use cyclotomic_rings::rings::SuitableRing;
use lattirust_ring::OverField;

use crate::{
    arith::{Witness, CCCS, CCS, LCCCS},
    commitment::AjtaiCommitmentScheme,
    decomposition_parameters::DecompositionParams,
    transcript::TranscriptWithShortChallenges,
};
use decomposition::*;
use error::LatticefoldError;
use folding::{FoldingProof, FoldingProver, FoldingVerifier, LFFoldingProver, LFFoldingVerifier};
use linearization::{
    LFLinearizationProver, LFLinearizationVerifier, LinearizationProof, LinearizationProver,
    LinearizationVerifier,
};
pub mod decomposition;
pub mod error;
pub mod folding;
pub mod linearization;
mod mle_helpers;

/// `C` is the length of Ajtai commitment vectors.
/// `NTT` is a cyclotomic ring in the NTT form.
#[derive(Clone)]
pub struct LFProof<const C: usize, NTT: OverField> {
    pub linearization_proof: LinearizationProof<NTT>,
    pub decomposition_proof_l: DecompositionProof<C, NTT>,
    pub decomposition_proof_r: DecompositionProof<C, NTT>,
    pub folding_proof: FoldingProof<NTT>,
}

/// `C` is the length of commitment vectors or, equivalently, the number of rows of the Ajtai matrix.
/// `W` is the length of witness vectors or, equivalently, the number of columns of the Ajtai matrix.
/// `NTT` is a suitable cyclotomic ring.
/// `P` is the decomposition parameters.
/// `T` is the FS-transform transcript.
pub struct NIFSProver<const C: usize, const W: usize, NTT, P, T> {
    _r: PhantomData<NTT>,
    _p: PhantomData<P>,
    _t: PhantomData<T>,
}

impl<
        const C: usize,
        const W: usize,
        NTT: SuitableRing,
        P: DecompositionParams,
        T: TranscriptWithShortChallenges<NTT>,
    > NIFSProver<C, W, NTT, P, T>
{
    pub fn prove(
        r_l: Vec<NTT>,
        acc: &LCCCS<C, NTT>,
        w_acc: &Witness<NTT>,
        cm_i: &CCCS<C, NTT>,
        w_i: &Witness<NTT>,
        transcript: &mut impl TranscriptWithShortChallenges<NTT>,
        ccs: &CCS<NTT>,
        scheme: &AjtaiCommitmentScheme<C, W, NTT>,
    ) -> Result<(Vec<NTT>, LCCCS<C, NTT>, Witness<NTT>, LFProof<C, NTT>), LatticefoldError<NTT>>
    {
        let (r_r, linearized_cm_i, linearization_proof) =
            LFLinearizationProver::<_, T>::prove(cm_i, w_i, transcript, ccs)?;
        let (r_l, decomposed_lcccs_l, decomposed_wit_l, decomposition_proof_l) =
            LFDecompositionProver::<_, T>::prove::<W, C, P>(
                r_l, acc, w_acc, transcript, ccs, scheme,
            )?;
        let (r_r, decomposed_lcccs_r, decomposed_wit_r, decomposition_proof_r) =
            LFDecompositionProver::<_, T>::prove::<W, C, P>(
                r_r,
                &linearized_cm_i,
                w_i,
                transcript,
                ccs,
                scheme,
            )?;

        let (lcccs, wit_s) = {
            let mut lcccs = decomposed_lcccs_l;
            let mut lcccs_r = decomposed_lcccs_r;
            lcccs.append(&mut lcccs_r);

            let mut wit_s = decomposed_wit_l;
            let mut wit_s_r = decomposed_wit_r;
            wit_s.append(&mut wit_s_r);

            (lcccs, wit_s)
        };

        let (r_0, folded_lcccs, wit, folding_proof) =
            LFFoldingProver::<_, T>::prove::<C, P>((r_l, r_r), &lcccs, &wit_s, transcript, ccs)?;

        Ok((
            r_0,
            folded_lcccs,
            wit,
            LFProof {
                linearization_proof,
                decomposition_proof_l,
                decomposition_proof_r,
                folding_proof,
            },
        ))
    }
}

/// `C` is the length of commitment vectors or, equivalently, the number of rows of the Ajtai matrix.
/// `W` is the length of witness vectors or, equivalently, the number of columns of the Ajtai matrix.
/// `NTT` is a suitable cyclotomic ring.
/// `P` is the decomposition parameters.
/// `T` is the FS-transform transcript.
pub struct NIFSVerifier<const C: usize, NTT, P, T> {
    _r: PhantomData<NTT>,
    _p: PhantomData<P>,
    _t: PhantomData<T>,
}

impl<
        const C: usize,
        NTT: SuitableRing,
        P: DecompositionParams,
        T: TranscriptWithShortChallenges<NTT>,
    > NIFSVerifier<C, NTT, P, T>
{
    pub fn verify(
        r_l: Vec<NTT>,
        acc: &LCCCS<C, NTT>,
        cm_i: &CCCS<C, NTT>,
        proof: &LFProof<C, NTT>,
        transcript: &mut impl TranscriptWithShortChallenges<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(Vec<NTT>, LCCCS<C, NTT>), LatticefoldError<NTT>> {
        let (r_r, linearized_cm_i) = LFLinearizationVerifier::<_, T>::verify(
            cm_i,
            &proof.linearization_proof,
            transcript,
            ccs,
        )?;
        let (r_l, decomposed_acc) = LFDecompositionVerifier::<_, T>::verify::<C, P>(
            r_l,
            acc,
            &proof.decomposition_proof_l,
            transcript,
            ccs,
        )?;
        let (r_r, decomposed_cm_i) = LFDecompositionVerifier::<_, T>::verify::<C, P>(
            r_r,
            &linearized_cm_i,
            &proof.decomposition_proof_r,
            transcript,
            ccs,
        )?;

        let lcccs_s = {
            let mut decomposed_acc = decomposed_acc;
            let mut decomposed_cm_i = decomposed_cm_i;

            decomposed_acc.append(&mut decomposed_cm_i);

            decomposed_acc
        };

        Ok(LFFoldingVerifier::<NTT, T>::verify::<C, P>(
            (r_l, r_r),
            &lcccs_s,
            &proof.folding_proof,
            transcript,
            ccs,
        )?)
    }
}
