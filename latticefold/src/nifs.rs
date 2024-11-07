use ark_std::marker::PhantomData;

use cyclotomic_rings::SuitableRing;
use lattirust_ring::OverField;

use crate::nifs::structs::LatticefoldState;
use crate::{
    arith::{Witness, CCCS, CCS, LCCCS},
    commitment::AjtaiCommitmentScheme,
    decomposition_parameters::DecompositionParams,
    transcript::TranscriptWithSmallChallenges,
};
use decomposition::{
    DecompositionProof, DecompositionProver, DecompositionVerifier, LFDecompositionProver,
    LFDecompositionVerifier,
};
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
pub mod structs;

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
        T: TranscriptWithSmallChallenges<NTT>,
    > NIFSProver<C, W, NTT, P, T>
{
    pub fn prove(
        acc: &LCCCS<C, NTT>,
        w_acc: &Witness<NTT>,
        cm_i: &CCCS<C, NTT>,
        w_i: &Witness<NTT>,
        transcript: &mut impl TranscriptWithSmallChallenges<NTT>,
        ccs: &CCS<NTT>,
        scheme: &AjtaiCommitmentScheme<C, W, NTT>,
    ) -> Result<(LCCCS<C, NTT>, Witness<NTT>, LFProof<C, NTT>), LatticefoldError<NTT>> {
        let mut latticefold_state = LatticefoldState::<C, NTT>::default();
        let linearization_proof = LFLinearizationProver::<_, T>::prove(
            cm_i,
            w_i,
            transcript,
            ccs,
            &mut latticefold_state,
        )?;
        let decomposition_proof_l = LFDecompositionProver::<_, T>::prove::<W, C, P>(
            acc,
            w_acc,
            transcript,
            ccs,
            scheme,
            &mut latticefold_state,
        )?;
        // Work around for the mutable and inmmutable borrow at the same time, needs a better solution
        let lcccs_for_decomposition = latticefold_state.lcccs.clone();
        let decomposition_proof_r = LFDecompositionProver::<_, T>::prove::<W, C, P>(
            &lcccs_for_decomposition,
            w_i,
            transcript,
            ccs,
            scheme,
            &mut latticefold_state,
        )?;

        let (folded_lcccs, wit, folding_proof) =
            LFFoldingProver::<_, T>::prove::<C, P>(&lcccs, &wit_s, transcript, ccs)?;

        Ok((
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
        T: TranscriptWithSmallChallenges<NTT>,
    > NIFSVerifier<C, NTT, P, T>
{
    pub fn verify(
        acc: &LCCCS<C, NTT>,
        cm_i: &CCCS<C, NTT>,
        proof: &LFProof<C, NTT>,
        transcript: &mut impl TranscriptWithSmallChallenges<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<LCCCS<C, NTT>, LatticefoldError<NTT>> {
        let linearized_cm_i = LFLinearizationVerifier::<_, T>::verify(
            cm_i,
            &proof.linearization_proof,
            transcript,
            ccs,
        )?;
        let decomposed_acc = LFDecompositionVerifier::<_, T>::verify::<C, P>(
            acc,
            &proof.decomposition_proof_l,
            transcript,
            ccs,
        )?;
        let decomposed_cm_i = LFDecompositionVerifier::<_, T>::verify::<C, P>(
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
            &lcccs_s,
            &proof.folding_proof,
            transcript,
            ccs,
        )?)
    }
}
