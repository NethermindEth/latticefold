use lattirust_arithmetic::{challenge_set::latticefold_challenge_set::OverField, ring::PolyRing};
use std::marker::PhantomData;

use crate::{
    arith::{Witness, CCCS, CCS, LCCCS},
    parameters::DecompositionParams,
    transcript::Transcript,
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
/// `CR` is the type parameter for the coefficient representation of the ring.
/// `NTT` is the NTT representation of the same ring.
/// `P` is the decomposition parameters.
/// `T` is the FS-transform transcript.
pub struct NIFSProver<const C: usize, const W: usize, CR, NTT, P, T> {
    _cr: PhantomData<CR>,
    _r: PhantomData<NTT>,
    _p: PhantomData<P>,
    _t: PhantomData<T>,
}

impl<
        const C: usize,
        const W: usize,
        CR: PolyRing,
        NTT: OverField,
        P: DecompositionParams,
        T: Transcript<NTT>,
    > NIFSProver<C, W, CR, NTT, P, T>
{
    pub fn prove(
        acc: &LCCCS<C, NTT>,
        w_acc: &Witness<NTT>,
        cm_i: &CCCS<C, NTT>,
        w_i: &Witness<NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(LCCCS<C, NTT>, Witness<NTT>, LFProof<C, NTT>), LatticefoldError<NTT>> {
        let (linearized_cm_i, linearization_proof) =
            LFLinearizationProver::<_, T>::prove(cm_i, w_i, transcript, ccs)?;
        let (decomposed_lcccs_l, decomposed_wit_l, decomposition_proof_l) =
            LFDecompositionProver::<_, T>::prove::<C, CR, P>(acc, w_acc, transcript, ccs)?;
        let (decomposed_lcccs_r, decomposed_wit_r, decomposition_proof_r) =
            LFDecompositionProver::<_, T>::prove::<C, CR, P>(
                &linearized_cm_i,
                w_i,
                transcript,
                ccs,
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

        let (folded_lcccs, wit, folding_proof) =
            LFFoldingProver::<_, T>::prove::<C, CR, P>(&lcccs, &wit_s, transcript, ccs)?;

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
/// `CR` is the type parameter for the coefficient representation of the ring.
/// `NTT` is the NTT representation of the same ring.
/// `P` is the decomposition parameters.
/// `T` is the FS-transform transcript.
pub struct NIFSVerifier<const C: usize, CR, NTT, P, T> {
    _cr: PhantomData<CR>,
    _r: PhantomData<NTT>,
    _p: PhantomData<P>,
    _t: PhantomData<T>,
}

impl<const C: usize, CR: PolyRing, NTT: OverField, P: DecompositionParams, T: Transcript<NTT>>
    NIFSVerifier<C, CR, NTT, P, T>
{
    pub fn verify(
        acc: &LCCCS<C, NTT>,
        cm_i: &CCCS<C, NTT>,
        proof: &LFProof<C, NTT>,
        transcript: &mut impl Transcript<NTT>,
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
