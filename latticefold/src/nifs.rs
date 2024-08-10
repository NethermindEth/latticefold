use lattirust_arithmetic::{challenge_set::latticefold_challenge_set::OverField, ring::PolyRing};
use std::marker::PhantomData;

use crate::{
    arith::{Witness, CCCS, CCS, LCCCS},
    commitment::AjtaiParams,
    transcript::Transcript,
};
use decomposition::{DecompositionProver, DecompositionVerifier};
use error::LatticefoldError;
use folding::{FoldingProver, FoldingVerifier};
use linearization::{LinearizationProver, LinearizationVerifier};

pub mod decomposition;
pub mod error;
pub mod folding;
pub mod linearization;

/// `CR` is the type parameter for the coefficient representation of the ring
/// `NTT` is the NTT representation of the same ring.
/// `P` is the Ajtai commitment parameters.
/// `T` is the FS-transform transcript.
#[derive(Debug, Clone)]
pub struct ComposedProof<
    CR: PolyRing,
    NTT: OverField,
    P: AjtaiParams,
    T: Transcript<NTT>,
    L: LinearizationProver<NTT, P, T>,
    D: DecompositionProver<CR, NTT, P, T>,
    FD: FoldingProver<CR, NTT, P, T>,
> {
    pub linearization_proof: L::Proof,
    pub decomposition_proof_l: D::Proof,
    pub decomposition_proof_r: D::Proof,
    pub folding_proof: FD::Proof,
}

type LatticefoldProof<CR, NTT, P, T> = ComposedProof<
    CR,
    NTT,
    P,
    T,
    NIFSProver<CR, NTT, P, T>,
    NIFSProver<CR, NTT, P, T>,
    NIFSProver<CR, NTT, P, T>,
>;

/// `CR` is the type parameter for the coefficient representation of the ring
/// `NTT` is the NTT representation of the same ring.
/// `P` is the Ajtai commitment parameters.
/// `T` is the FS-transform transcript.
pub struct NIFSProver<CR: PolyRing, NTT: OverField, P: AjtaiParams, T: Transcript<NTT>> {
    _cr: PhantomData<CR>,
    _r: PhantomData<NTT>,
    _p: PhantomData<P>,
    _t: PhantomData<T>,
}

impl<CR: PolyRing, NTT: OverField, P: AjtaiParams, T: Transcript<NTT>> NIFSProver<CR, NTT, P, T> {
    pub fn prove(
        acc: &LCCCS<NTT, P>,
        w_acc: &Witness<NTT>,
        cm_i: &CCCS<NTT, P>,
        w_i: &Witness<NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(LCCCS<NTT, P>, Witness<NTT>, LatticefoldProof<CR, NTT, P, T>), LatticefoldError<NTT>>
    {
        Self::prove_aux(acc, w_acc, cm_i, w_i, transcript, ccs)
    }

    fn prove_aux<
        L: LinearizationProver<NTT, P, T>,
        D: DecompositionProver<CR, NTT, P, T>,
        FP: FoldingProver<CR, NTT, P, T>,
        E: From<L::Error> + From<D::Error> + From<FP::Error>,
    >(
        acc: &LCCCS<NTT, P>,
        w_acc: &Witness<NTT>,
        cm_i: &CCCS<NTT, P>,
        w_i: &Witness<NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<
        (
            LCCCS<NTT, P>,
            Witness<NTT>,
            ComposedProof<CR, NTT, P, T, L, D, FP>,
        ),
        E,
    > {
        let (linearized_cm_i, linearization_proof) = L::prove(cm_i, w_i, transcript, ccs)?;
        let (decomposed_lcccs_l, decomposed_wit_l, decomposition_proof_l) =
            D::prove(acc, w_acc, transcript, ccs)?;
        let (decomposed_lcccs_r, decomposed_wit_r, decomposition_proof_r) =
            D::prove(&linearized_cm_i, w_i, transcript, ccs)?;

        let (lcccs, wit_s) = {
            let mut lcccs = decomposed_lcccs_l;
            let mut lcccs_r = decomposed_lcccs_r;
            lcccs.append(&mut lcccs_r);

            let mut wit_s = decomposed_wit_l;
            let mut wit_s_r = decomposed_wit_r;
            wit_s.append(&mut wit_s_r);

            (lcccs, wit_s)
        };

        let (folded_lcccs, wit, folding_proof) = FP::prove(&lcccs, &wit_s, transcript, ccs)?;

        Ok((
            folded_lcccs,
            wit,
            ComposedProof {
                linearization_proof,
                decomposition_proof_l,
                decomposition_proof_r,
                folding_proof,
            },
        ))
    }
}

/// `CR` is the type parameter for the coefficient representation of the ring
/// `NTT` is the NTT representation of the same ring.
/// `P` is the Ajtai commitment parameters.
/// `T` is the FS-transform transcript.
pub struct NIFSVerifier<CR: PolyRing, NTT: OverField, P: AjtaiParams, T: Transcript<NTT>> {
    _cr: PhantomData<CR>,
    _r: PhantomData<NTT>,
    _p: PhantomData<P>,
    _t: PhantomData<T>,
}

impl<CR: PolyRing, NTT: OverField, P: AjtaiParams, T: Transcript<NTT>> NIFSVerifier<CR, NTT, P, T> {
    pub fn verify(
        acc: &LCCCS<NTT, P>,
        cm_i: &CCCS<NTT, P>,
        proof: &LatticefoldProof<CR, NTT, P, T>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<LCCCS<NTT, P>, LatticefoldError<NTT>> {
        Self::verify_aux::<
            NIFSVerifier<CR, NTT, P, T>,
            NIFSVerifier<CR, NTT, P, T>,
            NIFSVerifier<CR, NTT, P, T>,
            LatticefoldError<NTT>,
        >(acc, cm_i, proof, transcript, ccs)
    }

    fn verify_aux<
        L: LinearizationVerifier<NTT, P, T>,
        D: DecompositionVerifier<CR, NTT, P, T>,
        FV: FoldingVerifier<CR, NTT, P, T>,
        E: From<L::Error> + From<D::Error> + From<FV::Error>,
    >(
        acc: &LCCCS<NTT, P>,
        cm_i: &CCCS<NTT, P>,
        proof: &ComposedProof<CR, NTT, P, T, L::Prover, D::Prover, FV::Prover>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<LCCCS<NTT, P>, E> {
        let linearized_cm_i = L::verify(cm_i, &proof.linearization_proof, transcript, ccs)?;
        let decomposed_acc = D::verify(acc, &proof.decomposition_proof_l, transcript, ccs)?;
        let decomposed_cm_i = D::verify(
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

        Ok(FV::verify(&lcccs_s, &proof.folding_proof, transcript, ccs)?)
    }
}
