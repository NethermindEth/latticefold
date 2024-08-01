pub mod decomposition;
pub mod error;
pub mod folding;
#[allow(non_snake_case)]
pub mod linearization;

use std::marker::PhantomData;

use decomposition::{DecompositionProver, DecompositionVerifier};
use error::LatticefoldError;
use folding::{FoldingProver, FoldingVerifier};
use lattirust_arithmetic::challenge_set::latticefold_challenge_set::OverField;
use linearization::{LinearizationProver, LinearizationVerifier};

use crate::arith::{Witness, CCS, LCCCS};
use crate::{arith::CCCS, transcript::Transcript};

#[derive(Debug, Clone)]
pub struct ComposedProof<
    R: OverField,
    T: Transcript<R>,
    L: LinearizationProver<R, T>,
    D: DecompositionProver<R, T>,
    FD: FoldingProver<R, T>,
> {
    pub linearization_proof: L::Proof,
    pub decomposition_proof_l: D::Proof,
    pub decomposition_proof_r: D::Proof,
    pub folding_proof: FD::Proof,
}

type LatticefoldProof<R, T> =
    ComposedProof<R, T, NIFSProver<R, T>, NIFSProver<R, T>, NIFSProver<R, T>>;

pub struct NIFSProver<R: OverField, T: Transcript<R>> {
    _r: PhantomData<R>,
    _t: PhantomData<T>,
}

impl<R: OverField, T: Transcript<R>> NIFSProver<R, T> {
    pub fn prove(
        acc: &LCCCS<R>,
        w_acc: &Witness<R>,
        cm_i: &CCCS<R>,
        w_i: &Witness<R>,
        transcript: &mut impl Transcript<R>,
        ccs: &CCS<R>,
    ) -> Result<(LCCCS<R>, Witness<R>, LatticefoldProof<R, T>), LatticefoldError<R>> {
        Self::prove_aux(acc, w_acc, cm_i, w_i, transcript, ccs)
    }

    fn prove_aux<
        L: LinearizationProver<R, T>,
        D: DecompositionProver<R, T>,
        FP: FoldingProver<R, T>,
        E: From<L::Error> + From<D::Error> + From<FP::Error>,
    >(
        acc: &LCCCS<R>,
        w_acc: &Witness<R>,
        cm_i: &CCCS<R>,
        w_i: &Witness<R>,
        transcript: &mut impl Transcript<R>,
        ccs: &CCS<R>,
    ) -> Result<(LCCCS<R>, Witness<R>, ComposedProof<R, T, L, D, FP>), E> {
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

pub struct NIFSVerifier<R: OverField, T: Transcript<R>> {
    _r: PhantomData<R>,
    _t: PhantomData<T>,
}

impl<R: OverField, T: Transcript<R>> NIFSVerifier<R, T> {
    pub fn verify(
        acc: &LCCCS<R>,
        cm_i: &CCCS<R>,
        proof: &LatticefoldProof<R, T>,
        transcript: &mut impl Transcript<R>,
        ccs: &CCS<R>,
    ) -> Result<LCCCS<R>, LatticefoldError<R>> {
        Self::verify_aux::<
            NIFSVerifier<R, T>,
            NIFSVerifier<R, T>,
            NIFSVerifier<R, T>,
            LatticefoldError<R>,
        >(acc, cm_i, proof, transcript, ccs)
    }

    fn verify_aux<
        L: LinearizationVerifier<R, T>,
        D: DecompositionVerifier<R, T>,
        FV: FoldingVerifier<R, T>,
        E: From<L::Error> + From<D::Error> + From<FV::Error>,
    >(
        acc: &LCCCS<R>,
        cm_i: &CCCS<R>,
        proof: &ComposedProof<R, T, L::Prover, D::Prover, FV::Prover>,
        transcript: &mut impl Transcript<R>,
        ccs: &CCS<R>,
    ) -> Result<LCCCS<R>, E> {
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
