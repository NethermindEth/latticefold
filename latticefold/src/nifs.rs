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
use lattirust_arithmetic::ring::ConvertibleRing;
use linearization::{LinearizationProver, LinearizationVerifier};

use crate::arith::{Witness, CCS, LCCCS};
use crate::commitment::AjtaiParams;
use crate::{arith::CCCS, transcript::Transcript};

#[derive(Debug, Clone)]
pub struct ComposedProof<
    CR: ConvertibleRing,
    R: OverField,
    P: AjtaiParams<CR>,
    T: Transcript<R>,
    L: LinearizationProver<CR, R, P, T>,
    D: DecompositionProver<CR, R, P, T>,
    FD: FoldingProver<CR, R, P, T>,
> {
    pub linearization_proof: L::Proof,
    pub decomposition_proof_l: D::Proof,
    pub decomposition_proof_r: D::Proof,
    pub folding_proof: FD::Proof,
}

type LatticefoldProof<CR, R, P, T> = ComposedProof<
    CR,
    R,
    P,
    T,
    NIFSProver<CR, R, P, T>,
    NIFSProver<CR, R, P, T>,
    NIFSProver<CR, R, P, T>,
>;

pub struct NIFSProver<CR: ConvertibleRing, R: OverField, P: AjtaiParams<CR>, T: Transcript<R>> {
    _cr: PhantomData<CR>,
    _r: PhantomData<R>,
    _p: PhantomData<P>,
    _t: PhantomData<T>,
}

impl<CR: ConvertibleRing, R: OverField, P: AjtaiParams<CR>, T: Transcript<R>>
    NIFSProver<CR, R, P, T>
{
    pub fn prove(
        acc: &LCCCS<CR, R, P>,
        w_acc: &Witness<CR, R>,
        cm_i: &CCCS<CR, R, P>,
        w_i: &Witness<CR, R>,
        transcript: &mut impl Transcript<R>,
        ccs: &CCS<R>,
    ) -> Result<
        (
            LCCCS<CR, R, P>,
            Witness<CR, R>,
            LatticefoldProof<CR, R, P, T>,
        ),
        LatticefoldError<R>,
    > {
        Self::prove_aux(acc, w_acc, cm_i, w_i, transcript, ccs)
    }

    fn prove_aux<
        L: LinearizationProver<CR, R, P, T>,
        D: DecompositionProver<CR, R, P, T>,
        FP: FoldingProver<CR, R, P, T>,
        E: From<L::Error> + From<D::Error> + From<FP::Error>,
    >(
        acc: &LCCCS<CR, R, P>,
        w_acc: &Witness<CR, R>,
        cm_i: &CCCS<CR, R, P>,
        w_i: &Witness<CR, R>,
        transcript: &mut impl Transcript<R>,
        ccs: &CCS<R>,
    ) -> Result<
        (
            LCCCS<CR, R, P>,
            Witness<CR, R>,
            ComposedProof<CR, R, P, T, L, D, FP>,
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

pub struct NIFSVerifier<CR: ConvertibleRing, R: OverField, P: AjtaiParams<CR>, T: Transcript<R>> {
    _cr: PhantomData<CR>,
    _r: PhantomData<R>,
    _p: PhantomData<P>,
    _t: PhantomData<T>,
}

impl<CR: ConvertibleRing, R: OverField, P: AjtaiParams<CR>, T: Transcript<R>>
    NIFSVerifier<CR, R, P, T>
{
    pub fn verify(
        acc: &LCCCS<CR, R, P>,
        cm_i: &CCCS<CR, R, P>,
        proof: &LatticefoldProof<CR, R, P, T>,
        transcript: &mut impl Transcript<R>,
        ccs: &CCS<R>,
    ) -> Result<LCCCS<CR, R, P>, LatticefoldError<R>> {
        Self::verify_aux::<
            NIFSVerifier<CR, R, P, T>,
            NIFSVerifier<CR, R, P, T>,
            NIFSVerifier<CR, R, P, T>,
            LatticefoldError<R>,
        >(acc, cm_i, proof, transcript, ccs)
    }

    fn verify_aux<
        L: LinearizationVerifier<CR, R, P, T>,
        D: DecompositionVerifier<CR, R, P, T>,
        FV: FoldingVerifier<CR, R, P, T>,
        E: From<L::Error> + From<D::Error> + From<FV::Error>,
    >(
        acc: &LCCCS<CR, R, P>,
        cm_i: &CCCS<CR, R, P>,
        proof: &ComposedProof<CR, R, P, T, L::Prover, D::Prover, FV::Prover>,
        transcript: &mut impl Transcript<R>,
        ccs: &CCS<R>,
    ) -> Result<LCCCS<CR, R, P>, E> {
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
