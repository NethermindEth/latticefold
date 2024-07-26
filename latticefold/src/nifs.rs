pub mod decomposition;
pub mod folding;
pub mod linearization;

use std::marker::PhantomData;

use ark_crypto_primitives::sponge::Absorb;
use ark_ff::Field;
use decomposition::{DecompositionError, DecompositionProver};
use folding::{FoldingError, FoldingProver};
use lattirust_arithmetic::challenge_set::latticefold_challenge_set::OverField;
use lattirust_arithmetic::ring::Ring;
use linearization::{LinearizationError, LinearizationProver};
use thiserror::Error;

use crate::arith::{Witness, LCCCS};
use crate::{arith::CCCS, transcript::Transcript};

pub struct NIFSProver<F: Field, R: OverField<F>, T: Transcript<F, R>> {
    _f: PhantomData<F>,
    _r: PhantomData<R>,
    _t: PhantomData<T>,
}

pub struct NIFSVerifier<F: Field, R: OverField<F>, T: Transcript<F, R>> {
    _f: PhantomData<F>,
    _r: PhantomData<R>,
    _t: PhantomData<T>,
}

#[derive(Debug, Clone)]
pub struct ComposedProof<
    F: Field,
    R: OverField<F>,
    T: Transcript<F, R>,
    L: LinearizationProver<F, R, T>,
    D: DecompositionProver<F, R, T>,
    FD: FoldingProver<F, R, T>,
> {
    pub linearization_proof: L::Proof,
    pub decomposition_proof_l: D::Proof,
    pub decomposition_proof_r: D::Proof,
    pub folding_proof: FD::Proof,
}

type LatticefoldProof<F, R, T> =
    ComposedProof<F, R, T, NIFSProver<F, R, T>, NIFSProver<F, R, T>, NIFSProver<F, R, T>>;

#[derive(Debug, Error)]
pub enum LatticefoldError<R: Ring> {
    #[error("linearization failed: {0}")]
    LinearizationError(#[from] LinearizationError<R>),
    #[error("decomposition failed: {0}")]
    DecompositionError(#[from] DecompositionError<R>),
    #[error("folding failed: {0}")]
    FoldingError(#[from] FoldingError<R>),
}

impl<F: Field + Absorb, R: OverField<F>, T: Transcript<F, R>> NIFSProver<F, R, T> {
    fn prove_aux<
        L: LinearizationProver<F, R, T>,
        D: DecompositionProver<F, R, T>,
        FP: FoldingProver<F, R, T>,
        E: From<L::Error> + From<D::Error> + From<FP::Error>,
    >(
        acc: &LCCCS<R>,
        w_acc: &Witness<R>,
        cm_i: &CCCS<R>,
        w_i: &Witness<R>,
        transcript: &mut impl Transcript<F, R>,
    ) -> Result<(LCCCS<R>, Witness<R>, ComposedProof<F, R, T, L, D, FP>), E> {
        let (linearized_cm_i, linearization_proof) = L::prove(cm_i, w_i, transcript)?;
        let (decomposed_lcccs_l, decomposed_wit_l, decomposition_proof_l) =
            D::prove(acc, w_acc, transcript)?;
        let (decomposed_lcccs_r, decomposed_wit_r, decomposition_proof_r) =
            D::prove(&linearized_cm_i, w_i, transcript)?;

        let (lcccs, wit_s) = {
            let mut lcccs = decomposed_lcccs_l;
            let mut lcccs_r = decomposed_lcccs_r;
            lcccs.append(&mut lcccs_r);

            let mut wit_s = decomposed_wit_l;
            let mut wit_s_r = decomposed_wit_r;
            wit_s.append(&mut wit_s_r);

            (lcccs, wit_s)
        };

        let (folded_lcccs, wit, folding_proof) = FP::prove(&lcccs, &wit_s, transcript)?;

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

    pub fn prove(
        acc: &LCCCS<R>,
        w_acc: &Witness<R>,
        cm_i: &CCCS<R>,
        w_i: &Witness<R>,
        transcript: &mut impl Transcript<F, R>,
    ) -> Result<(LCCCS<R>, Witness<R>, LatticefoldProof<F, R, T>), LatticefoldError<R>> {
        Self::prove_aux(acc, w_acc, cm_i, w_i, transcript)
    }
}
