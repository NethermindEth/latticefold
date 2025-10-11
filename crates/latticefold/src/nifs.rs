//! The NIFS module defines the behaviour of the [LatticeFold](https://eprint.iacr.org/2024/257.pdf) protocol
//!
//! NIFS = Non Interactive Folding Scheme

use ark_ff::{Field, PrimeField};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{marker::PhantomData, vec::Vec};
use cyclotomic_rings::rings::SuitableRing;
use stark_rings::OverField;

use self::{decomposition::*, error::LatticefoldError, folding::*, linearization::*};
use crate::{
    arith::{error::CSError, Witness, CCCS, CCS, LCCCS},
    commitment::AjtaiCommitmentScheme,
    decomposition_parameters::DecompositionParams,
    transcript::{Transcript, TranscriptWithShortChallenges},
};

pub mod decomposition;
pub mod error;
pub mod folding;
pub mod linearization;

#[cfg(test)]
mod tests;

/// `NTT` is a cyclotomic ring in the NTT form.
#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct LFProof<NTT: OverField> {
    pub linearization_proof: LinearizationProof<NTT>,
    pub decomposition_proof_l: DecompositionProof<NTT>,
    pub decomposition_proof_r: DecompositionProof<NTT>,
    pub folding_proof: FoldingProof<NTT>,
}

/// `NTT` is a suitable cyclotomic ring.
/// `T` is the FS-transform transcript.
#[derive(Clone)]
pub struct NIFSProver<NTT, T> {
    _r: PhantomData<NTT>,
    dparams: DecompositionParams,
    transcript: T,
}

impl<NTT: SuitableRing, T: TranscriptWithShortChallenges<NTT> + Sync> NIFSProver<NTT, T> {
    pub fn new(dparams: DecompositionParams, transcript: T) -> Self {
        Self {
            _r: Default::default(),
            dparams,
            transcript,
        }
    }
}

impl<NTT: SuitableRing, T: TranscriptWithShortChallenges<NTT> + Sync> NIFSProver<NTT, T> {
    pub fn prove(
        &mut self,
        acc: &LCCCS<NTT>,
        w_acc: &Witness<NTT>,
        cm_i: &CCCS<NTT>,
        w_i: &Witness<NTT>,
        ccs: &CCS<NTT>,
        scheme: &AjtaiCommitmentScheme<NTT>,
    ) -> Result<(LCCCS<NTT>, Witness<NTT>, LFProof<NTT>), LatticefoldError<NTT>> {
        sanity_check::<NTT>(ccs, self.dparams.l)?;

        absorb_public_input::<NTT>(acc, cm_i, &mut self.transcript);

        let mut linearization_prover = LFLinearizationProver::new(&mut self.transcript);
        let (linearized_cm_i, linearization_proof) = linearization_prover.prove(cm_i, w_i, ccs)?;

        let mut decomposition_prover =
            LFDecompositionProver::new(self.dparams.clone(), &mut self.transcript);
        let (mz_mles_l, decomposed_lcccs_l, decomposed_wit_l, decomposition_proof_l) =
            decomposition_prover.prove(acc, w_acc, ccs, scheme)?;
        let (mz_mles_r, decomposed_lcccs_r, decomposed_wit_r, decomposition_proof_r) =
            decomposition_prover.prove(&linearized_cm_i, w_i, ccs, scheme)?;

        let (mz_mles, lcccs, wit_s) = {
            let mut lcccs = decomposed_lcccs_l;
            let mut lcccs_r = decomposed_lcccs_r;
            lcccs.append(&mut lcccs_r);

            let mut wit_s = decomposed_wit_l;
            let mut wit_s_r = decomposed_wit_r;
            wit_s.append(&mut wit_s_r);

            let mut mz_mles = mz_mles_l;
            let mut mz_mles_r = mz_mles_r;
            mz_mles.append(&mut mz_mles_r);
            (mz_mles, lcccs, wit_s)
        };

        let mut folding_prover = LFFoldingProver::new(self.dparams.clone(), &mut self.transcript);
        let (folded_lcccs, wit, folding_proof) =
            folding_prover.prove(&lcccs, wit_s, ccs, &mz_mles)?;

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

/// `NTT` is a suitable cyclotomic ring.
/// `T` is the FS-transform transcript.
#[derive(Clone)]
pub struct NIFSVerifier<NTT, T> {
    _r: PhantomData<NTT>,
    dparams: DecompositionParams,
    transcript: T,
}

impl<NTT: SuitableRing, T: TranscriptWithShortChallenges<NTT> + Sync> NIFSVerifier<NTT, T> {
    pub fn new(dparams: DecompositionParams, transcript: T) -> Self {
        Self {
            _r: Default::default(),
            dparams,
            transcript,
        }
    }
}

impl<NTT: SuitableRing, T: TranscriptWithShortChallenges<NTT> + Sync> NIFSVerifier<NTT, T> {
    pub fn verify(
        &mut self,
        acc: &LCCCS<NTT>,
        cm_i: &CCCS<NTT>,
        proof: &LFProof<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<LCCCS<NTT>, LatticefoldError<NTT>> {
        sanity_check::<NTT>(ccs, self.dparams.l)?;

        absorb_public_input::<NTT>(acc, cm_i, &mut self.transcript);

        let mut linearization_verifier = LFLinearizationVerifier::new(&mut self.transcript);
        let linearized_cm_i =
            linearization_verifier.verify(cm_i, &proof.linearization_proof, ccs)?;

        let mut decomposition_verifier =
            LFDecompositionVerifier::new(self.dparams.clone(), &mut self.transcript);
        let decomposed_acc =
            decomposition_verifier.verify(acc, &proof.decomposition_proof_l, ccs)?;
        let decomposed_cm_i =
            decomposition_verifier.verify(&linearized_cm_i, &proof.decomposition_proof_r, ccs)?;

        let lcccs_s = {
            let mut decomposed_acc = decomposed_acc;
            let mut decomposed_cm_i = decomposed_cm_i;

            decomposed_acc.append(&mut decomposed_cm_i);

            decomposed_acc
        };

        let mut folding_verifier =
            LFFoldingVerifier::new(self.dparams.clone(), &mut self.transcript);
        Ok(folding_verifier.verify(&lcccs_s, &proof.folding_proof, ccs)?)
    }
}

fn sanity_check<NTT: SuitableRing>(ccs: &CCS<NTT>, l: usize) -> Result<(), LatticefoldError<NTT>> {
    if ccs.m != usize::max((ccs.n - ccs.l - 1) * l, ccs.m).next_power_of_two() {
        return Err(CSError::InvalidSizeBounds(ccs.m, ccs.n, l).into());
    }

    Ok(())
}

fn absorb_public_input<NTT: SuitableRing>(
    acc: &LCCCS<NTT>,
    cm_i: &CCCS<NTT>,
    transcript: &mut impl Transcript<NTT>,
) {
    transcript.absorb_field_element(&<NTT::BaseRing as Field>::from_base_prime_field(
        <NTT::BaseRing as Field>::BasePrimeField::from_be_bytes_mod_order(b"acc"),
    ));

    transcript.absorb_slice(&acc.r);
    transcript.absorb_slice(&acc.v);
    transcript.absorb_slice(acc.cm.as_ref());
    transcript.absorb_slice(&acc.u);
    transcript.absorb_slice(&acc.x_w);
    transcript.absorb(&acc.h);

    transcript.absorb_field_element(&<NTT::BaseRing as Field>::from_base_prime_field(
        <NTT::BaseRing as Field>::BasePrimeField::from_be_bytes_mod_order(b"cm_i"),
    ));

    transcript.absorb_slice(cm_i.cm.as_ref());
    transcript.absorb_slice(&cm_i.x_ccs);
}
