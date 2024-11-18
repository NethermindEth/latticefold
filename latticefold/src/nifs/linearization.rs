use ark_std::cfg_iter;
use cyclotomic_rings::rings::SuitableRing;
use lattirust_poly::{
    mle::DenseMultilinearExtension,
    polynomials::{eq_eval, VPAuxInfo},
};

use utils::{compute_u, prepare_lin_sumcheck_polynomial};

use super::error::LinearizationError;
use crate::{
    arith::{utils::mat_vec_mul, Instance, Witness, CCCS, CCS, LCCCS},
    transcript::Transcript,
    utils::sumcheck::{MLSumcheck, SumCheckError::SumCheckFailed},
};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::utils::sumcheck::Proof;
pub use structs::*;

pub(crate) mod challenge_generator;
mod structs;

#[cfg(test)]
mod tests;
mod utils;

impl<NTT: SuitableRing, T: Transcript<NTT>> LFLinearizationProver<NTT, T> {
    fn prepare_prover_state<const C: usize>(
        wit: &Witness<NTT>,
        cm_i: &CCCS<C, NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<ProverState<NTT>, LinearizationError<NTT>> {
        let beta_s = BetaChallengeGenerator::generate_challenges(transcript, ccs.s);
        let z_ccs = cm_i.get_z_vector(&wit.w_ccs);

        let Mz_mles = ccs
            .M
            .iter()
            .map(|M| {
                Ok(DenseMultilinearExtension::from_slice(
                    ccs.s,
                    &mat_vec_mul(M, &z_ccs)?,
                ))
            })
            .collect::<Result<_, LinearizationError<_>>>()?;

        Ok(ProverState {
            beta_s,
            #[cfg(test)]
            z_ccs,
            Mz_mles,
        })
    }

    fn generate_sumcheck_proof(
        state: &ProverState<NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(Proof<NTT>, Vec<NTT>), LinearizationError<NTT>> {
        let g =
            prepare_lin_sumcheck_polynomial(ccs.s, &ccs.c, &state.Mz_mles, &ccs.S, &state.beta_s)?;

        let (sum_check_proof, prover_state) = MLSumcheck::prove_as_subprotocol(transcript, &g);
        let point_r = prover_state
            .randomness
            .into_iter()
            .map(|x| x.into())
            .collect::<Vec<NTT>>();

        Ok((sum_check_proof, point_r))
    }
}

impl<NTT: SuitableRing, T: Transcript<NTT>> LinearizationProver<NTT, T>
    for LFLinearizationProver<NTT, T>
{
    fn prove<const C: usize>(
        cm_i: &CCCS<C, NTT>,
        wit: &Witness<NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(LCCCS<C, NTT>, LinearizationProof<NTT>), LinearizationError<NTT>> {
        let state = Self::prepare_prover_state(wit, cm_i, transcript, ccs)?;
        let (sum_check_proof, point_r) = Self::generate_sumcheck_proof(&state, transcript, ccs)?;

        let v: Vec<NTT> = cfg_iter!(wit.f_hat)
            .map(|f_hat_row| {
                DenseMultilinearExtension::from_slice(ccs.s, f_hat_row)
                    .evaluate(&point_r)
                    .expect("cannot end up here, because the sumcheck subroutine must yield a point of the length log m")
            })
            .collect();

        let u = compute_u(&state.Mz_mles, &point_r)?;

        transcript.absorb_slice(&v);
        transcript.absorb_slice(&u);

        let linearization_proof = LinearizationProof {
            linearization_sumcheck: sum_check_proof,
            v: v.clone(),
            u: u.clone(),
        };

        let lcccs = LCCCS {
            r: point_r,
            v,
            cm: cm_i.cm.clone(),
            u,
            x_w: cm_i.x_ccs.clone(),
            h: NTT::one(),
        };

        Ok((lcccs, linearization_proof))
    }
}

impl<NTT: SuitableRing, T: Transcript<NTT>> LFLinearizationVerifier<NTT, T> {
    fn prepare_verifier_state(
        transcript: &mut impl Transcript<NTT>,
        proof: &LinearizationProof<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<VerifierState<NTT>, LinearizationError<NTT>> {
        let beta_s = BetaChallengeGenerator::generate_challenges(transcript, ccs.s);
        let poly_info = VPAuxInfo::new(ccs.s, ccs.d + 1);

        let subclaim = MLSumcheck::verify_as_subprotocol(
            transcript,
            &poly_info,
            NTT::zero(),
            &proof.linearization_sumcheck,
        )?;

        Ok(VerifierState {
            beta_s,
            point_r: subclaim.point.into_iter().map(|x| x.into()).collect(),
            s: subclaim.expected_evaluation,
        })
    }

    fn verify_evaluation_claim(
        state: &VerifierState<NTT>,
        proof: &LinearizationProof<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(), LinearizationError<NTT>> {
        let e = eq_eval(&state.point_r, &state.beta_s)?;
        let should_equal_s = e * ccs
            .c
            .iter()
            .enumerate()
            .map(|(i, &c)| c * ccs.S[i].iter().map(|&j| proof.u[j]).product::<NTT>())
            .sum::<NTT>();

        if should_equal_s != state.s {
            return Err(LinearizationError::SumCheckError(SumCheckFailed(
                should_equal_s,
                state.s,
            )));
        }

        Ok(())
    }
}

impl<NTT: SuitableRing, T: Transcript<NTT>> LinearizationVerifier<NTT, T>
    for LFLinearizationVerifier<NTT, T>
{
    fn verify<const C: usize>(
        cm_i: &CCCS<C, NTT>,
        proof: &LinearizationProof<NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<LCCCS<C, NTT>, LinearizationError<NTT>> {
        let state = Self::prepare_verifier_state(transcript, proof, ccs)?;
        Self::verify_evaluation_claim(&state, proof, ccs)?;

        transcript.absorb_slice(&proof.v);
        transcript.absorb_slice(&proof.u);

        Ok(LCCCS::<C, NTT> {
            r: state.point_r,
            v: proof.v.clone(),
            cm: cm_i.cm.clone(),
            u: proof.u.clone(),
            x_w: cm_i.x_ccs.clone(),
            h: NTT::one(),
        })
    }
}
