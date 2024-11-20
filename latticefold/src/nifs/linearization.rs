use cyclotomic_rings::rings::SuitableRing;
use lattirust_poly::polynomials::VirtualPolynomial;
use lattirust_poly::{
    mle::DenseMultilinearExtension,
    polynomials::{eq_eval, VPAuxInfo},
};
use utils::{compute_u, prepare_lin_sumcheck_polynomial};

use super::error::LinearizationError;
use super::mle_helpers::{calculate_Mz_mles, calculate_f_hat_mle_evaluations, evaluate_mles};
use crate::ark_base::*;
use crate::{
    arith::{Witness, CCCS, CCS, LCCCS},
    transcript::Transcript,
    utils::sumcheck::{MLSumcheck, SumCheckError::SumCheckFailed},
};

use crate::arith::Instance;
use crate::nifs::linearization::utils::SqueezeBeta;
use crate::utils::sumcheck::Proof;
pub use structs::*;

mod structs;

#[cfg(test)]
mod tests;
pub mod utils;

impl<NTT: SuitableRing, T: Transcript<NTT>> LFLinearizationProver<NTT, T> {
    // Step 2 of Fig 5: Construct polynomial g and generate beta challenges
    fn construct_polynomial_g(
        z_ccs: &[NTT],
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<
        (VirtualPolynomial<NTT>, Vec<DenseMultilinearExtension<NTT>>),
        LinearizationError<NTT>,
    > {
        // Generate beta challenges from Step 1
        let beta_s =
            transcript.squeeze_beta_challenges(ccs.padded_length_of_decomposed_witness_log);

        // Prepare MLEs
        let Mz_mles = calculate_Mz_mles::<_, LinearizationError<_>>(ccs, z_ccs)?;

        // Construct the sumcheck polynomial g
        let g = prepare_lin_sumcheck_polynomial(
            ccs.padded_length_of_decomposed_witness_log,
            &ccs.c,
            &Mz_mles,
            &ccs.S,
            &beta_s,
        )?;

        Ok((g, Mz_mles))
    }

    // Step 2: Run sum-check protocol
    fn generate_sumcheck_proof(
        g: &VirtualPolynomial<NTT>,
        transcript: &mut impl Transcript<NTT>,
    ) -> Result<(Proof<NTT>, Vec<NTT>), LinearizationError<NTT>> {
        let (sum_check_proof, prover_state) = MLSumcheck::prove_as_subprotocol(transcript, g);
        let point_r = prover_state
            .randomness
            .into_iter()
            .map(|x| x.into())
            .collect::<Vec<NTT>>();

        Ok((sum_check_proof, point_r))
    }

    // Step 3: P sends V values
    fn compute_evaluation_vectors(
        ccs: &CCS<NTT>,
        wit: &Witness<NTT>,
        point_r: &[NTT],
        Mz_mles: &[DenseMultilinearExtension<NTT>],
    ) -> Result<(Vec<NTT>, Vec<NTT>, Vec<NTT>), LinearizationError<NTT>> {
        // Compute v
        let v: Vec<NTT> =
            calculate_f_hat_mle_evaluations::<_, LinearizationError<_>>(ccs, &wit.f_hat, point_r)?;

        // Compute u_j
        let u = compute_u(Mz_mles, point_r)?;

        Ok((point_r.to_vec(), v, u))
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
        // Step 1: Generate beta challenges (done in construct_polynomial_g because they are not needed
        // elsewhere.

        // Step 2: Sum check protocol.
        // z_ccs vector, i.e. concatenation x || 1 || w.
        let z_ccs = cm_i.get_z_vector(&wit.w_ccs);
        let (g, Mz_mles) = Self::construct_polynomial_g(&z_ccs, transcript, ccs)?;

        // Run sumcheck protocol.
        let (sumcheck_proof, point_r) = Self::generate_sumcheck_proof(&g, transcript)?;

        // Step 3: Compute v, u_vector.
        let (point_r, v, u) = Self::compute_evaluation_vectors(ccs, wit, &point_r, &Mz_mles)?;

        // Absorbing the prover's messages to the verifier.
        transcript.absorb_slice(&v);
        transcript.absorb_slice(&u);

        // Step 5: Output linearization_proof and lcccs
        let linearization_proof = LinearizationProof {
            linearization_sumcheck: sumcheck_proof,
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
    fn verify_sumcheck_proof(
        proof: &LinearizationProof<NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(Vec<NTT>, NTT), LinearizationError<NTT>> {
        // The polynomial has degree <= ccs.d + 1 and log_m (ccs.s) vars.
        let poly_info = VPAuxInfo::new(ccs.padded_length_of_decomposed_witness_log, ccs.d + 1);

        let subclaim = MLSumcheck::verify_as_subprotocol(
            transcript,
            &poly_info,
            NTT::zero(),
            &proof.linearization_sumcheck,
        )?;

        Ok((
            subclaim.point.into_iter().map(|x| x.into()).collect(),
            subclaim.expected_evaluation,
        ))
    }

    // Step 4: Verify that e·(Σ c_i·Π_{j∈S_i} u_j) = s
    fn verify_evaluation_claim(
        beta_s: &[NTT],
        point_r: &[NTT],
        s: NTT,
        proof: &LinearizationProof<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(), LinearizationError<NTT>> {
        let e = eq_eval(point_r, beta_s)?;
        let should_equal_s = e * ccs // e * (\sum c_i * \Pi_{j \in S_i} u_j)
            .c
            .iter()
            .enumerate()
            .map(|(i, &c)| c * ccs.S[i].iter().map(|&j| proof.u[j]).product::<NTT>()) // c_i * \Pi_{j \in S_i} u_j
            .sum::<NTT>(); // \sum c_i * \Pi_{j \in S_i} u_j

        if should_equal_s != s {
            return Err(LinearizationError::SumCheckError(SumCheckFailed(
                should_equal_s,
                s,
            )));
        }

        Ok(())
    }

    fn prepare_verifier_output<const C: usize>(
        cm_i: &CCCS<C, NTT>,
        point_r: Vec<NTT>,
        proof: &LinearizationProof<NTT>,
    ) -> LCCCS<C, NTT> {
        LCCCS {
            r: point_r,
            v: proof.v.clone(),
            cm: cm_i.cm.clone(),
            u: proof.u.clone(),
            x_w: cm_i.x_ccs.clone(),
            h: NTT::one(),
        }
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
        // Step 1: Generate the beta challenges.
        let beta_s =
            transcript.squeeze_beta_challenges(ccs.padded_length_of_decomposed_witness_log);

        //Step 2: The sumcheck.
        let (point_r, s) = Self::verify_sumcheck_proof(proof, transcript, ccs)?;

        Self::verify_evaluation_claim(&beta_s, &point_r, s, proof, ccs)?;

        // Absorbing the prover's messages to the verifier.
        transcript.absorb_slice(&proof.v);
        transcript.absorb_slice(&proof.u);

        // Step 5: Output z_o
        Ok(Self::prepare_verifier_output(cm_i, point_r, proof))
    }
}
