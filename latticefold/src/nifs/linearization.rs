use ark_ff::{Field, PrimeField};
use ark_std::cfg_iter;
use cyclotomic_rings::rings::SuitableRing;
use lattirust_poly::{
    mle::DenseMultilinearExtension,
    polynomials::{eq_eval, VPAuxInfo},
};

use super::error::LinearizationError;
use crate::{
    arith::{utils::mat_vec_mul, Instance, Witness, CCCS, CCS, LCCCS},
    nifs::common_helpers::to_mles_err,
    transcript::Transcript,
    utils::sumcheck::{MLSumcheck, SumCheckError::SumCheckFailed},
};
use utils::{compute_u, prepare_lin_sumcheck_polynomial};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

pub use structs::*;
mod structs;

#[cfg(test)]
mod tests;
mod utils;

impl<NTT: SuitableRing, T: Transcript<NTT>> LinearizationProver<NTT, T>
    for LFLinearizationProver<NTT, T>
{
    fn prove<const C: usize>(
        cm_i: &CCCS<C, NTT>,
        wit: &Witness<NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(LCCCS<C, NTT>, LinearizationProof<NTT>), LinearizationError<NTT>> {
        let log_m = ccs.s;
        // Step 1: Generate the beta challenges.
        transcript.absorb_field_element(&<NTT::BaseRing as Field>::from_base_prime_field(
            <NTT::BaseRing as Field>::BasePrimeField::from_be_bytes_mod_order(b"beta_s"),
        ));
        let beta_s: Vec<NTT> = transcript
            .get_challenges(log_m)
            .into_iter()
            .map(|x| x.into())
            .collect();
        // Step 2: Sum check protocol

        // z_ccs vector, i.e. concatenation x || 1 || w.
        let z_ccs: Vec<NTT> = cm_i.get_z_vector(&wit.w_ccs);

        let Mz_mles = to_mles_err::<_, _, LinearizationError<NTT>, _>(
            log_m,
            cfg_iter!(ccs.M).map(|M| mat_vec_mul(M, &z_ccs)),
        )?;
        // The sumcheck polynomial
        let g = prepare_lin_sumcheck_polynomial(log_m, &ccs.c, &Mz_mles, &ccs.S, &beta_s)?;

        // Run sum check prover
        let (sum_check_proof, prover_state) = MLSumcheck::prove_as_subprotocol(transcript, &g);

        // Extract the evaluation point
        let r = prover_state
            .randomness
            .into_iter()
            .map(|x| x.into())
            .collect::<Vec<NTT>>();

        // Step 3: Compute v, u_vector
        let v: Vec<NTT> = cfg_iter!(wit.f_hat).map(|f_hat_row| DenseMultilinearExtension::from_slice(log_m, f_hat_row).evaluate(&r).expect("cannot end up here, because the sumcheck subroutine must yield a point of the length log m")).collect();

        let u = compute_u(&Mz_mles, &r)?;

        // Absorbing the prover's messages to the verifier.
        transcript.absorb_slice(&v);
        transcript.absorb_slice(&u);

        // Step 5: Output linearization_proof and lcccs
        let linearization_proof = LinearizationProof {
            linearization_sumcheck: sum_check_proof,
            v: v.clone(),
            u: u.clone(),
        };
        let lcccs = LCCCS {
            r,
            v,
            cm: cm_i.cm.clone(),
            u,
            x_w: cm_i.x_ccs.clone(),
            h: NTT::one(),
        };
        Ok((lcccs, linearization_proof))
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
        let log_m = ccs.s;
        // Step 1: Generate the beta challenges.
        transcript.absorb_field_element(&<NTT::BaseRing as Field>::from_base_prime_field(
            <NTT::BaseRing as Field>::BasePrimeField::from_be_bytes_mod_order(b"beta_s"),
        ));
        let beta_s: Vec<NTT> = transcript
            .get_challenges(log_m)
            .into_iter()
            .map(|x| x.into())
            .collect();

        //Step 2: The sumcheck.
        // The polynomial has degree <= ccs.d + 1 and log_m vars.
        let poly_info = VPAuxInfo::new(log_m, ccs.d + 1);

        // Verify the sumcheck proof.
        let subclaim = MLSumcheck::verify_as_subprotocol(
            transcript,
            &poly_info,
            NTT::zero(),
            &proof.linearization_sumcheck,
        )?;

        // Absorbing the prover's messages to the verifier.
        transcript.absorb_slice(&proof.v);
        transcript.absorb_slice(&proof.u);

        // The final evaluation claim from the sumcheck.
        let s = subclaim.expected_evaluation;

        let point_r = subclaim
            .point
            .into_iter()
            .map(|x| x.into())
            .collect::<Vec<NTT>>();

        // Step 4: reshaping the evaluation claim.
        // eq(beta, r)
        let e = eq_eval(&point_r, &beta_s)?;
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

        Ok(LCCCS::<C, NTT> {
            r: point_r,
            v: proof.v.clone(),
            cm: cm_i.cm.clone(),
            u: proof.u.clone(),
            x_w: cm_i.x_ccs.clone(),
            h: NTT::one(),
        })
    }
}
