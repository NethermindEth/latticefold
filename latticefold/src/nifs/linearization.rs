use cyclotomic_rings::rings::SuitableRing;
use lattirust_poly::mle::DenseMultilinearExtension;
use lattirust_ring::OverField;
use utils::{compute_u, prepare_lin_sumcheck_polynomial};

use super::error::LinearizationError;
use crate::ark_base::*;
use crate::utils::mle_helpers::{calculate_Mz_mles, evaluate_mles};
use crate::{
    arith::{Witness, CCCS, CCS, LCCCS},
    transcript::Transcript,
    utils::sumcheck::{
        virtual_polynomial::{eq_eval, VPAuxInfo, VirtualPolynomial},
        MLSumcheck,
        SumCheckError::SumCheckFailed,
    },
};

use crate::arith::Instance;
use crate::nifs::linearization::utils::SqueezeBeta;
use crate::utils::sumcheck::Proof;
pub use structs::*;

mod structs;

#[cfg(test)]
mod tests;
pub mod utils;

/// Prover for the Linearization subprotocol
pub trait LinearizationProver<NTT: SuitableRing, T: Transcript<NTT>> {
    /// Generates a proof for the linearization subprotocol
    ///
    /// # Arguments
    ///
    /// * `cm_i` - A reference to a committed CCS statement to be linearized, i.e. a CCCS<C, NTT>.
    /// * `wit` - A reference to a CCS witness for the statement cm_i.
    /// * `transcript` - A mutable reference to a sponge for generating NI challenges.
    /// * `ccs` - A reference to a Customizable Constraint System circuit representation.
    ///
    /// # Returns
    ///
    /// On success, returns a tuple `(LCCCS<C, NTT>, LinearizationProof<NTT>)` where:
    ///   * `LCCCS<C, NTT>` is a linearized version of the CCS witness commitment.
    ///   * `LinearizationProof<NTT>` is a proof that the linearization subprotocol was executed correctly.
    ///
    /// # Errors
    ///
    /// Returns an error if asked to evaluate MLEs with incorrect number of variables
    ///
    fn prove<const C: usize>(
        cm_i: &CCCS<C, NTT>,
        wit: &Witness<NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(LCCCS<C, NTT>, LinearizationProof<NTT>), LinearizationError<NTT>>;
}

/// Verifier for the Linearization subprotocol.
pub trait LinearizationVerifier<NTT: OverField, T: Transcript<NTT>> {
    /// Verifies a proof for the linearization subprotocol.
    ///
    /// # Arguments
    ///
    /// * `cm_i` - A reference to a `CCCS<C, NTT>`, which represents a CCS statement and a commitment to a witness.
    /// * `proof` - A reference to a `LinearizationProof<NTT>` containing the linearization proof.
    /// * `transcript` - A mutable reference to a sponge for generating NI challenges.
    /// * `ccs` - A reference to a Customizable Constraint System instance used in the protocol.
    ///
    /// # Returns
    ///
    /// * `Ok(LCCCS<C, NTT>)` - On success, returns a linearized version of the CCS witness commitment.
    /// * `Err(LinearizationError<NTT>)` - If verification fails, returns a `LinearizationError<NTT>`.
    ///
    fn verify<const C: usize>(
        cm_i: &CCCS<C, NTT>,
        proof: &LinearizationProof<NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<LCCCS<C, NTT>, LinearizationError<NTT>>;
}

impl<NTT: SuitableRing, T: Transcript<NTT>> LFLinearizationProver<NTT, T> {
    /// Step 2 of Fig 5: Construct polynomial $g$ and generate $\beta$ challenges.
    fn construct_polynomial_g(
        z_ccs: &[NTT],
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<
        (VirtualPolynomial<NTT>, Vec<DenseMultilinearExtension<NTT>>),
        LinearizationError<NTT>,
    > {
        // Generate beta challenges from Step 1
        let beta_s = transcript.squeeze_beta_challenges(ccs.s);

        // Prepare MLEs
        let Mz_mles = calculate_Mz_mles::<NTT, LinearizationError<NTT>>(ccs, z_ccs)?;

        // Construct the sumcheck polynomial g
        let g = prepare_lin_sumcheck_polynomial(ccs.s, &ccs.c, &Mz_mles, &ccs.S, &beta_s)?;

        Ok((g, Mz_mles))
    }

    /// Step 2: Run linearization sum-check protocol.
    fn generate_sumcheck_proof(
        g: &VirtualPolynomial<NTT>,
        transcript: &mut impl Transcript<NTT>,
        comb_fn: impl Fn(&[NTT]) -> NTT + Sync + Send,
    ) -> Result<(Proof<NTT>, Vec<NTT>), LinearizationError<NTT>> {
        let (sum_check_proof, prover_state) =
            MLSumcheck::prove_as_subprotocol(transcript, g, comb_fn);
        let point_r = prover_state
            .randomness
            .into_iter()
            .map(|x| x.into())
            .collect::<Vec<NTT>>();

        Ok((sum_check_proof, point_r))
    }

    /// Step 3: the mle evaluations that the prover sends to the verifier.
    /// I.e. f-hat rows mle evaluations and Mz mle evaluations.
    fn compute_evaluation_vectors(
        wit: &Witness<NTT>,
        point_r: &[NTT],
        Mz_mles: &[DenseMultilinearExtension<NTT>],
    ) -> Result<(Vec<NTT>, Vec<NTT>, Vec<NTT>), LinearizationError<NTT>> {
        // Compute v

        let v: Vec<NTT> = evaluate_mles::<NTT, _, _, LinearizationError<NTT>>(&wit.f_hat, point_r)?;

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

        let comb_fn = |vals: &[NTT]| -> NTT {
            let mut result = NTT::zero();
            'outer: for (i, &c) in ccs.c.iter().enumerate() {
                if c.is_zero() {
                    continue;
                }
                let mut term = c;
                for &j in &ccs.S[i] {
                    if vals[j].is_zero() {
                        continue 'outer;
                    }
                    term *= vals[j];
                }
                result += term;
            }
            // eq() is the last term added
            result * vals[vals.len() - 1]
        };

        // Run sumcheck protocol.
        let (sumcheck_proof, point_r) = Self::generate_sumcheck_proof(&g, transcript, comb_fn)?;

        // Step 3: Compute v, u_vector.
        let (point_r, v, u) = Self::compute_evaluation_vectors(wit, &point_r, &Mz_mles)?;

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
        let poly_info = VPAuxInfo::new(ccs.s, ccs.d + 1);

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

    /// Step 4: Verify
    /// $$
    /// \mathbf{e} \cdot \left( \sum\_{i=1}^{n\_s} c_i \cdot \prod\_{j \in S\_i} \mathbf{u}\_j \right) \stackrel{?}{=} s.
    /// $$
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
        let beta_s = transcript.squeeze_beta_challenges(ccs.s);

        //Step 2: The sumcheck.
        let (point_r, s) = Self::verify_sumcheck_proof(proof, transcript, ccs)?;

        Self::verify_evaluation_claim(&beta_s, &point_r, s, proof, ccs)?;

        // Absorbing the prover's mmessages to the verifier.
        transcript.absorb_slice(&proof.v);
        transcript.absorb_slice(&proof.u);

        // Step 5: Output z_o
        Ok(Self::prepare_verifier_output(cm_i, point_r, proof))
    }
}
