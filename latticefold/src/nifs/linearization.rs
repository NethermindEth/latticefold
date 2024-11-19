use ark_ff::Field;
use ark_ff::PrimeField;
use ark_std::cfg_iter;
use cyclotomic_rings::rings::SuitableRing;
use lattirust_poly::polynomials::VirtualPolynomial;
use lattirust_poly::{
    mle::DenseMultilinearExtension,
    polynomials::{eq_eval, VPAuxInfo},
};
use lattirust_ring::OverField;
use utils::{compute_u, prepare_lin_sumcheck_polynomial};

use super::error::LinearizationError;
use crate::{
    arith::{utils::mat_vec_mul, Witness, CCCS, CCS, LCCCS},
    transcript::Transcript,
    utils::sumcheck::{MLSumcheck, SumCheckError::SumCheckFailed},
};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::utils::sumcheck::Proof;
pub use structs::*;

mod structs;

#[cfg(test)]
mod tests;
mod utils;

impl<NTT: SuitableRing, T: Transcript<NTT>> LinearizationProver<NTT, T>
    for LFLinearizationProver<NTT, T>
{
    fn compute_z_ccs<const C: usize>(
        wit: &Witness<NTT>,
        x_ccs: &[NTT],
    ) -> Result<Vec<NTT>, LinearizationError<NTT>> {
        let mut z_ccs = Vec::with_capacity(x_ccs.len() + 1 + wit.w_ccs.len());
        z_ccs.extend_from_slice(x_ccs);
        z_ccs.push(NTT::one());
        z_ccs.extend_from_slice(&wit.w_ccs);

        Ok(z_ccs)
    }

    fn construct_polynomial_g(
        z_ccs: &[NTT],
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(VirtualPolynomial<NTT>, Vec<NTT>), LinearizationError<NTT>> {
        let beta_s = BetaChallengeGenerator::generate_challenges(transcript, ccs.s);

        let Mz_mles = ccs
            .M
            .iter()
            .map(|M| {
                Ok(DenseMultilinearExtension::from_slice(
                    ccs.s,
                    &mat_vec_mul(M, z_ccs)?,
                ))
            })
            .collect::<Result<Vec<_>, LinearizationError<_>>>()?;

        let g = prepare_lin_sumcheck_polynomial(ccs.s, &ccs.c, &Mz_mles, &ccs.S, &beta_s)?;

        Ok((g, beta_s))
    }

    fn generate_sumcheck_proof(
        g: &VirtualPolynomial<NTT>,
        _beta_s: &[NTT],
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

    fn compute_evaluation_vectors(
        wit: &Witness<NTT>,
        point_r: &[NTT],
        ccs: &CCS<NTT>,
        z_ccs: &[NTT],
    ) -> Result<(Vec<NTT>, Vec<NTT>, Vec<NTT>), LinearizationError<NTT>> {
        let v: Vec<NTT> = cfg_iter!(wit.f_hat)
            .map(|f_hat_row| {
                DenseMultilinearExtension::from_slice(ccs.s, f_hat_row)
                    .evaluate(&point_r)
                    .expect("cannot end up here, because the sumcheck subroutine must yield a point of the length log m")
            })
            .collect();

        let Mz_mles = ccs
            .M
            .iter()
            .map(|M| {
                Ok(DenseMultilinearExtension::from_slice(
                    ccs.s,
                    &mat_vec_mul(M, z_ccs)?,
                ))
            })
            .collect::<Result<Vec<_>, LinearizationError<_>>>()?;

        let u = compute_u(&Mz_mles, point_r)?;

        Ok((point_r.to_vec(), v, u))
    }

    fn prove<const C: usize>(
        cm_i: &CCCS<C, NTT>,
        wit: &Witness<NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(LCCCS<C, NTT>, LinearizationProof<NTT>), LinearizationError<NTT>> {
        let z_state = Self::compute_z_ccs::<C>(wit, &cm_i.x_ccs)?;
        let (g, beta_s) = Self::construct_polynomial_g(&z_state, transcript, ccs)?;

        let (sumcheck_proof, point_r) = Self::generate_sumcheck_proof(&g, &beta_s, transcript)?;

        let (point_r, v, u) = Self::compute_evaluation_vectors(wit, &point_r, ccs, &z_state)?;

        transcript.absorb_slice(&v);
        transcript.absorb_slice(&u);

        let linearization_proof = LinearizationProof {
            linearization_sumcheck: sumcheck_proof,
            v: v.clone(),
            u: u.clone(),
        };

        let lcccs = LCCCS {
            r: point_r,
            v: v,
            cm: cm_i.cm.clone(),
            u: u,
            x_w: cm_i.x_ccs.clone(),
            h: NTT::one(),
        };

        Ok((lcccs, linearization_proof))
    }
}

impl<NTT: SuitableRing, T: Transcript<NTT>> LinearizationVerifier<NTT, T>
    for LFLinearizationVerifier<NTT, T>
{
    fn verify_sumcheck_proof(
        proof: &LinearizationProof<NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(Vec<NTT>, NTT), LinearizationError<NTT>> {
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

    fn verify_evaluation_claim(
        beta_s: &[NTT],
        point_r: &[NTT],
        s: NTT,
        proof: &LinearizationProof<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(), LinearizationError<NTT>> {
        let e = eq_eval(point_r, beta_s)?;
        let should_equal_s = e * ccs
            .c
            .iter()
            .enumerate()
            .map(|(i, &c)| c * ccs.S[i].iter().map(|&j| proof.u[j]).product::<NTT>())
            .sum::<NTT>();

        if should_equal_s != s {
            return Err(LinearizationError::SumCheckError(SumCheckFailed(
                should_equal_s,
                s,
            )));
        }

        Ok(())
    }

    fn prepare_final_state<const C: usize>(
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

    fn verify<const C: usize>(
        cm_i: &CCCS<C, NTT>,
        proof: &LinearizationProof<NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<LCCCS<C, NTT>, LinearizationError<NTT>> {
        let beta_s = BetaChallengeGenerator::generate_challenges(transcript, ccs.s);

        let (point_r, s) = Self::verify_sumcheck_proof(proof, transcript, ccs)?;

        Self::verify_evaluation_claim(&beta_s, &point_r, s, proof, ccs)?;

        transcript.absorb_slice(&proof.v);
        transcript.absorb_slice(&proof.u);

        Ok(Self::prepare_final_state(cm_i, point_r, proof))
    }
}

impl<NTT: OverField> ChallengeGenerator<NTT> for BetaChallengeGenerator<NTT> {
    fn generate_challenges(transcript: &mut impl Transcript<NTT>, log_m: usize) -> Vec<NTT> {
        transcript.absorb_field_element(&<NTT::BaseRing as Field>::from_base_prime_field(
            <NTT::BaseRing as Field>::BasePrimeField::from_be_bytes_mod_order(b"beta_s"),
        ));

        transcript
            .get_challenges(log_m)
            .into_iter()
            .map(|x| x.into())
            .collect()
    }
}
