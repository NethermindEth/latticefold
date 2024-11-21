#![allow(non_snake_case)]
use ark_std::cfg_iter;
use ark_std::iter::successors;
use ark_std::iterable::Iterable;
use cyclotomic_rings::rings::SuitableRing;
use lattirust_ring::cyclotomic_ring::CRT;
use utils::get_alphas_betas_zetas_mus;

use super::error::FoldingError;
use crate::transcript::TranscriptWithShortChallenges;
use crate::utils::sumcheck::{MLSumcheck, SumCheckError::SumCheckFailed};
use crate::{
    arith::{utils::mat_vec_mul, Instance, Witness, CCS, LCCCS},
    decomposition_parameters::DecompositionParams,
};

use lattirust_poly::{
    mle::DenseMultilinearExtension,
    polynomials::{eq_eval, VPAuxInfo},
};
use utils::*;

use crate::commitment::Commitment;
use crate::utils::sumcheck::prover::ProverState;
#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[cfg(test)]
mod tests;

mod utils;
pub use structs::*;
mod structs;

impl<NTT: SuitableRing, T: TranscriptWithShortChallenges<NTT>> LFFoldingProver<NTT, T> {
    fn setup_f_hat_mles(
        log_m: usize,
        w_s: &[Witness<NTT>],
    ) -> Vec<Vec<DenseMultilinearExtension<NTT>>> {
        w_s.iter()
            .map(|w| {
                w.f_hat
                    .iter()
                    .map(|f_hat_row| {
                        DenseMultilinearExtension::from_evaluations_slice(log_m, f_hat_row)
                    })
                    .collect()
            })
            .collect::<Vec<Vec<DenseMultilinearExtension<NTT>>>>()
    }

    fn get_zis_ris<const C: usize>(
        cm_i_s: &[LCCCS<C, NTT>],
        w_s: &[Witness<NTT>],
    ) -> (Vec<Vec<NTT>>, Vec<Vec<NTT>>) {
        let zis = cm_i_s
            .iter()
            .zip(w_s.iter())
            .map(|(cm_i, w_i)| cm_i.get_z_vector(&w_i.w_ccs))
            .collect::<Vec<_>>();

        let ris = cm_i_s.iter().map(|cm_i| cm_i.r.clone()).collect::<Vec<_>>();

        (zis, ris)
    }
    fn calculate_Mz_mles(
        ccs: &CCS<NTT>,
        zis: Vec<Vec<NTT>>,
    ) -> Result<Vec<Vec<DenseMultilinearExtension<NTT>>>, FoldingError<NTT>> {
        zis.iter()
            .map(|zi| {
                let Mz_mle = ccs
                    .M
                    .iter()
                    .map(|M| {
                        Ok(DenseMultilinearExtension::from_slice(
                            ccs.s,
                            &mat_vec_mul(M, zi)?,
                        ))
                    })
                    .collect::<Result<_, FoldingError<_>>>()?;
                Ok(Mz_mle)
            })
            .collect::<Result<_, FoldingError<_>>>()
    }

    fn sample_randomness(state: ProverState<NTT>) -> Vec<NTT> {
        state
            .randomness
            .into_iter()
            .map(|x| x.into())
            .collect::<Vec<NTT>>()
    }

    fn get_thetas_etas(
        f_hat_mles: Vec<Vec<DenseMultilinearExtension<NTT>>>,
        Mz_mles_vec: Vec<Vec<DenseMultilinearExtension<NTT>>>,
        r_0: &Vec<NTT>,
    ) -> Result<(Vec<Vec<NTT>>, Vec<Vec<NTT>>), FoldingError<NTT>> {
        let theta_s: Vec<Vec<NTT>> = cfg_iter!(f_hat_mles)
            .map(|f_hat_row| {
                f_hat_row
                    .iter()
                    .map(|f_hat_mle| {
                        f_hat_mle
                            .evaluate(r_0)
                            .ok_or(FoldingError::<NTT>::EvaluationError("f_hat".to_string()))
                    })
                    .collect::<Result<Vec<_>, _>>()
            })
            .collect::<Result<Vec<_>, _>>()?;

        let eta_s: Vec<Vec<NTT>> = cfg_iter!(Mz_mles_vec)
            .map(|Mz_mles| {
                Mz_mles
                    .iter()
                    .map(|mle| {
                        mle.evaluate(r_0.as_slice())
                            .ok_or(FoldingError::<NTT>::EvaluationError("Mz".to_string()))
                    })
                    .collect::<Result<Vec<_>, _>>()
            })
            .collect::<Result<Vec<_>, _>>()?;

        Ok((theta_s, eta_s))
    }

    fn prepare_lccs<const C: usize>(
        r_0: Vec<NTT>,
        v_0: Vec<NTT>,
        cm_0: Commitment<C, NTT>,
        u_0: Vec<NTT>,
        x_0: Vec<NTT>,
        h: NTT,
    ) -> LCCCS<C, NTT> {
        LCCCS {
            r: r_0,
            v: v_0,
            cm: cm_0,
            u: u_0,
            x_w: x_0[0..x_0.len() - 1].to_vec(),
            h,
        }
    }

    fn compute_f_0(rho_s: Vec<NTT::CoefficientRepresentation>, w_s: &[Witness<NTT>]) -> Vec<NTT> {
        rho_s
            .iter()
            .zip(w_s)
            .fold(vec![NTT::ZERO; w_s[0].f.len()], |acc, (&rho_i, w_i)| {
                let rho_i: NTT = rho_i.crt();

                acc.into_iter()
                    .zip(w_i.f.iter())
                    .map(|(acc_j, w_ij)| acc_j + rho_i * w_ij)
                    .collect()
            })
    }
}
impl<NTT: SuitableRing, T: TranscriptWithShortChallenges<NTT>> FoldingProver<NTT, T>
    for LFFoldingProver<NTT, T>
{
    fn prove<const C: usize, P: DecompositionParams>(
        cm_i_s: &[LCCCS<C, NTT>],
        w_s: &[Witness<NTT>],
        transcript: &mut impl TranscriptWithShortChallenges<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(LCCCS<C, NTT>, Witness<NTT>, FoldingProof<NTT>), FoldingError<NTT>> {
        if cm_i_s.len() != 2 * P::K {
            return Err(FoldingError::IncorrectLength);
        }

        let log_m = ccs.s;

        // Step 1: Generate alpha, zeta, mu, beta challenges
        let (alpha_s, beta_s, zeta_s, mu_s) =
            get_alphas_betas_zetas_mus::<_, _, P>(log_m, transcript);

        // Step 2: Compute g polynomial and sumcheck on it
        // Setup f_hat_mle for later evaluation of thetas
        let f_hat_mles = Self::setup_f_hat_mles(log_m, w_s);

        let (zis, ris) = Self::get_zis_ris(cm_i_s, w_s);

        let Mz_mles_vec: Vec<Vec<DenseMultilinearExtension<NTT>>> =
            Self::calculate_Mz_mles(ccs, zis)?;

        let g = create_sumcheck_polynomial::<_, P>(
            log_m,
            &f_hat_mles,
            &alpha_s,
            &Mz_mles_vec,
            &zeta_s,
            &ris,
            &beta_s,
            &mu_s,
        )?;

        // Step 5: Run sum check prover
        let (sum_check_proof, prover_state) = MLSumcheck::prove_as_subprotocol(transcript, &g);

        let r_0 = Self::sample_randomness(prover_state);

        // Step 3: Evaluate thetas and etas
        let (theta_s, eta_s) = Self::get_thetas_etas(f_hat_mles, Mz_mles_vec, &r_0)?;

        // Absorb them into the transcript
        theta_s
            .iter()
            .for_each(|thetas| transcript.absorb_slice(thetas));
        eta_s.iter().for_each(|etas| transcript.absorb_slice(etas));

        // Step 5 get rho challenges
        let rho_s = get_rhos::<_, _, P>(transcript);

        // Step 6 compute v0, u0, y0, x_w0
        let (v_0, cm_0, u_0, x_0) = compute_v0_u0_x0_cm_0(&rho_s, &theta_s, cm_i_s, &eta_s, ccs);

        // Step 7: Compute f0 and Witness_0
        let h = x_0.last().copied().ok_or(FoldingError::IncorrectLength)?;

        let lcccs = Self::prepare_lccs(r_0, v_0, cm_0, u_0, x_0, h);

        let f_0: Vec<NTT> = Self::compute_f_0(rho_s, w_s);

        let w_0 = Witness::from_f::<P>(f_0);

        let folding_proof = FoldingProof {
            pointshift_sumcheck_proof: sum_check_proof,
            theta_s,
            eta_s,
        };

        Ok((lcccs, w_0, folding_proof))
    }
}

impl<NTT: SuitableRing, T: TranscriptWithShortChallenges<NTT>> LFFoldingVerifier<NTT, T> {
    fn verify_evaluation<const C: usize, P: DecompositionParams>(
        alpha_s: &[NTT],
        beta_s: &[NTT],
        mu_s: &[NTT],
        zeta_s: &[NTT],
        r_0: &[NTT],
        expected_evaluation: NTT,
        proof: &FoldingProof<NTT>,
        cm_i_s: &[LCCCS<C, NTT>],
    ) -> Result<(), FoldingError<NTT>> {
        let ris = cm_i_s.iter().map(|cm_i| cm_i.r.clone()).collect::<Vec<_>>();

        let e_asterisk = eq_eval(beta_s, r_0)?;
        let e_s: Vec<NTT> = ris
            .iter()
            .map(|r_i: &Vec<NTT>| eq_eval(r_i, r_0))
            .collect::<Result<Vec<_>, _>>()?;

        let should_equal_s: NTT = compute_sumcheck_claim_expected_value::<NTT, P>(
            alpha_s,
            mu_s,
            &proof.theta_s,
            e_asterisk,
            &e_s,
            zeta_s,
            &proof.eta_s,
        );

        if should_equal_s != expected_evaluation {
            return Err(FoldingError::SumCheckError(SumCheckFailed(
                should_equal_s,
                expected_evaluation,
            )));
        }

        Ok(())
    }

    // Previously refactored methods stay the same
    fn get_alphas_betas_zetas_mus<const C: usize, P: DecompositionParams>(
        cm_i_s: &[LCCCS<C, NTT>],
        log_m: usize,
        transcript: &mut impl TranscriptWithShortChallenges<NTT>,
    ) -> Result<(Vec<NTT>, Vec<NTT>, Vec<NTT>, Vec<NTT>, VPAuxInfo<NTT>), FoldingError<NTT>> {
        // Validate input length
        if cm_i_s.len() != 2 * P::K {
            return Err(FoldingError::IncorrectLength);
        }

        let log_m = log_m;

        // Step 1: Generate alpha, zeta, mu, beta challenges
        let (alpha_s, beta_s, zeta_s, mu_s) =
            get_alphas_betas_zetas_mus::<_, _, P>(log_m, transcript);

        let poly_info = VPAuxInfo::new(log_m, 2 * P::B_SMALL);

        Ok((alpha_s, beta_s, zeta_s, mu_s, poly_info))
    }

    fn calculate_claims<const C: usize>(
        alpha_s: &[NTT],
        zeta_s: &[NTT],
        cm_i_s: &[LCCCS<C, NTT>],
    ) -> (NTT, NTT) {
        let vs = cm_i_s
            .iter()
            .map(|cm_i| cm_i.v.clone())
            .collect::<Vec<Vec<NTT>>>();
        let us = cm_i_s.iter().map(|cm_i| cm_i.u.clone()).collect::<Vec<_>>();

        // Calculate claim_g1
        let claim_g1: NTT = alpha_s
            .iter()
            .zip(vs.iter())
            .map(|(&alpha_i, v_i)| {
                successors(Some(alpha_i), |&alpha| Some(alpha * alpha_i))
                    .zip(v_i.iter())
                    .map(|(pow_of_alpha, v_i_j)| pow_of_alpha * v_i_j)
                    .sum::<NTT>()
            })
            .sum();

        // Calculate claim_g3
        let claim_g3: NTT = zeta_s
            .iter()
            .zip(us.iter())
            .map(|(&zeta_i, ui)| {
                successors(Some(zeta_i), |&zeta| Some(zeta * zeta_i))
                    .zip(ui.iter())
                    .map(|(pow_of_zeta, u_i_j)| pow_of_zeta * u_i_j)
                    .sum::<NTT>()
            })
            .sum();

        (claim_g1, claim_g3)
    }

    fn verify_sumcheck_proof(
        transcript: &mut impl TranscriptWithShortChallenges<NTT>,
        poly_info: &VPAuxInfo<NTT>,
        total_claim: NTT,
        proof: &FoldingProof<NTT>,
    ) -> Result<(Vec<NTT>, NTT), FoldingError<NTT>> {
        //Step 2: The sumcheck.
        // Verify the sumcheck proof.
        let sub_claim = MLSumcheck::verify_as_subprotocol(
            transcript,
            poly_info,
            total_claim,
            &proof.pointshift_sumcheck_proof,
        )?;

        let r_0 = sub_claim
            .point
            .into_iter()
            .map(|x| x.into())
            .collect::<Vec<NTT>>();

        Ok((r_0, sub_claim.expected_evaluation))
    }
}

impl<NTT: SuitableRing, T: TranscriptWithShortChallenges<NTT>> FoldingVerifier<NTT, T>
    for LFFoldingVerifier<NTT, T>
{
    fn verify<const C: usize, P: DecompositionParams>(
        cm_i_s: &[LCCCS<C, NTT>],
        proof: &FoldingProof<NTT>,
        transcript: &mut impl TranscriptWithShortChallenges<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<LCCCS<C, NTT>, FoldingError<NTT>> {
        // Step 1: Generate alpha, zeta, mu, beta challenges and validate input
        let (alpha_s, beta_s, zeta_s, mu_s, poly_info) =
            Self::get_alphas_betas_zetas_mus::<C, P>(cm_i_s, ccs.s, transcript)?;

        // Calculate claims for sumcheck verification
        let (claim_g1, claim_g3) = Self::calculate_claims(&alpha_s, &zeta_s, cm_i_s);

        //Step 2: The sumcheck.
        let (r_0, expected_evaluation) =
            Self::verify_sumcheck_proof(transcript, &poly_info, claim_g1 + claim_g3, proof)?;

        // Verify evaluation claim
        Self::verify_evaluation::<C, P>(
            &alpha_s,
            &beta_s,
            &mu_s,
            &zeta_s,
            &r_0,
            expected_evaluation,
            proof,
            cm_i_s,
        )?;

        // Step 5
        proof
            .theta_s
            .iter()
            .for_each(|thetas| transcript.absorb_slice(thetas));
        proof
            .eta_s
            .iter()
            .for_each(|etas| transcript.absorb_slice(etas));
        let rho_s = get_rhos::<_, _, P>(transcript);

        // Step 6
        let (v_0, cm_0, u_0, x_0) =
            compute_v0_u0_x0_cm_0(&rho_s, &proof.theta_s, cm_i_s, &proof.eta_s, ccs);

        // Step 7: Compute f0 and Witness_0
        let h = x_0.last().copied().ok_or(FoldingError::IncorrectLength)?;

        Ok(LCCCS {
            r: r_0,
            v: v_0,
            cm: cm_0,
            u: u_0,
            x_w: x_0[0..x_0.len() - 1].to_vec(),
            h,
        })
    }
}
