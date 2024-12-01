#![allow(non_snake_case)]

use ark_std::cfg_iter;
use ark_std::iter::successors;
use ark_std::iterable::Iterable;
use cyclotomic_rings::rings::SuitableRing;
use lattirust_ring::cyclotomic_ring::CRT;

use super::error::FoldingError;
use super::mle_helpers::{calculate_Mz_mles, evaluate_mles, to_mles};
use crate::ark_base::*;
use crate::transcript::TranscriptWithShortChallenges;
use crate::utils::sumcheck::{
    virtual_polynomial::{eq_eval, VPAuxInfo},
    MLSumcheck,
    SumCheckError::SumCheckFailed,
};
use crate::{
    arith::{error::CSError, Instance, Witness, CCS, LCCCS},
    decomposition_parameters::DecompositionParams,
};

use lattirust_poly::mle::DenseMultilinearExtension;
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

fn prepare_public_output<const C: usize, NTT: SuitableRing>(
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

impl<NTT: SuitableRing, T: TranscriptWithShortChallenges<NTT>> LFFoldingProver<NTT, T> {
    fn setup_f_hat_mles(
        log_m: usize,
        w_s: &[Witness<NTT>],
    ) -> Vec<Vec<DenseMultilinearExtension<NTT>>> {
        w_s.iter()
            .map(|w| to_mles::<_, _, FoldingError<NTT>>(log_m, &w.f_hat).unwrap())
            .collect::<Vec<Vec<DenseMultilinearExtension<NTT>>>>()
    }

    fn get_zis<const C: usize>(cm_i_s: &[LCCCS<C, NTT>], w_s: &[Witness<NTT>]) -> Vec<Vec<NTT>> {
        cm_i_s
            .iter()
            .zip(w_s.iter())
            .map(|(cm_i, w_i)| cm_i.get_z_vector(&w_i.w_ccs))
            .collect::<Vec<_>>()
    }

    fn get_ris<const C: usize>(cm_i_s: &[LCCCS<C, NTT>]) -> Vec<Vec<NTT>> {
        cm_i_s.iter().map(|cm_i| cm_i.r.clone()).collect::<Vec<_>>()
    }

    fn calculate_Mz_mles(
        ccs: &CCS<NTT>,
        zis: &Vec<Vec<NTT>>,
    ) -> Result<Vec<Vec<DenseMultilinearExtension<NTT>>>, FoldingError<NTT>> {
        cfg_iter!(zis)
            .map(|zi| calculate_Mz_mles::<NTT, FoldingError<NTT>>(ccs, zi))
            .collect::<Result<_, FoldingError<_>>>()
    }

    fn get_sumcheck_randomness(sumcheck_prover_state: ProverState<NTT>) -> Vec<NTT> {
        sumcheck_prover_state
            .randomness
            .into_iter()
            .map(|x| x.into())
            .collect::<Vec<NTT>>()
    }

    fn get_thetas(
        f_hat_mles: &Vec<Vec<DenseMultilinearExtension<NTT>>>,
        r_0: &[NTT],
    ) -> Result<Vec<Vec<NTT>>, FoldingError<NTT>> {
        let theta_s: Vec<Vec<NTT>> = cfg_iter!(f_hat_mles)
            .map(|f_hat_row| evaluate_mles::<_, _, _, FoldingError<NTT>>(f_hat_row, r_0))
            .collect::<Result<Vec<_>, _>>()?;

        Ok(theta_s)
    }

    fn get_etas(
        Mz_mles_vec: &Vec<Vec<DenseMultilinearExtension<NTT>>>,
        r_0: &[NTT],
    ) -> Result<Vec<Vec<NTT>>, FoldingError<NTT>> {
        let eta_s: Vec<Vec<NTT>> = cfg_iter!(Mz_mles_vec)
            .map(|Mz_mles| evaluate_mles::<_, _, _, FoldingError<NTT>>(Mz_mles, r_0))
            .collect::<Result<Vec<_>, _>>()?;

        Ok(eta_s)
    }

    fn compute_f_0(rho_s: &Vec<NTT::CoefficientRepresentation>, w_s: &[Witness<NTT>]) -> Vec<NTT> {
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
        sanity_check::<NTT, P>(ccs)?;

        if cm_i_s.len() != 2 * P::K {
            return Err(FoldingError::IncorrectLength);
        }

        let log_m = ccs.s;

        // Step 1: Generate alpha, zeta, mu, beta challenges
        let (alpha_s, beta_s, zeta_s, mu_s) = transcript.squeeze_alpha_beta_zeta_mu::<P>(log_m);

        // Step 2: Compute g polynomial and sumcheck on it
        // Setup f_hat_mle for later evaluation of thetas
        let f_hat_mles = Self::setup_f_hat_mles(log_m, w_s);

        let zis = Self::get_zis(cm_i_s, w_s);
        let ris = Self::get_ris(cm_i_s);

        let Mz_mles_vec: Vec<Vec<DenseMultilinearExtension<NTT>>> =
            Self::calculate_Mz_mles(ccs, &zis)?;

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

        #[cfg(feature = "jolt-sumcheck")]
        let comb_fn = |_: &ProverState<NTT>, vals: &[NTT]| -> NTT {
            let mut result = vals[0] * vals[1]
                + vals[P::K * P::B_SMALL as usize + 3] * vals[P::K * P::B_SMALL as usize + 4];
            result += &vals[2..P::K * P::B_SMALL as usize + 2].iter().product();
            result += &vals[P::K * P::B_SMALL as usize + 5..2 * P::K * P::B_SMALL as usize + 5]
                .iter()
                .product();
            result
        };

        // Step 5: Run sum check prover
        let (sum_check_proof, prover_state) = MLSumcheck::prove_as_subprotocol(
            transcript,
            &g,
            #[cfg(feature = "jolt-sumcheck")]
            comb_fn,
        );

        let r_0 = Self::get_sumcheck_randomness(prover_state);

        // Step 3: Evaluate thetas and etas
        let theta_s = Self::get_thetas(&f_hat_mles, &r_0)?;
        let eta_s = Self::get_etas(&Mz_mles_vec, &r_0)?;

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

        let lcccs = prepare_public_output(r_0, v_0, cm_0, u_0, x_0, h);

        let f_0: Vec<NTT> = Self::compute_f_0(&rho_s, w_s);

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
    #[allow(clippy::too_many_arguments)]
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
        sanity_check::<NTT, P>(ccs)?;

        // Step 1: Generate alpha, zeta, mu, beta challenges and validate input
        let (alpha_s, beta_s, zeta_s, mu_s) = transcript.squeeze_alpha_beta_zeta_mu::<P>(ccs.s);

        // Calculate claims for sumcheck verification
        let (claim_g1, claim_g3) = Self::calculate_claims(&alpha_s, &zeta_s, cm_i_s);

        let poly_info = VPAuxInfo::new(ccs.s, 2 * P::B_SMALL);

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
        Ok(prepare_public_output(r_0, v_0, cm_0, u_0, x_0, h))
    }
}

fn sanity_check<NTT: SuitableRing, DP: DecompositionParams>(
    ccs: &CCS<NTT>,
) -> Result<(), FoldingError<NTT>> {
    if ccs.m != usize::max((ccs.n - ccs.l - 1) * DP::L, ccs.m).next_power_of_two() {
        return Err(CSError::InvalidSizeBounds(ccs.m, ccs.n, DP::L).into());
    }

    Ok(())
}
