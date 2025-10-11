//! The folding module defines the behaviour for the folding subprotocol
//! of the LatticeFold protocol.
//!
//! <https://eprint.iacr.org/2024/257.pdf#page=45>

#![allow(non_snake_case)]

use ark_ff::Zero;
use ark_std::{cfg_iter, iter::successors, iterable::Iterable};
use cyclotomic_rings::rings::SuitableRing;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use stark_rings_poly::mle::DenseMultilinearExtension;

use self::utils::*;
use super::error::FoldingError;
use crate::{
    arith::{error::CSError, Witness, CCS, LCCCS},
    ark_base::*,
    commitment::Commitment,
    transcript::TranscriptWithShortChallenges,
    utils::{
        mle_helpers::evaluate_mles,
        sumcheck::{
            prover::ProverState,
            utils::{build_eq_x_r, eq_eval},
            MLSumcheck,
            SumCheckError::SumCheckFailed,
        },
    },
};

#[cfg(test)]
mod tests;

mod utils;
pub use structs::*;

mod structs;

impl<'t, NTT: SuitableRing, T: TranscriptWithShortChallenges<NTT>> FoldingProver<NTT, T>
    for LFFoldingProver<'t, NTT, T>
{
    fn prove(
        &mut self,
        cm_i_s: &[LCCCS<NTT>],
        mut w_s: Vec<Witness<NTT>>,
        ccs: &CCS<NTT>,
        mz_mles: &[Vec<DenseMultilinearExtension<NTT>>],
    ) -> Result<(LCCCS<NTT>, Witness<NTT>, FoldingProof<NTT>), FoldingError<NTT>> {
        sanity_check::<NTT>(ccs, self.dparams.l)?;

        // Free some unneeded vars/memory
        w_s.iter_mut().for_each(|w_i| {
            w_i.f_coeff.clear();
            w_i.f_coeff.shrink_to_fit();
            w_i.w_ccs.clear();
            w_i.w_ccs.shrink_to_fit();
        });

        if cm_i_s.len() != 2 * self.dparams.k {
            return Err(FoldingError::IncorrectLength);
        }

        let log_m = ccs.s;

        // Step 1: Generate alpha, zeta, mu, beta challenges
        let (alpha_s, beta_s, zeta_s, mu_s) = self
            .transcript
            .squeeze_alpha_beta_zeta_mu(log_m, self.dparams.k);

        // Step 2: Compute g polynomial and sumcheck on it
        // Setup f_hat_mle for later evaluation of thetas
        let f_hat_mles = Self::setup_f_hat_mles(&mut w_s);

        let ris = Self::get_ris(cm_i_s);

        let prechallenged_Ms_1 = Self::calculate_challenged_mz_mle(
            &mz_mles[0..self.dparams.k],
            &zeta_s[0..self.dparams.k],
        )?;
        let prechallenged_Ms_2 = Self::calculate_challenged_mz_mle(
            &mz_mles[self.dparams.k..2 * self.dparams.k],
            &zeta_s[self.dparams.k..2 * self.dparams.k],
        )?;
        let (g_mles, g_degree) = self.create_sumcheck_polynomial(
            log_m,
            f_hat_mles.clone(),
            &alpha_s,
            &prechallenged_Ms_1,
            &prechallenged_Ms_2,
            &ris,
            &beta_s,
            &mu_s,
        )?;

        let comb_fn = |vals: &[NTT]| -> NTT {
            sumcheck_polynomial_comb_fn::<NTT>(vals, &mu_s, self.dparams.b)
        };

        // Step 5: Run sum check prover
        let (sum_check_proof, prover_state) =
            MLSumcheck::prove_as_subprotocol(self.transcript, g_mles, log_m, g_degree, comb_fn);

        let r_0 = Self::get_sumcheck_randomness(prover_state);

        // Step 3: Evaluate thetas and etas
        let theta_s = Self::get_thetas(&f_hat_mles, &r_0)?;
        let eta_s = Self::get_etas(mz_mles, &r_0)?;

        // Absorb them into the transcript
        theta_s
            .iter()
            .for_each(|thetas| self.transcript.absorb_slice(thetas));
        eta_s
            .iter()
            .for_each(|etas| self.transcript.absorb_slice(etas));

        // Step 5 get rho challenges
        let (rho_s_coeff, rho_s) = get_rhos::<_, _>(self.dparams.k, self.transcript);

        let f_0: Vec<NTT> = Self::compute_f_0(&rho_s, &w_s);

        // Step 6 compute v0, u0, y0, x_w0
        let (v_0, cm_0, u_0, x_0) =
            compute_v0_u0_x0_cm_0(&rho_s_coeff, &rho_s, &theta_s, cm_i_s, &eta_s, ccs);

        // Step 7: Compute f0 and Witness_0
        let h = x_0.last().copied().ok_or(FoldingError::IncorrectLength)?;

        let lcccs = prepare_public_output(r_0, v_0, cm_0, u_0, x_0, h);

        let w_0 = Witness::from_f(f_0, self.dparams.B, self.dparams.l);

        let folding_proof = FoldingProof {
            pointshift_sumcheck_proof: sum_check_proof,
            theta_s,
            eta_s,
        };

        Ok((lcccs, w_0, folding_proof))
    }
}

impl<'t, NTT: SuitableRing, T: TranscriptWithShortChallenges<NTT>> FoldingVerifier<NTT, T>
    for LFFoldingVerifier<'t, NTT, T>
{
    fn verify(
        &mut self,
        cm_i_s: &[LCCCS<NTT>],
        proof: &FoldingProof<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<LCCCS<NTT>, FoldingError<NTT>> {
        sanity_check::<NTT>(ccs, self.dparams.l)?;

        // Step 1: Generate alpha, zeta, mu, beta challenges and validate input
        let (alpha_s, beta_s, zeta_s, mu_s) = self
            .transcript
            .squeeze_alpha_beta_zeta_mu(ccs.s, self.dparams.k);

        // Calculate claims for sumcheck verification
        let (claim_g1, claim_g3) = Self::calculate_claims(&alpha_s, &zeta_s, cm_i_s);

        let nvars = ccs.s;
        let degree = 2 * self.dparams.b;

        //Step 2: The sumcheck.
        let (r_0, expected_evaluation) =
            self.verify_sumcheck_proof(nvars, degree, claim_g1 + claim_g3, proof)?;

        // Verify evaluation claim
        self.verify_evaluation(
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
            .for_each(|thetas| self.transcript.absorb_slice(thetas));
        proof
            .eta_s
            .iter()
            .for_each(|etas| self.transcript.absorb_slice(etas));
        let (rho_s_coeff, rho_s) = get_rhos::<_, _>(self.dparams.k, self.transcript);

        // Step 6
        let (v_0, cm_0, u_0, x_0) = compute_v0_u0_x0_cm_0(
            &rho_s_coeff,
            &rho_s,
            &proof.theta_s,
            cm_i_s,
            &proof.eta_s,
            ccs,
        );

        // Step 7: Compute f0 and Witness_0

        let h = x_0.last().copied().ok_or(FoldingError::IncorrectLength)?;
        Ok(prepare_public_output(r_0, v_0, cm_0, u_0, x_0, h))
    }
}

impl<'t, NTT: SuitableRing, T> LFFoldingProver<'t, NTT, T> {
    /// Creates sumcheck polynomial
    ///
    /// $$
    /// g(\vec{x}) := \sum_{i=1}^{2k} \left[\alpha_i g_{1,i}(\vec{x}) + \mu_i g_{2,i}(\vec{x}) + \zeta_i g_{3,i}(\vec{x})\right]
    /// $$
    ///
    /// where, for all $i \in \[2k\]$,
    ///
    /// $$
    /// g_{1,i}(\vec{x}) := \sum\_{j=0}^{\tau - 1} \alpha_i^j \cdot \left( eq(\vec{r}_i, \vec{x}) \cdot \mathrm{mle} \[\hat{f}\_{ij}\](\vec{x}) \right),
    /// $$
    ///
    /// $$
    /// g_{2,i,j}(\vec{x}) := \sum\_{j=0}^{\tau - 1} \mu_i^j \left( eq(\vec{\beta}, \vec{x}) \cdot
    /// \prod_{j=-(b-1)}^{b-1} \( \mathrm{mle} \[\hat{f}\_{ij}\](\vec{x}) - j \)\right),
    /// $$
    ///
    /// $$
    /// g_{3,i}(\vec{x}) := \sum\_{j=0}^{t-1} \zeta_i^j \cdot \left(eq(\vec{r}_i, \vec{x}) \cdot
    /// \left(
    /// \sum\_{
    /// \vec{b} \in \\{0,1\\}^\{log\(n + n\_{in}\)\}
    /// }
    /// \text{mle}\[M_j\]\(\vec{x}, \vec{b}\) \cdot \text{mle}\[z_i\]\(\vec{b}\)
    /// \right)
    /// \right).
    /// $$
    ///
    /// # Arguments
    ///
    /// - `log_m: usize`  
    ///   The number of variables in the final polynomial.
    ///
    /// - `f_hat_mles: &[Vec<DenseMultilinearExtension<NTT>>]`  
    ///   A reference to the multilinear extension of the decomposed NTT witnesses
    ///
    /// - `alpha_s: &[NTT]`  
    ///   A slice containing the $\alpha$ challenges.
    ///
    /// - `challenged_Ms_1: &DenseMultilinearExtension<NTT>`  
    ///   A reference to the M matrices multiplied by the first $k$ decomposed vectors, and then taken a linear combination of.
    ///
    /// - `challenged_Ms_2: &DenseMultilinearExtension<NTT>`  
    ///   A reference to the M matrices multiplied by the second $k$ decomposed vectors, and then taken a linear combination of.
    ///
    /// - `r_s: &[Vec<NTT>]`  
    ///   The linearization challenge vectors
    ///
    /// - `beta_s: &[NTT]`  
    ///   The $\beta$ challenges
    ///
    /// - `mu_s: &[NTT]`  
    ///   The $\mu$ challenges
    ///
    /// # Returns
    ///
    /// - `Result<(Vec<RefCounter<DenseMultilinearExtension<NTT>>>, usize), FoldingError<NTT>>`  
    ///   - On success, returns a tuple containing:
    ///     - A `Vec<RefCounter<DenseMultilinearExtension<NTT>>>`, the MLEs that make up the polynomial.
    ///     - A `usize` of the degree of the final polynomial.
    ///
    /// # Errors
    ///
    /// This function will return a `FoldingError<NTT>` if any of the multilinear extensions or vectors are of the wrong size.
    ///
    /// $$
    #[allow(clippy::too_many_arguments)]
    pub(super) fn create_sumcheck_polynomial(
        &self,
        log_m: usize,
        f_hat_mles: Vec<Vec<DenseMultilinearExtension<NTT>>>,
        alpha_s: &[NTT],
        challenged_Ms_1: &DenseMultilinearExtension<NTT>,
        challenged_Ms_2: &DenseMultilinearExtension<NTT>,
        r_s: &[Vec<NTT>],
        beta_s: &[NTT],
        mu_s: &[NTT],
    ) -> Result<(Vec<DenseMultilinearExtension<NTT>>, usize), FoldingError<NTT>> {
        let k = self.dparams.k;

        if alpha_s.len() != 2 * k
            || f_hat_mles.len() != 2 * k
            || r_s.len() != 2 * k
            || beta_s.len() != log_m
            || mu_s.len() != 2 * k
        {
            return Err(FoldingError::IncorrectLength);
        }

        #[cfg(test)]
        {
            if r_s[..k].iter().any(|r| r != &r_s[0]) || r_s[k..].iter().any(|r| r != &r_s[k]) {
                return Err(FoldingError::SumcheckChallengeError);
            }
        }

        let len = 2 + 2 + // g1 + g3
            1 + f_hat_mles.len() * f_hat_mles[0].len(); // g2
        let mut mles = Vec::with_capacity(len);

        // We assume here that decomposition subprotocol puts the same r challenge point
        // into all decomposed linearized commitments
        let r_i_eq = build_eq_x_r(&r_s[0])?;
        prepare_g1_and_3_k_mles_list(
            &mut mles,
            r_i_eq.clone(),
            &f_hat_mles[0..k],
            &alpha_s[0..k],
            challenged_Ms_1,
        );

        let r_i_eq = build_eq_x_r(&r_s[k])?;
        prepare_g1_and_3_k_mles_list(
            &mut mles,
            r_i_eq,
            &f_hat_mles[k..2 * k],
            &alpha_s[k..2 * k],
            challenged_Ms_2,
        );

        // g2
        let beta_eq_x = build_eq_x_r(beta_s)?;
        prepare_g2_i_mle_list(&mut mles, beta_eq_x, f_hat_mles);

        let degree = 2 * self.dparams.b;

        Ok((mles, degree))
    }

    fn setup_f_hat_mles(w_s: &mut [Witness<NTT>]) -> Vec<Vec<DenseMultilinearExtension<NTT>>> {
        cfg_iter_mut!(w_s)
            .map(|w| w.take_f_hat())
            .collect::<Vec<Vec<DenseMultilinearExtension<NTT>>>>()
    }

    fn get_ris(cm_i_s: &[LCCCS<NTT>]) -> Vec<Vec<NTT>> {
        cm_i_s.iter().map(|cm_i| cm_i.r.clone()).collect::<Vec<_>>()
    }

    fn calculate_challenged_mz_mle(
        Mz_mles_vec: &[Vec<DenseMultilinearExtension<NTT>>],
        zeta_s: &[NTT],
    ) -> Result<DenseMultilinearExtension<NTT>, FoldingError<NTT>> {
        let mut combined_mle: DenseMultilinearExtension<NTT> = DenseMultilinearExtension::zero();

        zeta_s
            .iter()
            .zip(Mz_mles_vec)
            .for_each(|(zeta_i, Mz_mles)| {
                let mut mle: DenseMultilinearExtension<NTT> = DenseMultilinearExtension::zero();
                for M in Mz_mles.iter().rev() {
                    mle += M;
                    mle *= *zeta_i;
                }
                combined_mle += mle;
            });
        Ok(combined_mle)
    }

    fn get_sumcheck_randomness(sumcheck_prover_state: ProverState<NTT>) -> Vec<NTT> {
        sumcheck_prover_state
            .randomness
            .into_iter()
            .map(|x| x.into())
            .collect::<Vec<NTT>>()
    }

    fn get_thetas(
        f_hat_mles: &[Vec<DenseMultilinearExtension<NTT>>],
        r_0: &[NTT],
    ) -> Result<Vec<Vec<NTT>>, FoldingError<NTT>> {
        let theta_s: Vec<Vec<NTT>> = cfg_iter!(f_hat_mles)
            .map(|f_hat_row| evaluate_mles::<_, _, _, FoldingError<NTT>>(f_hat_row, r_0))
            .collect::<Result<Vec<_>, _>>()?;

        Ok(theta_s)
    }

    fn get_etas(
        Mz_mles_vec: &[Vec<DenseMultilinearExtension<NTT>>],
        r_0: &[NTT],
    ) -> Result<Vec<Vec<NTT>>, FoldingError<NTT>> {
        let eta_s: Vec<Vec<NTT>> = cfg_iter!(Mz_mles_vec)
            .map(|Mz_mles| evaluate_mles::<_, _, _, FoldingError<NTT>>(Mz_mles, r_0))
            .collect::<Result<Vec<_>, _>>()?;

        Ok(eta_s)
    }

    fn compute_f_0(rho_s: &[NTT], w_s: &[Witness<NTT>]) -> Vec<NTT> {
        rho_s
            .iter()
            .zip(w_s)
            .fold(vec![NTT::ZERO; w_s[0].f.len()], |acc, (&rho_i, w_i)| {
                acc.into_iter()
                    .zip(w_i.f.iter())
                    .map(|(acc_j, w_ij)| acc_j + rho_i * w_ij)
                    .collect()
            })
    }
}

impl<'t, NTT: SuitableRing, T: TranscriptWithShortChallenges<NTT>> LFFoldingVerifier<'t, NTT, T> {
    #[allow(clippy::too_many_arguments)]
    fn verify_evaluation(
        &self,
        alpha_s: &[NTT],
        beta_s: &[NTT],
        mu_s: &[NTT],
        zeta_s: &[NTT],
        r_0: &[NTT],
        expected_evaluation: NTT,
        proof: &FoldingProof<NTT>,
        cm_i_s: &[LCCCS<NTT>],
    ) -> Result<(), FoldingError<NTT>> {
        let ris = cm_i_s.iter().map(|cm_i| cm_i.r.clone()).collect::<Vec<_>>();

        let e_asterisk = eq_eval(beta_s, r_0)?;
        let e_s: Vec<NTT> = ris
            .iter()
            .map(|r_i: &Vec<NTT>| eq_eval(r_i, r_0))
            .collect::<Result<Vec<_>, _>>()?;

        let should_equal_s: NTT = self.compute_sumcheck_claim_expected_value(
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

    fn calculate_claims(alpha_s: &[NTT], zeta_s: &[NTT], cm_i_s: &[LCCCS<NTT>]) -> (NTT, NTT) {
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

    /// Computes the grand sum from point 4 of the Latticefold folding protocol.
    ///
    /// # Arguments
    ///
    /// - `alpha_s: &[NTT]`  
    ///   A slice containing the $\alpha$ challenges.
    ///
    /// - `mu_s: &[NTT]`  
    ///   A slice containing the $\mu$ challenges.
    ///
    /// - `theta_s: &[Vec<NTT>]`  
    ///   $$
    ///   \left[\theta\_{i} := \text{mle}\[\hat{f}\_i\](\vec{r}_o) \right]\_{i=1}^{2k},
    ///   $$
    ///
    /// - `e_asterisk: NTT`  
    ///   $$
    ///   \mathbf{e}^* := eq(\boldsymbol{\beta}, \mathbf{r}_o)
    ///   $$
    ///
    /// - `e_s: &[NTT]`  
    ///   $$
    ///   \left[ e_i := eq(\vec{r}\_i, \vec{r}\_o) \right]\_{i=1}^{2k}
    ///   $$
    /// - `zeta_s: &[NTT]`  
    ///
    ///     A slice containing the $\zeta$ challenges.
    ///
    /// - `eta_s: &[Vec<NTT>]`  
    ///   $$
    ///   \eta[i] :=
    ///   \sum\_{
    ///   \vec{b} \in \\{0,1\\}^\{log\(n + n\_{in}\)\}
    ///   }
    ///   \text{mle}\[M_1\]\(\vec{r}\_o, \vec{b}\) \cdot \text{mle}\[z_i\]\(\vec{b}\)
    ///   $$
    ///
    /// # Returns
    ///
    /// - `NTT`  
    ///   Returns the expected value of the sumcheck claim.
    ///
    #[allow(clippy::too_many_arguments)]
    pub(super) fn compute_sumcheck_claim_expected_value(
        &self,
        alpha_s: &[NTT],
        mu_s: &[NTT],
        theta_s: &[Vec<NTT>],
        e_asterisk: NTT,
        e_s: &[NTT],
        zeta_s: &[NTT],
        eta_s: &[Vec<NTT>],
    ) -> NTT {
        (0..(2 * self.dparams.k))
            .map(|i| {
                // Evaluation claims about f hats.
                let mut s_summand: NTT = successors(Some(alpha_s[i]), |alpha_power| {
                    Some(alpha_s[i] * alpha_power)
                })
                .zip(theta_s[i].iter())
                .map(|(pow_of_alpha_i, theta)| pow_of_alpha_i * e_s[i] * theta) // Might need to change e_s[i] double check
                .sum();

                // norm range check contribution
                s_summand += e_asterisk
                    * successors(Some(mu_s[i]), |mu_power| Some(mu_s[i] * mu_power))
                        .zip(theta_s[i].iter())
                        .map(|(mu_power, &theta)| {
                            mu_power
                                * theta
                                * (1..self.dparams.b)
                                    .map(|x| NTT::from(x as u128))
                                    .map(|j_hat| (theta - j_hat) * (theta + j_hat))
                                    .product::<NTT>()
                        })
                        .sum::<NTT>();

                // linearisation claims contribuition
                s_summand += e_s[i]
                    * successors(Some(zeta_s[i]), |&zeta| Some(zeta * zeta_s[i]))
                        .zip(eta_s[i].iter())
                        .map(|(pow_of_zeta, eta_i_j)| pow_of_zeta * eta_i_j)
                        .sum::<NTT>();

                s_summand
            })
            .sum()
    }

    fn verify_sumcheck_proof(
        &mut self,
        nvars: usize,
        degree: usize,
        total_claim: NTT,
        proof: &FoldingProof<NTT>,
    ) -> Result<(Vec<NTT>, NTT), FoldingError<NTT>> {
        //Step 2: The sumcheck.
        // Verify the sumcheck proof.
        let sub_claim = MLSumcheck::verify_as_subprotocol(
            self.transcript,
            nvars,
            degree,
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

fn sanity_check<NTT: SuitableRing>(ccs: &CCS<NTT>, l: usize) -> Result<(), FoldingError<NTT>> {
    if ccs.m != usize::max((ccs.n - ccs.l - 1) * l, ccs.m).next_power_of_two() {
        return Err(CSError::InvalidSizeBounds(ccs.m, ccs.n, l).into());
    }

    Ok(())
}

fn prepare_public_output<NTT: SuitableRing>(
    r_0: Vec<NTT>,
    v_0: Vec<NTT>,
    cm_0: Commitment<NTT>,
    u_0: Vec<NTT>,
    x_0: Vec<NTT>,
    h: NTT,
) -> LCCCS<NTT> {
    LCCCS {
        r: r_0,
        v: v_0,
        cm: cm_0,
        u: u_0,
        x_w: x_0[0..x_0.len() - 1].to_vec(),
        h,
    }
}
