#![allow(non_snake_case)]

use ark_ff::Zero;
use ark_ff::{Field, PrimeField};
use ark_std::iter::successors;
use ark_std::iterable::Iterable;

// use ark_std::sync::Arc;
use cyclotomic_rings::{rings::SuitableRing, rotation::rot_lin_combination};
use lattirust_poly::polynomials::{ArithErrors, RefCounter};
use lattirust_ring::{cyclotomic_ring::CRT, Ring};

use crate::ark_base::*;
use crate::commitment::Commitment;
use crate::nifs::error::FoldingError;
use crate::transcript::TranscriptWithShortChallenges;
use crate::utils::sumcheck::dense_polynomial::{build_eq_x_r, DensePolynomial};
use crate::{
    arith::{CCS, LCCCS},
    decomposition_parameters::DecompositionParams,
    transcript::Transcript,
};
use lattirust_poly::mle::DenseMultilinearExtension;
use lattirust_ring::OverField;

pub(crate) trait SqueezeAlphaBetaZetaMu<NTT: SuitableRing> {
    fn squeeze_alpha_beta_zeta_mu<P: DecompositionParams>(
        &mut self,
        log_m: usize,
    ) -> (Vec<NTT>, Vec<NTT>, Vec<NTT>, Vec<NTT>);
}

impl<NTT: SuitableRing, T: Transcript<NTT>> SqueezeAlphaBetaZetaMu<NTT> for T {
    fn squeeze_alpha_beta_zeta_mu<P: DecompositionParams>(
        &mut self,
        log_m: usize,
    ) -> (Vec<NTT>, Vec<NTT>, Vec<NTT>, Vec<NTT>) {
        self.absorb_field_element(&<NTT::BaseRing as Field>::from_base_prime_field(
            <NTT::BaseRing as Field>::BasePrimeField::from_be_bytes_mod_order(b"alpha_s"),
        ));
        let alpha_s = self
            .get_challenges(2 * P::K)
            .into_iter()
            .map(|x| NTT::from(x))
            .collect::<Vec<_>>();

        self.absorb_field_element(&<NTT::BaseRing as Field>::from_base_prime_field(
            <NTT::BaseRing as Field>::BasePrimeField::from_be_bytes_mod_order(b"zeta_s"),
        ));
        let zeta_s = self
            .get_challenges(2 * P::K)
            .into_iter()
            .map(|x| NTT::from(x))
            .collect::<Vec<_>>();

        self.absorb_field_element(&<NTT::BaseRing as Field>::from_base_prime_field(
            <NTT::BaseRing as Field>::BasePrimeField::from_be_bytes_mod_order(b"mu_s"),
        ));
        let mut mu_s = self
            .get_challenges((2 * P::K) - 1)
            .into_iter()
            .map(|x| NTT::from(x))
            .collect::<Vec<_>>(); // Note is one challenge less

        mu_s.push(NTT::ONE);

        self.absorb_field_element(&<NTT::BaseRing as Field>::from_base_prime_field(
            <NTT::BaseRing as Field>::BasePrimeField::from_be_bytes_mod_order(b"beta_s"),
        ));
        let beta_s = self
            .get_challenges(log_m)
            .into_iter()
            .map(|x| NTT::from(x))
            .collect::<Vec<_>>();

        (alpha_s, beta_s, zeta_s, mu_s)
    }
}

pub(super) fn get_rhos<
    R: SuitableRing,
    T: TranscriptWithShortChallenges<R>,
    P: DecompositionParams,
>(
    transcript: &mut T,
) -> Vec<R::CoefficientRepresentation> {
    transcript.absorb_field_element(&<R::BaseRing as Field>::from_base_prime_field(
        <R::BaseRing as Field>::BasePrimeField::from_be_bytes_mod_order(b"rho_s"),
    ));

    let mut rhos = transcript.get_small_challenges((2 * P::K) - 1); // Note that we are missing the first element
    rhos.push(R::CoefficientRepresentation::ONE);

    rhos
}

#[allow(clippy::too_many_arguments)]
pub(super) fn create_sumcheck_polynomial<NTT: OverField, DP: DecompositionParams>(
    log_m: usize,
    f_hat_mles: &[Vec<DenseMultilinearExtension<NTT>>],
    alpha_s: &[NTT],
    challenged_Ms_1: &DenseMultilinearExtension<NTT>,
    challenged_Ms_2: &DenseMultilinearExtension<NTT>,
    r_s: &[Vec<NTT>],
    beta_s: &[NTT],
    mu_s: &[NTT],
) -> Result<DensePolynomial<NTT>, FoldingError<NTT>> {
    if alpha_s.len() != 2 * DP::K
        || f_hat_mles.len() != 2 * DP::K
        || r_s.len() != 2 * DP::K
        || beta_s.len() != log_m
        || mu_s.len() != 2 * DP::K
    {
        return Err(FoldingError::IncorrectLength);
    }

    #[cfg(test)]
    {
        if r_s[..DP::K].iter().any(|r| r != &r_s[0])
            || r_s[DP::K..].iter().any(|r| r != &r_s[DP::K])
        {
            return Err(FoldingError::SumcheckChallengeError);
        }
    }

    let mut g = DensePolynomial::<NTT>::new(log_m);

    let beta_eq_x = build_eq_x_r(beta_s)?;

    let f_hat_mles: Vec<Vec<RefCounter<DenseMultilinearExtension<NTT>>>> = f_hat_mles
        .iter()
        .map(|f_hat_mles_i| {
            f_hat_mles_i
                .clone()
                .into_iter()
                .map(RefCounter::from)
                .collect::<Vec<_>>()
        })
        .collect();

    // We assume here that decomposition subprotocol puts the same r challenge point
    // into all decomposed linearized commitments
    let r_i_eq = build_eq_x_r(&r_s[0])?;
    prepare_g1_and_3_k_mles_list(
        &mut g,
        r_i_eq.clone(),
        &f_hat_mles[0..DP::K],
        &alpha_s[0..DP::K],
        challenged_Ms_1,
    )?;

    let r_i_eq = build_eq_x_r(&r_s[DP::K])?;
    prepare_g1_and_3_k_mles_list(
        &mut g,
        r_i_eq,
        &f_hat_mles[DP::K..2 * DP::K],
        &alpha_s[DP::K..2 * DP::K],
        challenged_Ms_2,
    )?;

    // g2
    g.mul_mle(beta_eq_x)?;

    for mles in f_hat_mles.iter().take(DP::K) {
        prepare_g2_i_mle_list(&mut g, mles)?;
    }

    for mles in f_hat_mles.iter().take(2 * DP::K).skip(DP::K) {
        prepare_g2_i_mle_list(&mut g, mles)?;
    }

    g.aux_info.max_degree = 2 * DP::B_SMALL;

    Ok(g)
}

/// The grand sum from point 4 of the Latticefold folding protocol.
pub(super) fn compute_sumcheck_claim_expected_value<NTT: Ring, P: DecompositionParams>(
    alpha_s: &[NTT],
    mu_s: &[NTT],
    theta_s: &[Vec<NTT>],
    e_asterisk: NTT,
    e_s: &[NTT],
    zeta_s: &[NTT],
    eta_s: &[Vec<NTT>],
) -> NTT {
    (0..(2 * P::K))
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
                            * (1..P::B_SMALL)
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

pub(super) fn compute_v0_u0_x0_cm_0<const C: usize, NTT: SuitableRing>(
    rho_s: &[NTT::CoefficientRepresentation],
    theta_s: &[Vec<NTT>],
    cm_i_s: &[LCCCS<C, NTT>],
    eta_s: &[Vec<NTT>],
    ccs: &CCS<NTT>,
) -> (Vec<NTT>, Commitment<C, NTT>, Vec<NTT>, Vec<NTT>) {
    let v_0: Vec<NTT> = rot_lin_combination(rho_s, theta_s);

    let cm_0: Commitment<C, NTT> = rho_s
        .iter()
        .zip(cm_i_s.iter())
        .map(|(&rho_i, cm_i)| cm_i.cm.clone() * rho_i.crt())
        .sum();

    let u_0: Vec<NTT> = rho_s
        .iter()
        .zip(eta_s.iter())
        .map(|(&rho_i, etas_i)| {
            etas_i
                .iter()
                .map(|etas_i_j| rho_i.crt() * etas_i_j)
                .collect::<Vec<NTT>>()
        })
        .fold(vec![NTT::zero(); ccs.l], |mut acc, rho_i_times_etas_i| {
            acc.iter_mut()
                .zip(rho_i_times_etas_i)
                .for_each(|(acc_j, rho_i_times_etas_i_j)| {
                    *acc_j += rho_i_times_etas_i_j;
                });

            acc
        });

    let x_0: Vec<NTT> = rho_s
        .iter()
        .zip(cm_i_s.iter())
        .map(|(&rho_i, cm_i)| {
            cm_i.x_w
                .iter()
                .map(|x_w_i| rho_i.crt() * x_w_i)
                .collect::<Vec<NTT>>()
        })
        .fold(vec![NTT::zero(); ccs.n], |mut acc, rho_i_times_x_w_i| {
            acc.iter_mut()
                .zip(rho_i_times_x_w_i)
                .for_each(|(acc_j, rho_i_times_x_w_i)| {
                    *acc_j += rho_i_times_x_w_i;
                });

            acc
        });

    (v_0, cm_0, u_0, x_0)
}

fn prepare_g1_and_3_k_mles_list<NTT: OverField>(
    g: &mut DensePolynomial<NTT>,
    r_i_eq: RefCounter<DenseMultilinearExtension<NTT>>,
    f_hat_mle_s: &[Vec<RefCounter<DenseMultilinearExtension<NTT>>>],
    alpha_s: &[NTT],
    challenged_Ms: &DenseMultilinearExtension<NTT>,
) -> Result<(), ArithErrors> {
    let mut combined_mle: DenseMultilinearExtension<NTT> = DenseMultilinearExtension::zero();

    for (fi_hat_mle_s, alpha_i) in f_hat_mle_s.iter().zip(alpha_s.iter()) {
        let mut mle = DenseMultilinearExtension::zero();
        for fi_hat_mle in fi_hat_mle_s.iter().rev() {
            mle += fi_hat_mle.as_ref();
            mle *= *alpha_i;
        }
        combined_mle += mle;
    }

    combined_mle += challenged_Ms;

    g.add_mles([r_i_eq, RefCounter::from(combined_mle)])?;

    Ok(())
}

fn prepare_g2_i_mle_list<NTT: OverField>(
    g: &mut DensePolynomial<NTT>,
    fi_hat_mle_s: &[RefCounter<DenseMultilinearExtension<NTT>>],
) -> Result<(), ArithErrors> {
    for fi_hat_mle in fi_hat_mle_s.iter() {
        g.add_mles([fi_hat_mle.clone()])?;
    }

    Ok(())
}
