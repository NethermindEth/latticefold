#![allow(non_snake_case)]
use ark_ff::{Field, PrimeField};
use ark_std::iter::successors;
use ark_std::iterable::Iterable;
use ark_std::sync::Arc;
use cyclotomic_rings::SuitableRing;
use lattirust_poly::polynomials::ArithErrors;
use lattirust_ring::Ring;

use crate::commitment::Commitment;
use crate::nifs::error::FoldingError;
use crate::transcript::TranscriptWithSmallChallenges;
use crate::{
    arith::{CCS, LCCCS},
    decomposition_parameters::DecompositionParams,
    transcript::Transcript,
};
use cyclotomic_rings::rot_sum;
use lattirust_poly::{
    mle::DenseMultilinearExtension,
    polynomials::{build_eq_x_r, VirtualPolynomial},
};
use lattirust_ring::OverField;

pub(super) fn get_alphas_betas_zetas_mus<
    NTT: OverField,
    T: Transcript<NTT>,
    P: DecompositionParams,
>(
    log_m: usize,
    transcript: &mut T,
) -> (Vec<NTT>, Vec<NTT>, Vec<NTT>, Vec<NTT>) {
    transcript.absorb_field_element(
        &<NTT::BaseRing as Field>::from_base_prime_field_elems(&[
            <NTT::BaseRing as Field>::BasePrimeField::from_be_bytes_mod_order(b"alpha_s"),
        ])
        .unwrap(),
    );
    let alpha_s = transcript
        .get_challenges(2 * P::K)
        .into_iter()
        .map(|x| NTT::from(x))
        .collect::<Vec<_>>();

    transcript.absorb_field_element(
        &<NTT::BaseRing as Field>::from_base_prime_field_elems(&[
            <NTT::BaseRing as Field>::BasePrimeField::from_be_bytes_mod_order(b"zeta_s"),
        ])
        .unwrap(),
    );
    let zeta_s = transcript
        .get_challenges(2 * P::K)
        .into_iter()
        .map(|x| NTT::from(x))
        .collect::<Vec<_>>();

    transcript.absorb_field_element(
        &<NTT::BaseRing as Field>::from_base_prime_field_elems(&[
            <NTT::BaseRing as Field>::BasePrimeField::from_be_bytes_mod_order(b"mu_s"),
        ])
        .unwrap(),
    );
    let mut mu_s = transcript
        .get_challenges((2 * P::K) - 1)
        .into_iter()
        .map(|x| NTT::from(x))
        .collect::<Vec<_>>(); // Note is one challenge less

    mu_s.push(NTT::ONE);

    transcript.absorb_field_element(
        &<NTT::BaseRing as Field>::from_base_prime_field_elems(&[
            <NTT::BaseRing as Field>::BasePrimeField::from_be_bytes_mod_order(b"beta_s"),
        ])
        .unwrap(),
    );
    let beta_s = transcript
        .get_challenges(log_m)
        .into_iter()
        .map(|x| NTT::from(x))
        .collect::<Vec<_>>();

    (alpha_s, beta_s, zeta_s, mu_s)
}

pub(super) fn get_rhos<
    R: SuitableRing,
    T: TranscriptWithSmallChallenges<R>,
    P: DecompositionParams,
>(
    transcript: &mut T,
) -> Vec<R::CoefficientRepresentation> {
    transcript.absorb_field_element(
        &<R::BaseRing as Field>::from_base_prime_field_elems(&[
            <R::BaseRing as Field>::BasePrimeField::from_be_bytes_mod_order(b"rho_s"),
        ])
        .unwrap(),
    );

    let mut rhos = transcript.get_small_challenges((2 * P::K) - 1); // Note that we are missing the first element
    rhos.push(R::CoefficientRepresentation::ONE);

    rhos
}

#[allow(clippy::too_many_arguments)]
pub(super) fn create_sumcheck_polynomial<NTT: OverField, DP: DecompositionParams>(
    log_m: usize,
    f_hat_mles: &[Vec<DenseMultilinearExtension<NTT>>],
    alpha_s: &[NTT],
    Mz_mles: &[Vec<DenseMultilinearExtension<NTT>>],
    zeta_s: &[NTT],
    r_s: &[Vec<NTT>],
    beta_s: &[NTT],
    mu_s: &[NTT],
) -> Result<VirtualPolynomial<NTT>, FoldingError<NTT>> {
    if alpha_s.len() != 2 * DP::K
        || f_hat_mles.len() != 2 * DP::K
        || Mz_mles.len() != 2 * DP::K
        || zeta_s.len() != 2 * DP::K
        || r_s.len() != 2 * DP::K
        || beta_s.len() != log_m
        || mu_s.len() != 2 * DP::K
    {
        return Err(FoldingError::IncorrectLength);
    }

    let mut g = VirtualPolynomial::<NTT>::new(log_m);

    let beta_eq_x = build_eq_x_r(beta_s)?;

    for i in 0..2 * DP::K {
        let r_i_eq = build_eq_x_r(&r_s[i])?;

        prepare_g1_i_mle_list(&mut g, f_hat_mles[i].clone(), r_i_eq.clone(), alpha_s[i])?;
        prepare_g2_i_mle_list(
            &mut g,
            f_hat_mles[i].clone(),
            DP::B_SMALL,
            mu_s[i],
            beta_eq_x.clone(),
        )?;
        prepare_g3_i_mle_list(&mut g, &Mz_mles[i], zeta_s[i], r_i_eq)?;
    }

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
    let transposed_theta_s = (0..theta_s[0].len())
        .map(|i| theta_s.iter().map(|row| row[i]).collect::<Vec<NTT>>())
        .collect::<Vec<Vec<NTT>>>();

    let v_0: Vec<NTT> = rho_s
        .iter()
        .zip(transposed_theta_s.iter())
        .map(|(&rho_i, theta_j)| {
            theta_j
                .iter()
                .map(|theta_j_i| NTT::from(rot_sum::<NTT>(rho_i, theta_j_i.coeffs()))) // TODO: Double check
                .sum()
        })
        .collect();

    let cm_0: Commitment<C, NTT> = rho_s
        .iter()
        .zip(cm_i_s.iter())
        .map(|(&rho_i, cm_i)| cm_i.cm.clone() * NTT::from(rho_i))
        .sum();

    let u_0: Vec<NTT> = rho_s
        .iter()
        .zip(eta_s.iter())
        .map(|(&rho_i, etas_i)| {
            etas_i
                .iter()
                .map(|etas_i_j| NTT::from(rho_i) * etas_i_j)
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
                .map(|x_w_i| NTT::from(rho_i) * x_w_i)
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

fn prepare_g1_i_mle_list<NTT: OverField>(
    g: &mut VirtualPolynomial<NTT>,
    fi_hat_mle_s: Vec<DenseMultilinearExtension<NTT>>,
    r_i_eq: Arc<DenseMultilinearExtension<NTT>>,
    alpha_i: NTT,
) -> Result<(), ArithErrors> {
    for (alpha, fi_hat_mle) in successors(Some(alpha_i), |alpha_power| Some(alpha_i * alpha_power))
        .zip(fi_hat_mle_s.iter())
    {
        g.add_mle_list(vec![r_i_eq.clone(), Arc::from(fi_hat_mle.clone())], alpha)?;
    }

    Ok(())
}

fn prepare_g2_i_mle_list<NTT: OverField>(
    g: &mut VirtualPolynomial<NTT>,
    fi_hat_mle_s: Vec<DenseMultilinearExtension<NTT>>,
    b: usize,
    mu_i: NTT,
    beta_eq_x: Arc<DenseMultilinearExtension<NTT>>,
) -> Result<(), ArithErrors> {
    for (mu, fi_hat_mle) in
        successors(Some(mu_i), |mu_power| Some(mu_i * mu_power)).zip(fi_hat_mle_s.iter())
    {
        let mut mle_list: Vec<Arc<DenseMultilinearExtension<NTT>>> = Vec::new();

        for i in 1..b {
            let i_hat = NTT::from(i as u128);

            mle_list.push(Arc::from(fi_hat_mle.clone() - i_hat));
            mle_list.push(Arc::from(fi_hat_mle.clone() + i_hat));
        }

        mle_list.push(Arc::from(fi_hat_mle.clone()));

        mle_list.push(beta_eq_x.clone());

        let _ = g.add_mle_list(mle_list, mu);
    }
    Ok(())
}

fn prepare_g3_i_mle_list<NTT: OverField>(
    g: &mut VirtualPolynomial<NTT>,
    Mz_mles: &[DenseMultilinearExtension<NTT>],
    zeta_i: NTT,
    r_i_eq: Arc<DenseMultilinearExtension<NTT>>,
) -> Result<(), ArithErrors> {
    for (zeta, M) in successors(Some(zeta_i), |x| Some(zeta_i * x)).zip(Mz_mles.iter()) {
        g.add_mle_list(vec![Arc::from(M.clone()), r_i_eq.clone()], zeta)?;
    }
    Ok(())
}
