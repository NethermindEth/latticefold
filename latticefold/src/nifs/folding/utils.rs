#![allow(non_snake_case)]
use std::iter::successors;

use ark_ff::{Field, PrimeField};
use ark_std::iterable::Iterable;
use ark_std::log2;
use ark_std::marker::PhantomData;
use ark_std::sync::Arc;
use cyclotomic_rings::SuitableRing;
use lattirust_poly::mle::SparseMultilinearExtension;
use lattirust_poly::polynomials::ArithErrors;
use lattirust_ring::Ring;

use crate::nifs::error::FoldingError;
use crate::transcript::TranscriptWithSmallChallenges;
use crate::utils::sumcheck::{MLSumcheck, SumCheckError::SumCheckFailed};
use crate::{
    arith::{utils::mat_vec_mul, Instance, Witness, CCS, LCCCS},
    parameters::DecompositionParams,
    transcript::Transcript,
    utils::{mle::dense_vec_to_dense_mle, sumcheck},
};

use lattirust_poly::{
    mle::DenseMultilinearExtension,
    polynomials::{build_eq_x_r, eq_eval, VPAuxInfo, VirtualPolynomial},
};
use lattirust_ring::{OverField, PolyRing};

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

pub(super) fn create_sumcheck_polynomial<NTT: OverField, DP: DecompositionParams>(
    log_m: usize,
    f_hat_mles: &[DenseMultilinearExtension<NTT>],
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
        || beta_s.len() != 2 * DP::K
        || mu_s.len() != 2 * DP::K
    {
        return Err(FoldingError::IncorrectLength);
    }

    let mut g = VirtualPolynomial::<NTT>::new(log_m);

    let beta_eq_x = build_eq_x_r(beta_s)?;

    for i in 0..2 * DP::K {
        let r_i_eq = build_eq_x_r(&r_s[i])?;

        prepare_g1_i_mle_list(
            &mut g,
            log_m,
            f_hat_mles[i].clone(),
            r_i_eq.clone(),
            alpha_s[i],
        )?;
        prepare_g2_i_mle_list(
            &mut g,
            log_m,
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
    theta_s: &[NTT],
    e_asterisk: NTT,
    e_s: &[NTT],
    zeta_s: &[NTT],
    eta_s: &[Vec<NTT>],
) -> NTT {
    (0..(2 * P::K))
        .map(|i| {
            // Evaluation claims about f hats.
            let mut s_summand = alpha_s[i] * e_s[i] * theta_s[i];

            // norm range check contribution
            s_summand += e_asterisk
                * mu_s[i]
                * theta_s[i]
                * (1..P::B_SMALL)
                    .map(NTT::from)
                    .map(|j_hat| (theta_s[i] - j_hat) * (theta_s[i] + j_hat))
                    .product::<NTT>();

            // linearisation claims contribuition
            s_summand += zeta_s
                .iter()
                .zip(eta_s.iter())
                .map(|(&zeta_i, eta_i_s)| {
                    successors(Some(zeta_i), |&zeta| Some(zeta * zeta_i))
                        .zip(eta_i_s.iter())
                        .map(|(pow_of_zeta, eta_i_j)| pow_of_zeta * eta_i_j)
                        .sum::<NTT>()
                })
                .sum::<NTT>();

            s_summand
        })
        .sum()
}

fn prepare_g1_i_mle_list<NTT: OverField>(
    g: &mut VirtualPolynomial<NTT>,
    log_m: usize,
    fi_mle: DenseMultilinearExtension<NTT>,
    r_i_eq: Arc<DenseMultilinearExtension<NTT>>,
    alpha_i: NTT,
) -> Result<(), ArithErrors> {
    g.add_mle_list(vec![r_i_eq, Arc::from(fi_mle)], alpha_i)
}

fn prepare_g2_i_mle_list<NTT: OverField>(
    g: &mut VirtualPolynomial<NTT>,
    log_m: usize,
    fi_mle: DenseMultilinearExtension<NTT>,
    b: u128,
    mu_i: NTT,
    beta_eq_x: Arc<DenseMultilinearExtension<NTT>>,
) -> Result<(), ArithErrors> {
    let mut mle_list: Vec<Arc<DenseMultilinearExtension<NTT>>> = Vec::new();

    for i in 1..b {
        let i_hat = NTT::from(i);

        mle_list.push(Arc::from(fi_mle.clone() - i_hat));
        mle_list.push(Arc::from(fi_mle.clone() + i_hat));
    }

    mle_list.push(Arc::from(fi_mle));

    g.add_mle_list(mle_list, mu_i)
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
