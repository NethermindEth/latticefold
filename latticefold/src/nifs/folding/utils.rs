use std::iter::successors;

use ark_std::iterable::Iterable;
use ark_std::log2;
use ark_std::marker::PhantomData;
use ark_std::sync::Arc;
use cyclotomic_rings::SuitableRing;
use lattirust_poly::mle::SparseMultilinearExtension;

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

pub(super) fn get_alphas_betas_zetas_mus<R: OverField, T: Transcript<R>, P: DecompositionParams>(
    log_m: usize,
    transcript: &mut T,
) -> (Vec<R>, Vec<R>, Vec<R>, Vec<R>) {
    let alpha_s = transcript
        .get_challenges(2 * P::K)
        .into_iter()
        .map(|x| R::from(x))
        .collect::<Vec<_>>();
    let zeta_s = transcript
        .get_challenges(2 * P::K)
        .into_iter()
        .map(|x| R::from(x))
        .collect::<Vec<_>>();
    let mu_s = transcript
        .get_challenges((2 * P::K) - 1)
        .into_iter()
        .map(|x| R::from(x))
        .collect::<Vec<_>>(); // Note is one challenge less
    let beta_s = transcript
        .get_challenges(log_m)
        .into_iter()
        .map(|x| R::from(x))
        .collect::<Vec<_>>();

    (alpha_s, beta_s, zeta_s, mu_s)
}

pub(super) fn create_sumcheck_polynomial<NTT: OverField, DP: DecompositionParams>(
    log_m: usize,
    f_hat_mles: &Vec<DenseMultilinearExtension<NTT>>,
    alpha_is: &Vec<NTT>,
    matrix_mles: &Vec<Vec<DenseMultilinearExtension<NTT>>>,
    zeta_is: &Vec<NTT>,
    ris: Vec<Vec<NTT>>,
    beta_s: &Vec<NTT>,
    mus: &Vec<NTT>,
) -> Result<VirtualPolynomial<NTT>, FoldingError<NTT>> {
    let k = DP::K;
    let mut g = VirtualPolynomial::new(2 * k);
    let mut g1_plus_g3 = VirtualPolynomial::new(2 * k);

    for i in 0..2 * k {
        let gi_1 = prepare_g1_i_mle(&f_hat_mles[i], alpha_is[i]);
        let gi_3 = prepare_g3_i_mle(&matrix_mles[i], zeta_is[i]);
        let gi_1_plus_gi_3 = Arc::from(gi_1 + gi_3);

        let mut g1_and_g3_virtual = VirtualPolynomial::new_from_mle(&gi_1_plus_gi_3, NTT::ONE);
        let eq_r_i = build_eq_x_r(ris[i].as_slice()).unwrap();
        g1_and_g3_virtual.mul_by_mle(eq_r_i, NTT::one());
        g1_plus_g3 = &g1_plus_g3 + &g1_and_g3_virtual;
    }

    let b = DP::B_SMALL; // Get this from the decomposition step
    let mut g2 = prepare_g2_i_mle(log_m, &f_hat_mles[0], b, mus[0]);
    for i in 1..2 * k - 1 {
        let gi_2 = prepare_g2_i_mle(log_m, &f_hat_mles[i], b, mus[i]);
        g2 = &g2 + &gi_2;
    }
    // Recall that the last mu is 1
    let gi_2 = prepare_g2_i_mle(log_m, &f_hat_mles[2 * k - 1], b, NTT::one());
    g2 = &g2 + &gi_2;

    let eq_beta = build_eq_x_r::<NTT>(beta_s.as_slice()).unwrap();
    g2.mul_by_mle(eq_beta, NTT::one());
    g = &g1_plus_g3 + &g2;

    Ok(g)
}

fn prepare_g1_i_mle<NTT: OverField>(
    log_m: usize,
    fi_mle: &DenseMultilinearExtension<NTT>,
    r_i_eq: Arc<DenseMultilinearExtension<NTT>>,
    alpha_i: NTT,
) -> Result<VirtualPolynomial<NTT>, FoldingError<NTT>> {
    let mut mle = fi_mle.clone(); // remove clone

    let mut g1 = VirtualPolynomial::new(log_m);

    g1.add_mle_list(vec![Arc::from(fi_mle.clone())], NTT::ONE)?;

    g1.mul_by_mle(r_i_eq, alpha_i)?;

    Ok(g1)
}

fn prepare_g2_i_mle<NTT: OverField>(
    log_m: usize,
    fi_mle: &DenseMultilinearExtension<NTT>,
    b: u128,
    beta_s: &Vec<NTT>,
    mu_i: NTT,
    beta_eq_r: Arc<DenseMultilinearExtension<NTT>>,
) -> Result<VirtualPolynomial<NTT>, FoldingError<NTT>> {
    let mut mle_list: Vec<Arc<DenseMultilinearExtension<NTT>>> = Vec::new();

    mle_list.push(Arc::from(fi_mle.clone()));

    for i in 0..b {
        let i_hat = NTT::from(i);

        mle_list.push(Arc::from(fi_mle.clone() - i_hat));
        mle_list.push(Arc::from(fi_mle.clone() + i_hat));
    }

    let mut gi_2 = VirtualPolynomial::new(log_m);
    gi_2.add_mle_list(mle_list, NTT::ONE)?;
    gi_2.mul_by_mle(beta_eq_r, mu_i)?;

    Ok(gi_2)
}

fn prepare_g3_i_mle<NTT: OverField>(
    Mz_mles: &Vec<DenseMultilinearExtension<NTT>>,
    zeta_i: NTT,
    r_i_eq: Arc<DenseMultilinearExtension<NTT>>,
) -> Result<VirtualPolynomial<NTT>, FoldingError<NTT>> {
    let mut gi_3 = VirtualPolynomial::new(log_m);

    for (zeta, M) in successors(Some(NTT::ONE), |x| Some(zeta_i * x)).zip(Mz_mles.iter()) {
        gi_3.add_mle_list(vec![Arc::from(M.clone()), r_i_eq], zeta)?;
    }

    Ok(gi_3)
}

pub(super) fn mle_val_from_vector<NTT: OverField>(vector: &Vec<NTT>, values: &Vec<NTT>) -> NTT {
    assert_eq!(values.len(), log2(vector.len()) as usize);
    let mle = DenseMultilinearExtension::from_evaluations_vec(values.len(), vector.clone());
    mle.evaluate(values.as_slice()).unwrap()
}

// Convert a bivariate MLE to a univariate MLE by evaluating the second vector
pub(super) fn mle_matrix_to_val_eval_second<NTT: OverField>(
    matrix: &Vec<Vec<NTT>>,
    values_y: &Vec<NTT>,
) -> Vec<NTT> {
    assert_eq!(values_y.len(), log2(matrix.len()) as usize);
    (0..matrix[0].len())
        .into_iter()
        .map(|i| mle_val_from_vector(&matrix.iter().map(|col| col[i]).collect(), values_y))
        .collect()
}

pub(super) fn usize_to_binary_vector<NTT: OverField>(n: usize, length: usize) -> Vec<NTT> {
    let mut bits = Vec::new();
    let mut current = n;

    // Extract bits from the number
    while current > 0 {
        bits.push((current & 1) as u8);
        current >>= 1;
    }

    // Reverse to get the bits in correct order
    bits.reverse();

    // Pad with leading zeros if necessary
    if bits.len() < length {
        let padding = length - bits.len();
        bits.splice(0..0, std::iter::repeat(0).take(padding));
    }

    // Convert to the target field elements
    bits.into_iter()
        .map(|bit| if bit == 1 { NTT::one() } else { NTT::zero() })
        .collect()
}
