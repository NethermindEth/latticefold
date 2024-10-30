use ark_std::sync::Arc;

use lattirust_poly::{
    mle::DenseMultilinearExtension,
    polynomials::{build_eq_x_r, VirtualPolynomial},
};
use lattirust_ring::OverField;

use crate::nifs::error::LinearizationError;

/// Batch compute the values of mles at the point r.
pub fn compute_u<NTT: OverField>(
    Mz_mles: &[DenseMultilinearExtension<NTT>],
    r: &[NTT],
) -> Result<Vec<NTT>, LinearizationError<NTT>> {
    Mz_mles
        .iter()
        .map(|M_i_mle| {
            M_i_mle
                .evaluate(r)
                .ok_or(LinearizationError::ParametersError(format!(
                    "one of the CCS matrices has an incorrect length {}, expected {}",
                    M_i_mle.evaluations.len(),
                    1 << r.len(),
                )))
        })
        .collect()
}

/// Prepare the main linearization polynomial.
pub fn prepare_lin_sumcheck_polynomial<NTT: OverField>(
    log_m: usize,
    c: &[NTT],
    M_mles: &[DenseMultilinearExtension<NTT>],
    S: &[Vec<usize>],
    beta_s: &[NTT],
) -> Result<VirtualPolynomial<NTT>, LinearizationError<NTT>> {
    let mut g = VirtualPolynomial::new(log_m);

    for (i, coefficient) in c.iter().enumerate().filter(|(_, c)| !c.is_zero()) {
        let mut mle_list: Vec<Arc<DenseMultilinearExtension<NTT>>> = Vec::with_capacity(S[i].len());

        for &j in &S[i] {
            mle_list.push(Arc::new(M_mles[j].clone()));
        }

        g.add_mle_list(mle_list, *coefficient)?;
    }

    g.mul_by_mle(build_eq_x_r(beta_s)?, NTT::one())?;

    Ok(g)
}
