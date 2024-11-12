use crate::{ark_base::Vec, nifs::mle_helpers::evaluate_mles};
use ark_ff::PrimeField;

use lattirust_poly::{
    mle::DenseMultilinearExtension,
    polynomials::{build_eq_x_r, RefCounter, VirtualPolynomial},
};
use lattirust_ring::OverField;

use crate::nifs::error::LinearizationError;
use crate::transcript::Transcript;
use ark_ff::Field;
use cyclotomic_rings::rings::SuitableRing;

/// Batch compute the values of mles at the point r.
pub fn compute_u<NTT: OverField>(
    Mz_mles: &[DenseMultilinearExtension<NTT>],
    r: &[NTT],
) -> Result<Vec<NTT>, LinearizationError<NTT>> {
    evaluate_mles::<NTT, _, _, LinearizationError<NTT>>(Mz_mles, r)
}

/// Prepare the main linearization polynomial.
///
/// # Parameters:
///
/// * `log_m` (`usize`): The number of variables in the polynomial
///
/// * `c` (`&[NTT]`): Polynomial is a linear comibination lists of MLEs, c is the coefficients of the lists
///
/// * `M_mles` (`&[DenseMultilinearExtension<NTT>]`): MLEs that the polynomial is constructed from
///
/// * `S` (`&[Vec<usize>]`): ] indices for the MLE lists
///
/// * `beta_s` (`&[NTT]`): Randomness
///
/// # Returns:
///
/// * `VirtualPolynomial<NTT>`: The linearization sumcheck polynomial
///
/// # Errors:
/// * Will return an error if na of the MLEs are of the wrong size
///
pub fn prepare_lin_sumcheck_polynomial<NTT: OverField>(
    log_m: usize,
    c: &[NTT],
    M_mles: &[DenseMultilinearExtension<NTT>],
    S: &[Vec<usize>],
    beta_s: &[NTT],
) -> Result<VirtualPolynomial<NTT>, LinearizationError<NTT>> {
    let mut g = VirtualPolynomial::new(log_m);

    for (i, coefficient) in c.iter().enumerate().filter(|(_, c)| !c.is_zero()) {
        let mut mle_list: Vec<RefCounter<DenseMultilinearExtension<NTT>>> =
            Vec::with_capacity(S[i].len());

        for &j in &S[i] {
            mle_list.push(RefCounter::new(M_mles[j].clone()));
        }

        g.add_mle_list(mle_list, *coefficient)?;
    }

    g.mul_by_mle(build_eq_x_r(beta_s)?, NTT::one())?;

    Ok(g)
}

pub(crate) trait SqueezeBeta<NTT: SuitableRing> {
    fn squeeze_beta_challenges(&mut self, n: usize) -> Vec<NTT>;
}

impl<R: SuitableRing, T: Transcript<R>> SqueezeBeta<R> for T {
    fn squeeze_beta_challenges(&mut self, n: usize) -> Vec<R> {
        self.absorb_field_element(&<R::BaseRing as Field>::from_base_prime_field(
            <R::BaseRing as Field>::BasePrimeField::from_be_bytes_mod_order(b"beta_s"),
        ));

        self.get_challenges(n)
            .into_iter()
            .map(|x| x.into())
            .collect()
    }
}
