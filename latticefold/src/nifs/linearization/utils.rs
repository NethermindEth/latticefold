use ark_ff::{Field, PrimeField};
use cyclotomic_rings::rings::SuitableRing;
use stark_rings::OverField;
use stark_rings_poly::mle::DenseMultilinearExtension;

use crate::{
    ark_base::Vec,
    nifs::{error::LinearizationError, CCS},
    transcript::Transcript,
    utils::{mle_helpers::evaluate_mles, sumcheck::utils::build_eq_x_r},
};

/// Computes the evaluation of the MLEs of $\\{ M_j \mathbf{z} \mid j = 1, 2, \dots, t \\}$ at the sumcheck challenge point.
///
/// # Parameters
///
/// * `Mz_mles` (`&[DenseMultilinearExtension<NTT>]`): The MLEs of $\\{ M_j \mathbf{z} \mid j = 1, 2, \dots, t \\}$
/// * `r` (`&[NTT]`): The sumcheck challenge point
///
/// # Errors
///
/// This function may return a `LinearizationError` if there are inconsistencies in the dimensions of `Mz_mles` and `r`.
///
pub fn compute_u<NTT: OverField>(
    Mz_mles: &[DenseMultilinearExtension<NTT>],
    r: &[NTT],
) -> Result<Vec<NTT>, LinearizationError<NTT>> {
    evaluate_mles::<NTT, _, _, LinearizationError<NTT>>(Mz_mles, r)
}

/// Prepare the main linearization polynomial.
///
/// $$ g(\vec{\mathbf{x}}) := eq(\vec{\beta}, \vec{\mathbf{x}}) \cdot
/// \left(
/// \sum\_{i=1}^{n\_s} c\_i \cdot
/// \left[
/// \prod\_{j \in S\_i}
/// \left(
/// \sum\_{\vec{\mathbf{b}} \in \\{0,1\\}^{\log n\_c}}
/// \text{mle}[M\_j](\vec{\mathbf{x}}, \vec{\mathbf{b}}) \cdot \text{mle}\[\mathbf{z}\_{ccs}\](\vec{b})
/// \right)
/// \right]
/// \right) $$  
///  
/// # Parameters:
///
/// * `c` (`&[NTT]`): The second multiplicand of the polynomial is a linear combination of products of lists of MLEs, c is the coefficients of the lists
///
/// * `M_mles` (`&[DenseMultilinearExtension<NTT>]`): MLEs that the polynomial is constructed from
///
/// * `S` (`&[Vec<usize>]`): ] indices for the MLE lists
///
/// * `beta_s` (`&[NTT]`): Randomness
///
/// # Returns:
///
/// * The MLEs which form the polynomial
/// * The max degree of the polynomial
///
/// # Errors:
/// * Will return an error if any of the MLEs are of the wrong size
///
pub fn prepare_lin_sumcheck_polynomial<NTT: OverField>(
    c: &[NTT],
    d: &usize,
    M_mles: &[DenseMultilinearExtension<NTT>],
    S: &[Vec<usize>],
    beta_s: &[NTT],
) -> Result<(Vec<DenseMultilinearExtension<NTT>>, usize), LinearizationError<NTT>> {
    let len = 1 + c
        .iter()
        .enumerate()
        .filter(|(_, c)| !c.is_zero())
        .map(|(i, _)| S[i].len())
        .sum::<usize>();

    let mut mles = Vec::with_capacity(len);

    for (i, _) in c.iter().enumerate().filter(|(_, c)| !c.is_zero()) {
        for &j in &S[i] {
            mles.push(M_mles[j].clone());
        }
    }

    mles.push(build_eq_x_r(beta_s)?);

    Ok((mles, d + 1))
}

pub(crate) fn sumcheck_polynomial_comb_fn<NTT: SuitableRing>(vals: &[NTT], ccs: &CCS<NTT>) -> NTT {
    let mut result = NTT::zero();
    'outer: for (i, &c) in ccs.c.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let mut term = c;
        for &j in &ccs.S[i] {
            if vals[j].is_zero() {
                continue 'outer;
            }
            term *= vals[j];
        }
        result += term;
    }
    // eq() is the last term added
    result * vals[vals.len() - 1]
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
