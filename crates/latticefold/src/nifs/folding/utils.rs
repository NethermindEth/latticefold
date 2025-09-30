#![allow(non_snake_case)]

use ark_ff::{Field, PrimeField, Zero};
use ark_std::{
    iter::{self, successors},
    iterable::Iterable,
};
use cyclotomic_rings::{rings::SuitableRing, rotation::rot_lin_combination};
use stark_rings::{cyclotomic_ring::CRT, OverField, PolyRing, Ring};
use stark_rings_poly::mle::DenseMultilinearExtension;

use crate::{
    arith::{CCS, LCCCS},
    ark_base::*,
    commitment::Commitment,
    decomposition_parameters::DecompositionParams,
    nifs::error::FoldingError,
    transcript::{Transcript, TranscriptWithShortChallenges},
    utils::sumcheck::utils::build_eq_x_r,
};

/// A trait for squeezing challenges (`alpha`, `beta`, `zeta`, `mu`) from a cryptographic sponge.
///
///
/// # Type Parameters
/// - `NTT`: A type that implements the `SuitableRing` trait, representing a ring that can be used in the
///   LatticeFold protocol.
///
pub(crate) trait SqueezeAlphaBetaZetaMu<NTT: SuitableRing> {
    /// Extracts the cryptographic challenge vectors of provided length
    ///
    /// ### Arguments
    /// - `log_m`: The length of the $\beta$ challenge vector.
    ///
    /// ### Type Parameters
    /// - `P`: The decomposition parameters of the protocol.
    ///
    /// ### Returns
    /// - `(Vec<NTT>, Vec<NTT>, Vec<NTT>, Vec<NTT>)`: A tuple containing four challenge vectors:
    ///   - `alpha`: A challenge vector of length $2 \cdot k$, where $k$ is defined in the decomposition parameters.
    ///   - `beta`: A challenge vector of length `log_m`.
    ///   - `zeta`: A challenge vector of length $2 \cdot k$, where $k$ is defined in the decomposition parameters.
    ///   - `mu`: A challenge vector of length $2 \cdot k$, where $k$ is defined in the decomposition parameters.
    ///
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

/// Generates `rho` values based on the provided transcript and decomposition parameters.
///
/// This function is used within the module to extract or compute values required for further
/// operations, based on the interaction with a transcript that supports short challenges.
///
/// # Type Parameters
/// - `R`: A ring suitable to be used in the LatticeFold protocol.
/// - `T`: A type implementing a cryptographic sponge construction.
/// - `P`: The decomposition parameters of the protocol.
///
/// # Arguments
/// - `transcript`: A mutable reference to the transcript `T` from which we squeeze the challenges.
///
/// # Returns
/// - `(Vec<R::CoefficientRepresentation>, Vec<R>)`:
///   - The first element is a vector of challenges in coefficient form.
///   - The second element is the same vector of challenges in NTT form.
///
pub(super) fn get_rhos<
    R: SuitableRing,
    T: TranscriptWithShortChallenges<R>,
    P: DecompositionParams,
>(
    transcript: &mut T,
) -> (Vec<R::CoefficientRepresentation>, Vec<R>) {
    transcript.absorb_field_element(&<R::BaseRing as Field>::from_base_prime_field(
        <R::BaseRing as Field>::BasePrimeField::from_be_bytes_mod_order(b"rho_s"),
    ));

    let mut rhos_coeff = transcript.get_small_challenges((2 * P::K) - 1); // Note that we are missing the first element
    rhos_coeff.push(R::CoefficientRepresentation::ONE);
    let rhos = CRT::elementwise_crt(rhos_coeff.clone());
    (rhos_coeff, rhos)
}

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
pub(super) fn create_sumcheck_polynomial<NTT: OverField, DP: DecompositionParams>(
    log_m: usize,
    f_hat_mles: Vec<Vec<DenseMultilinearExtension<NTT>>>,
    alpha_s: &[NTT],
    challenged_Ms_1: &DenseMultilinearExtension<NTT>,
    challenged_Ms_2: &DenseMultilinearExtension<NTT>,
    r_s: &[Vec<NTT>],
    beta_s: &[NTT],
    mu_s: &[NTT],
) -> Result<(Vec<DenseMultilinearExtension<NTT>>, usize), FoldingError<NTT>> {
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

    let len = 2 + 2 + // g1 + g3
        1 + f_hat_mles.len() * f_hat_mles[0].len(); // g2
    let mut mles = Vec::with_capacity(len);

    // We assume here that decomposition subprotocol puts the same r challenge point
    // into all decomposed linearized commitments
    let r_i_eq = build_eq_x_r(&r_s[0])?;
    prepare_g1_and_3_k_mles_list(
        &mut mles,
        r_i_eq.clone(),
        &f_hat_mles[0..DP::K],
        &alpha_s[0..DP::K],
        challenged_Ms_1,
    );

    let r_i_eq = build_eq_x_r(&r_s[DP::K])?;
    prepare_g1_and_3_k_mles_list(
        &mut mles,
        r_i_eq,
        &f_hat_mles[DP::K..2 * DP::K],
        &alpha_s[DP::K..2 * DP::K],
        challenged_Ms_2,
    );

    // g2
    let beta_eq_x = build_eq_x_r(beta_s)?;
    prepare_g2_i_mle_list(&mut mles, beta_eq_x, f_hat_mles);

    let degree = 2 * DP::B_SMALL;

    Ok((mles, degree))
}

/// Combines evaluations of MLE into evaluation of folding sumcheck polynomial
///
/// # Arguments
///
/// - `vals: &[NTT]`:
///   The evaluations of the multilinear extensions produced by the `create_sumcheck_polynomial` function
/// - `mu_s: &[NTT]`
///   The $\mu$ challenges
///
///  # Returns
///  - NTT:
///    The value of the same evaluation point evaluated by the folding sumcheck polynomial
pub(crate) fn sumcheck_polynomial_comb_fn<NTT: SuitableRing, P: DecompositionParams>(
    vals: &[NTT],
    mu_s: &[NTT],
) -> NTT {
    let extension_degree = NTT::CoefficientRepresentation::dimension() / <NTT>::dimension();

    // Add eq_r * g1 * g3 for first k
    let mut result = vals[0] * vals[1];

    // Add eq_r * g1 * g3 for second k
    result += vals[2] * vals[3];

    // We have k * extension degree mles of b
    // each one consists of (2 * small_b) -1 extensions
    // We start at index 5
    // Multiply each group of (2 * small_b) -1 extensions
    // Then multiply by the eq_beta evaluation at index 4
    for (k, mu) in mu_s.iter().enumerate() {
        let mut inter_result = NTT::zero();
        for d in (0..extension_degree).rev() {
            let i = k * extension_degree + d;

            let f_i = vals[5 + i];

            if f_i.is_zero() {
                if !inter_result.is_zero() {
                    inter_result *= mu;
                }
                continue;
            }

            // start with eq_b
            let mut eval = vals[4];

            let f_i_squared = f_i * f_i;

            for b in 1..<P>::B_SMALL {
                let multiplicand = f_i_squared - NTT::from(b as u128 * b as u128);
                if multiplicand.is_zero() {
                    eval = NTT::zero();
                    break;
                }
                eval *= multiplicand
            }
            eval *= f_i;
            inter_result += eval;
            inter_result *= mu
        }
        result += inter_result;
    }

    result
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

/// Computes `v0`, `u0`, `x0`, and `cm_0` as folding subprotocol.
///
/// # Type Parameters
///
/// - `NTT`: A ring suitable to be used in the LatticeFold protocol.
///
/// # Arguments
///
/// - `rho_s: &[NTT::CoefficientRepresentation]`  
///
///     $\rho$ challenges
///
/// - `theta_s: &[Vec<NTT>]`
///   $$
///   \left[\theta\_{i} := \text{mle}\[\hat{f}\_i\](\vec{r}_o) \right]\_{i=1}^{2k},
///   $$
/// - `cm_i_s: &[LCCCS<NTT>]`
///
///     Decomposed linearized commitments
///
/// - `eta_s: &[Vec<NTT>]`  
///
///     $$
///   \eta[i] :=
///   \sum\_{
///   \vec{b} \in \\{0,1\\}^\{log\(n + n\_{in}\)\}
///   }
///   \text{mle}\[M_1\]\(\vec{r}\_o, \vec{b}\) \cdot \text{mle}\[z_i\]\(\vec{b}\)
///   $$
///
/// - `ccs: &CCS<NTT>`  
///
///     A reference to a Customizable Constraint System instance used in the protocol.
///
/// # Returns
///
/// - `(Vec<NTT>, Commitment<NTT>, Vec<NTT>, Vec<NTT>)`  
///   A tuple containing:
///   - `v0: Vec<NTT>`  
///     Evaluation of linearized folded witness at $\vec{r}\_o$
///   - `u_0: Commitment<NTT>`
///     A linear combination of $\left[ eta_s[i] \right]\_{i=1}^{2k}$
///   - `x0: Vec<NTT>`
///     Folded CCS statement
///   - `cm_0: Vec<NTT>`
///     Folded commitment
pub(super) fn compute_v0_u0_x0_cm_0<NTT: SuitableRing>(
    rho_s_coeff: &[NTT::CoefficientRepresentation],
    rho_s: &[NTT],
    theta_s: &[Vec<NTT>],
    cm_i_s: &[LCCCS<NTT>],
    eta_s: &[Vec<NTT>],
    ccs: &CCS<NTT>,
) -> (Vec<NTT>, Commitment<NTT>, Vec<NTT>, Vec<NTT>) {
    let v_0: Vec<NTT> = rot_lin_combination(rho_s_coeff, theta_s);

    let empty_commitment = Commitment::zeroed(cm_i_s.first().map(|cm_i| cm_i.cm.len()).unwrap_or(0));
    let cm_0: Commitment<NTT> = rho_s
        .iter()
        .zip(cm_i_s.iter())
        .map(|(&rho_i, cm_i)| cm_i.cm.clone() * rho_i)
        .fold(empty_commitment, |acc, c| acc + c);

    let u_0: Vec<NTT> = rho_s
        .iter()
        .zip(eta_s.iter())
        .map(|(&rho_i, etas_i)| {
            etas_i
                .iter()
                .map(|etas_i_j| rho_i * etas_i_j)
                .collect::<Vec<NTT>>()
        })
        .fold(vec![NTT::zero(); ccs.t], |mut acc, rho_i_times_etas_i| {
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
                .chain(iter::once(&cm_i.h))
                .map(|x_w_i| rho_i * x_w_i)
                .collect::<Vec<NTT>>()
        })
        .fold(
            vec![NTT::zero(); ccs.l + 1],
            |mut acc, rho_i_times_x_w_i| {
                acc.iter_mut()
                    .zip(rho_i_times_x_w_i)
                    .for_each(|(acc_j, rho_i_times_x_w_i)| {
                        *acc_j += rho_i_times_x_w_i;
                    });

                acc
            },
        );

    (v_0, cm_0, u_0, x_0)
}

/// Get the MLEs needed for $k$ g1 and g3 components of the sumcheck polynomial
fn prepare_g1_and_3_k_mles_list<NTT: OverField>(
    mles: &mut Vec<DenseMultilinearExtension<NTT>>,
    r_i_eq: DenseMultilinearExtension<NTT>,
    f_hat_mle_s: &[Vec<DenseMultilinearExtension<NTT>>],
    alpha_s: &[NTT],
    challenged_Ms: &DenseMultilinearExtension<NTT>,
) {
    let mut combined_mle: DenseMultilinearExtension<NTT> = DenseMultilinearExtension::zero();

    for (fi_hat_mle_s, alpha_i) in f_hat_mle_s.iter().zip(alpha_s.iter()) {
        let mut mle = DenseMultilinearExtension::zero();
        for fi_hat_mle in fi_hat_mle_s.iter().rev() {
            mle += fi_hat_mle;
            mle *= *alpha_i;
        }
        combined_mle += mle;
    }

    combined_mle += challenged_Ms;

    mles.push(r_i_eq);
    mles.push(combined_mle);
}
/// Get the MLEs needed for one g2 component of the sumcheck polynomial
fn prepare_g2_i_mle_list<NTT: OverField>(
    mles: &mut Vec<DenseMultilinearExtension<NTT>>,
    beta_eq_x: DenseMultilinearExtension<NTT>,
    f_hat_mles: Vec<Vec<DenseMultilinearExtension<NTT>>>,
) {
    mles.push(beta_eq_x);
    f_hat_mles
        .into_iter()
        .for_each(|mut fhms| mles.append(&mut fhms))
}
