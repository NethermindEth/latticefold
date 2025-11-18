#![allow(non_snake_case)]

use ark_ff::{Field, PrimeField, Zero};
use ark_std::{iter, iterable::Iterable};
use cyclotomic_rings::{rings::SuitableRing, rotation::rot_lin_combination};
use stark_rings::{cyclotomic_ring::CRT, OverField, PolyRing, Ring};
use stark_rings_poly::mle::DenseMultilinearExtension;

use crate::{
    arith::{CCS, LCCCS},
    ark_base::*,
    commitment::Commitment,
    transcript::{Transcript, TranscriptWithShortChallenges},
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
    fn squeeze_alpha_beta_zeta_mu(
        &mut self,
        log_m: usize,
        k: usize,
    ) -> (Vec<NTT>, Vec<NTT>, Vec<NTT>, Vec<NTT>);
}

impl<NTT: SuitableRing, T: Transcript<NTT>> SqueezeAlphaBetaZetaMu<NTT> for T {
    fn squeeze_alpha_beta_zeta_mu(
        &mut self,
        log_m: usize,
        k: usize,
    ) -> (Vec<NTT>, Vec<NTT>, Vec<NTT>, Vec<NTT>) {
        self.absorb_field_element(&<NTT::BaseRing as Field>::from_base_prime_field(
            <NTT::BaseRing as Field>::BasePrimeField::from_be_bytes_mod_order(b"alpha_s"),
        ));
        let alpha_s = self
            .get_challenges(2 * k)
            .into_iter()
            .map(|x| NTT::from(x))
            .collect::<Vec<_>>();

        self.absorb_field_element(&<NTT::BaseRing as Field>::from_base_prime_field(
            <NTT::BaseRing as Field>::BasePrimeField::from_be_bytes_mod_order(b"zeta_s"),
        ));
        let zeta_s = self
            .get_challenges(2 * k)
            .into_iter()
            .map(|x| NTT::from(x))
            .collect::<Vec<_>>();

        self.absorb_field_element(&<NTT::BaseRing as Field>::from_base_prime_field(
            <NTT::BaseRing as Field>::BasePrimeField::from_be_bytes_mod_order(b"mu_s"),
        ));
        let mut mu_s = self
            .get_challenges((2 * k) - 1)
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
pub(super) fn get_rhos<R: SuitableRing, T: TranscriptWithShortChallenges<R>>(
    k: usize,
    transcript: &mut T,
) -> (Vec<R::CoefficientRepresentation>, Vec<R>) {
    transcript.absorb_field_element(&<R::BaseRing as Field>::from_base_prime_field(
        <R::BaseRing as Field>::BasePrimeField::from_be_bytes_mod_order(b"rho_s"),
    ));

    let mut rhos_coeff = transcript.get_small_challenges((2 * k) - 1); // Note that we are missing the first element
    rhos_coeff.push(R::CoefficientRepresentation::ONE);
    let rhos = CRT::elementwise_crt(rhos_coeff.clone());
    (rhos_coeff, rhos)
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
pub(crate) fn sumcheck_polynomial_comb_fn<NTT: SuitableRing>(
    vals: &[NTT],
    mu_s: &[NTT],
    b: usize,
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

            for b in 1..b {
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

    let empty_commitment =
        Commitment::zeroed(cm_i_s.first().map(|cm_i| cm_i.cm.len()).unwrap_or(0));
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
pub fn prepare_g1_and_3_k_mles_list<NTT: OverField>(
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
pub fn prepare_g2_i_mle_list<NTT: OverField>(
    mles: &mut Vec<DenseMultilinearExtension<NTT>>,
    beta_eq_x: DenseMultilinearExtension<NTT>,
    f_hat_mles: Vec<Vec<DenseMultilinearExtension<NTT>>>,
) {
    mles.push(beta_eq_x);
    f_hat_mles
        .into_iter()
        .for_each(|mut fhms| mles.append(&mut fhms))
}
