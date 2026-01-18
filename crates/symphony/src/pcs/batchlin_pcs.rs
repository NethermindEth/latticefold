//! Batchlin PCS integration for Symphony.
//!
//! This module provides the bridge between Symphony's batchlin relation (evaluation at `r'`)
//! and the Cini ℓ=2 folding PCS (evaluation at `(x0, x1, x2)`).
//!
//! ## Architecture
//!
//! - Π_mon prover messages (`cm_g`) are transcript-bound before Π_mon coins are sampled (via a
//!   compressed aggregate commitment `cm_g_agg` that is itself bound in-gate).
//! - The folded batchlin object is bound/verified via the batched scalar PCS (PCS#2) in the WE gate.
//!
//! ## Point Conversion
//!
//! Symphony's batchlin evaluation point is `r' ∈ F^{log2(n)}` (MLE format).
//! The Cini PCS expects `(x0, x1, x2)` where each `xk ∈ F^r` is a tensor weight table.
//!
//! For n = 2^24 with r = 256 = 2^8:
//! - Split `r'` into three blocks of 8 coordinates: `r' = (r0 || r1 || r2)`
//! - Convert each block to tensor weights: `xk[a] = eq(rk, bits(a))` for `a ∈ [0..r)`
//!
//! Then `MLE(g)(r') = ⟨x0 ⊗ x1 ⊗ x2, g⟩` which is what the PCS verifier computes.

use ark_ff::PrimeField;
use ark_std::vec::Vec;
use rayon::prelude::*;

use super::folding_pcs_l2::FoldingPcsL2ProofCore;
use super::folding_pcs_l2::{DenseMatrix, FoldingPcsL2Params};
use rand::{RngCore, SeedableRng};
use rand_chacha::ChaCha20Rng;
use sha2::{Digest, Sha256};
use crate::symphony_coins::ts_weights;

/// Configuration for batchlin PCS.
#[derive(Clone, Debug)]
pub struct BatchlinPcsConfig {
    /// Log2 of the vector length n = m*d.
    pub log_n: usize,
    /// Branching factor r for ℓ=2 folding.
    ///
    /// We allow padding: we choose `r` as a power-of-two such that `r^3 >= 2^{log_n}`.
    /// This removes the hard `log_n % 3 == 0` restriction while keeping the ℓ=2 PCS core unchanged.
    pub r: usize,
    /// Security parameter κ.
    pub kappa: usize,
    /// Number of digits k_g.
    pub k_g: usize,
}

impl BatchlinPcsConfig {
    /// Create config for standard SP1 parameters.
    pub fn for_sp1(log_n: usize, kappa: usize, k_g: usize) -> Result<Self, String> {
        // For ℓ=2, the PCS core expects a cubic domain of size r^3, with x0,x1,x2 each length r.
        // We choose the smallest power-of-two r such that r^3 >= 2^{log_n} and treat the message
        // as zero-padded on the larger domain.
        let log_r = (log_n + 2) / 3; // ceil(log_n/3)
        let r = 1usize << log_r;
        Ok(Self { log_n, r, kappa, k_g })
    }

    /// Vector length n = 2^{log_n}.
    pub fn n(&self) -> usize {
        1 << self.log_n
    }

    /// Number of coordinates in each block (log_r = log_n / 3).
    pub fn log_r(&self) -> usize {
        (self.log_n + 2) / 3
    }
}

/// Convert MLE evaluation point `r' ∈ F^{log_n}` to tensor format `(x0, x1, x2)`.
///
/// Each output `xk` has length `r` and represents the eq-weights for the k-th block.
pub fn mle_point_to_tensor<F: PrimeField>(
    r_prime: &[F],
    config: &BatchlinPcsConfig,
) -> Result<(Vec<F>, Vec<F>, Vec<F>), String> {
    if r_prime.len() != config.log_n {
        return Err(format!(
            "r_prime length {} != log_n {}",
            r_prime.len(),
            config.log_n
        ));
    }

    let log_r = config.log_r();
    let r = config.r;

    // If log_n is not divisible by 3, we conceptually extend the MLE to log_n' = 3*log_r variables
    // by padding the evaluation table with zeros on the extra dimensions, and we evaluate at
    // extra coordinates fixed to 0.
    //
    // Equivalently: embed f: {0,1}^{log_n} -> F into f': {0,1}^{log_n'} -> F where
    // f'(x, y) = f(x) if y == 0...0 else 0, and evaluate at (r', 0...0).
    let need = 3 * log_r;
    let mut r_ext: Vec<F> = Vec::with_capacity(need);
    r_ext.extend_from_slice(r_prime);
    r_ext.resize(need, F::ZERO);

    // Split r' into three blocks (in protocol order) and convert using the same
    // tensor-weight convention as the rest of this codebase (`ts_weights`).
    let r0 = &r_ext[0..log_r];
    let r1 = &r_ext[log_r..2 * log_r];
    let r2 = &r_ext[2 * log_r..3 * log_r];

    let x0 = ts_weights(r0);
    let x1 = ts_weights(r1);
    let x2 = ts_weights(r2);

    debug_assert_eq!(x0.len(), r);
    debug_assert_eq!(x1.len(), r);
    debug_assert_eq!(x2.len(), r);

    Ok((x0, x1, x2))
}

/// Batchlin PCS commitment for k_g digits.
#[derive(Clone, Debug)]
pub struct BatchlinPcsCommitment<F: PrimeField> {
    /// PCS commitment roots for each digit, t^{(dig)} of length κ*n_inner.
    /// Here n_inner is the internal PCS dimension (typically smaller due to structure).
    pub t: Vec<Vec<F>>,
}

/// Batchlin PCS opening proof for k_g digits.
#[derive(Clone, Debug)]
pub struct BatchlinPcsProof<F: PrimeField> {
    /// Individual proofs per digit.
    pub proofs: Vec<FoldingPcsL2ProofCore<F>>,
    /// Batching challenge γ (derived from transcript, stored for verification).
    /// Not strictly needed if we re-derive, but useful for debugging.
    pub gamma: Option<F>,
}

/// Compute the batched commitment and claimed value for k_g digits.
///
/// Given:
/// - commitments `t[dig]` for each digit
/// - claimed values `u[dig]` for each digit
/// - batching challenge `γ`
///
/// Returns:
/// - `t_batch = Σ_dig γ^dig * t[dig]`
/// - `u_batch = Σ_dig γ^dig * u[dig]`
pub fn batch_commitments_and_values<F: PrimeField>(
    t: &[Vec<F>],
    u: &[F],
    gamma: F,
) -> (Vec<F>, F) {
    let k_g = t.len();
    assert_eq!(u.len(), k_g);

    if k_g == 0 {
        return (vec![], F::ZERO);
    }

    let t_len = t[0].len();

    // Compute batched commitment: t_batch = Σ_dig γ^dig * t[dig]
    let t_batch: Vec<F> = (0..t_len)
        .into_par_iter()
        .map(|j| {
            let mut gamma_pow = F::ONE;
            let mut sum = F::ZERO;
            for dig in 0..k_g {
                sum += gamma_pow * t[dig][j];
                gamma_pow *= gamma;
            }
            sum
        })
        .collect();

    // Compute batched value: u_batch = Σ_dig γ^dig * u[dig]
    let mut gamma_pow = F::ONE;
    let mut u_batch = F::ZERO;
    for dig in 0..k_g {
        u_batch += gamma_pow * u[dig];
        gamma_pow *= gamma;
    }

    (t_batch, u_batch)
}

/// Domain separator for batchlin PCS in transcript.
pub const BATCHLIN_PCS_DOMAIN_SEP: u128 = 0x4241_5443_484C_494E_5F50_4353_5F56_3100; // "BATCHLIN_PCS_V1\0"

/// Deterministic PCS params for **batched scalar** batchlin PCS.
///
/// Shape:
/// - `kappa=1`, `n=1` ⇒ commitment length 1 and scalar evaluation
/// - `r = 2^{log_n/3}` with `r^3 = 2^{log_n}`
///
/// Notes:
/// - This is a *benchmark/prototype* parameter generator: it derives A from a fixed seed + log_n.
/// - `delta,alpha` are chosen so `delta^alpha >= modulus` for 64-bit-ish prime fields (Frog).
pub fn batchlin_scalar_pcs_params<F: PrimeField>(log_n: usize) -> Result<FoldingPcsL2Params<F>, String> {
    let log_r = (log_n + 2) / 3; // ceil(log_n/3)
    let r = 1usize << log_r;
    let kappa = 1usize;
    let n = 1usize;

    // Chosen to satisfy delta^alpha >= modulus for ~64-bit moduli:
    // (2^32)^2 = 2^64.
    let delta = 1u64 << 32;
    let alpha = 2usize;

    // Generous coefficient bounds for the small vectors (prototype).
    let beta0 = 1u64 << 63;
    let beta1 = beta0;
    let beta2 = beta0;

    // A is 1 × (r*n*alpha).
    let cols = r * n * alpha;
    let mut h = Sha256::new();
    h.update(b"BATCHLIN_SCALAR_PCS_A_V1");
    h.update(&(log_n as u64).to_le_bytes());
    let seed: [u8; 32] = h.finalize().into();
    let mut rng = ChaCha20Rng::from_seed(seed);
    let mut a_data = Vec::with_capacity(cols);
    for _ in 0..cols {
        a_data.push(F::from(rng.next_u64()));
    }
    let a = DenseMatrix::new(n, cols, a_data);

    Ok(FoldingPcsL2Params {
        r,
        kappa,
        n,
        delta,
        alpha,
        beta0,
        beta1,
        beta2,
        a,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::Field;
    use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as R;
    use stark_rings::PolyRing;

    // Use Frog's base prime field (already in symphony deps)
    type F = <<R as PolyRing>::BaseRing as Field>::BasePrimeField;

    #[test]
    fn test_ts_weights_sum_to_one() {
        // `ts_weights` is the canonical tensor-weight convention used across this codebase.
        // It should always sum to 1 (partition of unity).
        let point = vec![F::from(2u64), F::from(3u64)]; // arbitrary non-binary values
        let weights = ts_weights(&point);
        assert_eq!(weights.len(), 4);
        let sum: F = weights.iter().sum();
        assert_eq!(sum, F::ONE);
    }

    #[test]
    fn test_ts_weights_at_binary_point() {
        // At a binary point, `ts_weights(point)` should be 1 at the corresponding index and 0 elsewhere.
        // Note: `ts_weights` follows this repo's MLE variable ordering (see `symphony_coins::ts_weights`).
        let point = vec![F::ONE, F::ZERO];
        let weights = ts_weights(&point);
        assert_eq!(weights.len(), 4);

        // Determine the selected index by brute force (avoid hardcoding bit-order assumptions).
        let ones = weights.iter().filter(|&&w| w == F::ONE).count();
        assert_eq!(ones, 1);
        let idx = weights.iter().position(|&w| w == F::ONE).unwrap();
        for (i, &w) in weights.iter().enumerate() {
            if i == idx {
                assert_eq!(w, F::ONE);
            } else {
                assert_eq!(w, F::ZERO);
            }
        }
    }

    #[test]
    fn test_mle_point_to_tensor() {
        let config = BatchlinPcsConfig::for_sp1(24, 4, 6).unwrap();

        // Create a test point of length 24
        let r_prime: Vec<F> = (0..24).map(|i| F::from(i as u64)).collect();

        let (x0, x1, x2) = mle_point_to_tensor(&r_prime, &config).unwrap();

        // Each xk should have length r = 256
        assert_eq!(x0.len(), 256);
        assert_eq!(x1.len(), 256);
        assert_eq!(x2.len(), 256);

        // Each xk should sum to 1
        let sum0: F = x0.iter().sum();
        let sum1: F = x1.iter().sum();
        let sum2: F = x2.iter().sum();
        assert_eq!(sum0, F::ONE);
        assert_eq!(sum1, F::ONE);
        assert_eq!(sum2, F::ONE);
    }

    #[test]
    fn test_batch_commitments() {
        let t = vec![
            vec![F::from(1u64), F::from(2u64)],
            vec![F::from(3u64), F::from(4u64)],
        ];
        let u = vec![F::from(10u64), F::from(20u64)];
        let gamma = F::from(5u64);

        let (t_batch, u_batch) = batch_commitments_and_values(&t, &u, gamma);

        // t_batch[j] = t[0][j] + γ * t[1][j]
        assert_eq!(t_batch[0], F::from(1u64) + F::from(5u64) * F::from(3u64)); // 1 + 15 = 16
        assert_eq!(t_batch[1], F::from(2u64) + F::from(5u64) * F::from(4u64)); // 2 + 20 = 22

        // u_batch = u[0] + γ * u[1] = 10 + 5*20 = 110
        assert_eq!(u_batch, F::from(110u64));
    }
}
