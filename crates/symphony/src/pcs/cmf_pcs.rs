//! PCS-backed `cm_f` commitment surface (Option A).
//!
//! This module provides helper utilities to:
//! - choose deterministic FoldingPCS(ℓ=2) parameters for a given message length
//! - commit/open a **padded** flattened witness vector
//! - pack/unpack the PCS commitment surface into ring elements for existing Π_fold plumbing

use ark_ff::PrimeField;
use rand::{RngCore, SeedableRng};
use sha2::{Digest, Sha256};

use crate::pcs::folding_pcs_l2::{DenseMatrix, FoldingPcsL2Params};

/// Domain separator for PCS#1 (cm_f PCS coin phase).
pub const CMF_PCS_DOMAIN_SEP: u128 = 0x434d465f5043535f5631000000000000; // "CMF_PCS_V1\0..."

/// Bytes needed for Fig.5 `C1/C2` coins:
/// each C is (r*κ)×κ bits, so total bytes = ceil(2*r*κ^2 / 8) = ceil(r*κ^2/4).
pub fn cmf_pcs_coin_bytes_len<F: PrimeField>(p: &FoldingPcsL2Params<F>) -> usize {
    let bits = 2usize
        .saturating_mul(p.r)
        .saturating_mul(p.kappa)
        .saturating_mul(p.kappa);
    (bits + 7) / 8
}

/// Deterministically build FoldingPCS(ℓ=2) params for committing to a (possibly padded) message.
///
/// We use the simplest scalar-output shape:
/// - `n = 1`
/// - `kappa = kappa_commit`  (commitment length)
/// - choose `r` as a power-of-two so `f_len = r^3 * kappa >= flat_len`
///
/// This avoids requiring `log2(flat_len) % 3 == 0` by allowing padding.
pub fn cmf_pcs_params_for_flat_len<F: PrimeField>(
    flat_len: usize,
    kappa_commit: usize,
) -> Result<FoldingPcsL2Params<F>, String> {
    if kappa_commit == 0 {
        return Err("cmf_pcs_params_for_flat_len: kappa_commit must be > 0".to_string());
    }
    if flat_len == 0 {
        return Err("cmf_pcs_params_for_flat_len: flat_len must be > 0".to_string());
    }

    // Choose alpha/delta to satisfy delta^alpha >= modulus for ~64-bit prime fields.
    let delta = 1u64 << 32;
    let alpha = 2usize;

    // Use generous bounds; this is a commitment surface, not a tight range-proof.
    let beta0 = 1u64 << 63;
    let beta1 = beta0;
    let beta2 = beta0;

    let n = 1usize;
    let blocks = (flat_len + kappa_commit - 1) / kappa_commit; // ceil(flat_len / kappa)
    let blocks_pow2 = blocks.next_power_of_two();

    // Find smallest power-of-two r with r^3 >= blocks_pow2.
    let mut r = 1usize;
    while r
        .checked_mul(r)
        .and_then(|x| x.checked_mul(r))
        .ok_or_else(|| "cmf_pcs_params_for_flat_len: r^3 overflow".to_string())?
        < blocks_pow2
    {
        r <<= 1;
    }

    // A is (n × (r*n*alpha)) = (1 × (r*alpha)).
    let cols = r * n * alpha;
    let mut h = Sha256::new();
    h.update(b"CMF_PCS_A_V1");
    h.update(&(flat_len as u64).to_le_bytes());
    h.update(&(kappa_commit as u64).to_le_bytes());
    h.update(&(r as u64).to_le_bytes());
    let seed: [u8; 32] = h.finalize().into();
    let mut rng = rand_chacha::ChaCha20Rng::from_seed(seed);
    let mut a_data = Vec::with_capacity(n * cols);
    for _ in 0..(n * cols) {
        a_data.push(F::from(rng.next_u64()));
    }
    let a = DenseMatrix::new(n, cols, a_data);

    Ok(FoldingPcsL2Params {
        r,
        kappa: kappa_commit,
        n,
        delta,
        alpha,
        beta0,
        beta1,
        beta2,
        a,
    })
}

/// Pad a flattened message to `p.f_len()` (zero-extended).
pub fn pad_flat_message<F: PrimeField>(p: &FoldingPcsL2Params<F>, flat: &[F]) -> Vec<F> {
    let mut out = vec![F::ZERO; p.f_len()];
    let take = core::cmp::min(out.len(), flat.len());
    out[..take].copy_from_slice(&flat[..take]);
    out
}

/// Pack a base-field commitment surface `t` into ring elements by filling coefficients.
///
/// This is a compatibility shim so existing Π_fold APIs that expect `Vec<R>` can absorb `cm_f`
/// without changing types.
pub fn pack_t_as_ring<R>(t: &[R::BaseRing]) -> Vec<R>
where
    R: stark_rings::Ring + stark_rings::PolyRing,
    R::BaseRing: PrimeField,
{
    let d = R::dimension();
    let mut out = Vec::with_capacity((t.len() + d - 1) / d);
    for chunk in t.chunks(d) {
        let mut re = R::ZERO;
        for (i, &x) in chunk.iter().enumerate() {
            re.coeffs_mut()[i] = x;
        }
        out.push(re);
    }
    out
}

