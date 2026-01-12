//! Helpers to bind WE locks to a statement (proof-agnostic arming).
//!
//! Following Symphony’s model, the armer derives per-lock coins from a **statement hash** that binds:
//! - a program/verifier id (`vk_hash`)
//! - an instance digest (`r1cs_digest`)
//! - a gate digest (`gate_digest`, i.e. exact WE gate relation version)
//! - and the **public inputs**
//!
//! Importantly, this intentionally does **not** include prover-chosen proof artifacts (e.g. Ajtai
//! commitments to witness), since arming must be possible without observing a specific proving run.

use ark_ff::{BigInteger, Field, PrimeField};
use sha2::{Digest, Sha256};
use stark_rings::OverField;

/// Fixed parameters that must be statement-bound for WE/DPP.
///
/// These are implicit in native Rust (loop bounds, type params), but must be explicit in an
/// arithmetized WE gate to prevent “reinterpretation under different sizes”.
#[derive(Clone, Debug)]
pub struct WeParams {
    pub nvars_setchk: u64,
    pub degree_setchk: u64,
    pub nvars_cm: u64,
    pub degree_cm: u64,
    pub kappa: u64,
    pub ring_dim_d: u64,
    /// Balanced decomposition base `b` used in LF+ range-check (`RgChk`) digit decomposition.
    ///
    /// For SP1 BabyBear-in-Frog boundedness, this should be `2^16`.
    pub decomp_b: u64,
    pub k: u64,
    pub l: u64,
    pub mlen: u64,
}

impl WeParams {
    pub fn to_field_vec<BF: PrimeField>(&self) -> Vec<BF> {
        vec![
            BF::from(self.nvars_setchk),
            BF::from(self.degree_setchk),
            BF::from(self.nvars_cm),
            BF::from(self.degree_cm),
            BF::from(self.kappa),
            BF::from(self.ring_dim_d),
            BF::from(self.decomp_b),
            BF::from(self.k),
            BF::from(self.l),
            BF::from(self.mlen),
        ]
    }
}

/// SHA256 hash of the **public** statement for `R_WE` (LF+).
///
/// This intentionally excludes proof artifacts (e.g. witness commitments).
pub fn we_statement_hash_lf_plus<R: OverField>(
    vk_hash: [u8; 32],
    r1cs_digest: [u8; 32],
    gate_digest: [u8; 32],
    public_inputs: &[R::BaseRing],
) -> [u8; 32]
where
    R::BaseRing: Field,
{
    let mut h = Sha256::new();
    h.update(b"LATTICEFOLD_PLUS_WE_STATEMENT_V1");
    h.update(&vk_hash);
    h.update(&r1cs_digest);
    h.update(&gate_digest);

    h.update(&(public_inputs.len() as u64).to_le_bytes());
    for x in public_inputs {
        for fp in x.to_base_prime_field_elements() {
            h.update(fp.into_bigint().to_bytes_le());
        }
    }

    h.finalize().into()
}

/// Canonical encoding for a 32-byte digest as a field element.
///
/// Interpret the digest as a **little-endian** integer and reduce mod p.
pub fn digest32_to_field<BF: PrimeField>(digest: [u8; 32]) -> BF {
    BF::from_le_bytes_mod_order(&digest)
}

/// Collision-robust encoding of a 32-byte digest as **256 boolean field elements** (bits).
///
/// - Output length is exactly 256.
/// - Bit order is little-endian within each byte (LSB-first), matching typical `from_le_bytes` conventions.
///
/// This is intended for WE/DPP public inputs when we want:
/// - statement binding without relying on `mod p` reduction uniqueness, and
/// - DPP parameters that assume small / boolean public inputs.
pub fn digest32_to_bits_field<BF: PrimeField>(digest: [u8; 32]) -> Vec<BF> {
    let mut out = Vec::with_capacity(256);
    for &b in &digest {
        for i in 0..8 {
            let bit = (b >> i) & 1;
            out.push(BF::from(bit as u64));
        }
    }
    debug_assert_eq!(out.len(), 256);
    out
}


/// Gate digest for the LF+ WE gate relation.
///
/// In the SP1/WE arming model this should be a **precomputed constant** identifying the exact
/// WE gate version (code + constraint system generation), not something derived per-proof.
///
/// TODO: fill with the production gate digest for the SP1 WE gate.
pub const LFP_WE_GATE_DIGEST_V1: [u8; 32] = [0u8; 32];

/// Deterministic statement encoding for WE/DPP.
///
/// Current encoding:
/// - `x[0] = 1` (shared constant slot convention)
/// - fixed params (see `WeParams`)
/// - optional extra statement elements (e.g. commitment surface limbs, transcript-bound absorbs)
///
/// NOTE: The exact set of extra statement elements is decided by the WE arithmetizer; this module
/// just provides a stable prefix layout for params.
pub fn encode_public_x<BF: PrimeField>(params: &WeParams, extra: &[BF]) -> Vec<BF> {
    let mut out = Vec::with_capacity(1 + 10 + extra.len());
    out.push(BF::ONE);
    out.extend(params.to_field_vec::<BF>());
    out.extend_from_slice(extra);
    out
}

