//! High-level “Rev2 pipeline” builder helpers.
//!
//! The large-field DPP path in Rev2 (Section 5) is typically used as:
//!   (1) Ensure Boolean proofs (Claim 5.3) via binary decomposition
//!   (2) Embed to a bounded promise FLPCP (Theorem 5.6 / Construction 5.13)
//!   (3) Pack k queries into 1 query (Theorem 5.7 / Construction 5.21)
//!
//! This module provides a helper to compose these wrappers into a single DPP verifier object.

use ark_ff::PrimeField;

use crate::{
    boolean_proof::{BooleanProofFlpcp, BooleanProofFlpcpSparse},
    embedding::{EmbeddedFlpcp, EmbeddedFlpcpSparse, EmbeddingParams},
    packing::{
        BoundedFlpcp, BoundedFlpcpSparse, DppFromBoundedFlpcp, DppFromBoundedFlpcpSparse, PackedDppParams,
        PackingError,
    },
};

/// Build the Rev2 “Booleanize → Embed → Pack” DPP verifier.
pub fn build_rev2_dpp<FSmall: PrimeField, FLarge: PrimeField, V: BoundedFlpcp<FSmall>>(
    flpcp: V,
    embed: EmbeddingParams,
    pack: PackedDppParams,
) -> DppFromBoundedFlpcp<FLarge, EmbeddedFlpcp<FSmall, FLarge, BooleanProofFlpcp<FSmall, V>>> {
    let boolized = BooleanProofFlpcp::<FSmall, V>::new(flpcp);
    let embedded = EmbeddedFlpcp::<FSmall, FLarge, _>::new(boolized, embed);
    DppFromBoundedFlpcp::<FLarge, _>::new(embedded, pack)
}

/// Build the Rev2 “Embed → Pack” DPP verifier in sparse-query form (no Booleanization).
pub fn build_rev2_dpp_sparse<FSmall: PrimeField, FLarge: PrimeField, V: BoundedFlpcpSparse<FSmall>>(
    flpcp: V,
    embed: EmbeddingParams,
    pack: PackedDppParams,
) -> DppFromBoundedFlpcpSparse<FLarge, EmbeddedFlpcpSparse<FSmall, FLarge, V>> {
    let embedded = EmbeddedFlpcpSparse::<FSmall, FLarge, _>::new(flpcp, embed);
    DppFromBoundedFlpcpSparse::<FLarge, _>::new(embedded, pack)
}

/// Build the Rev2 “Booleanize → Embed → Pack” DPP verifier in sparse-query form.
pub fn build_rev2_dpp_sparse_boolean<
    FSmall: PrimeField,
    FLarge: PrimeField,
    V: BoundedFlpcpSparse<FSmall>,
>(
    flpcp: V,
    embed: EmbeddingParams,
    pack: PackedDppParams,
) -> DppFromBoundedFlpcpSparse<FLarge, EmbeddedFlpcpSparse<FSmall, FLarge, BooleanProofFlpcpSparse<FSmall, V>>> {
    let boolized = BooleanProofFlpcpSparse::<FSmall, V>::new(flpcp);
    let embedded = EmbeddedFlpcpSparse::<FSmall, FLarge, _>::new(boolized, embed);
    DppFromBoundedFlpcpSparse::<FLarge, _>::new(embedded, pack)
}

/// Build the Rev2 “Booleanize → Embed → Pack” DPP verifier, choosing `PackedDppParams`
/// automatically from the resulting bounded FLPCP’s bounds and the target field modulus.
///
/// This removes ad-hoc `ell` tuning and uses the deterministic safe upper bound derived
/// from Construction 5.21’s weight ranges (see `PackedDppParams::from_bounds_max`).
pub fn build_rev2_dpp_auto<FSmall: PrimeField, FLarge: PrimeField, V: BoundedFlpcp<FSmall>>(
    flpcp: V,
    embed: EmbeddingParams,
) -> Result<
    DppFromBoundedFlpcp<FLarge, EmbeddedFlpcp<FSmall, FLarge, BooleanProofFlpcp<FSmall, V>>>,
    PackingError,
> {
    let boolized = BooleanProofFlpcp::<FSmall, V>::new(flpcp);
    let embedded = EmbeddedFlpcp::<FSmall, FLarge, _>::new(boolized, embed);
    let pack = PackedDppParams::from_bounds_max::<FLarge>(&embedded.bounds_b())?;
    Ok(DppFromBoundedFlpcp::<FLarge, _>::new(embedded, pack))
}

/// Sparse-query version of `build_rev2_dpp_auto` (no Booleanization).
pub fn build_rev2_dpp_sparse_auto<FSmall: PrimeField, FLarge: PrimeField, V: BoundedFlpcpSparse<FSmall>>(
    flpcp: V,
    embed: EmbeddingParams,
) -> Result<DppFromBoundedFlpcpSparse<FLarge, EmbeddedFlpcpSparse<FSmall, FLarge, V>>, PackingError> {
    let embedded = EmbeddedFlpcpSparse::<FSmall, FLarge, _>::new(flpcp, embed);
    let pack = PackedDppParams::from_bounds_max::<FLarge>(&embedded.bounds_b())?;
    Ok(DppFromBoundedFlpcpSparse::<FLarge, _>::new(embedded, pack))
}

/// Sparse-query version of `build_rev2_dpp_auto` with Booleanization.
pub fn build_rev2_dpp_sparse_boolean_auto<
    FSmall: PrimeField,
    FLarge: PrimeField,
    V: BoundedFlpcpSparse<FSmall>,
>(
    flpcp: V,
    embed: EmbeddingParams,
) -> Result<
    DppFromBoundedFlpcpSparse<FLarge, EmbeddedFlpcpSparse<FSmall, FLarge, BooleanProofFlpcpSparse<FSmall, V>>>,
    PackingError,
> {
    let boolized = BooleanProofFlpcpSparse::<FSmall, V>::new(flpcp);
    let embedded = EmbeddedFlpcpSparse::<FSmall, FLarge, _>::new(boolized, embed);
    let pack = PackedDppParams::from_bounds_max::<FLarge>(&embedded.bounds_b())?;
    Ok(DppFromBoundedFlpcpSparse::<FLarge, _>::new(embedded, pack))
}

/// Helper: choose k' such that (0.51)^{k'} <= 2^{-security_bits}.
pub fn k_prime_for_128bit() -> usize {
    k_prime_for_security_bits(128)
}

pub fn k_prime_for_security_bits(security_bits: usize) -> usize {
    // k' >= security_bits / log2(1/0.51)
    let denom = (1.0f64 / 0.51f64).log2();
    ((security_bits as f64) / denom).ceil() as usize
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::packing::{FlpcpPredicate, PackingError};
    use ark_ff::{Fp64, MontBackend, MontConfig};
    use num_bigint::BigInt;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[derive(MontConfig)]
    #[modulus = "101"]
    #[generator = "2"]
    pub struct F101Config;
    type FSmall = Fp64<MontBackend<F101Config, 1>>;

    #[derive(MontConfig)]
    // Goldilocks prime (2^64 - 2^32 + 1): large 2-adicity and plenty of headroom for tests.
    #[modulus = "18446744069414584321"]
    #[generator = "7"]
    pub struct FBigConfig;
    type FLarge = Fp64<MontBackend<FBigConfig, 1>>;

    #[derive(Clone, Debug)]
    struct ToyMulEq;

    impl crate::packing::BoundedFlpcp<FSmall> for ToyMulEq {
        fn n(&self) -> usize { 1 }
        fn m(&self) -> usize { 0 }
        fn k(&self) -> usize { 3 }
        fn bounds_b(&self) -> Vec<BigInt> { vec![BigInt::from(2), BigInt::from(2), BigInt::from(2)] }
        fn sample_queries_and_predicate(
            &self,
            _rng: &mut dyn rand::RngCore,
            _x: &[FSmall],
        ) -> Result<(Vec<Vec<FSmall>>, FlpcpPredicate<FSmall>), String> {
            // Answers are (1,1,1) regardless of x (toy).
            Ok((vec![vec![FSmall::from(1u64)], vec![FSmall::from(1u64)], vec![FSmall::from(1u64)]], FlpcpPredicate::MulEq))
        }
    }

    #[test]
    fn test_rev2_pipeline_toy_accepts() {
        let flpcp = ToyMulEq;
        let mut rng = ChaCha20Rng::seed_from_u64(123);

        // Keep parameters small so the packing “no wrap” condition (Claim 5.22) holds.
        let dpp = build_rev2_dpp::<FSmall, FLarge, _>(
            flpcp,
            EmbeddingParams { gamma: 2, assume_boolean_proof: true, k_prime: 0 },
            PackedDppParams { ell: 2 },
        );

        let x = vec![FLarge::from(1u64)];
        let ok = dpp.verify(&mut rng, &x, &[]).expect("verify result");
        assert!(ok);
    }

    #[test]
    fn test_rev2_pipeline_rejects_bad_length() {
        let flpcp = ToyMulEq;
        let mut rng = ChaCha20Rng::seed_from_u64(124);
        let dpp = build_rev2_dpp::<FSmall, FLarge, _>(
            flpcp,
            EmbeddingParams { gamma: 3, assume_boolean_proof: true, k_prime: 0 },
            PackedDppParams { ell: 16 },
        );
        let x = vec![FLarge::from(1u64), FLarge::from(2u64)];
        let res = dpp.verify(&mut rng, &x, &[]);
        assert!(matches!(res, Err(PackingError::InvalidParams)));
    }
}

