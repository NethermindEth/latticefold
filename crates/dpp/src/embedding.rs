//! Bounded embedding transformation (Construction 5.13 / Theorem 5.6) — prototype.
//!
//! This module implements the *mechanical* embedding wrapper that:
//! - embeds k FLPCP query vectors from a smaller prime field F_p into a larger prime field F_{p'}
//! - adds a random “bound test” component `p * u'` to the first query on the proof coordinates
//! - enforces boundedness of answers in integer representatives
//! - checks the original predicate **mod p** (via `MulEqModP`)
//!
//! Important:
//! - This is only the code-level transformation. Correctness/soundness depends on parameter regimes
//!   from Theorem 5.6, which must be instantiated carefully.

use ark_ff::{BigInteger, PrimeField};
use num_bigint::BigInt;
use num_traits::{One, Zero};
use rand::RngCore;
use thiserror::Error;

use crate::packing::{BoundedFlpcp, BoundedFlpcpSparse, FlpcpPredicate};
use crate::sparse::SparseVec;

#[derive(Debug, Error)]
pub enum EmbeddingError {
    #[error("invalid parameters")]
    InvalidParams,
}

/// Wrapper parameters for the embedding.
#[derive(Clone, Debug)]
pub struct EmbeddingParams {
    /// γ: randomness range for the bound test (`u ← F_γ` in the paper).
    pub gamma: u64,
    /// If true, assume the proof has been Booleanized (Claim 5.3), so inner-product
    /// answers scale linearly in `p` (not quadratically in `p`).
    ///
    /// This materially affects the boundedness parameters (and thus the packing modulus condition).
    pub assume_boolean_proof: bool,
    /// If nonzero, enable the Claim 5.28 strong-soundness variant by adding `k_prime`
    /// additional random bound-test queries.
    pub k_prime: usize,
}

/// Embed an FLPCP over `FSmall` into a bounded FLPCP over `FLarge`.
#[derive(Clone, Debug)]
pub struct EmbeddedFlpcp<FSmall: PrimeField, FLarge: PrimeField, V: BoundedFlpcp<FSmall>> {
    pub inner: V,
    pub params: EmbeddingParams,
    _pd: core::marker::PhantomData<(FSmall, FLarge)>,
}

impl<FSmall: PrimeField, FLarge: PrimeField, V: BoundedFlpcp<FSmall>> EmbeddedFlpcp<FSmall, FLarge, V> {
    pub fn new(inner: V, params: EmbeddingParams) -> Self {
        Self { inner, params, _pd: core::marker::PhantomData }
    }
}

impl<FSmall: PrimeField, FLarge: PrimeField, V: BoundedFlpcp<FSmall>> BoundedFlpcp<FLarge>
    for EmbeddedFlpcp<FSmall, FLarge, V>
{
    fn n(&self) -> usize {
        self.inner.n()
    }

    fn m(&self) -> usize {
        self.inner.m()
    }

    fn k(&self) -> usize {
        self.inner.k() + self.params.k_prime
    }

    fn bounds_b(&self) -> Vec<BigInt> {
        let k_inner = self.inner.k();
        let k_prime = self.params.k_prime;
        let k_total = k_inner + k_prime;
        let mut out = Vec::with_capacity(k_total);

        if k_total == 0 {
            return out;
        }

        // Two modes:
        //
        // - Standard (Construction 5.13 / Theorem 5.6 mechanics):
        //     b1 = μ̂ (pγ-1)/2, b2 = μ̂ (p-1)/2 for remaining inner queries.
        //
        // - Strong-soundness variant (Claim 5.28):
        //     b = (μ̂ p γ / 2, ..., μ̂ p γ / 2) for all k+k' queries.
        //
        // Here μ̂ is the total length of (x||π).
        //
        // Two modes:
        // - If `assume_boolean_proof`, we use the (much tighter) linear-in-p bounds that match the
        //   Rev2 “Booleanize → Embed → Pack” pipeline (Claim 5.3 + Construction 5.13).
        // - Otherwise, we use conservative bounds that are valid for general proofs, where each term
        //   contributes up to half^2 (quadratic-in-p).
        let mu_hat = BigInt::from((self.n() + self.m()) as u64);
        let p_small = modulus_bigint::<FSmall>();
        let gamma = BigInt::from(self.params.gamma);
        if self.params.assume_boolean_proof {
            // Linear-in-p bounds (matches earlier prototype and Claim 5.28 / Construction 5.13 analysis).
            if k_prime > 0 {
                let b = &mu_hat * (&p_small * gamma) / BigInt::from(2u64);
                out.resize(k_total, b);
                return out;
            }
            let b1 = &mu_hat * ((&p_small * gamma - BigInt::one()) / BigInt::from(2u64));
            let b2 = &mu_hat * ((&p_small - BigInt::one()) / BigInt::from(2u64));
            out.push(b1);
            for _ in 1..k_inner {
                out.push(b2.clone());
            }
            return out;
        }

        // General-proof (quadratic-in-p) bounds.
        let half = (&p_small - BigInt::one()) / BigInt::from(2u64);
        let half2 = &half * &half;
        let p_gamma_half = (&p_small * gamma * &half) / BigInt::from(2u64);

        let b1 = &mu_hat * (&half2 + &p_gamma_half);
        let b2 = &mu_hat * &half2;

        if k_prime > 0 {
            out.resize(k_total, b1);
            return out;
        }

        out.push(b1);
        for _ in 1..k_inner {
            out.push(b2.clone());
        }
        out
    }

    fn sample_queries_and_predicate(
        &self,
        rng: &mut dyn RngCore,
        x: &[FLarge],
    ) -> Result<(Vec<Vec<FLarge>>, FlpcpPredicate<FLarge>), String> {
        if x.len() != self.n() {
            return Err("bad input".to_string());
        }

        // Downcast x into FSmall assuming x is in {0,1} or otherwise within FSmall embed.
        // (For general use, inputs should be provided already as small-field elements.)
        let x_small = x
            .iter()
            .map(|xi| {
                let bi = BigInt::from_bytes_le(num_bigint::Sign::Plus, &xi.into_bigint().to_bytes_le());
                let p = modulus_bigint::<FSmall>();
                let r = bi % &p;
                FSmall::from_le_bytes_mod_order(&r.to_bytes_le().1)
            })
            .collect::<Vec<_>>();

        let (q_small, pred_small) = self.inner.sample_queries_and_predicate(rng, &x_small)?;
        let k_inner = q_small.len();
        let k_prime = self.params.k_prime;

        // Embed small-field queries into large field.
        let mut q_large = q_small
            .iter()
            .map(|row| row.iter().map(|c| embed_small_to_large::<FSmall, FLarge>(c)).collect::<Vec<_>>())
            .collect::<Vec<_>>();

        // Modify first query: q1' = q1 + p * u' where u' = (0^n || u) and u is uniform in F_gamma.
        if k_inner > 0 && self.m() > 0 {
            let p_small = modulus_bigint::<FSmall>();
            let p_as_large = FLarge::from_le_bytes_mod_order(&p_small.to_bytes_le().1);

            for j in 0..self.m() {
                let u = sample_centered_int(rng, self.params.gamma);
                let u_large = FLarge::from_le_bytes_mod_order(&u.to_bytes_le().1);
                let idx = self.n() + j;
                q_large[0][idx] += p_as_large * u_large;
            }
        }

        // Strong-soundness variant (Claim 5.28): add k' extra bound-test queries
        //   q_{k+i} = (0 || p u_i),  u_i ← F_gamma^m
        if k_prime > 0 && self.m() > 0 {
            let p_small = modulus_bigint::<FSmall>();
            let p_as_large = FLarge::from_le_bytes_mod_order(&p_small.to_bytes_le().1);
            for _ in 0..k_prime {
                let mut row = vec![FLarge::ZERO; self.n() + self.m()];
                for j in 0..self.m() {
                    let u = sample_centered_int(rng, self.params.gamma);
                    let u_large = FLarge::from_le_bytes_mod_order(&u.to_bytes_le().1);
                    row[self.n() + j] = p_as_large * u_large;
                }
                q_large.push(row);
            }
        } else if k_prime > 0 && self.m() == 0 {
            // Still add the correct number of queries (they will be zero rows).
            for _ in 0..k_prime {
                q_large.push(vec![FLarge::ZERO; self.n() + self.m()]);
            }
        }

        // Predicate: check underlying predicate mod p_small when needed, and ignore the tail
        // answers (which are only required to be bounded, enforced by decoding).
        let base_pred = match pred_small {
            FlpcpPredicate::MulEq => FlpcpPredicate::MulEqModP { p_small: modulus_bigint::<FSmall>() },
            FlpcpPredicate::MulEqAndZero => FlpcpPredicate::Unsupported(core::marker::PhantomData),
            FlpcpPredicate::AllZero => FlpcpPredicate::AllZero,
            FlpcpPredicate::MulEqModP { .. } => FlpcpPredicate::Unsupported(core::marker::PhantomData),
            FlpcpPredicate::Unsupported(_) => FlpcpPredicate::Unsupported(core::marker::PhantomData),
            FlpcpPredicate::IgnoreTail { .. } => FlpcpPredicate::Unsupported(core::marker::PhantomData),
        };
        let pred = if k_prime > 0 {
            FlpcpPredicate::IgnoreTail { inner: Box::new(base_pred), tail_len: k_prime }
        } else {
            base_pred
        };

        Ok((q_large, pred))
    }
}

/// Embed a sparse-query FLPCP over `FSmall` into a bounded sparse-query FLPCP over `FLarge`.
#[derive(Clone, Debug)]
pub struct EmbeddedFlpcpSparse<FSmall: PrimeField, FLarge: PrimeField, V: BoundedFlpcpSparse<FSmall>> {
    pub inner: V,
    pub params: EmbeddingParams,
    _pd: core::marker::PhantomData<(FSmall, FLarge)>,
}

impl<FSmall: PrimeField, FLarge: PrimeField, V: BoundedFlpcpSparse<FSmall>>
    EmbeddedFlpcpSparse<FSmall, FLarge, V>
{
    pub fn new(inner: V, params: EmbeddingParams) -> Self {
        Self { inner, params, _pd: core::marker::PhantomData }
    }
}

impl<FSmall: PrimeField, FLarge: PrimeField, V: BoundedFlpcpSparse<FSmall>> BoundedFlpcpSparse<FLarge>
    for EmbeddedFlpcpSparse<FSmall, FLarge, V>
{
    fn n(&self) -> usize {
        self.inner.n()
    }

    fn m(&self) -> usize {
        self.inner.m()
    }

    fn k(&self) -> usize {
        self.inner.k() + self.params.k_prime
    }

    fn bounds_b(&self) -> Vec<BigInt> {
        // Same logic as the dense wrapper.
        let k_inner = self.inner.k();
        let k_prime = self.params.k_prime;
        let k_total = k_inner + k_prime;
        let mut out = Vec::with_capacity(k_total);

        if k_total == 0 {
            return out;
        }

        let mu_hat = BigInt::from((self.n() + self.m()) as u64);
        let p_small = modulus_bigint::<FSmall>();
        let gamma = BigInt::from(self.params.gamma);
        if self.params.assume_boolean_proof {
            if k_prime > 0 {
                let b = &mu_hat * (&p_small * gamma) / BigInt::from(2u64);
                out.resize(k_total, b);
                return out;
            }
            let b1 = &mu_hat * ((&p_small * gamma - BigInt::one()) / BigInt::from(2u64));
            let b2 = &mu_hat * ((&p_small - BigInt::one()) / BigInt::from(2u64));
            out.push(b1);
            for _ in 1..k_inner {
                out.push(b2.clone());
            }
            return out;
        }

        let half = (&p_small - BigInt::one()) / BigInt::from(2u64);
        let half2 = &half * &half;
        let p_gamma_half = (&p_small * gamma * &half) / BigInt::from(2u64);

        let b1 = &mu_hat * (&half2 + &p_gamma_half);
        let b2 = &mu_hat * &half2;

        if k_prime > 0 {
            out.resize(k_total, b1);
            return out;
        }

        out.push(b1);
        for _ in 1..k_inner {
            out.push(b2.clone());
        }
        out
    }

    fn sample_queries_and_predicate_sparse(
        &self,
        rng: &mut dyn RngCore,
        x: &[FLarge],
    ) -> Result<(Vec<SparseVec<FLarge>>, FlpcpPredicate<FLarge>), String> {
        if x.len() != self.n() {
            return Err("bad input".to_string());
        }

        // Downcast x into FSmall assuming x is in the embedded range (typically {0,1} after Booleanization).
        let x_small = x
            .iter()
            .map(|xi| {
                let bi = BigInt::from_bytes_le(num_bigint::Sign::Plus, &xi.into_bigint().to_bytes_le());
                let p = modulus_bigint::<FSmall>();
                let r = bi % &p;
                FSmall::from_le_bytes_mod_order(&r.to_bytes_le().1)
            })
            .collect::<Vec<_>>();

        let (q_small, pred_small) = self
            .inner
            .sample_queries_and_predicate_sparse(rng, &x_small)?;
        let k_inner = q_small.len();
        let k_prime = self.params.k_prime;

        // Embed small-field sparse queries into large field.
        let mut q_large = q_small
            .iter()
            .map(|row| {
                SparseVec::new(
                    row.terms
                        .iter()
                        .map(|(c, idx)| (embed_small_to_large::<FSmall, FLarge>(c), *idx))
                        .collect::<Vec<_>>(),
                )
            })
            .collect::<Vec<_>>();

        // Modify first query: q1' = q1 + p * u' where u' = (0^n || u).
        if k_inner > 0 && self.m() > 0 {
            let p_small = modulus_bigint::<FSmall>();
            let p_as_large = FLarge::from_le_bytes_mod_order(&p_small.to_bytes_le().1);

            // Accumulate into a map to merge indices if needed.
            let mut acc: std::collections::BTreeMap<usize, FLarge> = std::collections::BTreeMap::new();
            for (c, idx) in q_large[0].terms.iter() {
                *acc.entry(*idx).or_insert(FLarge::ZERO) += *c;
            }
            for j in 0..self.m() {
                let u = sample_centered_int(rng, self.params.gamma);
                let u_large = FLarge::from_le_bytes_mod_order(&u.to_bytes_le().1);
                let idx = self.n() + j;
                *acc.entry(idx).or_insert(FLarge::ZERO) += p_as_large * u_large;
            }
            q_large[0] = SparseVec::new(
                acc.into_iter()
                    .filter_map(|(idx, c)| if c.is_zero() { None } else { Some((c, idx)) })
                    .collect::<Vec<_>>(),
            );
        }

        // Strong-soundness tail (Claim 5.28): add k' extra bound-test queries q=(0||p u).
        if k_prime > 0 {
            let p_small = modulus_bigint::<FSmall>();
            let p_as_large = FLarge::from_le_bytes_mod_order(&p_small.to_bytes_le().1);
            for _ in 0..k_prime {
                let mut terms: Vec<(FLarge, usize)> = Vec::new();
                if self.m() > 0 {
                    for j in 0..self.m() {
                        let u = sample_centered_int(rng, self.params.gamma);
                        let u_large = FLarge::from_le_bytes_mod_order(&u.to_bytes_le().1);
                        let c = p_as_large * u_large;
                        if !c.is_zero() {
                            terms.push((c, self.n() + j));
                        }
                    }
                }
                q_large.push(SparseVec::new(terms));
            }
        }

        // Predicate same as dense wrapper.
        let base_pred = match pred_small {
            FlpcpPredicate::MulEq => FlpcpPredicate::MulEqModP { p_small: modulus_bigint::<FSmall>() },
            FlpcpPredicate::MulEqAndZero => FlpcpPredicate::Unsupported(core::marker::PhantomData),
            FlpcpPredicate::AllZero => FlpcpPredicate::AllZero,
            FlpcpPredicate::MulEqModP { .. } => FlpcpPredicate::Unsupported(core::marker::PhantomData),
            FlpcpPredicate::Unsupported(_) => FlpcpPredicate::Unsupported(core::marker::PhantomData),
            FlpcpPredicate::IgnoreTail { .. } => FlpcpPredicate::Unsupported(core::marker::PhantomData),
        };
        let pred = if k_prime > 0 {
            FlpcpPredicate::IgnoreTail { inner: Box::new(base_pred), tail_len: k_prime }
        } else {
            base_pred
        };

        Ok((q_large, pred))
    }
}

fn modulus_bigint<F: PrimeField>() -> BigInt {
    BigInt::from_bytes_le(num_bigint::Sign::Plus, &F::MODULUS.to_bytes_le())
}

fn embed_small_to_large<FSmall: PrimeField, FLarge: PrimeField>(x: &FSmall) -> FLarge {
    let bi = BigInt::from_bytes_le(num_bigint::Sign::Plus, &x.into_bigint().to_bytes_le());
    FLarge::from_le_bytes_mod_order(&bi.to_bytes_le().1)
}

/// Sample a centered integer representative from F_gamma.
fn sample_centered_int(rng: &mut dyn RngCore, gamma: u64) -> BigInt {
    // Pick uniform in [0, gamma-1], then center to [-(gamma-1)/2, (gamma-1)/2].
    //
    // Important: the embedding analysis (Theorem 5.6 / Claim 5.28) assumes **uniform**
    // sampling in F_γ. Do not use `mod gamma` reduction here (bias).
    if gamma == 0 {
        return BigInt::zero();
    }
    let r = sample_uniform_u64_below(rng, gamma) as i128;
    let half = ((gamma - 1) / 2) as i128;
    let centered = if r > half { r - gamma as i128 } else { r };
    BigInt::from(centered)
}

fn sample_uniform_u64_below(rng: &mut dyn RngCore, n: u64) -> u64 {
    debug_assert!(n > 0);
    // Rejection sample from u64 to avoid modulo bias.
    //
    // Let limit be the largest multiple of n below 2^64, then accept r < limit.
    let limit = u64::MAX - (u64::MAX % n);
    loop {
        let r = rng.next_u64();
        if r < limit {
            return r % n;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::packing::{BoundedFlpcp, BoundedFlpcpSparse, DppFromBoundedFlpcp, DppFromBoundedFlpcpSparse, PackedDppParams};
    use ark_ff::{Fp64, MontBackend, MontConfig};
    use ark_ff::Field;
    use num_bigint::BigInt;
    use num_traits::ToPrimitive;
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

    impl BoundedFlpcp<FSmall> for ToyMulEq {
        fn n(&self) -> usize { 1 }
        fn m(&self) -> usize { 0 }
        fn k(&self) -> usize { 3 }
        fn bounds_b(&self) -> Vec<BigInt> { vec![BigInt::from(2), BigInt::from(2), BigInt::from(2)] }
        fn sample_queries_and_predicate(
            &self,
            _rng: &mut dyn RngCore,
            _x: &[FSmall],
        ) -> Result<(Vec<Vec<FSmall>>, FlpcpPredicate<FSmall>), String> {
            // Answers are (1,1,1) regardless of x (toy).
            Ok((vec![vec![FSmall::from(1u64)], vec![FSmall::from(1u64)], vec![FSmall::from(1u64)]], FlpcpPredicate::MulEq))
        }
    }

    #[test]
    fn test_embedding_then_packing_roundtrip_toy() {
        let inner = ToyMulEq;
        let emb = EmbeddedFlpcp::<FSmall, FLarge, _>::new(
            inner,
            EmbeddingParams { gamma: 3, assume_boolean_proof: true, k_prime: 0 },
        );

        // Choose ell per Construction 5.21: ell <= p / (2b)^k.
        let b = emb.bounds_b();
        let mut b_max = BigInt::zero();
        for bi in &b { if bi > &b_max { b_max = bi.clone(); } }
        let two_b = BigInt::from(2u64) * &b_max;
        let denom = two_b.pow(emb.k() as u32);
        let p = BigInt::from_bytes_le(num_bigint::Sign::Plus, &FLarge::MODULUS.to_bytes_le());
        let ell = (&p / denom).to_u64().unwrap_or(1).max(2);

        let dpp = DppFromBoundedFlpcp::<FLarge, _>::new(emb, PackedDppParams { ell });
        let mut rng = ChaCha20Rng::seed_from_u64(42);
        let x = vec![FLarge::from(1u64)];
        let ok = dpp.verify(&mut rng, &x, &[]).unwrap();
        assert!(ok);
    }

    #[test]
    fn test_strong_embedding_adds_tail_bound_tests_but_ignores_them_in_predicate() {
        // Underlying verifier accepts iff a0 == 0.
        #[derive(Clone, Debug)]
        struct OneQueryZero;
        impl BoundedFlpcp<FSmall> for OneQueryZero {
            fn n(&self) -> usize { 0 }
            fn m(&self) -> usize { 1 }
            fn k(&self) -> usize { 1 }
            fn bounds_b(&self) -> Vec<BigInt> { vec![BigInt::from(10)] }
            fn sample_queries_and_predicate(
                &self,
                _rng: &mut dyn RngCore,
                _x: &[FSmall],
            ) -> Result<(Vec<Vec<FSmall>>, FlpcpPredicate<FSmall>), String> {
                Ok((vec![vec![FSmall::from(1u64)]], FlpcpPredicate::AllZero))
            }
        }

        let emb = EmbeddedFlpcp::<FSmall, FLarge, _>::new(
            OneQueryZero,
            EmbeddingParams { gamma: 7, assume_boolean_proof: true, k_prime: 5 },
        );
        assert_eq!(emb.k(), 6);
        assert_eq!(emb.bounds_b().len(), 6);

        let mut rng = ChaCha20Rng::seed_from_u64(7);
        let (q, pred) = emb.sample_queries_and_predicate(&mut rng, &[]).unwrap();
        assert_eq!(q.len(), 6);

        // First answer 0 makes the base predicate pass; tail answers are ignored by predicate,
        // but would still be required to be bounded by decoding in the packed verifier.
        let answers = vec![FLarge::ZERO; 6];
        assert!(pred.check(&answers));
    }

    #[test]
    fn test_sparse_embedding_then_sparse_packing_roundtrip_toy() {
        #[derive(Clone, Debug)]
        struct ToyMulEqSparse;

        impl BoundedFlpcpSparse<FSmall> for ToyMulEqSparse {
            fn n(&self) -> usize { 1 }
            fn m(&self) -> usize { 0 }
            fn k(&self) -> usize { 3 }
            fn bounds_b(&self) -> Vec<BigInt> { vec![BigInt::from(2), BigInt::from(2), BigInt::from(2)] }

            fn sample_queries_and_predicate_sparse(
                &self,
                _rng: &mut dyn RngCore,
                _x: &[FSmall],
            ) -> Result<(Vec<SparseVec<FSmall>>, FlpcpPredicate<FSmall>), String> {
                // Same trivial toy schedule: all queries read x[0].
                Ok((
                    vec![
                        SparseVec::new(vec![(FSmall::from(1u64), 0)]),
                        SparseVec::new(vec![(FSmall::from(1u64), 0)]),
                        SparseVec::new(vec![(FSmall::from(1u64), 0)]),
                    ],
                    FlpcpPredicate::MulEq,
                ))
            }
        }

        let inner = ToyMulEqSparse;
        let emb = EmbeddedFlpcpSparse::<FSmall, FLarge, _>::new(
            inner,
            EmbeddingParams { gamma: 2, assume_boolean_proof: true, k_prime: 0 },
        );

        // Pick small ell so modulus condition holds (toy).
        let dpp = DppFromBoundedFlpcpSparse::<FLarge, _>::new(emb, PackedDppParams { ell: 2 });
        let mut rng = ChaCha20Rng::seed_from_u64(99);
        let x = vec![FLarge::from(1u64)];
        let q = dpp.sample_query(&mut rng, &x).unwrap();
        let ok = dpp.verify_with_query(&x, &[], &q).unwrap();
        assert!(ok);
    }
}

