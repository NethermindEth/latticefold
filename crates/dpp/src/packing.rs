//! Query packing transformation for bounded FLPCPs (Construction 5.21 / Theorem 5.7).
//!
//! This compiles a bounded k-query (F)LPCP into a 1-query “dot-product proof” verifier,
//! by sampling a structured weight vector `w` and packing the k queries `Q` into one query
//! `q = Q^T w`.

use ark_ff::{BigInteger, PrimeField};
use num_bigint::{BigInt, BigUint};
use num_traits::{One, Signed, ToPrimitive, Zero};
use rand::RngCore;
use std::collections::{BTreeMap, HashMap};
use thiserror::Error;

use rayon::prelude::*;

use crate::subset_sum::{decode_bounded_subset_sum, SubsetSumError};
use crate::sparse::SparseVec;

const PAR_SPARSE_PACKING_MIN_TOTAL_TERMS: usize = 4096;
const PAR_SPARSE_PACKING_MIN_K: usize = 8;
const PAR_DENSE_DOT_MIN_LEN: usize = 8192;

#[derive(Debug, Error)]
pub enum PackingError {
    #[error("invalid parameters")]
    InvalidParams,
    #[error("packed field modulus too small for bounded decoding")]
    ModulusTooSmall,
    #[error("subset-sum decode failed: {0}")]
    Decode(#[from] SubsetSumError),
}

/// Minimal interface of a bounded FLPCP verifier needed for Construction 5.21.
///
/// We model the verifier output as:
/// - `Q`: k linear query vectors (each length n+m) over `F`
/// - `accept`: a deterministic predicate on the k answers (as field elements)
pub trait BoundedFlpcp<F: PrimeField> {
    /// Input length n (field elements).
    fn n(&self) -> usize;
    /// Proof length m (field elements).
    fn m(&self) -> usize;
    /// Query count k.
    fn k(&self) -> usize;

    /// Bounds `b_i` such that honest answers satisfy |[a_i]_p| <= b_i - 1 (integer rep).
    fn bounds_b(&self) -> Vec<BigInt>;

    /// Generate k query vectors `Q` and an accept predicate for a given instance `x`.
    ///
    /// The accept predicate is returned as an enum tag + parameters (to keep this trait object-safe).
    fn sample_queries_and_predicate(
        &self,
        rng: &mut dyn RngCore,
        x: &[F],
    ) -> Result<(Vec<Vec<F>>, FlpcpPredicate<F>), String>;
}

/// Sparse-query interface for a bounded FLPCP verifier.
///
/// Same as `BoundedFlpcp`, but query vectors are represented sparsely as `(coeff, index)` lists.
pub trait BoundedFlpcpSparse<F: PrimeField> {
    /// Input length n (field elements).
    fn n(&self) -> usize;
    /// Proof length m (field elements).
    fn m(&self) -> usize;
    /// Query count k.
    fn k(&self) -> usize;
    /// Bounds `b_i` such that honest answers satisfy |[a_i]_p| <= b_i - 1 (integer rep).
    fn bounds_b(&self) -> Vec<BigInt>;

    fn sample_queries_and_predicate_sparse(
        &self,
        rng: &mut dyn RngCore,
        x: &[F],
    ) -> Result<(Vec<SparseVec<F>>, FlpcpPredicate<F>), String>;
}

/// A small set of built-in predicates for tests and for the dR1CS FLPCP.
#[derive(Clone, Debug)]
pub enum FlpcpPredicate<F: PrimeField> {
    /// Accept iff a0 * a1 == a2.
    MulEq,
    /// Accept iff (a0 * a1 - a2) == 0 mod p_small, where p_small is an integer modulus.
    ///
    /// This is used by the bounded-embedding transformation (Section 5.1), where verification
    /// is performed on answers reduced mod the original smaller field.
    MulEqModP { p_small: BigInt },
    /// Accept iff a0 * a1 == a2 and a3 == 0 (used for 4-query variants).
    MulEqAndZero,
    /// Accept iff all answers equal 0 (useful for toy tests).
    AllZero,
    /// Custom predicate (not supported in packing verifier).
    #[allow(dead_code)]
    Unsupported(core::marker::PhantomData<F>),
    /// Ignore the last `tail_len` answers and apply `inner` to the prefix.
    ///
    /// This is useful for strong-soundness embedding (Claim 5.28), where extra “bound-test”
    /// queries contribute answers that must only be bounded (enforced by decoding), with no
    /// additional predicate constraints.
    IgnoreTail { inner: Box<FlpcpPredicate<F>>, tail_len: usize },
}

impl<F: PrimeField> FlpcpPredicate<F> {
    pub fn check(&self, answers: &[F]) -> bool {
        match self {
            FlpcpPredicate::MulEq => {
                if answers.len() != 3 {
                    return false;
                }
                answers[0] * answers[1] == answers[2]
            }
            FlpcpPredicate::MulEqModP { p_small } => {
                if answers.len() != 3 {
                    return false;
                }
                let a0 = BigInt::from_bytes_le(num_bigint::Sign::Plus, &answers[0].into_bigint().to_bytes_le()) % p_small;
                let a1 = BigInt::from_bytes_le(num_bigint::Sign::Plus, &answers[1].into_bigint().to_bytes_le()) % p_small;
                let a2 = BigInt::from_bytes_le(num_bigint::Sign::Plus, &answers[2].into_bigint().to_bytes_le()) % p_small;
                ((a0 * a1) - a2) % p_small == BigInt::zero()
            }
            FlpcpPredicate::MulEqAndZero => {
                if answers.len() != 4 {
                    return false;
                }
                answers[0] * answers[1] == answers[2] && answers[3].is_zero()
            }
            FlpcpPredicate::AllZero => answers.iter().all(|a| a.is_zero()),
            FlpcpPredicate::Unsupported(_) => false,
            FlpcpPredicate::IgnoreTail { inner, tail_len } => {
                if *tail_len > answers.len() {
                    return false;
                }
                inner.check(&answers[..(answers.len() - *tail_len)])
            }
        }
    }
}

/// Parameters for the packed DPP (Construction 5.21).
#[derive(Clone, Debug)]
pub struct PackedDppParams {
    /// ℓ (denoted `ell` in paper): the per-coordinate sampling range size.
    pub ell: u64,
}

impl PackedDppParams {
    /// Compute a **deterministic safe upper bound** on ℓ that guarantees the packing
    /// “no wrap” condition (Claim 5.22 style) for **all** possible sampled weights `w`
    /// from Construction 5.21.
    ///
    /// In the construction, weights satisfy:
    /// - Let `B = max_i b_i` and `two_b = 2B`.
    /// - For i=0..k-1, sampled `w_i ∈ [ (two_b^i - 1) * ℓ + 1 , two_b^i * ℓ ]`.
    ///
    /// Therefore, for any sample, we have the worst-case bound:
    ///   <b-1, w> ≤ ℓ * Σ_i (b_i - 1) * two_b^i.
    ///
    /// To guarantee `p > 2 * <b-1, w>` for all samples, it suffices to choose:
    ///   ℓ ≤ floor( (p-1) / ( 2 * Σ_i (b_i - 1) * two_b^i ) ).
    ///
    /// Returns `None` if parameters are invalid or no positive ℓ exists.
    pub fn max_ell_for_bounds<F: PrimeField>(b: &[BigInt]) -> Option<u64> {
        if b.is_empty() {
            return None;
        }

        // B = max_i b_i
        let mut b_max = BigInt::zero();
        for bi in b {
            if bi > &b_max {
                b_max = bi.clone();
            }
        }
        if b_max <= BigInt::one() {
            return None;
        }
        let two_b = BigInt::from(2u64) * b_max;

        // denom = 2 * Σ_i (b_i - 1) * (two_b^i)
        let mut denom = BigInt::zero();
        let mut pow = BigInt::one(); // two_b^0
        for bi in b {
            let bound_i = bi - BigInt::one();
            if bound_i.is_negative() {
                return None;
            }
            denom += &bound_i * &pow;
            pow *= &two_b;
        }
        denom *= BigInt::from(2u64);
        if denom.is_zero() {
            return None;
        }

        // ell_max = floor((p-1)/denom)
        let p = modulus_bigint::<F>();
        let num = p - BigInt::one();
        let ell_max = &num / denom;
        if ell_max <= BigInt::zero() {
            return None;
        }

        // Clamp to u64 (the current `PackedDppParams` type).
        //
        // IMPORTANT:
        // If the computed safe upper bound exceeds `u64::MAX`, that does *not* mean the modulus
        // is too small — it means the safe range is “effectively unbounded” for our `u64`-typed ℓ.
        // In that case we should saturate, not error (otherwise callers see `ModulusTooSmall`).
        ell_max.to_u64().or(Some(u64::MAX))
    }

    /// Construct params picking the maximum safe ℓ (clamped to u64).
    pub fn from_bounds_max<F: PrimeField>(b: &[BigInt]) -> Result<Self, PackingError> {
        let ell = Self::max_ell_for_bounds::<F>(b).ok_or(PackingError::ModulusTooSmall)?;
        Ok(Self { ell })
    }
}

/// Sample packing weights `w` as in Construction 5.21, and enforce the Claim 5.22 “no wrap”
/// condition `p > 2 * <b-1, w>`.
///
/// This is the heavy-hitter metadata needed for decoding, but does **not** require building the
/// packed query vector `q = Q^T w`.
pub fn sample_packing_weights<F: PrimeField>(
    rng: &mut dyn RngCore,
    ell: u64,
    b: &[BigInt],
) -> Result<Vec<BigInt>, PackingError> {
    let k = b.len();
    if k == 0 {
        return Err(PackingError::InvalidParams);
    }
    let mut b_max = BigInt::zero();
    for bi in b {
        if bi > &b_max {
            b_max = bi.clone();
        }
    }
    if b_max <= BigInt::one() {
        return Err(PackingError::InvalidParams);
    }

    let ell = BigInt::from(ell);
    if ell.is_zero() {
        return Err(PackingError::InvalidParams);
    }
    let two_b = BigInt::from(2u64) * &b_max;
    let mut w: Vec<BigInt> = Vec::with_capacity(k);
    let mut pow = BigInt::one();
    for _i in 0..k {
        let low = (&pow - BigInt::one()) * &ell + BigInt::one();
        let high = &pow * &ell;
        let span = (&high - &low) + BigInt::one();
        let span_u = big_int_to_big_uint(&span).ok_or(PackingError::InvalidParams)?;
        let r_u = sample_uniform_below(rng, &span_u);
        let samp = &low + BigInt::from(r_u);
        w.push(samp);
        pow *= &two_b;
    }

    let p_mod = modulus_bigint::<F>();
    let mut bw = BigInt::zero();
    for i in 0..k {
        let bound_i = &b[i] - BigInt::one();
        if bound_i.is_negative() {
            return Err(PackingError::InvalidParams);
        }
        bw += bound_i * &w[i];
    }
    if p_mod <= BigInt::from(2u64) * bw {
        return Err(PackingError::ModulusTooSmall);
    }

    Ok(w)
}

/// A packed DPP verifier built from a bounded FLPCP via Construction 5.21.
#[derive(Clone, Debug)]
pub struct DppFromBoundedFlpcp<F: PrimeField, V: BoundedFlpcp<F>> {
    pub flpcp: V,
    pub params: PackedDppParams,
    _pd: core::marker::PhantomData<F>,
}

/// A sampled packed DPP query (Construction 5.21), plus decoding metadata.
///
/// This is the artifact a WE “armer” wants: a single linear query vector `q` over the
/// concatenated witness vector `v=(x||π)` such that the verifier accepts based on the
/// packed answer `a = <q, v>`.
#[derive(Clone, Debug)]
pub struct PackedDppQuery<F: PrimeField> {
    /// Packed query vector (length n+m).
    pub q: Vec<F>,
    /// Packing weights `w` (integers, length k).
    pub w: Vec<BigInt>,
    /// Bounds `b` (integers, length k).
    pub b: Vec<BigInt>,
    /// Underlying k-answer predicate (on decoded answers).
    pub pred: FlpcpPredicate<F>,
}

/// A sampled packed DPP query in sparse form.
#[derive(Clone, Debug)]
pub struct PackedDppQuerySparse<F: PrimeField> {
    /// Packed query vector (sparse) over v=(x||π), indices in [0..n+m).
    pub q: SparseVec<F>,
    /// Packing weights `w` (integers, length k).
    pub w: Vec<BigInt>,
    /// Bounds `b` (integers, length k).
    pub b: Vec<BigInt>,
    /// Underlying k-answer predicate (on decoded answers).
    pub pred: FlpcpPredicate<F>,
}

impl<F: PrimeField, V: BoundedFlpcp<F>> DppFromBoundedFlpcp<F, V> {
    pub fn new(flpcp: V, params: PackedDppParams) -> Self {
        Self { flpcp, params, _pd: core::marker::PhantomData }
    }

    /// Sample the packed query vector `q` (length n+m) and the packing weights `w`.
    pub fn sample_packed_query(
        &self,
        rng: &mut dyn RngCore,
        x: &[F],
    ) -> Result<(Vec<F>, Vec<BigInt>, FlpcpPredicate<F>, Vec<Vec<F>>), PackingError> {
        let k = self.flpcp.k();
        let b = self.flpcp.bounds_b();
        if k == 0 || b.len() != k || x.len() != self.flpcp.n() {
            return Err(PackingError::InvalidParams);
        }
        let (q_mat, pred) = self
            .flpcp
            .sample_queries_and_predicate(rng, x)
            .map_err(|_| PackingError::InvalidParams)?;
        if q_mat.len() != k || q_mat.iter().any(|row| row.len() != self.flpcp.n() + self.flpcp.m()) {
            return Err(PackingError::InvalidParams);
        }

        // Choose scalar b_max = max_i b_i (as BigInt).
        let mut b_max = BigInt::zero();
        for bi in &b {
            if bi > &b_max {
                b_max = bi.clone();
            }
        }
        if b_max <= BigInt::one() {
            return Err(PackingError::InvalidParams);
        }

        // Sample w as in Construction 5.21 (integer ranges).
        let ell = BigInt::from(self.params.ell);
        if ell.is_zero() {
            return Err(PackingError::InvalidParams);
        }

        let two_b = BigInt::from(2u64) * &b_max;
        let mut w: Vec<BigInt> = Vec::with_capacity(k);
        let mut pow = BigInt::one(); // (2b)^(i-1)
        for _i in 0..k {
            // low = ((pow - 1) * ell + 1), high = pow * ell
            let low = (&pow - BigInt::one()) * &ell + BigInt::one();
            let high = &pow * &ell;
            // Sample uniformly in [low, high] (inclusive).
            let span = (&high - &low) + BigInt::one();
            let span_u = big_int_to_big_uint(&span).ok_or(PackingError::InvalidParams)?;
            let r_u = sample_uniform_below(rng, &span_u);
            let samp = &low + BigInt::from(r_u);
            w.push(samp);
            pow *= &two_b;
        }

        // Ensure the field modulus is large enough so that bounded answers do not wrap mod p.
        // This is the condition in Claim 5.22: p > 2 * <b-1, w>.
        let p_mod = modulus_bigint::<F>();
        let mut bw = BigInt::zero();
        for i in 0..k {
            let bound_i = &b[i] - BigInt::one();
            if bound_i.is_negative() {
                return Err(PackingError::InvalidParams);
            }
            bw += bound_i * &w[i];
        }
        if p_mod <= BigInt::from(2u64) * bw {
            return Err(PackingError::ModulusTooSmall);
        }

        // Compute packed query coefficients: q = Q^T w over integers, then reduce mod p.
        let len = self.flpcp.n() + self.flpcp.m();
        // This is typically the dominant cost for large instances (len = n+m can be huge).
        // Parallelize over output coordinates (each j is independent).
        let p_mod = modulus_bigint::<F>();
        let half = (&p_mod - BigInt::one()) / BigInt::from(2u64);
        let field_to_centered_bigint_fast = |x: &F| -> BigInt {
            let mut t = BigInt::from_bytes_le(num_bigint::Sign::Plus, &x.into_bigint().to_bytes_le());
            if t > half {
                t -= &p_mod;
            }
            t
        };

        let q: Vec<F> = (0..len)
            .into_par_iter()
            .map(|j| {
                let mut acc = BigInt::zero();
                for i in 0..k {
                    let qij = field_to_centered_bigint_fast(&q_mat[i][j]);
                    acc += &w[i] * qij;
                }
                centered_bigint_to_field::<F>(&acc)
            })
            .collect();

        Ok((q, w, pred, q_mat))
    }

    /// Sample a reusable packed DPP query artifact for this instance `x`.
    pub fn sample_query(&self, rng: &mut dyn RngCore, x: &[F]) -> Result<PackedDppQuery<F>, PackingError> {
        let b = self.flpcp.bounds_b();
        let (q, w, pred, _q_mat) = self.sample_packed_query(rng, x)?;
        Ok(PackedDppQuery { q, w, b, pred })
    }

    /// Verify a packed answer `a = <q, (x||π)>` using the sampled query metadata.
    ///
    /// This is the verifier’s “accepting set” check in algorithmic form:
    /// decode `a` into k bounded answers, then apply the underlying predicate.
    pub fn verify_packed_answer(&self, a: &F, query: &PackedDppQuery<F>) -> Result<bool, PackingError> {
        let a_int = field_to_centered_bigint::<F>(a);
        let ans_int = decode_bounded_subset_sum(&a_int, &query.w, &query.b)?;
        let ans_field = ans_int
            .iter()
            .map(|z| centered_bigint_to_field::<F>(z))
            .collect::<Vec<_>>();
        Ok(query.pred.check(&ans_field))
    }

    /// Verify using a caller-supplied packed query (no resampling).
    pub fn verify_with_query(
        &self,
        x: &[F],
        pi: &[F],
        query: &PackedDppQuery<F>,
    ) -> Result<bool, PackingError> {
        if x.len() != self.flpcp.n() || pi.len() != self.flpcp.m() {
            return Err(PackingError::InvalidParams);
        }
        if query.q.len() != self.flpcp.n() + self.flpcp.m() {
            return Err(PackingError::InvalidParams);
        }
        // Avoid allocating/copying a huge `(x||pi)` vector; compute the dot directly.
        let a = dot_two_slices::<F>(&query.q, x, pi);
        self.verify_packed_answer(&a, query)
    }

    /// Verify by computing the packed dot product and decoding it back into k bounded answers.
    ///
    /// This corresponds to running `V'` from Construction 5.21 on the packed answer.
    pub fn verify(
        &self,
        rng: &mut dyn RngCore,
        x: &[F],
        pi: &[F],
    ) -> Result<bool, PackingError> {
        if x.len() != self.flpcp.n() || pi.len() != self.flpcp.m() {
            return Err(PackingError::InvalidParams);
        }
        let (q, w, pred, _q_mat) = self.sample_packed_query(rng, x)?;
        // Avoid allocating/copying a huge `(x||pi)` vector; compute the dot directly.
        let a = dot_two_slices::<F>(&q, x, pi);
        let a_int = field_to_centered_bigint::<F>(&a);

        // Decode a_int as <ans, w> with bounds b.
        let b = self.flpcp.bounds_b();
        let ans_int = decode_bounded_subset_sum(&a_int, &w, &b)?;
        let ans_field = ans_int
            .iter()
            .map(|z| centered_bigint_to_field::<F>(z))
            .collect::<Vec<_>>();

        Ok(pred.check(&ans_field))
    }
}

/// Sparse packed DPP verifier built from a sparse bounded FLPCP.
#[derive(Clone, Debug)]
pub struct DppFromBoundedFlpcpSparse<F: PrimeField, V: BoundedFlpcpSparse<F>> {
    pub flpcp: V,
    pub params: PackedDppParams,
    _pd: core::marker::PhantomData<F>,
}

impl<F: PrimeField, V: BoundedFlpcpSparse<F>> DppFromBoundedFlpcpSparse<F, V> {
    pub fn new(flpcp: V, params: PackedDppParams) -> Self {
        Self { flpcp, params, _pd: core::marker::PhantomData }
    }

    pub fn sample_query(
        &self,
        rng: &mut dyn RngCore,
        x: &[F],
    ) -> Result<PackedDppQuerySparse<F>, PackingError> {
        let k = self.flpcp.k();
        let b = self.flpcp.bounds_b();
        if k == 0 || b.len() != k || x.len() != self.flpcp.n() {
            return Err(PackingError::InvalidParams);
        }
        let (q_mat, pred) = self
            .flpcp
            .sample_queries_and_predicate_sparse(rng, x)
            .map_err(|_| PackingError::InvalidParams)?;
        if q_mat.len() != k {
            return Err(PackingError::InvalidParams);
        }

        // Choose scalar b_max = max_i b_i (as BigInt).
        let mut b_max = BigInt::zero();
        for bi in &b {
            if bi > &b_max {
                b_max = bi.clone();
            }
        }
        if b_max <= BigInt::one() {
            return Err(PackingError::InvalidParams);
        }

        // Sample w as in Construction 5.21 (integer ranges).
        let ell = BigInt::from(self.params.ell);
        if ell.is_zero() {
            return Err(PackingError::InvalidParams);
        }
        let two_b = BigInt::from(2u64) * &b_max;
        let mut w: Vec<BigInt> = Vec::with_capacity(k);
        let mut pow = BigInt::one(); // (2b)^(i-1)
        for _i in 0..k {
            let low = (&pow - BigInt::one()) * &ell + BigInt::one();
            let high = &pow * &ell;
            let span = (&high - &low) + BigInt::one();
            let span_u = big_int_to_big_uint(&span).ok_or(PackingError::InvalidParams)?;
            let r_u = sample_uniform_below(rng, &span_u);
            let samp = &low + BigInt::from(r_u);
            w.push(samp);
            pow *= &two_b;
        }

        // Ensure the field modulus is large enough so that bounded answers do not wrap mod p.
        let p_mod = modulus_bigint::<F>();
        let mut bw = BigInt::zero();
        for i in 0..k {
            let bound_i = &b[i] - BigInt::one();
            if bound_i.is_negative() {
                return Err(PackingError::InvalidParams);
            }
            bw += bound_i * &w[i];
        }
        if p_mod <= BigInt::from(2u64) * bw {
            return Err(PackingError::ModulusTooSmall);
        }

        // Compute packed query coefficients q = Σ_i w_i * q_i over integers, then reduce mod p.
        let p_mod = modulus_bigint::<F>();
        let half = (&p_mod - BigInt::one()) / BigInt::from(2u64);
        let field_to_centered_bigint_fast = |x: &F| -> BigInt {
            let mut t = BigInt::from_bytes_le(num_bigint::Sign::Plus, &x.into_bigint().to_bytes_le());
            if t > half {
                t -= &p_mod;
            }
            t
        };

        let total_terms: usize = q_mat.iter().map(|row| row.terms.len()).sum();
        let use_par = k >= PAR_SPARSE_PACKING_MIN_K && total_terms >= PAR_SPARSE_PACKING_MIN_TOTAL_TERMS;

        let mut acc_items: Vec<(usize, BigInt)> = if use_par {
            let merged: HashMap<usize, BigInt> = q_mat
                .par_iter()
                .enumerate()
                .fold(
                    || HashMap::<usize, BigInt>::new(),
                    |mut local, (i, row)| {
                        let wi = &w[i];
                        for (coeff, idx) in &row.terms {
                            let cij = field_to_centered_bigint_fast(coeff);
                            let entry = local.entry(*idx).or_insert_with(BigInt::zero);
                            *entry += wi * cij;
                        }
                        local
                    },
                )
                .reduce(
                    || HashMap::<usize, BigInt>::new(),
                    |mut a, b| {
                        for (idx, z) in b {
                            let entry = a.entry(idx).or_insert_with(BigInt::zero);
                            *entry += z;
                        }
                        a
                    },
                );
            merged.into_iter().collect()
        } else {
            // Small instances: avoid Rayon/HashMap overhead.
            let mut acc: BTreeMap<usize, BigInt> = BTreeMap::new();
            for i in 0..k {
                let wi = &w[i];
                for (coeff, idx) in &q_mat[i].terms {
                    let cij = field_to_centered_bigint_fast(coeff);
                    let entry = acc.entry(*idx).or_insert_with(BigInt::zero);
                    *entry += wi * cij;
                }
            }
            acc.into_iter().collect()
        };

        // Make output stable/deterministic regardless of parallel merging order.
        acc_items.sort_by_key(|(idx, _)| *idx);
        let mut terms: Vec<(F, usize)> = Vec::with_capacity(acc_items.len());
        for (idx, z) in acc_items {
            let f = centered_bigint_to_field::<F>(&z);
            if !f.is_zero() {
                terms.push((f, idx));
            }
        }

        Ok(PackedDppQuerySparse { q: SparseVec::new(terms), w, b, pred })
    }

    /// Sample packing metadata (q_i queries, weights w, bounds b, predicate), but **do not**
    /// build the packed query vector `q = Q^T w`.
    ///
    /// This is useful for benchmarking/decapper-side evaluation where we want the packed answer
    /// `a = <q, v>` but do not want to allocate a gigantic packed query vector.
    pub fn sample_query_components(
        &self,
        rng: &mut dyn RngCore,
        x: &[F],
    ) -> Result<(Vec<SparseVec<F>>, Vec<BigInt>, Vec<BigInt>, FlpcpPredicate<F>), PackingError> {
        let k = self.flpcp.k();
        let b = self.flpcp.bounds_b();
        if k == 0 || b.len() != k || x.len() != self.flpcp.n() {
            return Err(PackingError::InvalidParams);
        }
        let (q_mat, pred) = self
            .flpcp
            .sample_queries_and_predicate_sparse(rng, x)
            .map_err(|_| PackingError::InvalidParams)?;
        if q_mat.len() != k {
            return Err(PackingError::InvalidParams);
        }

        // Choose scalar b_max = max_i b_i (as BigInt).
        let mut b_max = BigInt::zero();
        for bi in &b {
            if bi > &b_max {
                b_max = bi.clone();
            }
        }
        if b_max <= BigInt::one() {
            return Err(PackingError::InvalidParams);
        }

        // Sample w as in Construction 5.21 (integer ranges).
        let ell = BigInt::from(self.params.ell);
        if ell.is_zero() {
            return Err(PackingError::InvalidParams);
        }
        let two_b = BigInt::from(2u64) * &b_max;
        let mut w: Vec<BigInt> = Vec::with_capacity(k);
        let mut pow = BigInt::one(); // (2b)^(i-1)
        for _i in 0..k {
            let low = (&pow - BigInt::one()) * &ell + BigInt::one();
            let high = &pow * &ell;
            let span = (&high - &low) + BigInt::one();
            let span_u = big_int_to_big_uint(&span).ok_or(PackingError::InvalidParams)?;
            let r_u = sample_uniform_below(rng, &span_u);
            let samp = &low + BigInt::from(r_u);
            w.push(samp);
            pow *= &two_b;
        }

        // Ensure the field modulus is large enough so that bounded answers do not wrap mod p.
        let p_mod = modulus_bigint::<F>();
        let mut bw = BigInt::zero();
        for i in 0..k {
            let bound_i = &b[i] - BigInt::one();
            if bound_i.is_negative() {
                return Err(PackingError::InvalidParams);
            }
            bw += bound_i * &w[i];
        }
        if p_mod <= BigInt::from(2u64) * bw {
            return Err(PackingError::ModulusTooSmall);
        }

        Ok((q_mat, w, b, pred))
    }

    pub fn verify_packed_answer(
        &self,
        a: &F,
        query: &PackedDppQuerySparse<F>,
    ) -> Result<bool, PackingError> {
        let a_int = field_to_centered_bigint::<F>(a);
        let ans_int = decode_bounded_subset_sum(&a_int, &query.w, &query.b)?;
        let ans_field = ans_int
            .iter()
            .map(|z| centered_bigint_to_field::<F>(z))
            .collect::<Vec<_>>();
        Ok(query.pred.check(&ans_field))
    }

    pub fn verify_with_query(
        &self,
        x: &[F],
        pi: &[F],
        query: &PackedDppQuerySparse<F>,
    ) -> Result<bool, PackingError> {
        if x.len() != self.flpcp.n() || pi.len() != self.flpcp.m() {
            return Err(PackingError::InvalidParams);
        }
        // Avoid allocating/copying a huge `(x||pi)` vector; compute the sparse dot directly.
        let a = query.q.dot_two_slices(x, pi);
        self.verify_packed_answer(&a, query)
    }
}

impl<F: PrimeField> SparseVec<F> {
    /// Dot-product against a concatenation `(x || pi)` without allocating/copying it.
    #[inline]
    pub fn dot_two_slices(&self, x: &[F], pi: &[F]) -> F {
        const PAR_SPARSE_DOT_MIN_TERMS: usize = 2048;
        let n = x.len();
        let m = pi.len();
        if self.terms.iter().any(|(_, idx)| *idx >= n + m) {
            // Keep this as a debug assertion to avoid overhead in optimized builds.
            debug_assert!(false, "SparseVec::dot_two_slices: index out of range");
        }
        // Fast path for Boolean witnesses: if v ∈ {0,1}, then c*v is either 0 or c.
        // This is a big win for large `FLarge` (e.g. 384-bit) fields.
        if self.terms.len() >= PAR_SPARSE_DOT_MIN_TERMS {
            self.terms
                .par_iter()
                .map(|(c, idx)| {
                    let v = if *idx < n { x[*idx] } else { pi[*idx - n] };
                    if v.is_zero() {
                        F::ZERO
                    } else if v == F::ONE {
                        *c
                    } else {
                        *c * v
                    }
                })
                .reduce(|| F::ZERO, |acc, t| acc + t)
        } else {
            self.terms.iter().fold(F::ZERO, |acc, (c, idx)| {
                let v = if *idx < n { x[*idx] } else { pi[*idx - n] };
                if v.is_zero() {
                    acc
                } else if v == F::ONE {
                    acc + *c
                } else {
                    acc + (*c * v)
                }
            })
        }
    }
}

fn dot_two_slices<F: PrimeField>(q: &[F], x: &[F], pi: &[F]) -> F {
    debug_assert_eq!(q.len(), x.len() + pi.len());
    let n = x.len();
    if q.len() >= PAR_DENSE_DOT_MIN_LEN {
        q.par_iter()
            .enumerate()
            .map(|(i, qi)| {
                let vi = if i < n { x[i] } else { pi[i - n] };
                *qi * vi
            })
            .reduce(|| F::ZERO, |acc, t| acc + t)
    } else {
        q.iter().enumerate().fold(F::ZERO, |acc, (i, qi)| {
            let vi = if i < n { x[i] } else { pi[i - n] };
            acc + (*qi * vi)
        })
    }
}

fn big_int_to_big_uint(x: &BigInt) -> Option<BigUint> {
    if x.is_negative() {
        None
    } else {
        Some(BigUint::from_bytes_le(&x.to_bytes_le().1))
    }
}

/// Unbiased sampling of a uniform integer in [0, n).
fn sample_uniform_below(rng: &mut dyn RngCore, n: &BigUint) -> BigUint {
    if n.is_zero() {
        return BigUint::ZERO;
    }
    let bits = n.bits();
    let nbytes = ((bits + 7) / 8) as usize;
    let max = BigUint::one() << (nbytes * 8);
    let limit = &max - (&max % n); // largest multiple of n below 2^(8*nbytes)

    loop {
        let mut buf = vec![0u8; nbytes];
        rng.fill_bytes(&mut buf);
        let r = BigUint::from_bytes_le(&buf);
        if r < limit {
            return r % n;
        }
    }
}

// NOTE: keep the old dense dot helper around for tests/bench experiments.
// (Currently unused in production paths; we use `dot_two_slices` to avoid allocating `(x||pi)`.)
#[allow(dead_code)]
fn dot<F: PrimeField>(a: &[F], b: &[F]) -> F {
    debug_assert_eq!(a.len(), b.len());
    if a.len() >= PAR_DENSE_DOT_MIN_LEN {
        a.par_iter()
            .zip(b.par_iter())
            .map(|(x, y)| *x * *y)
            .reduce(|| F::ZERO, |acc, t| acc + t)
    } else {
        a.iter().zip(b.iter()).fold(F::ZERO, |acc, (x, y)| acc + (*x * *y))
    }
}

/// Map a field element to its centered integer representative in [-(p-1)/2, (p-1)/2].
pub fn field_to_centered_bigint<F: PrimeField>(x: &F) -> BigInt {
    let p = modulus_bigint::<F>();
    let half = (&p - BigInt::one()) / BigInt::from(2u64);
    let mut t = BigInt::from_bytes_le(num_bigint::Sign::Plus, &x.into_bigint().to_bytes_le());
    // If t > half, interpret as negative: t - p
    if t > half {
        t -= p;
    }
    t
}

/// Reduce an integer (centered or not) into the field.
pub fn centered_bigint_to_field<F: PrimeField>(z: &BigInt) -> F {
    // Convert via modular reduction into [0,p).
    let p = modulus_bigint::<F>();
    let mut t = z.mod_floor(&p);
    if t.is_negative() {
        t += &p;
    }
    let bytes = t.to_bytes_le().1;
    F::from_le_bytes_mod_order(&bytes)
}

fn modulus_bigint<F: PrimeField>() -> BigInt {
    BigInt::from_bytes_le(num_bigint::Sign::Plus, &F::MODULUS.to_bytes_le())
}

trait ModFloor {
    fn mod_floor(&self, m: &BigInt) -> BigInt;
}

impl ModFloor for BigInt {
    fn mod_floor(&self, m: &BigInt) -> BigInt {
        let r = self % m;
        if r.is_negative() { r + m } else { r }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::Field;
    use ark_ff::{Fp64, MontBackend, MontConfig};
    use rand_chacha::ChaCha20Rng;
    use rand::SeedableRng;
    use num_traits::ToPrimitive;

    // Small prime field for tests.
    #[derive(MontConfig)]
    #[modulus = "10007"]
    #[generator = "5"]
    pub struct F97Config;
    type F = Fp64<MontBackend<F97Config, 1>>;

    #[derive(Clone, Debug)]
    struct Toy3Query;

    impl BoundedFlpcp<F> for Toy3Query {
        fn n(&self) -> usize { 1 }
        fn m(&self) -> usize { 0 }
        fn k(&self) -> usize { 3 }
        fn bounds_b(&self) -> Vec<BigInt> { vec![BigInt::from(2), BigInt::from(2), BigInt::from(2)] } // |a_i|<=1
        fn sample_queries_and_predicate(
            &self,
            _rng: &mut dyn RngCore,
            _x: &[F],
        ) -> Result<(Vec<Vec<F>>, FlpcpPredicate<F>), String> {
            // Q answers are simply a0=x, a1=1, a2=x
            Ok((
                vec![vec![F::one()], vec![F::one()], vec![F::one()]],
                FlpcpPredicate::MulEq,
            ))
        }
    }

    #[test]
    fn test_packing_roundtrip_toy() {
        let flpcp = Toy3Query;
        let dpp = DppFromBoundedFlpcp::<F, _>::new(flpcp, PackedDppParams { ell: 32 });
        let mut rng = ChaCha20Rng::seed_from_u64(1);
        let x = vec![F::from(1u64)];
        let ok = dpp.verify(&mut rng, &x, &[]).unwrap();
        assert!(ok);
    }

    #[test]
    fn test_packing_query_object_roundtrip_toy() {
        let flpcp = Toy3Query;
        let dpp = DppFromBoundedFlpcp::<F, _>::new(flpcp, PackedDppParams { ell: 32 });
        let mut rng = ChaCha20Rng::seed_from_u64(2);
        let x = vec![F::from(1u64)];
        let q = dpp.sample_query(&mut rng, &x).unwrap();
        let ok = dpp.verify_with_query(&x, &[], &q).unwrap();
        assert!(ok);
    }

    #[derive(Clone, Debug)]
    struct Toy3QuerySparse;

    impl BoundedFlpcpSparse<F> for Toy3QuerySparse {
        fn n(&self) -> usize { 1 }
        fn m(&self) -> usize { 0 }
        fn k(&self) -> usize { 3 }
        fn bounds_b(&self) -> Vec<BigInt> { vec![BigInt::from(2), BigInt::from(2), BigInt::from(2)] } // |a_i|<=1

        fn sample_queries_and_predicate_sparse(
            &self,
            _rng: &mut dyn RngCore,
            _x: &[F],
        ) -> Result<(Vec<SparseVec<F>>, FlpcpPredicate<F>), String> {
            Ok((
                vec![
                    SparseVec::new(vec![(F::one(), 0)]),
                    SparseVec::new(vec![(F::one(), 0)]),
                    SparseVec::new(vec![(F::one(), 0)]),
                ],
                FlpcpPredicate::MulEq,
            ))
        }
    }

    #[test]
    fn test_packing_sparse_roundtrip_toy() {
        let flpcp = Toy3QuerySparse;
        let dpp = DppFromBoundedFlpcpSparse::<F, _>::new(flpcp, PackedDppParams { ell: 32 });
        let mut rng = ChaCha20Rng::seed_from_u64(3);
        let x = vec![F::from(1u64)];
        let q = dpp.sample_query(&mut rng, &x).unwrap();
        let ok = dpp.verify_with_query(&x, &[], &q).unwrap();
        assert!(ok);
    }

    // Large prime field to keep the modulus condition (Claim 5.22) satisfied even for k>=8.
    #[derive(MontConfig)]
    #[modulus = "18446744069414584321"]
    #[generator = "7"]
    pub struct FBigConfig;
    type FBig = Fp64<MontBackend<FBigConfig, 1>>;

    #[derive(Clone, Debug)]
    struct BigSparseFlpcp {
        n: usize,
        k: usize,
    }

    impl BoundedFlpcpSparse<FBig> for BigSparseFlpcp {
        fn n(&self) -> usize {
            self.n
        }
        fn m(&self) -> usize {
            0
        }
        fn k(&self) -> usize {
            self.k
        }
        fn bounds_b(&self) -> Vec<BigInt> {
            // Tight bounds to keep packing weights small.
            vec![BigInt::from(2u64); self.k]
        }

        fn sample_queries_and_predicate_sparse(
            &self,
            _rng: &mut dyn RngCore,
            _x: &[FBig],
        ) -> Result<(Vec<SparseVec<FBig>>, FlpcpPredicate<FBig>), String> {
            // Deterministic, dense-ish sparse queries:
            // every query row includes every index with coefficient 1.
            let mut rows = Vec::with_capacity(self.k);
            let terms = (0..self.n).map(|idx| (FBig::ONE, idx)).collect::<Vec<_>>();
            for _ in 0..self.k {
                rows.push(SparseVec::new(terms.clone()));
            }
            Ok((rows, FlpcpPredicate::AllZero))
        }
    }

    #[test]
    fn test_sparse_packed_query_matches_manual_reference_on_large_terms() {
        // Force the parallel sparse accumulation path by using many terms.
        let flpcp = BigSparseFlpcp { n: 8192, k: 8 };
        let dpp = DppFromBoundedFlpcpSparse::<FBig, _>::new(flpcp.clone(), PackedDppParams { ell: 2 });
        let mut rng = ChaCha20Rng::seed_from_u64(4242);

        let x: Vec<FBig> = vec![FBig::ONE; flpcp.n()];
        let query = dpp.sample_query(&mut rng, &x).unwrap();

        // Recompute packed q sequentially from the (deterministic) q_mat and sampled w.
        let (q_mat, _pred) = flpcp
            .sample_queries_and_predicate_sparse(&mut rng, &x)
            .expect("deterministic q_mat");
        let k = flpcp.k();

        let mut acc: BTreeMap<usize, BigInt> = BTreeMap::new();
        for i in 0..k {
            let wi = &query.w[i];
            for (coeff, idx) in &q_mat[i].terms {
                let cij = field_to_centered_bigint::<FBig>(coeff);
                let entry = acc.entry(*idx).or_insert_with(BigInt::zero);
                *entry += wi * cij;
            }
        }
        let mut exp_terms: Vec<(FBig, usize)> = acc
            .into_iter()
            .filter_map(|(idx, z)| {
                let f = centered_bigint_to_field::<FBig>(&z);
                if f.is_zero() {
                    None
                } else {
                    Some((f, idx))
                }
            })
            .collect();
        exp_terms.sort_by_key(|(_, idx)| *idx);

        let mut got_terms = query.q.terms.clone();
        got_terms.sort_by_key(|(_, idx)| *idx);

        assert_eq!(got_terms, exp_terms);
    }

    #[test]
    fn test_max_ell_for_bounds_is_sound_for_worst_case_weights() {
        // Use a large 64-bit-ish prime field so ell_max is nontrivial.
        type F = FBig;

        // Example bounds vector (k = 8), with small b to keep ell reasonably large.
        let k = 8usize;
        let b = vec![BigInt::from(3u64); k]; // |ans_i| <= 2

        let ell = PackedDppParams::max_ell_for_bounds::<F>(&b).expect("ell exists");
        assert!(ell >= 1);

        // Verify the sufficient condition derived in max_ell_for_bounds:
        // p > 2 * ell * Σ_i (b_i-1) * (2B)^i
        let p = modulus_bigint::<F>();
        let b_max = BigInt::from(3u64);
        let two_b = BigInt::from(2u64) * b_max;
        let mut sum = BigInt::zero();
        let mut pow = BigInt::one();
        for _ in 0..k {
            sum += BigInt::from(2u64) * &pow; // (b_i-1)=2
            pow *= &two_b;
        }
        let rhs = BigInt::from(2u64) * BigInt::from(ell) * sum;
        assert!(
            p > rhs,
            "p={} should be > rhs={}",
            p.to_u128().unwrap_or(0),
            rhs.to_u128().unwrap_or(0)
        );
    }
}

