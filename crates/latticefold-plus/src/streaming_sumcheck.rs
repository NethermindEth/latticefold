//! Streaming sumcheck prover for LF+ (memory-friendly).
//!
//! Produces the *same* `latticefold::utils::sumcheck::Proof<R>` format as the dense prover,
//! so the existing verifier (`MLSumcheck::verify_as_subprotocol`) remains unchanged.

use std::sync::Arc;

use latticefold::transcript::Transcript;
use latticefold::utils::sumcheck::prover::ProverMsg;
use latticefold::utils::sumcheck::Proof;
use stark_rings::{OverField, PolyRing, Ring};
use stark_rings_linalg::{Matrix, SparseMatrix};
use crate::setchk::DigitsMatrix;
use crate::utils::maybe_print_rss;
use core::mem::MaybeUninit;

// A small per-thread direct-mapped cache for streaming-h evaluations `h[idx]` when
// `h` is represented as `HFromMfDigitsConstCol0`. This targets repeated column indices
// across sparse mat-vec scans without materializing full `h`.
//
// NOTE: This cache is performance-critical when streaming-h is active and `h` is dense in the ring.
// The access pattern can have many conflict misses with a naive `idx & (size-1)` mapping, so we:
// - use a somewhat larger table, and
// - apply a multiplicative hash before masking.
//
// IMPORTANT (d-scaling):
// `vals` stores full ring elements, so its footprint scales linearly with `R::dimension()`.
// We size it to target ~4 MiB/thread for the *values* across different rings (e.g. Frog d=16 vs d=64)
// to avoid a 4x memory jump when switching to d=64.
const CM_HFROM_TARGET_VAL_BYTES: usize = 4 * 1024 * 1024;

#[inline]
fn cm_hfrom_cache_size_for<R: OverField + PolyRing>() -> usize {
    // Size by the stored ring element type (compile-time constant).
    let elem = core::mem::size_of::<R>().max(1);
    // Clamp to keep the table reasonably sized even for very small/large rings.
    let raw = (CM_HFROM_TARGET_VAL_BYTES / elem).clamp(1024, 1 << 15);
    raw.next_power_of_two()
}

pub(crate) struct HFromIndexCache<R: OverField + PolyRing>
where
    R::BaseRing: Ring,
{
    id: usize,          // identifies the underlying precomps Arc
    mask: usize,        // size - 1
    keys: Vec<usize>,   // cached index
    vals: Vec<R>,       // cached h[idx]
}

impl<R: OverField + PolyRing> HFromIndexCache<R>
where
    R::BaseRing: Ring,
{
    #[inline]
    fn new(id: usize) -> Self {
        let size = cm_hfrom_cache_size_for::<R>();
        Self {
            id,
            mask: size - 1,
            keys: vec![usize::MAX; size],
            vals: vec![R::ZERO; size],
        }
    }

    #[inline]
    fn get_or_compute(
        &mut self,
        idx: usize,
        compute: impl FnOnce() -> R,
    ) -> R {
        // Multiplicative hashing to reduce conflict misses for structured index patterns.
        let slot = idx.wrapping_mul(11400714819323198485usize) & self.mask;
        if self.keys[slot] == idx {
            return self.vals[slot];
        }
        let v = compute();
        self.keys[slot] = idx;
        self.vals[slot] = v;
        v
    }
}

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// Precomputed pieces for `h[row]` when `M_f` matrices are `DigitsBacking::ConstCol0`.
#[derive(Clone)]
pub struct HCol0Precomp<Rr: PolyRing> {
    pub col0: Arc<Vec<u16>>,
    pub zero_idx: u16,
    pub term0_tab: Arc<Vec<Rr>>,
    pub term_rest: Rr,
}

/// Grouped lookup for fast streamed-`h` evaluation.
///
/// For a group of `g` digit streams (each digit in `[0, base)`), we precompute a table
/// of length `base^g` containing the sum of the corresponding `term0_tab[d]` ring elements.
#[derive(Clone)]
pub struct HFromGroup<Rr: PolyRing> {
    pub base: usize,
    /// Minimum length across `cols` (used to skip per-column bounds checks on the hot path).
    pub len: usize,
    pub cols: Arc<Vec<Arc<Vec<u16>>>>,
    pub zero_idx: Arc<Vec<u16>>,
    pub table: Arc<Vec<Rr>>,
}

impl<Rr: PolyRing> HFromGroup<Rr> {
    #[inline]
    fn packed_at(&self, idx: usize) -> usize {
        let base = self.base;
        let mut packed: usize = 0;
        if idx < self.len {
            // Fast path: `idx` is in-bounds for all `cols` by construction (`len` is the min length).
            for col0 in self.cols.iter() {
                let d = unsafe { *col0.get_unchecked(idx) as usize };
                packed = packed * base + d;
            }
        } else {
            // Fallback: preserve the old "short table => implicit zero_idx padding" semantics.
            for (col0, &z) in self.cols.iter().zip(self.zero_idx.iter()) {
                let d = col0.get(idx).copied().unwrap_or(z) as usize;
                packed = packed * base + d;
            }
        }
        packed
    }
}

#[cfg(feature = "parallel")]
#[inline]
fn alloc_init_par<T: Send + Sync>(len: usize, f: impl Fn(usize) -> T + Sync) -> Vec<T> {
    // Allocate without zeroing, then fill in parallel.
    let mut v: Vec<MaybeUninit<T>> = Vec::with_capacity(len);
    unsafe { v.set_len(len) };
    v.par_iter_mut()
        .enumerate()
        .for_each(|(i, slot)| {
            slot.write(f(i));
        });
    // SAFETY: all elements were written exactly once.
    unsafe { core::mem::transmute::<Vec<MaybeUninit<T>>, Vec<T>>(v) }
}

/// A structured multilinear function that supports:
/// - evaluating at a hypercube vertex index
/// - fixing variables LSB-first (same schedule as LF sumcheck)
#[derive(Clone)]
pub enum StreamingMleEnum<R: OverField + PolyRing>
where
    R::BaseRing: Ring,
{
    DenseOwned { evals: Vec<R>, num_vars: usize },
    DenseArc { evals: Arc<Vec<R>>, num_vars: usize },
    /// Base-ring scalars (constant-coeff ring elements), stored in base ring for memory.
    BaseScalarOwned {
        evals: Vec<R::BaseRing>,
        num_vars: usize,
    },
    /// Arc-wrapped base-scalar table (to share between m and m^{∘2}).
    BaseScalarArc {
        evals: Arc<Vec<R::BaseRing>>,
        num_vars: usize,
        /// If true, interpret this MLE as **vertex-wise squares** (see Symphony notes).
        square: bool,
    },
    /// Constant base-scalar MLE (same value at every vertex), optionally vertex-squared.
    ///
    /// Useful when a whole column is constant by construction (e.g. DigitsMatrix::ConstCol0 where
    /// all col>0 entries are the same "zero digit").
    BaseScalarConst {
        value: R::BaseRing,
        num_vars: usize,
        square: bool,
    },
    /// On-demand column evaluation from a dense matrix:
    /// evals[row] = ev(mat[row][col], beta) in the base ring, optionally vertex-squared.
    ///
    /// This avoids materializing the full length-2^n table up front. On the first `fix_variable`,
    /// we materialize the half-sized base-scalar table (and then proceed as `BaseScalarOwned`).
    DenseMatrixColEv {
        mat: Arc<Matrix<R>>,
        col: usize,
        beta_pows: Arc<Vec<R::BaseRing>>,
        num_vars: usize,
        square: bool,
    },
    /// On-demand column evaluation from a compact digit-backed monomial matrix:
    /// evals[row] = ev(exp(digit[row,col]), beta) in the base ring, optionally vertex-squared.
    ///
    /// On the first `fix_variable`, materializes into a half-sized `BaseScalarOwned` table.
    DigitsMatrixColEv {
        mat: Arc<DigitsMatrix<R>>,
        col: usize,
        beta_pows: Arc<Vec<R::BaseRing>>,
        num_vars: usize,
        square: bool,
    },
    /// eq(bits(index), r) in the base ring, then lifted to `R`.
    EqBase {
        scale: R::BaseRing,
        r: Vec<R::BaseRing>,
        one_minus_r: Vec<R::BaseRing>,
    },
    /// y[row] = (M * w)[row], computed from sparse rows (only before the first fix).
    SparseMatVec {
        matrix: Arc<SparseMatrix<R>>,
        witness: Arc<Vec<R>>,
        num_vars: usize,
    },
    /// y[row] = (M * w)[row] **in the base ring**, then lifted to `R` as a constant-coeff ring element.
    ///
    /// This is a fast path for the case where:
    /// - all `matrix` coefficients are constant-coeff ring elements, and
    /// - the witness vector is also constant-coeff (represented here by its base scalar table).
    ///
    /// IMPORTANT: This variant is only correct under those conditions.
    SparseMatVecConstCoeff {
        matrix: Arc<SparseMatrix<R>>,
        witness0: Arc<Vec<R::BaseRing>>,
        num_vars: usize,
    },
    /// y[row] = (M * w)[row] **in the base ring**, where `M` is stored directly as base-ring scalars.
    ///
    /// This is the natural representation for SP1 R1LF chunks (const-coeff by construction) and
    /// avoids materializing full ring elements per nonzero.
    SparseMatVecConstCoeffBase {
        matrix: Arc<SparseMatrix<R::BaseRing>>,
        witness0: Arc<Vec<R::BaseRing>>,
        num_vars: usize,
    },
    /// y[row] = (M0 * w)[row] where `M0` has **base-ring coefficients** and `w` is ring-valued.
    ///
    /// This avoids storing constant-coeff matrices as full ring elements (huge win for d64).
    SparseMatVecBaseCoeffRing {
        matrix0: Arc<SparseMatrix<R::BaseRing>>,
        witness: Arc<Vec<R>>,
        num_vars: usize,
    },
    /// Monomial vector given by small digits into `exp_table` (e.g. m_tau).
    MonomialDigitsArc {
        digits: Arc<Vec<u16>>,
        exp_table: Arc<Vec<R>>,
        num_vars: usize,
    },
    /// Sparse mat-vec where witness is monomial digits (e.g. M * m_tau).
    SparseMatVecMonomialDigits {
        matrix: Arc<SparseMatrix<R>>,
        digits: Arc<Vec<u16>>,
        exp_table: Arc<Vec<R>>,
        num_vars: usize,
    },
    /// Sparse mat-vec where matrix is base-ring coeffs but witness is monomial digits.
    SparseMatVecBaseCoeffMonomialDigits {
        matrix0: Arc<SparseMatrix<R::BaseRing>>,
        digits: Arc<Vec<u16>>,
        exp_table: Arc<Vec<R>>,
        num_vars: usize,
    },
    /// y[row] = (M * w)[row] where witness is another MLE (on-demand).
    ///
    /// Used to avoid materializing huge derived witness tables (e.g. `h`) while still supporting
    /// downstream mat-vec MLEs.
    SparseMatVecFromMle {
        matrix: Arc<SparseMatrix<R>>,
        witness_mle: Arc<StreamingMleEnum<R>>,
        num_vars: usize,
    },
    /// Same as `SparseMatVecFromMle`, but with base-ring matrix coefficients.
    SparseMatVecBaseCoeffFromMle {
        matrix0: Arc<SparseMatrix<R::BaseRing>>,
        witness_mle: Arc<StreamingMleEnum<R>>,
        num_vars: usize,
    },
    /// h[row] computed on demand from digit-backed `M_f` (ConstCol0) and `s'`.
    ///
    /// This avoids materializing the full length-2^n `h: Vec<R>` at all.
    HFromMfDigitsConstCol0 {
        /// Grouped digit streams and precomputed group tables for fast lookup.
        ///
        /// This turns `h[idx] = rest_sum + Σ_i term0_tab_i[d_i(idx)]` into:
        /// `h[idx] = rest_sum + Σ_groups table_group[pack_digits(idx)]`,
        /// reducing the number of dense ring additions per evaluation from `O(#precomps)` to `O(#groups)`.
        groups: Arc<Vec<HFromGroup<R>>>,
        /// Precomputed sum of `term_rest` across all `precomps` (dense ring element in SP1).
        ///
        /// This is a *huge* win because CM sumcheck evaluates `h[cj]` extremely often; without
        /// this, we'd re-sum dense `term_rest` values at every evaluation site.
        rest_sum: R,
        num_vars: usize,
    },
    /// A padded 4-way tensor-product table:
    /// t = t1 ⊗ t2 ⊗ t3 ⊗ t4, then padded with zeros up to 2^num_vars.
    ///
    /// Indexing matches `utils::tensor_product` nesting:
    /// for a in t1:
    ///   for b in t2:
    ///     for c in t3:
    ///       for d in t4:
    ///         push(a*b*c*d)
    Tensor4Padded {
        t1: Arc<Vec<R>>,
        t2: Arc<Vec<R>>,
        t3: Arc<Vec<R>>,
        t4: Arc<Vec<R>>,
        tensor_len: usize,
        num_vars: usize,
    },
    /// Lazily apply the first `max_lazy` variable fixes without materializing a dense half-table.
    ///
    /// This is used to **avoid the huge RSS jump** at the first fix in CM sumcheck by delaying
    /// densification to a later round (so the materialized table is smaller by 2^max_lazy).
    ///
    /// Semantics: this represents the function
    ///   g(x_high) = Σ_{b in {0,1}^k} inner( (x_high << k) | b ) * w[b]
    /// where `k = fixed.len()` and `w` is the corresponding eq-weight table for the fixed bits.
    ///
    /// Once `fixed.len() == max_lazy`, the *next* fix triggers materialization into `DenseOwned`.
    LazyFixed {
        inner: Box<StreamingMleEnum<R>>,
        /// Remaining variables in the *outer* function (excluding fixed prefix).
        num_vars: usize,
        /// Fixed prefix bits (LSB-first), stored as base-ring scalars.
        fixed: Vec<R::BaseRing>,
        /// Precomputed eq weights for the fixed prefix bits (length 2^fixed.len()).
        weights: Vec<R::BaseRing>,
        /// Maximum number of lazy fixed bits before materializing.
        max_lazy: usize,
    },

    /// CM optimization: four sparse mat-vecs that share the same sparse matrix row structure.
    ///
    /// This is used to compute the four CM mat-vec MLEs per external matrix \(M_i\):
    /// `M_i * tau`, `M_i * m_tau`, `M_i * f`, `M_i * h`.
    ///
    /// The **semantics** of each part are identical to the corresponding standalone MLE; the
    /// only difference is that the sumcheck prover can cache the shared row scan and reuse it
    /// across the four parts when evaluating at the same row index.
    CmMatVec4Part {
        shared: Arc<CmMatVec4Shared<R>>,
        which: u8, // 0..4
        num_vars: usize,
    },
}

/// Shared data for [`StreamingMleEnum::CmMatVec4Part`].
#[derive(Clone)]
pub struct CmMatVec4Shared<R: OverField + PolyRing>
where
    R::BaseRing: Ring,
{
    pub matrix0: Arc<SparseMatrix<R::BaseRing>>,
    pub w0: CmMatVecWitness<R>,
    pub w1: CmMatVecWitness<R>,
    pub w2: CmMatVecWitness<R>,
    pub w3: CmMatVecWitness<R>,
}

/// Witness sources for CM fused sparse mat-vecs.
#[derive(Clone)]
pub enum CmMatVecWitness<R: OverField + PolyRing>
where
    R::BaseRing: Ring,
{
    Base(Arc<Vec<R::BaseRing>>),
    Ring(Arc<Vec<R>>),
    MonomialDigits {
        digits: Arc<Vec<u16>>,
        exp_table: Arc<Vec<R>>,
    },
    /// Specialized monomial-digit witness where every `exp_table[d]` is monomial-like (<=1 nonzero coeff).
    ///
    /// This allows O(1) update of the accumulator coefficient instead of doing a full ring multiply
    /// by a monomial element inside the sparse row scan (important for d64).
    MonomialDigitsMonomial {
        digits: Arc<Vec<u16>>,
        mono_idx: Arc<Vec<u16>>,
        mono_coeff: Arc<Vec<R::BaseRing>>,
    },
    Mle(Arc<StreamingMleEnum<R>>),
}

impl<R: OverField + PolyRing> CmMatVec4Shared<R>
where
    R::BaseRing: Ring,
{
    #[inline(always)]
    fn add_scaled(acc: &mut R, v: &R, c0: R::BaseRing) {
        // Fused multiply-add into the accumulator to avoid creating temporaries like `v * c0`.
        // This is hot for dense `h` and ring-valued `m_tau`.
        let ac = acc.coeffs_mut();
        let vc = v.coeffs();
        debug_assert_eq!(ac.len(), vc.len());
        for i in 0..ac.len() {
            ac[i] += vc[i] * c0;
        }
    }

    #[inline]
    fn eval_witness(w: &CmMatVecWitness<R>, col: usize) -> Option<R> {
        match w {
            CmMatVecWitness::Base(v) => v.get(col).copied().map(R::from),
            CmMatVecWitness::Ring(v) => v.get(col).copied(),
            CmMatVecWitness::MonomialDigits { digits, exp_table } => {
                let cj = col;
                if cj < digits.len() {
                    Some(exp_table[digits[cj] as usize])
                } else {
                    None
                }
            }
            CmMatVecWitness::MonomialDigitsMonomial { digits, mono_idx, mono_coeff } => {
                let cj = col;
                if cj < digits.len() {
                    let d = digits[cj] as usize;
                    let idx = mono_idx[d] as usize;
                    // Return ring monomial with coefficient `mono_coeff[d]`.
                    let mut out = R::ZERO;
                    out.coeffs_mut()[idx] = mono_coeff[d];
                    Some(out)
                } else {
                    None
                }
            }
            CmMatVecWitness::Mle(m) => Some(m.eval_at_index(col)),
        }
    }

    #[inline]
    fn witness_ref<'a>(&'a self, which: u8) -> Option<&'a CmMatVecWitness<R>> {
        match which {
            0 => Some(&self.w0),
            1 => Some(&self.w1),
            2 => Some(&self.w2),
            3 => Some(&self.w3),
            _ => None,
        }
    }

    /// Evaluate all four mat-vec outputs at a given **row** index.
    #[inline]
    pub(crate) fn eval4_at_row(
        &self,
        row: usize,
        hfrom_cache: &mut Option<HFromIndexCache<R>>,
    ) -> [R; 4] {
        if row >= self.matrix0.coeffs.len() {
            return [R::ZERO; 4];
        }
        // Specialize common CM patterns to avoid per-nonzero enum branching.
        match (&self.w0, &self.w1, &self.w2, &self.w3) {
            (
                CmMatVecWitness::Base(w0),
                CmMatVecWitness::MonomialDigitsMonomial { digits, mono_idx, mono_coeff },
                CmMatVecWitness::Base(w2),
                CmMatVecWitness::Mle(w3),
            ) => {
                let w0s: &[R::BaseRing] = w0.as_ref();
                let w2s: &[R::BaseRing] = w2.as_ref();
                let digs: &[u16] = digits.as_ref();
                let mono_idx: &[u16] = mono_idx.as_ref();
                let mono_coeff: &[R::BaseRing] = mono_coeff.as_ref();
                // Peel `LazyFixed` layers with empty `fixed` (identity) so we can hit the dense/HFrom fast paths.
                let mut w3m: &StreamingMleEnum<R> = w3.as_ref();
                loop {
                    match w3m {
                        StreamingMleEnum::LazyFixed { inner, fixed, .. } if fixed.is_empty() => {
                            w3m = inner.as_ref();
                        }
                        _ => break,
                    }
                }
                enum W3Fast<'a, Rr: OverField + PolyRing>
                where
                    Rr::BaseRing: Ring,
                {
                    Dense(&'a [Rr]),
                    HFrom {
                        rest_sum: &'a Rr,
                        groups: &'a [HFromGroup<Rr>],
                        h_id: usize,
                    },
                    Other(&'a StreamingMleEnum<Rr>),
                }
                let w3fast = match w3m {
                    StreamingMleEnum::DenseArc { evals, .. } => W3Fast::Dense(evals.as_ref()),
                    StreamingMleEnum::DenseOwned { evals, .. } => W3Fast::Dense(evals.as_ref()),
                    StreamingMleEnum::HFromMfDigitsConstCol0 { groups, rest_sum, .. } => {
                        W3Fast::HFrom {
                            rest_sum,
                            groups: groups.as_ref(),
                            h_id: Arc::as_ptr(groups) as usize,
                        }
                    }
                    _ => W3Fast::Other(w3m),
                };

                let mut acc0 = R::BaseRing::ZERO;
                let mut acc2 = R::BaseRing::ZERO;
                let mut acc1 = R::ZERO;
                let mut acc3 = R::ZERO;
                let acc1_coeffs = acc1.coeffs_mut();
                // If `w3` is streamed-h, prepare the per-thread cache once per row (avoid per-nonzero checks).
                let mut hcache: Option<&mut HFromIndexCache<R>> = None;
                if let W3Fast::HFrom { h_id, .. } = &w3fast {
                    if hfrom_cache.as_ref().map(|c| c.id) != Some(*h_id) {
                        *hfrom_cache = Some(HFromIndexCache::<R>::new(*h_id));
                    }
                    hcache = Some(hfrom_cache.as_mut().unwrap());
                }
                // Collapse 3 bounds checks (w0/w2/digs) into one fast-path where possible.
                let len0 = w0s.len();
                let len2 = w2s.len();
                let len_d = digs.len();
                let len_fast = len0.min(len2).min(len_d);
                for (coeff0, col_idx) in &self.matrix0.coeffs[row] {
                    let c0 = *coeff0;
                    let cj = *col_idx;
                    if cj < len_fast {
                        // Safe: `cj < min(len0,len2,len_d)`.
                        unsafe {
                            acc0 += c0 * *w0s.get_unchecked(cj);
                            acc2 += c0 * *w2s.get_unchecked(cj);
                            let d = *digs.get_unchecked(cj) as usize;
                            // exp_table[d] is monomial-like: add into that coefficient directly.
                            acc1_coeffs[*mono_idx.get_unchecked(d) as usize] +=
                                *mono_coeff.get_unchecked(d) * c0;
                        }
                    } else {
                        if cj < len0 {
                            acc0 += c0 * w0s[cj];
                        }
                        if cj < len2 {
                            acc2 += c0 * w2s[cj];
                        }
                        if cj < len_d {
                            let d = digs[cj] as usize;
                            // exp_table[d] is monomial-like: add into that coefficient directly.
                            acc1_coeffs[mono_idx[d] as usize] += mono_coeff[d] * c0;
                        }
                    }
                    // `w3` is an MLE; avoid creating temporaries `hv * c0` in the hot loop.
                    match &w3fast {
                        W3Fast::Dense(v) => {
                            if let Some(hv) = v.get(cj) {
                                Self::add_scaled(&mut acc3, hv, c0);
                            }
                        }
                        W3Fast::HFrom { groups, rest_sum, .. } => {
                            let cache = hcache.as_mut().unwrap();
                            let hv = cache.get_or_compute(cj, || {
                                let mut hv = (*rest_sum).clone();
                                for g in groups.iter() {
                                    hv += &g.table[g.packed_at(cj)];
                                }
                                hv
                            });
                            Self::add_scaled(&mut acc3, &hv, c0);
                        }
                        W3Fast::Other(m) => {
                            let hv = m.eval_at_index(cj);
                            Self::add_scaled(&mut acc3, &hv, c0);
                        }
                    }
                }
                [R::from(acc0), acc1, R::from(acc2), acc3]
            }
            (
                CmMatVecWitness::Base(w0),
                CmMatVecWitness::Ring(w1),
                CmMatVecWitness::Base(w2),
                CmMatVecWitness::Mle(w3),
            ) => {
                let w0s: &[R::BaseRing] = w0.as_ref();
                let w1s: &[R] = w1.as_ref();
                let w2s: &[R::BaseRing] = w2.as_ref();
                // Peel `LazyFixed` layers with empty `fixed` (identity) so we can hit the dense/HFrom fast paths.
                let mut w3m: &StreamingMleEnum<R> = w3.as_ref();
                loop {
                    match w3m {
                        StreamingMleEnum::LazyFixed { inner, fixed, .. } if fixed.is_empty() => {
                            w3m = inner.as_ref();
                        }
                        _ => break,
                    }
                }
                enum W3Fast<'a, Rr: OverField + PolyRing>
                where
                    Rr::BaseRing: Ring,
                {
                    Dense(&'a [Rr]),
                    HFrom {
                        rest_sum: &'a Rr,
                        groups: &'a [HFromGroup<Rr>],
                        h_id: usize,
                    },
                    Other(&'a StreamingMleEnum<Rr>),
                }
                let w3fast = match w3m {
                    StreamingMleEnum::DenseArc { evals, .. } => W3Fast::Dense(evals.as_ref()),
                    StreamingMleEnum::DenseOwned { evals, .. } => W3Fast::Dense(evals.as_ref()),
                    StreamingMleEnum::HFromMfDigitsConstCol0 { groups, rest_sum, .. } => {
                        W3Fast::HFrom {
                            rest_sum,
                            groups: groups.as_ref(),
                            h_id: Arc::as_ptr(groups) as usize,
                        }
                    }
                    _ => W3Fast::Other(w3m),
                };

                let mut acc0 = R::BaseRing::ZERO;
                let mut acc2 = R::BaseRing::ZERO;
                let mut acc1 = R::ZERO;
                let mut acc3 = R::ZERO;
                // If `w3` is streamed-h, prepare the per-thread cache once per row (avoid per-nonzero checks).
                let mut hcache: Option<&mut HFromIndexCache<R>> = None;
                if let W3Fast::HFrom { h_id, .. } = &w3fast {
                    if hfrom_cache.as_ref().map(|c| c.id) != Some(*h_id) {
                        *hfrom_cache = Some(HFromIndexCache::<R>::new(*h_id));
                    }
                    hcache = Some(hfrom_cache.as_mut().unwrap());
                }
                // Collapse 2 bounds checks (w0/w2) into one fast-path where possible.
                let len0 = w0s.len();
                let len2 = w2s.len();
                let len_fast = len0.min(len2);
                for (coeff0, col_idx) in &self.matrix0.coeffs[row] {
                    let c0 = *coeff0;
                    let cj = *col_idx;
                    if cj < len_fast {
                        unsafe {
                            acc0 += c0 * *w0s.get_unchecked(cj);
                            acc2 += c0 * *w2s.get_unchecked(cj);
                        }
                    } else {
                        if cj < len0 {
                            acc0 += c0 * w0s[cj];
                        }
                        if cj < len2 {
                            acc2 += c0 * w2s[cj];
                        }
                    }
                    if let Some(v1) = w1s.get(cj) {
                        Self::add_scaled(&mut acc1, v1, c0);
                    }
                    match &w3fast {
                        W3Fast::Dense(v) => {
                            if let Some(hv) = v.get(cj) {
                                Self::add_scaled(&mut acc3, hv, c0);
                            }
                        }
                        W3Fast::HFrom { groups, rest_sum, .. } => {
                            let cache = hcache.as_mut().unwrap();
                            let hv = cache.get_or_compute(cj, || {
                                let mut hv = (*rest_sum).clone();
                                for g in groups.iter() {
                                    hv += &g.table[g.packed_at(cj)];
                                }
                                hv
                            });
                            Self::add_scaled(&mut acc3, &hv, c0);
                        }
                        W3Fast::Other(m) => {
                            let hv = m.eval_at_index(cj);
                            Self::add_scaled(&mut acc3, &hv, c0);
                        }
                    }
                }
                [R::from(acc0), acc1, R::from(acc2), acc3]
            }
            (CmMatVecWitness::Base(w0), CmMatVecWitness::Base(w1), CmMatVecWitness::Base(w2), CmMatVecWitness::Mle(w3)) => {
                let w0s: &[R::BaseRing] = w0.as_ref();
                let w1s: &[R::BaseRing] = w1.as_ref();
                let w2s: &[R::BaseRing] = w2.as_ref();
                // Peel `LazyFixed` layers with empty `fixed` (identity) so we can hit the dense/HFrom fast paths.
                let mut w3m: &StreamingMleEnum<R> = w3.as_ref();
                loop {
                    match w3m {
                        StreamingMleEnum::LazyFixed { inner, fixed, .. } if fixed.is_empty() => {
                            w3m = inner.as_ref();
                        }
                        _ => break,
                    }
                }
                enum W3Fast<'a, Rr: OverField + PolyRing>
                where
                    Rr::BaseRing: Ring,
                {
                    Dense(&'a [Rr]),
                    HFrom {
                        rest_sum: &'a Rr,
                        groups: &'a [HFromGroup<Rr>],
                        h_id: usize,
                    },
                    Other(&'a StreamingMleEnum<Rr>),
                }
                let w3fast = match w3m {
                    StreamingMleEnum::DenseArc { evals, .. } => W3Fast::Dense(evals.as_ref()),
                    StreamingMleEnum::DenseOwned { evals, .. } => W3Fast::Dense(evals.as_ref()),
                    StreamingMleEnum::HFromMfDigitsConstCol0 { groups, rest_sum, .. } => {
                        W3Fast::HFrom {
                            rest_sum,
                            groups: groups.as_ref(),
                            h_id: Arc::as_ptr(groups) as usize,
                        }
                    }
                    _ => W3Fast::Other(w3m),
                };

                let mut acc0 = R::BaseRing::ZERO;
                let mut acc1 = R::BaseRing::ZERO;
                let mut acc2 = R::BaseRing::ZERO;
                let mut acc3 = R::ZERO;
                // If `w3` is streamed-h, prepare the per-thread cache once per row (avoid per-nonzero checks).
                let mut hcache: Option<&mut HFromIndexCache<R>> = None;
                if let W3Fast::HFrom { h_id, .. } = &w3fast {
                    if hfrom_cache.as_ref().map(|c| c.id) != Some(*h_id) {
                        *hfrom_cache = Some(HFromIndexCache::<R>::new(*h_id));
                    }
                    hcache = Some(hfrom_cache.as_mut().unwrap());
                }
                // Collapse 3 bounds checks (w0/w1/w2) into one fast-path where possible.
                let len0 = w0s.len();
                let len1 = w1s.len();
                let len2 = w2s.len();
                let len_fast = len0.min(len1).min(len2);
                for (coeff0, col_idx) in &self.matrix0.coeffs[row] {
                    let c0 = *coeff0;
                    let cj = *col_idx;
                    if cj < len_fast {
                        unsafe {
                            acc0 += c0 * *w0s.get_unchecked(cj);
                            acc1 += c0 * *w1s.get_unchecked(cj);
                            acc2 += c0 * *w2s.get_unchecked(cj);
                        }
                    } else {
                        if cj < len0 {
                            acc0 += c0 * w0s[cj];
                        }
                        if cj < len1 {
                            acc1 += c0 * w1s[cj];
                        }
                        if cj < len2 {
                            acc2 += c0 * w2s[cj];
                        }
                    }
                    match &w3fast {
                        W3Fast::Dense(v) => {
                            if let Some(hv) = v.get(cj) {
                                Self::add_scaled(&mut acc3, hv, c0);
                            }
                        }
                        W3Fast::HFrom { groups, rest_sum, .. } => {
                            let cache = hcache.as_mut().unwrap();
                            let hv = cache.get_or_compute(cj, || {
                                let mut hv = (*rest_sum).clone();
                                for g in groups.iter() {
                                    hv += &g.table[g.packed_at(cj)];
                                }
                                hv
                            });
                            Self::add_scaled(&mut acc3, &hv, c0);
                        }
                        W3Fast::Other(m) => {
                            let hv = m.eval_at_index(cj);
                            Self::add_scaled(&mut acc3, &hv, c0);
                        }
                    }
                }
                [R::from(acc0), R::from(acc1), R::from(acc2), acc3]
            }
            _ => {
                // Generic fallback (slower): still keeps base witnesses in base ring.
                // Critical: if the witness is base-scalar, stay in the base ring for accumulation and lift once.
                let (mut b0, mut r0) = if matches!(&self.w0, CmMatVecWitness::Base(_)) {
                    (Some(R::BaseRing::ZERO), R::ZERO)
                } else {
                    (None, R::ZERO)
                };
                let (mut b1, mut r1) = if matches!(&self.w1, CmMatVecWitness::Base(_)) {
                    (Some(R::BaseRing::ZERO), R::ZERO)
                } else {
                    (None, R::ZERO)
                };
                let (mut b2, mut r2) = if matches!(&self.w2, CmMatVecWitness::Base(_)) {
                    (Some(R::BaseRing::ZERO), R::ZERO)
                } else {
                    (None, R::ZERO)
                };
                let (mut b3, mut r3) = if matches!(&self.w3, CmMatVecWitness::Base(_)) {
                    (Some(R::BaseRing::ZERO), R::ZERO)
                } else {
                    (None, R::ZERO)
                };

                for (coeff0, col_idx) in &self.matrix0.coeffs[row] {
                    let c0 = *coeff0;
                    let cj = *col_idx;

                    if let Some(acc) = b0.as_mut() {
                        if let CmMatVecWitness::Base(w) = &self.w0 {
                            if cj < w.len() {
                                *acc += c0 * w[cj];
                            }
                        }
                    } else if let Some(v) = Self::eval_witness(&self.w0, cj) {
                        r0 += v * c0;
                    }

                    if let Some(acc) = b1.as_mut() {
                        if let CmMatVecWitness::Base(w) = &self.w1 {
                            if cj < w.len() {
                                *acc += c0 * w[cj];
                            }
                        }
                    } else if let Some(v) = Self::eval_witness(&self.w1, cj) {
                        r1 += v * c0;
                    }

                    if let Some(acc) = b2.as_mut() {
                        if let CmMatVecWitness::Base(w) = &self.w2 {
                            if cj < w.len() {
                                *acc += c0 * w[cj];
                            }
                        }
                    } else if let Some(v) = Self::eval_witness(&self.w2, cj) {
                        r2 += v * c0;
                    }

                    if let Some(acc) = b3.as_mut() {
                        if let CmMatVecWitness::Base(w) = &self.w3 {
                            if cj < w.len() {
                                *acc += c0 * w[cj];
                            }
                        }
                    } else if let Some(v) = Self::eval_witness(&self.w3, cj) {
                        r3 += v * c0;
                    }
                }

                [
                    b0.map(R::from).unwrap_or(r0),
                    b1.map(R::from).unwrap_or(r1),
                    b2.map(R::from).unwrap_or(r2),
                    b3.map(R::from).unwrap_or(r3),
                ]
            }
        }
    }

    /// Evaluate a single mat-vec output at a given **row** index.
    ///
    /// This is used as a correctness-only fallback when the sumcheck prover's cache is not active
    /// (e.g. inside `LazyFixed` evaluation with non-empty fixed bits). It avoids the 4× overhead
    /// of computing all parts when only one is needed.
    #[inline]
    pub fn eval_part_at_row(&self, which: u8, row: usize) -> R {
        if row >= self.matrix0.coeffs.len() {
            return R::ZERO;
        }
        let Some(w) = self.witness_ref(which) else {
            return R::ZERO;
        };

        // If base-scalar, stay in base ring for the accumulation.
        if let CmMatVecWitness::Base(w0) = w {
            let mut sum0 = R::BaseRing::ZERO;
            for (coeff0, col_idx) in &self.matrix0.coeffs[row] {
                if *col_idx < w0.len() {
                    sum0 += *coeff0 * w0[*col_idx];
                }
            }
            return R::from(sum0);
        }

        let mut acc = R::ZERO;
        for (coeff0, col_idx) in &self.matrix0.coeffs[row] {
            let c0 = *coeff0;
            let cj = *col_idx;
            if let Some(v) = Self::eval_witness(w, cj) {
                acc += v * c0;
            }
        }
        acc
    }

    /// If the selected witness is base-scalar, evaluate the mat-vec and return the base scalar.
    #[inline]
    pub fn eval0_at_row_if_base(&self, which: u8, row: usize) -> Option<R::BaseRing> {
        if row >= self.matrix0.coeffs.len() {
            return Some(R::BaseRing::ZERO);
        }
        let w = match which {
            0 => &self.w0,
            1 => &self.w1,
            2 => &self.w2,
            3 => &self.w3,
            _ => return None,
        };
        let CmMatVecWitness::Base(w0) = w else {
            return None;
        };
        let mut sum0 = R::BaseRing::ZERO;
        for (coeff0, col_idx) in &self.matrix0.coeffs[row] {
            if *col_idx < w0.len() {
                sum0 += *coeff0 * w0[*col_idx];
            }
        }
        Some(sum0)
    }
}

impl<R: OverField + PolyRing> StreamingMleEnum<R>
where
    R::BaseRing: Ring,
{
    #[inline]
    pub fn num_vars(&self) -> usize {
        match self {
            StreamingMleEnum::DenseOwned { num_vars, .. } => *num_vars,
            StreamingMleEnum::DenseArc { num_vars, .. } => *num_vars,
            StreamingMleEnum::BaseScalarOwned { num_vars, .. } => *num_vars,
            StreamingMleEnum::BaseScalarArc { num_vars, .. } => *num_vars,
            StreamingMleEnum::BaseScalarConst { num_vars, .. } => *num_vars,
            StreamingMleEnum::DenseMatrixColEv { num_vars, .. } => *num_vars,
            StreamingMleEnum::DigitsMatrixColEv { num_vars, .. } => *num_vars,
            StreamingMleEnum::EqBase { r, .. } => r.len(),
            StreamingMleEnum::SparseMatVec { num_vars, .. } => *num_vars,
            StreamingMleEnum::SparseMatVecConstCoeff { num_vars, .. } => *num_vars,
            StreamingMleEnum::SparseMatVecConstCoeffBase { num_vars, .. } => *num_vars,
            StreamingMleEnum::SparseMatVecBaseCoeffRing { num_vars, .. } => *num_vars,
            StreamingMleEnum::MonomialDigitsArc { num_vars, .. } => *num_vars,
            StreamingMleEnum::SparseMatVecMonomialDigits { num_vars, .. } => *num_vars,
            StreamingMleEnum::SparseMatVecBaseCoeffMonomialDigits { num_vars, .. } => *num_vars,
            StreamingMleEnum::SparseMatVecFromMle { num_vars, .. } => *num_vars,
            StreamingMleEnum::SparseMatVecBaseCoeffFromMle { num_vars, .. } => *num_vars,
            StreamingMleEnum::HFromMfDigitsConstCol0 { num_vars, .. } => *num_vars,
            StreamingMleEnum::Tensor4Padded { num_vars, .. } => *num_vars,
            StreamingMleEnum::LazyFixed { num_vars, .. } => *num_vars,
            StreamingMleEnum::CmMatVec4Part { num_vars, .. } => *num_vars,
        }
    }

    #[inline]
    fn ev_fast_from_beta_pows(x: &R, beta_pows: &[R::BaseRing]) -> R::BaseRing {
        let coeffs = x.coeffs();
        debug_assert_eq!(coeffs.len(), beta_pows.len());

        // Fast monomial check: <=1 nonzero coefficient.
        let mut idx: Option<usize> = None;
        let mut c: R::BaseRing = R::BaseRing::ZERO;
        for (i, &ci) in coeffs.iter().enumerate() {
            if ci != R::BaseRing::ZERO {
                if idx.is_some() {
                    // fallback full dot
                    let mut acc = R::BaseRing::ZERO;
                    for (cj, pj) in coeffs.iter().zip(beta_pows.iter()) {
                        if *cj != R::BaseRing::ZERO {
                            acc += *cj * *pj;
                        }
                    }
                    return acc;
                }
                idx = Some(i);
                c = ci;
            }
        }
        match idx {
            None => R::BaseRing::ZERO,
            Some(i) => c * beta_pows[i],
        }
    }

    /// Evaluate the MLE at a vertex index and return the **base-ring scalar** (constant term).
    ///
    /// This is primarily for Symphony-style constant-coeff fast paths.
    #[inline]
    pub fn eval0_at_index(&self, index: usize) -> R::BaseRing {
        match self {
            StreamingMleEnum::BaseScalarOwned { evals, .. } => evals.get(index).copied().unwrap_or(R::BaseRing::ZERO),
            StreamingMleEnum::BaseScalarArc { evals, square, .. } => {
                let v = evals.get(index).copied().unwrap_or(R::BaseRing::ZERO);
                if *square { v * v } else { v }
            }
            StreamingMleEnum::BaseScalarConst { value, square, .. } => {
                if *square { *value * *value } else { *value }
            }
            StreamingMleEnum::DenseMatrixColEv {
                mat,
                col,
                beta_pows,
                square,
                ..
            } => {
                if index >= mat.nrows {
                    return R::BaseRing::ZERO;
                }
                let v0 = Self::ev_fast_from_beta_pows(&mat.vals[index][*col], beta_pows);
                if *square { v0 * v0 } else { v0 }
            }
            StreamingMleEnum::DigitsMatrixColEv {
                mat,
                col,
                beta_pows,
                square,
                ..
            } => {
                if index >= mat.nrows {
                    return R::BaseRing::ZERO;
                }
                let x = mat.get(index, *col);
                let v0 = Self::ev_fast_from_beta_pows(&x, beta_pows);
                if *square { v0 * v0 } else { v0 }
            }
            StreamingMleEnum::EqBase {
                scale,
                r,
                one_minus_r,
            } => {
                let mut prod = R::BaseRing::ONE;
                for i in 0..r.len() {
                    let bit = ((index >> i) & 1) == 1;
                    prod *= if bit { r[i] } else { one_minus_r[i] };
                }
                *scale * prod
            }
            StreamingMleEnum::SparseMatVecConstCoeff {
                matrix,
                witness0,
                ..
            } => {
                if index >= matrix.coeffs.len() {
                    return R::BaseRing::ZERO;
                }
                let mut sum0 = R::BaseRing::ZERO;
                for (coeff, col_idx) in &matrix.coeffs[index] {
                    if *col_idx < witness0.len() {
                        // Constant-coeff assumption: use only the constant term of `coeff`.
                        sum0 += coeff.coeffs()[0] * witness0[*col_idx];
                    }
                }
                sum0
            }
            StreamingMleEnum::SparseMatVecConstCoeffBase {
                matrix,
                witness0,
                ..
            } => eval0_sparse_matvec_const_coeff_base::<R>(matrix, witness0, index),
            StreamingMleEnum::SparseMatVecBaseCoeffRing { .. } => self.eval_at_index(index).coeffs()[0],
            StreamingMleEnum::SparseMatVecFromMle { .. } => self.eval_at_index(index).coeffs()[0],
            StreamingMleEnum::SparseMatVecBaseCoeffFromMle { .. } => self.eval_at_index(index).coeffs()[0],
            StreamingMleEnum::HFromMfDigitsConstCol0 { .. } => self.eval_at_index(index).coeffs()[0],
            StreamingMleEnum::LazyFixed { .. } => self.eval_at_index(index).coeffs()[0],
            StreamingMleEnum::CmMatVec4Part { .. } => self.eval_at_index(index).coeffs()[0],
            // Fallback: compute full ring value then project constant term.
            _ => self.eval_at_index(index).coeffs()[0],
        }
    }

    #[inline]
    pub fn eval_at_index(&self, index: usize) -> R {
        match self {
            // IMPORTANT: allow implicit zero-padding when the backing table is shorter than 2^num_vars.
            // This matches existing LF usage patterns where callers sometimes pass `nvars` for a
            // larger padded domain, while some intermediates live on a smaller row domain.
            StreamingMleEnum::DenseOwned { evals, .. } => evals.get(index).copied().unwrap_or(R::ZERO),
            StreamingMleEnum::DenseArc { evals, .. } => evals.get(index).copied().unwrap_or(R::ZERO),
            StreamingMleEnum::BaseScalarOwned { evals, .. } => {
                evals.get(index).copied().map(R::from).unwrap_or(R::ZERO)
            }
            StreamingMleEnum::BaseScalarArc {
                evals,
                square,
                ..
            } => {
                let v = evals.get(index).copied().unwrap_or(R::BaseRing::ZERO);
                let v = if *square { v * v } else { v };
                R::from(v)
            }
            StreamingMleEnum::BaseScalarConst { value, square, .. } => {
                let v = if *square { *value * *value } else { *value };
                R::from(v)
            }
            StreamingMleEnum::DenseMatrixColEv {
                mat,
                col,
                beta_pows,
                square,
                ..
            } => {
                if index >= mat.nrows {
                    return R::ZERO;
                }
                let v0 = Self::ev_fast_from_beta_pows(&mat.vals[index][*col], beta_pows);
                let v0 = if *square { v0 * v0 } else { v0 };
                R::from(v0)
            }
            StreamingMleEnum::DigitsMatrixColEv {
                mat,
                col,
                beta_pows,
                square,
                ..
            } => {
                if index >= mat.nrows {
                    return R::ZERO;
                }
                let x = mat.get(index, *col);
                let v0 = Self::ev_fast_from_beta_pows(&x, beta_pows);
                let v0 = if *square { v0 * v0 } else { v0 };
                R::from(v0)
            }
            StreamingMleEnum::EqBase {
                scale,
                r,
                one_minus_r,
            } => {
                let mut prod = R::BaseRing::ONE;
                for i in 0..r.len() {
                    let bit = ((index >> i) & 1) == 1;
                    prod *= if bit { r[i] } else { one_minus_r[i] };
                }
                R::from(*scale * prod)
            }
            StreamingMleEnum::SparseMatVec {
                matrix, witness, ..
            } => {
                if index >= matrix.coeffs.len() {
                    return R::ZERO;
                }
                let mut sum = R::ZERO;
                for (coeff, col_idx) in &matrix.coeffs[index] {
                    if *col_idx < witness.len() {
                        sum += *coeff * witness[*col_idx];
                    }
                }
                sum
            }
            StreamingMleEnum::SparseMatVecConstCoeff { .. } => R::from(self.eval0_at_index(index)),
            StreamingMleEnum::SparseMatVecConstCoeffBase { .. } => R::from(self.eval0_at_index(index)),
            StreamingMleEnum::SparseMatVecBaseCoeffRing {
                matrix0,
                witness,
                ..
            } => {
                if index >= matrix0.coeffs.len() {
                    return R::ZERO;
                }
                let mut sum = R::ZERO;
                for (coeff0, col_idx) in &matrix0.coeffs[index] {
                    if *col_idx < witness.len() {
                        sum += witness[*col_idx] * *coeff0;
                    }
                }
                sum
            }
            StreamingMleEnum::MonomialDigitsArc { digits, exp_table, .. } => {
                let di = digits.get(index).copied().unwrap_or(0) as usize;
                exp_table[di]
            }
            StreamingMleEnum::SparseMatVecMonomialDigits { matrix, digits, exp_table, .. } => {
                if index >= matrix.coeffs.len() {
                    return R::ZERO;
                }
                let mut sum = R::ZERO;
                for (coeff, col_idx) in &matrix.coeffs[index] {
                    let cj = *col_idx;
                    if cj < digits.len() {
                        sum += *coeff * exp_table[digits[cj] as usize];
                    }
                }
                sum
            }
            StreamingMleEnum::SparseMatVecBaseCoeffMonomialDigits {
                matrix0,
                digits,
                exp_table,
                ..
            } => {
                if index >= matrix0.coeffs.len() {
                    return R::ZERO;
                }
                let mut sum = R::ZERO;
                for (coeff0, col_idx) in &matrix0.coeffs[index] {
                    let cj = *col_idx;
                    if cj < digits.len() {
                        sum += exp_table[digits[cj] as usize] * *coeff0;
                    }
                }
                sum
            }
            StreamingMleEnum::SparseMatVecFromMle {
                matrix,
                witness_mle,
                ..
            } => {
                if index >= matrix.coeffs.len() {
                    return R::ZERO;
                }
                let mut sum = R::ZERO;
                for (coeff, col_idx) in &matrix.coeffs[index] {
                    sum += *coeff * witness_mle.eval_at_index(*col_idx);
                }
                sum
            }
            StreamingMleEnum::SparseMatVecBaseCoeffFromMle {
                matrix0,
                witness_mle,
                ..
            } => {
                if index >= matrix0.coeffs.len() {
                    return R::ZERO;
                }
                let mut sum = R::ZERO;
                for (coeff0, col_idx) in &matrix0.coeffs[index] {
                    sum += witness_mle.eval_at_index(*col_idx) * R::from(*coeff0);
                }
                sum
            }
            StreamingMleEnum::HFromMfDigitsConstCol0 { groups, rest_sum, .. } => {
                let mut acc = rest_sum.clone();
                for g in groups.iter() {
                    acc += &g.table[g.packed_at(index)];
                }
                acc
            }
            StreamingMleEnum::Tensor4Padded {
                t1,
                t2,
                t3,
                t4,
                tensor_len,
                ..
            } => {
                if index >= *tensor_len {
                    return R::ZERO;
                }
                // Index decomposition for the nested-loop order:
                // i = (((i1 * |t2| + i2) * |t3| + i3) * |t4| + i4)
                let n4 = t4.len();
                let n3 = t3.len();
                let n2 = t2.len();
                let i4 = index % n4;
                let q = index / n4;
                let i3 = q % n3;
                let q = q / n3;
                let i2 = q % n2;
                let i1 = q / n2;
                t1[i1] * t2[i2] * t3[i3] * t4[i4]
            }
            StreamingMleEnum::LazyFixed {
                inner,
                fixed,
                weights,
                ..
            } => {
                // Combine along the fixed low bits.
                let k = fixed.len();
                if k == 0 {
                    return inner.eval_at_index(index);
                }
                let mut acc = R::ZERO;
                // index refers to the remaining high bits.
                let base = index << k;
                for (b, &w) in weights.iter().enumerate() {
                    if w == R::BaseRing::ZERO {
                        continue;
                    }
                    acc += inner.eval_at_index(base | b) * w;
                }
                acc
            }
            StreamingMleEnum::CmMatVec4Part { shared, which, .. } => {
                debug_assert!((*which as usize) < 4);
                // NOTE: in the sumcheck hot loop we use a cache to compute all 4 at once and reuse.
                // This fallback must be efficient when the cache is not active (e.g. inside LazyFixed),
                // so it computes only the requested part.
                shared.eval_part_at_row(*which, index)
            }
        }
    }

    #[inline]
    pub fn fix_variable_in_place_base(&mut self, r0: R::BaseRing) {
        let nv = self.num_vars();
        assert!(nv > 0);
        let half_dom = 1usize << (nv - 1);
        let one_minus0 = R::BaseRing::ONE - r0;
        let r_ring = R::from(r0);
        match self {
            StreamingMleEnum::DenseOwned { evals, num_vars } => {
                // Allow implicit zero-padding: if `evals.len() < 2^nv`, only the prefix is stored.
                // After fixing one variable, the stored support shrinks to `ceil(len/2)`.
                let cur_len = evals.len();
                let new_len = ((cur_len + 1) >> 1).min(half_dom);
                let one_minus = R::ONE - r_ring;
                for i in 0..new_len {
                    // Allow implicit zero-padding (table shorter than 2^num_vars).
                    let a = evals.get(i << 1).copied().unwrap_or(R::ZERO);
                    let b = evals.get((i << 1) | 1).copied().unwrap_or(R::ZERO);
                    evals[i] = one_minus * a + r_ring * b;
                }
                evals.truncate(new_len);
                *num_vars -= 1;
            }
            StreamingMleEnum::DenseArc { evals, num_vars } => {
                // Avoid allocating a brand new half-sized Vec<R> (which is enormous for n=2^27):
                // - if the Arc is uniquely owned, take it and fix in-place
                // - otherwise, DO NOT clone the full table; compute the half-sized fixed table directly.
                let arc = std::mem::take(evals);
                match Arc::try_unwrap(arc) {
                    Ok(mut owned) => {
                        let cur_len = owned.len();
                        let new_len = ((cur_len + 1) >> 1).min(half_dom);
                        let one_minus = R::ONE - r_ring;
                        for i in 0..new_len {
                            let a = owned.get(i << 1).copied().unwrap_or(R::ZERO);
                            let b = owned.get((i << 1) | 1).copied().unwrap_or(R::ZERO);
                            owned[i] = one_minus * a + r_ring * b;
                        }
                        owned.truncate(new_len);
                        *self = StreamingMleEnum::DenseOwned {
                            evals: owned,
                            num_vars: *num_vars - 1,
                        };
                    }
                    Err(a) => {
                        // Shared table: allocate only the half-sized result.
                        let src: &[R] = a.as_ref();
                        let cur_len = src.len();
                        let new_len = ((cur_len + 1) >> 1).min(half_dom);
                        let one_minus = R::ONE - r_ring;
                        #[cfg(feature = "parallel")]
                        {
                            let out = alloc_init_par(new_len, |i| {
                                let aa = src.get(i << 1).copied().unwrap_or(R::ZERO);
                                let bb = src.get((i << 1) | 1).copied().unwrap_or(R::ZERO);
                                one_minus * aa + r_ring * bb
                            });
                            *self = StreamingMleEnum::DenseOwned {
                                evals: out,
                                num_vars: *num_vars - 1,
                            };
                            return;
                        }
                        #[cfg(not(feature = "parallel"))]
                        {
                            let mut out = vec![R::ZERO; new_len];
                            for i in 0..new_len {
                                let aa = src.get(i << 1).copied().unwrap_or(R::ZERO);
                                let bb = src.get((i << 1) | 1).copied().unwrap_or(R::ZERO);
                                out[i] = one_minus * aa + r_ring * bb;
                            }
                            *self = StreamingMleEnum::DenseOwned {
                                evals: out,
                                num_vars: *num_vars - 1,
                            };
                            return;
                        }
                    }
                }
            }
            StreamingMleEnum::BaseScalarOwned { evals, num_vars } => {
                let cur_len = evals.len();
                let new_len = ((cur_len + 1) >> 1).min(half_dom);
                for i in 0..new_len {
                    // Allow implicit zero-padding (table shorter than 2^num_vars).
                    let a = evals.get(i << 1).copied().unwrap_or(R::BaseRing::ZERO);
                    let b = evals
                        .get((i << 1) | 1)
                        .copied()
                        .unwrap_or(R::BaseRing::ZERO);
                    evals[i] = one_minus0 * a + r0 * b;
                }
                evals.truncate(new_len);
                *num_vars -= 1;
            }
            StreamingMleEnum::BaseScalarArc {
                evals,
                num_vars,
                square,
            } => {
                // Take ownership of the Arc if possible; otherwise clone.
                let arc = std::mem::take(evals);
                match Arc::try_unwrap(arc) {
                    Ok(mut owned) => {
                        let cur_len = owned.len();
                        let new_len = ((cur_len + 1) >> 1).min(half_dom);
                        if *square {
                            // Vertex-wise squares: square BEFORE combining.
                            for i in 0..new_len {
                                // Allow implicit zero-padding (table shorter than 2^num_vars).
                                let mut a = owned.get(i << 1).copied().unwrap_or(R::BaseRing::ZERO);
                                let mut b = owned
                                    .get((i << 1) | 1)
                                    .copied()
                                    .unwrap_or(R::BaseRing::ZERO);
                                a *= a;
                                b *= b;
                                owned[i] = one_minus0 * a + r0 * b;
                            }
                        } else {
                            for i in 0..new_len {
                                // Allow implicit zero-padding (table shorter than 2^num_vars).
                                let a = owned.get(i << 1).copied().unwrap_or(R::BaseRing::ZERO);
                                let b = owned
                                    .get((i << 1) | 1)
                                    .copied()
                                    .unwrap_or(R::BaseRing::ZERO);
                                owned[i] = one_minus0 * a + r0 * b;
                            }
                        }
                        owned.truncate(new_len);
                        // After fixing, the table is now the correct MLE values (square semantics consumed).
                        *self = StreamingMleEnum::BaseScalarOwned {
                            evals: owned,
                            num_vars: *num_vars - 1,
                        };
                    }
                    Err(a) => {
                        // Shared table: allocate only the half-sized result.
                        let src: &[R::BaseRing] = a.as_ref();
                        let cur_len = src.len();
                        let new_len = ((cur_len + 1) >> 1).min(half_dom);
                        if *square {
                            #[cfg(feature = "parallel")]
                            {
                                let out = alloc_init_par(new_len, |i| {
                                    let mut aa =
                                        src.get(i << 1).copied().unwrap_or(R::BaseRing::ZERO);
                                    let mut bb = src
                                        .get((i << 1) | 1)
                                        .copied()
                                        .unwrap_or(R::BaseRing::ZERO);
                                    aa *= aa;
                                    bb *= bb;
                                    one_minus0 * aa + r0 * bb
                                });
                                *self = StreamingMleEnum::BaseScalarOwned {
                                    evals: out,
                                    num_vars: *num_vars - 1,
                                };
                                return;
                            }
                            #[cfg(not(feature = "parallel"))]
                            {
                                let mut out = vec![R::BaseRing::ZERO; new_len];
                                for i in 0..new_len {
                                    let mut aa =
                                        src.get(i << 1).copied().unwrap_or(R::BaseRing::ZERO);
                                    let mut bb = src
                                        .get((i << 1) | 1)
                                        .copied()
                                        .unwrap_or(R::BaseRing::ZERO);
                                    aa *= aa;
                                    bb *= bb;
                                    out[i] = one_minus0 * aa + r0 * bb;
                                }
                                *self = StreamingMleEnum::BaseScalarOwned {
                                    evals: out,
                                    num_vars: *num_vars - 1,
                                };
                                return;
                            }
                        } else {
                            #[cfg(feature = "parallel")]
                            {
                                let out = alloc_init_par(new_len, |i| {
                                    let aa =
                                        src.get(i << 1).copied().unwrap_or(R::BaseRing::ZERO);
                                    let bb = src
                                        .get((i << 1) | 1)
                                        .copied()
                                        .unwrap_or(R::BaseRing::ZERO);
                                    one_minus0 * aa + r0 * bb
                                });
                                *self = StreamingMleEnum::BaseScalarOwned {
                                    evals: out,
                                    num_vars: *num_vars - 1,
                                };
                                return;
                            }
                            #[cfg(not(feature = "parallel"))]
                            {
                                let mut out = vec![R::BaseRing::ZERO; new_len];
                                for i in 0..new_len {
                                    let aa =
                                        src.get(i << 1).copied().unwrap_or(R::BaseRing::ZERO);
                                    let bb = src
                                        .get((i << 1) | 1)
                                        .copied()
                                        .unwrap_or(R::BaseRing::ZERO);
                                    out[i] = one_minus0 * aa + r0 * bb;
                                }
                                *self = StreamingMleEnum::BaseScalarOwned {
                                    evals: out,
                                    num_vars: *num_vars - 1,
                                };
                                return;
                            }
                        }
                    }
                }
            }
            StreamingMleEnum::BaseScalarConst { num_vars, .. } => {
                // Constant function stays constant after fixing; just decrement dimension.
                *num_vars -= 1;
            }
            StreamingMleEnum::DenseMatrixColEv {
                mat,
                col,
                beta_pows,
                num_vars,
                square,
            } => {
                // Materialize after the first fix into base-scalar owned table (half size).
                let half = 1usize << (*num_vars - 1);
                let one_minus0 = R::BaseRing::ONE - r0;
                let mut out = vec![R::BaseRing::ZERO; half];
                for i in 0..half {
                    let idx0 = i << 1;
                    let idx1 = (i << 1) | 1;
                    let a0 = if idx0 < mat.nrows {
                        Self::ev_fast_from_beta_pows(&mat.vals[idx0][*col], beta_pows)
                    } else {
                        R::BaseRing::ZERO
                    };
                    let b0 = if idx1 < mat.nrows {
                        Self::ev_fast_from_beta_pows(&mat.vals[idx1][*col], beta_pows)
                    } else {
                        R::BaseRing::ZERO
                    };
                    let (a0, b0) = if *square { (a0 * a0, b0 * b0) } else { (a0, b0) };
                    out[i] = one_minus0 * a0 + r0 * b0;
                }
                *self = StreamingMleEnum::BaseScalarOwned {
                    evals: out,
                    num_vars: *num_vars - 1,
                };
            }
            StreamingMleEnum::DigitsMatrixColEv {
                mat,
                col,
                beta_pows,
                num_vars,
                square,
            } => {
                let half = 1usize << (*num_vars - 1);
                let one_minus0 = R::BaseRing::ONE - r0;
                let mut out = vec![R::BaseRing::ZERO; half];
                for i in 0..half {
                    let idx0 = i << 1;
                    let idx1 = (i << 1) | 1;
                    let a0 = if idx0 < mat.nrows {
                        let x0 = mat.get(idx0, *col);
                        Self::ev_fast_from_beta_pows(&x0, beta_pows)
                    } else {
                        R::BaseRing::ZERO
                    };
                    let b0 = if idx1 < mat.nrows {
                        let x1 = mat.get(idx1, *col);
                        Self::ev_fast_from_beta_pows(&x1, beta_pows)
                    } else {
                        R::BaseRing::ZERO
                    };
                    let (a0, b0) = if *square { (a0 * a0, b0 * b0) } else { (a0, b0) };
                    out[i] = one_minus0 * a0 + r0 * b0;
                }
                *self = StreamingMleEnum::BaseScalarOwned {
                    evals: out,
                    num_vars: *num_vars - 1,
                };
            }
            StreamingMleEnum::EqBase {
                scale,
                r,
                one_minus_r,
            } => {
                let eq_factor = one_minus0 * one_minus_r[0] + r0 * r[0];
                *scale *= eq_factor;
                r.remove(0);
                one_minus_r.remove(0);
            }
            StreamingMleEnum::SparseMatVec { .. } => {
                let next = self.fix_variable(r_ring);
                *self = next;
            }
            StreamingMleEnum::MonomialDigitsArc { .. } => {
                let next = self.fix_variable(r_ring);
                *self = next;
            }
            StreamingMleEnum::SparseMatVecMonomialDigits { .. } => {
                let next = self.fix_variable(r_ring);
                *self = next;
            }
            StreamingMleEnum::SparseMatVecBaseCoeffRing { .. } => {
                let next = self.fix_variable(r_ring);
                *self = next;
            }
            StreamingMleEnum::SparseMatVecBaseCoeffMonomialDigits { .. } => {
                let next = self.fix_variable(r_ring);
                *self = next;
            }
            StreamingMleEnum::SparseMatVecFromMle { .. } => {
                let next = self.fix_variable(r_ring);
                *self = next;
            }
            StreamingMleEnum::SparseMatVecBaseCoeffFromMle { .. } => {
                let next = self.fix_variable(r_ring);
                *self = next;
            }
            StreamingMleEnum::HFromMfDigitsConstCol0 { .. } => {
                let next = self.fix_variable(r_ring);
                *self = next;
            }
            StreamingMleEnum::SparseMatVecConstCoeff {
                matrix,
                witness0,
                num_vars,
            } => {
                // Materialize after the first fix into a half-sized base-scalar table.
                //
                // IMPORTANT: avoid calling `self.eval0_at_index` here, since `self` is mutably borrowed
                // by this match arm (borrow checker).
                let nv0 = *num_vars;
                let half = 1usize << (nv0 - 1);
                let one_minus0 = R::BaseRing::ONE - r0;
                let m = matrix.clone();
                let w0 = witness0.clone();
                let mut out = vec![R::BaseRing::ZERO; half];
                for i in 0..half {
                    let idx0 = i << 1;
                    let idx1 = (i << 1) | 1;
                    let a0 = eval0_sparse_matvec_const_coeff::<R>(&m, &w0, idx0);
                    let b0 = eval0_sparse_matvec_const_coeff::<R>(&m, &w0, idx1);
                    out[i] = one_minus0 * a0 + r0 * b0;
                }
                *self = StreamingMleEnum::BaseScalarOwned {
                    evals: out,
                    num_vars: nv0 - 1,
                };
            }
            StreamingMleEnum::SparseMatVecConstCoeffBase {
                matrix,
                witness0,
                num_vars,
            } => {
                // Materialize after the first fix into a half-sized base-scalar table.
                let nv0 = *num_vars;
                let half = 1usize << (nv0 - 1);
                let one_minus0 = R::BaseRing::ONE - r0;
                let m = matrix.clone();
                let w0 = witness0.clone();
                let mut out = vec![R::BaseRing::ZERO; half];
                for i in 0..half {
                    let idx0 = i << 1;
                    let idx1 = (i << 1) | 1;
                    let a0 = eval0_sparse_matvec_const_coeff_base::<R>(&m, &w0, idx0);
                    let b0 = eval0_sparse_matvec_const_coeff_base::<R>(&m, &w0, idx1);
                    out[i] = one_minus0 * a0 + r0 * b0;
                }
                *self = StreamingMleEnum::BaseScalarOwned {
                    evals: out,
                    num_vars: nv0 - 1,
                };
            }
            StreamingMleEnum::Tensor4Padded { .. } => {
                // The tensor table is typically *tiny* (e.g. ~61k) and zero-padded to `2^num_vars`.
                // Materializing a full half-table of length `2^(nv-1)` would be catastrophic.
                //
                // After fixing one variable, the resulting table is still zero outside indices
                // i such that (2i) or (2i+1) was inside the original tensor_len. So the new support
                // length is `ceil(tensor_len / 2)`, which remains tiny.
                let nv0 = self.num_vars();
                debug_assert!(nv0 > 0);
                let half_dom = 1usize << (nv0 - 1);
                if let StreamingMleEnum::Tensor4Padded { tensor_len, .. } = self {
                    let new_len = ((*tensor_len) + 1) >> 1;
                    let new_len = new_len.min(half_dom);
                    let mut out = vec![R::ZERO; new_len];
                    let one_minus = R::ONE - r_ring;
                    for i in 0..new_len {
                        let a = self.eval_at_index(i << 1);
                        let b = self.eval_at_index((i << 1) | 1);
                        out[i] = one_minus * a + r_ring * b;
                    }
                    *self = StreamingMleEnum::DenseOwned {
                        evals: out,
                        num_vars: nv0 - 1,
                    };
                } else {
                    unreachable!();
                }
            }
            StreamingMleEnum::LazyFixed {
                inner,
                num_vars,
                fixed,
                weights,
                max_lazy,
            } => {
                // If we haven't reached the lazy threshold, just update the fixed-bit weights and shrink `num_vars`.
                if fixed.len() < *max_lazy {
                    fixed.push(r0);
                    // IMPORTANT: fixed bits are LSB-first in the sumcheck schedule.
                    // We maintain `weights[b] = Π_{j=0..k-1} (bit_j ? r_j : (1-r_j))` where `bit_j`
                    // is the j-th LSB of `b` (so the *first fixed variable* corresponds to bit 0).
                    //
                    // When appending a new fixed variable `r_k`, it becomes bit `k` (i.e. the new MSB
                    // among the fixed bits), so we update:
                    // - next[b]         = weights[b] * (1-r_k)
                    // - next[b + 2^k]   = weights[b] * r_k
                    let old_len = weights.len();
                    let mut next = vec![R::BaseRing::ZERO; old_len << 1];
                    let om = R::BaseRing::ONE - r0;
                    for b in 0..old_len {
                        let w = weights[b];
                        next[b] = w * om;
                        next[b + old_len] = w * r0;
                    }
                    *weights = next;
                    *num_vars -= 1;
                    return;
                }

                // Materialize the current outer function (with fixed bits applied) into a dense table,
                // then fall through by re-invoking the fix on the dense table.
                let cur_nv = *num_vars;
                let len = 1usize << cur_nv;
                let k = fixed.len();
                let wtab = weights.clone(); // small (<= 2^max_lazy)
                let inner_ref = inner.as_ref();

                #[cfg(feature = "parallel")]
                let dense = alloc_init_par(len, |i| {
                    let base = i << k;
                    let mut acc = R::ZERO;
                    for (b, &w) in wtab.iter().enumerate() {
                        if w == R::BaseRing::ZERO {
                            continue;
                        }
                        acc += inner_ref.eval_at_index(base | b) * w;
                    }
                    acc
                });
                #[cfg(not(feature = "parallel"))]
                let dense = {
                    let mut out = vec![R::ZERO; len];
                    for i in 0..len {
                        let base = i << k;
                        let mut acc = R::ZERO;
                        for (b, &w) in wtab.iter().enumerate() {
                            if w == R::BaseRing::ZERO {
                                continue;
                            }
                            acc += inner_ref.eval_at_index(base | b) * w;
                        }
                        out[i] = acc;
                    }
                    out
                };

                // Replace self with a dense table for the already-fixed function.
                *self = StreamingMleEnum::DenseOwned {
                    evals: dense,
                    num_vars: cur_nv,
                };
                // Now apply this fix to the dense table in-place.
                self.fix_variable_in_place_base(r0);
            }
            StreamingMleEnum::CmMatVec4Part { .. } => {
                // Critical: if this part is constant-coefficient (base-scalar), keep it base-scalar
                // after fixing (otherwise we allocate a gigantic Vec<R> and blow up RAM for d64).
                let nv0 = self.num_vars();
                let half = 1usize << (nv0 - 1);
                let one_minus0 = R::BaseRing::ONE - r0;

                // SAFETY: this arm is only entered when `self` is CmMatVec4Part.
                let (shared, which) = match self {
                    StreamingMleEnum::CmMatVec4Part { shared, which, .. } => (shared.clone(), *which),
                    _ => unreachable!(),
                };

                if let Some(_) = shared.eval0_at_row_if_base(which, 0) {
                    // Materialize into base scalars (same pattern as SparseMatVecConstCoeffBase).
                    let mut out = vec![R::BaseRing::ZERO; half];
                    for i in 0..half {
                        let a0 = shared.eval0_at_row_if_base(which, i << 1).unwrap();
                        let b0 = shared.eval0_at_row_if_base(which, (i << 1) | 1).unwrap();
                        out[i] = one_minus0 * a0 + r0 * b0;
                    }
                    *self = StreamingMleEnum::BaseScalarOwned {
                        evals: out,
                        num_vars: nv0 - 1,
                    };
                    return;
                }

                // Non-const-coeff: fall back (and rely on LazyFixed to avoid early materialization).
                let next = self.fix_variable(r_ring);
                *self = next;
            }
        }
    }

    pub fn fix_variable(&self, r: R) -> StreamingMleEnum<R> {
        let nv = self.num_vars();
        assert!(nv > 0);
        let half = 1usize << (nv - 1);
        match self {
            StreamingMleEnum::DenseOwned { evals, .. } => {
                let new_evals: Vec<R> = (0..half)
                    .map(|i| (R::ONE - r) * evals[i << 1] + r * evals[(i << 1) | 1])
                    .collect();
                StreamingMleEnum::DenseOwned {
                    evals: new_evals,
                    num_vars: nv - 1,
                }
            }
            StreamingMleEnum::DenseArc { evals, .. } => {
                let new_evals: Vec<R> = (0..half)
                    .map(|i| {
                        let a = evals.get(i << 1).copied().unwrap_or(R::ZERO);
                        let b = evals.get((i << 1) | 1).copied().unwrap_or(R::ZERO);
                        (R::ONE - r) * a + r * b
                    })
                    .collect();
                StreamingMleEnum::DenseOwned {
                    evals: new_evals,
                    num_vars: nv - 1,
                }
            }
            StreamingMleEnum::BaseScalarOwned { .. } => {
                // Keep base-scalar after fixing.
                let r0 = r.coeffs()[0];
                let one_minus0 = R::BaseRing::ONE - r0;
                let mut out = vec![R::BaseRing::ZERO; half];
                for i in 0..half {
                    let a = self.eval_at_index(i << 1).coeffs()[0];
                    let b = self.eval_at_index((i << 1) | 1).coeffs()[0];
                    out[i] = one_minus0 * a + r0 * b;
                }
                StreamingMleEnum::BaseScalarOwned {
                    evals: out,
                    num_vars: nv - 1,
                }
            }
            StreamingMleEnum::BaseScalarArc { .. } => {
                let mut c = self.clone();
                c.fix_variable_in_place_base(r.coeffs()[0]);
                c
            }
            StreamingMleEnum::BaseScalarConst { value, square, .. } => StreamingMleEnum::BaseScalarConst {
                value: *value,
                num_vars: nv - 1,
                square: *square,
            },
            StreamingMleEnum::DenseMatrixColEv { .. } => {
                let mut c = self.clone();
                c.fix_variable_in_place_base(r.coeffs()[0]);
                c
            }
            StreamingMleEnum::DigitsMatrixColEv { .. } => {
                let mut c = self.clone();
                c.fix_variable_in_place_base(r.coeffs()[0]);
                c
            }
            StreamingMleEnum::EqBase { .. } => {
                // Use in-place path by cloning and applying one fix.
                let mut c = self.clone();
                c.fix_variable_in_place_base(r.coeffs()[0]);
                c
            }
            StreamingMleEnum::SparseMatVec { .. } => {
                #[cfg(feature = "parallel")]
                let new_evals: Vec<R> = {
                    use rayon::prelude::*;
                    (0..half)
                        .into_par_iter()
                        .map(|i| {
                            let v0 = self.eval_at_index(i << 1);
                            let v1 = self.eval_at_index((i << 1) | 1);
                            (R::ONE - r) * v0 + r * v1
                        })
                        .collect()
                };
                #[cfg(not(feature = "parallel"))]
                let new_evals: Vec<R> = (0..half)
                    .map(|i| {
                        let v0 = self.eval_at_index(i << 1);
                        let v1 = self.eval_at_index((i << 1) | 1);
                        (R::ONE - r) * v0 + r * v1
                    })
                    .collect();
                StreamingMleEnum::DenseOwned {
                    evals: new_evals,
                    num_vars: nv - 1,
                }
            }
            StreamingMleEnum::SparseMatVecConstCoeff { .. } => {
                // Keep base-scalar after fixing.
                let r0 = r.coeffs()[0];
                let one_minus0 = R::BaseRing::ONE - r0;
                let mut out = vec![R::BaseRing::ZERO; half];
                for i in 0..half {
                    let a0 = self.eval0_at_index(i << 1);
                    let b0 = self.eval0_at_index((i << 1) | 1);
                    out[i] = one_minus0 * a0 + r0 * b0;
                }
                StreamingMleEnum::BaseScalarOwned {
                    evals: out,
                    num_vars: nv - 1,
                }
            }
            StreamingMleEnum::SparseMatVecConstCoeffBase { .. } => {
                // Keep base-scalar after fixing.
                let r0 = r.coeffs()[0];
                let one_minus0 = R::BaseRing::ONE - r0;
                let mut out = vec![R::BaseRing::ZERO; half];
                for i in 0..half {
                    let a0 = self.eval0_at_index(i << 1);
                    let b0 = self.eval0_at_index((i << 1) | 1);
                    out[i] = one_minus0 * a0 + r0 * b0;
                }
                StreamingMleEnum::BaseScalarOwned {
                    evals: out,
                    num_vars: nv - 1,
                }
            }
            StreamingMleEnum::MonomialDigitsArc { .. } => {
                // After fixing, the result is no longer guaranteed to be a pure monomial table,
                // so we fall back to dense materialization.
                #[cfg(feature = "parallel")]
                let new_evals: Vec<R> = {
                    use rayon::prelude::*;
                    (0..half)
                        .into_par_iter()
                        .map(|i| {
                            let v0 = self.eval_at_index(i << 1);
                            let v1 = self.eval_at_index((i << 1) | 1);
                            (R::ONE - r) * v0 + r * v1
                        })
                        .collect()
                };
                #[cfg(not(feature = "parallel"))]
                let new_evals: Vec<R> = (0..half)
                    .map(|i| {
                        let v0 = self.eval_at_index(i << 1);
                        let v1 = self.eval_at_index((i << 1) | 1);
                        (R::ONE - r) * v0 + r * v1
                    })
                    .collect();
                StreamingMleEnum::DenseOwned {
                    evals: new_evals,
                    num_vars: nv - 1,
                }
            }
            StreamingMleEnum::SparseMatVecMonomialDigits { .. } => {
                // After fixing, the result is a general dense table.
                #[cfg(feature = "parallel")]
                let new_evals: Vec<R> = {
                    use rayon::prelude::*;
                    (0..half)
                        .into_par_iter()
                        .map(|i| {
                            let v0 = self.eval_at_index(i << 1);
                            let v1 = self.eval_at_index((i << 1) | 1);
                            (R::ONE - r) * v0 + r * v1
                        })
                        .collect()
                };
                #[cfg(not(feature = "parallel"))]
                let new_evals: Vec<R> = (0..half)
                    .map(|i| {
                        let v0 = self.eval_at_index(i << 1);
                        let v1 = self.eval_at_index((i << 1) | 1);
                        (R::ONE - r) * v0 + r * v1
                    })
                    .collect();
                StreamingMleEnum::DenseOwned {
                    evals: new_evals,
                    num_vars: nv - 1,
                }
            }
            StreamingMleEnum::SparseMatVecBaseCoeffRing { .. } => {
                #[cfg(feature = "parallel")]
                let new_evals: Vec<R> = {
                    use rayon::prelude::*;
                    (0..half)
                        .into_par_iter()
                        .map(|i| {
                            let v0 = self.eval_at_index(i << 1);
                            let v1 = self.eval_at_index((i << 1) | 1);
                            (R::ONE - r) * v0 + r * v1
                        })
                        .collect()
                };
                #[cfg(not(feature = "parallel"))]
                let new_evals: Vec<R> = (0..half)
                    .map(|i| {
                        let v0 = self.eval_at_index(i << 1);
                        let v1 = self.eval_at_index((i << 1) | 1);
                        (R::ONE - r) * v0 + r * v1
                    })
                    .collect();
                StreamingMleEnum::DenseOwned {
                    evals: new_evals,
                    num_vars: nv - 1,
                }
            }
            StreamingMleEnum::SparseMatVecBaseCoeffMonomialDigits { .. } => {
                #[cfg(feature = "parallel")]
                let new_evals: Vec<R> = {
                    use rayon::prelude::*;
                    (0..half)
                        .into_par_iter()
                        .map(|i| {
                            let v0 = self.eval_at_index(i << 1);
                            let v1 = self.eval_at_index((i << 1) | 1);
                            (R::ONE - r) * v0 + r * v1
                        })
                        .collect()
                };
                #[cfg(not(feature = "parallel"))]
                let new_evals: Vec<R> = (0..half)
                    .map(|i| {
                        let v0 = self.eval_at_index(i << 1);
                        let v1 = self.eval_at_index((i << 1) | 1);
                        (R::ONE - r) * v0 + r * v1
                    })
                    .collect();
                StreamingMleEnum::DenseOwned {
                    evals: new_evals,
                    num_vars: nv - 1,
                }
            }
            StreamingMleEnum::SparseMatVecFromMle { .. } => {
                #[cfg(feature = "parallel")]
                let new_evals: Vec<R> = {
                    use rayon::prelude::*;
                    (0..half)
                        .into_par_iter()
                        .map(|i| {
                            let v0 = self.eval_at_index(i << 1);
                            let v1 = self.eval_at_index((i << 1) | 1);
                            (R::ONE - r) * v0 + r * v1
                        })
                        .collect()
                };
                #[cfg(not(feature = "parallel"))]
                let new_evals: Vec<R> = (0..half)
                    .map(|i| {
                        let v0 = self.eval_at_index(i << 1);
                        let v1 = self.eval_at_index((i << 1) | 1);
                        (R::ONE - r) * v0 + r * v1
                    })
                    .collect();
                StreamingMleEnum::DenseOwned {
                    evals: new_evals,
                    num_vars: nv - 1,
                }
            }
            StreamingMleEnum::SparseMatVecBaseCoeffFromMle { .. } => {
                #[cfg(feature = "parallel")]
                let new_evals: Vec<R> = {
                    use rayon::prelude::*;
                    (0..half)
                        .into_par_iter()
                        .map(|i| {
                            let v0 = self.eval_at_index(i << 1);
                            let v1 = self.eval_at_index((i << 1) | 1);
                            (R::ONE - r) * v0 + r * v1
                        })
                        .collect()
                };
                #[cfg(not(feature = "parallel"))]
                let new_evals: Vec<R> = (0..half)
                    .map(|i| {
                        let v0 = self.eval_at_index(i << 1);
                        let v1 = self.eval_at_index((i << 1) | 1);
                        (R::ONE - r) * v0 + r * v1
                    })
                    .collect();
                StreamingMleEnum::DenseOwned {
                    evals: new_evals,
                    num_vars: nv - 1,
                }
            }
            StreamingMleEnum::HFromMfDigitsConstCol0 { .. } => {
                #[cfg(feature = "parallel")]
                let new_evals: Vec<R> = {
                    use rayon::prelude::*;
                    (0..half)
                        .into_par_iter()
                        .map(|i| {
                            let v0 = self.eval_at_index(i << 1);
                            let v1 = self.eval_at_index((i << 1) | 1);
                            (R::ONE - r) * v0 + r * v1
                        })
                        .collect()
                };
                #[cfg(not(feature = "parallel"))]
                let new_evals: Vec<R> = (0..half)
                    .map(|i| {
                        let v0 = self.eval_at_index(i << 1);
                        let v1 = self.eval_at_index((i << 1) | 1);
                        (R::ONE - r) * v0 + r * v1
                    })
                    .collect();
                StreamingMleEnum::DenseOwned {
                    evals: new_evals,
                    num_vars: nv - 1,
                }
            }
            StreamingMleEnum::Tensor4Padded { .. } => {
                // Same idea as the in-place path: preserve sparsity of the padded tensor.
                let nv0 = nv;
                let half_dom = 1usize << (nv0 - 1);
                let tensor_len = match self {
                    StreamingMleEnum::Tensor4Padded { tensor_len, .. } => *tensor_len,
                    _ => unreachable!(),
                };
                let new_len = ((tensor_len + 1) >> 1).min(half_dom);
                #[cfg(feature = "parallel")]
                let new_evals: Vec<R> = {
                    use rayon::prelude::*;
                    (0..new_len)
                        .into_par_iter()
                        .map(|i| {
                            let v0 = self.eval_at_index(i << 1);
                            let v1 = self.eval_at_index((i << 1) | 1);
                            (R::ONE - r) * v0 + r * v1
                        })
                        .collect()
                };
                #[cfg(not(feature = "parallel"))]
                let new_evals: Vec<R> = (0..new_len)
                    .map(|i| {
                        let v0 = self.eval_at_index(i << 1);
                        let v1 = self.eval_at_index((i << 1) | 1);
                        (R::ONE - r) * v0 + r * v1
                    })
                    .collect();
                StreamingMleEnum::DenseOwned {
                    evals: new_evals,
                    num_vars: nv - 1,
                }
            }
            StreamingMleEnum::LazyFixed { .. } => {
                // Use clone + in-place to keep logic centralized.
                let mut c = self.clone();
                c.fix_variable_in_place_base(r.coeffs()[0]);
                c
            }
            StreamingMleEnum::CmMatVec4Part { shared, which, .. } => {
                // Mirror the in-place logic: if this is base-scalar, keep it base-scalar after fixing.
                let r0 = r.coeffs()[0];
                let one_minus0 = R::BaseRing::ONE - r0;
                if let Some(_) = shared.eval0_at_row_if_base(*which, 0) {
                    let mut out = vec![R::BaseRing::ZERO; half];
                    for i in 0..half {
                        let a0 = shared.eval0_at_row_if_base(*which, i << 1).unwrap();
                        let b0 = shared.eval0_at_row_if_base(*which, (i << 1) | 1).unwrap();
                        out[i] = one_minus0 * a0 + r0 * b0;
                    }
                    StreamingMleEnum::BaseScalarOwned {
                        evals: out,
                        num_vars: nv - 1,
                    }
                } else {
                    // General dense fallback.
                    #[cfg(feature = "parallel")]
                    let new_evals: Vec<R> = {
                        use rayon::prelude::*;
                        (0..half)
                            .into_par_iter()
                            .map(|i| {
                                let v0 = self.eval_at_index(i << 1);
                                let v1 = self.eval_at_index((i << 1) | 1);
                                (R::ONE - r) * v0 + r * v1
                            })
                            .collect()
                    };
                    #[cfg(not(feature = "parallel"))]
                    let new_evals: Vec<R> = (0..half)
                        .map(|i| {
                            let v0 = self.eval_at_index(i << 1);
                            let v1 = self.eval_at_index((i << 1) | 1);
                            (R::ONE - r) * v0 + r * v1
                        })
                        .collect();
                    StreamingMleEnum::DenseOwned {
                        evals: new_evals,
                        num_vars: nv - 1,
                    }
                }
            }
        }
    }
}

#[inline]
fn eval0_sparse_matvec_const_coeff<R: OverField + PolyRing>(
    matrix: &SparseMatrix<R>,
    witness0: &[R::BaseRing],
    row: usize,
) -> R::BaseRing
where
    R::BaseRing: Ring,
{
    if row >= matrix.coeffs.len() {
        return R::BaseRing::ZERO;
    }
    let mut sum0 = R::BaseRing::ZERO;
    for (coeff, col_idx) in &matrix.coeffs[row] {
        if *col_idx < witness0.len() {
            sum0 += coeff.coeffs()[0] * witness0[*col_idx];
        }
    }
    sum0
}

#[inline]
fn eval0_sparse_matvec_const_coeff_base<R: OverField + PolyRing>(
    matrix: &SparseMatrix<R::BaseRing>,
    witness0: &[R::BaseRing],
    row: usize,
) -> R::BaseRing
where
    R::BaseRing: Ring,
{
    if row >= matrix.coeffs.len() {
        return R::BaseRing::ZERO;
    }
    let mut sum0 = R::BaseRing::ZERO;
    for (coeff, col_idx) in &matrix.coeffs[row] {
        if *col_idx < witness0.len() {
            sum0 += *coeff * witness0[*col_idx];
        }
    }
    sum0
}

pub struct StreamingSumcheckState<R: OverField + PolyRing>
where
    R::BaseRing: Ring,
{
    pub mles: Vec<StreamingMleEnum<R>>,
    pub randomness: Vec<R::BaseRing>,
    pub num_vars: usize,
    pub max_degree: usize,
    pub round: usize,
}

impl<R: OverField + PolyRing> StreamingSumcheckState<R>
where
    R::BaseRing: Ring,
{
    #[inline]
    pub fn remaining_vars(&self) -> usize {
        self.mles[0].num_vars()
    }

    pub fn fix_last_variable(&mut self, r: R::BaseRing) {
        let nv = self.remaining_vars();
        assert!(nv == 1, "fix_last_variable expects 1 remaining var, got {nv}");
        #[cfg(feature = "parallel")]
        {
            self.mles
                .par_iter_mut()
                .for_each(|m| m.fix_variable_in_place_base(r));
        }
        #[cfg(not(feature = "parallel"))]
        {
            for m in self.mles.iter_mut() {
                m.fix_variable_in_place_base(r);
            }
        }
    }

    pub fn final_evals(&self) -> Vec<R> {
        let nv = self.remaining_vars();
        assert!(nv == 0, "final_evals expects 0 remaining vars, got {nv}");
        self.mles.iter().map(|m| m.eval_at_index(0)).collect()
    }
}

pub struct StreamingSumcheck;

impl StreamingSumcheck {
    pub fn prover_init<R: OverField + PolyRing>(
        mles: Vec<StreamingMleEnum<R>>,
        nvars: usize,
        degree: usize,
    ) -> StreamingSumcheckState<R>
    where
        R::BaseRing: Ring,
    {
        assert!(nvars > 0);
        assert!(!mles.is_empty());
        for m in &mles {
            assert_eq!(m.num_vars(), nvars);
        }
        StreamingSumcheckState {
            mles,
            randomness: Vec::with_capacity(nvars),
            num_vars: nvars,
            max_degree: degree,
            round: 0,
        }
    }

    pub fn prove_round<R: OverField + PolyRing>(
        state: &mut StreamingSumcheckState<R>,
        v_msg: Option<R::BaseRing>,
        comb_fn: &(dyn Fn(&[R]) -> R + Sync + Send),
    ) -> ProverMsg<R>
    where
        R::BaseRing: Ring,
    {
        if let Some(r) = v_msg {
            assert!(state.round > 0);
            state.randomness.push(r);
            // The first "fix" (right after round 1) is where the biggest allocations typically happen
            // (e.g. DenseArc -> DenseOwned half-table, mat-vec MLEs materializing, etc).
            if state.round == 1 {
                crate::utils::maybe_print_rss("streaming_sumcheck: fix(start)");
            }
            // This step is often O(total_table_size) and can dominate wall time if left serial,
            // especially when some MLE variants need to materialize on first fix.
            #[cfg(feature = "parallel")]
            {
                state
                    .mles
                    .par_iter_mut()
                    .for_each(|m| m.fix_variable_in_place_base(r));
            }
            #[cfg(not(feature = "parallel"))]
            {
                for m in state.mles.iter_mut() {
                    m.fix_variable_in_place_base(r);
                }
            }
            if state.round == 1 {
                crate::utils::maybe_print_rss("streaming_sumcheck: fix(done)");
            }
        } else {
            assert!(state.round == 0);
        }

        state.round += 1;
        assert!(state.round <= state.num_vars);

        let nv = state.mles[0].num_vars();
        let degree = state.max_degree;
        let domain_half = 1usize << (nv - 1);
        let num_polys = state.mles.len();

        struct Scratch<Rr: OverField + PolyRing>
        where
            Rr::BaseRing: Ring,
        {
            evals: Vec<Rr>,
            steps: Vec<Rr>,
            vals0: Vec<Rr>,
            vals1: Vec<Rr>,
            vals: Vec<Rr>,
            levals: Vec<Rr>,
            hfrom_cache: Option<HFromIndexCache<Rr>>,
            // CM optimization: cache one fused mat-vec row scan for the "even" and "odd" indices
            // within the current hypercube vertex pair.
            // NOTE: store cached `[R;4]` on the heap to avoid huge per-thread stack frames when
            // `R` is large (e.g. `FrogRing64`). Rayon's worker stack is limited and can overflow.
            cm_cache_even: Option<(usize, usize, Box<[Rr; 4]>)>, // (shared_id, row, vals)
            cm_cache_odd: Option<(usize, usize, Box<[Rr; 4]>)>,
            // CM optimization (LazyFixed): cache the *lazy-combined* value for a fused mat-vec
            // (summing over the fixed low bits). This prevents re-scanning the same 2^k rows for
            // each of the 4 parts.
            // Keyed only by (shared_id, base, k). For CM, the per-part LazyFixed weights are identical
            // across the 4 parts (they are fixed by the same challenges), so we can safely share.
            cm_lazy_cache_even: Option<(usize, usize, usize, Box<[Rr; 4]>)>, // (shared_id, base, k, vals)
            cm_lazy_cache_odd: Option<(usize, usize, usize, Box<[Rr; 4]>)>,
        }

        let scratch = || Scratch {
            evals: vec![R::ZERO; degree + 1],
            steps: vec![R::ZERO; num_polys],
            vals0: vec![R::ZERO; num_polys],
            vals1: vec![R::ZERO; num_polys],
            vals: vec![R::ZERO; num_polys],
            levals: vec![R::ZERO; degree + 1],
            hfrom_cache: None,
            cm_cache_even: None,
            cm_cache_odd: None,
            cm_lazy_cache_even: None,
            cm_lazy_cache_odd: None,
        };

        #[inline]
        fn eval_mle_with_cm_cache<R: OverField + PolyRing>(
            mle: &StreamingMleEnum<R>,
            index: usize,
            cache: &mut Option<(usize, usize, Box<[R; 4]>)>,
            lazy_cache: &mut Option<(usize, usize, usize, Box<[R; 4]>)>,
            hfrom_cache: &mut Option<HFromIndexCache<R>>,
        ) -> R
        where
            R::BaseRing: Ring,
        {
            // Peel/handle LazyFixed so CM caching is effective:
            // - if fixed.len()==0, just unwrap (identity)
            // - if inner is CmMatVec4Part and fixed.len()>0, compute the lazy-combined value for all 4 parts once
            if let StreamingMleEnum::LazyFixed { inner, fixed, weights, .. } = mle {
                if fixed.is_empty() {
                    return eval_mle_with_cm_cache::<R>(
                        inner.as_ref(),
                        index,
                        cache,
                        lazy_cache,
                        hfrom_cache,
                    );
                }
                if let StreamingMleEnum::CmMatVec4Part { shared, which, .. } = inner.as_ref() {
                    let k = fixed.len();
                    let wid = *which as usize;
                    debug_assert!(wid < 4);
                    let sid = Arc::as_ptr(shared) as usize;
                    let base = index << k;

                    if let Some((csid, cbase, ck, vals)) = lazy_cache.as_ref() {
                        if *csid == sid && *cbase == base && *ck == k {
                            return vals[wid];
                        }
                    }

                    // Compute all 4 parts at once: Σ_b weights[b] * shared.eval4_at_row(base|b).
                    let mut acc = [R::ZERO; 4];
                    for (b, &w) in weights.iter().enumerate() {
                        if w == R::BaseRing::ZERO {
                            continue;
                        }
                        let v = shared.eval4_at_row(base | b, hfrom_cache);
                        acc[0] += v[0] * w;
                        acc[1] += v[1] * w;
                        acc[2] += v[2] * w;
                        acc[3] += v[3] * w;
                    }
                    match lazy_cache.as_mut() {
                        Some((csid, cbase, ck, vals)) => {
                            *csid = sid;
                            *cbase = base;
                            *ck = k;
                            **vals = acc;
                        }
                        None => {
                            *lazy_cache = Some((sid, base, k, Box::new(acc)));
                        }
                    }
                    return acc[wid];
                }
            }
            if let StreamingMleEnum::CmMatVec4Part { shared, which, .. } = mle {
                let wid = *which as usize;
                debug_assert!(wid < 4);
                let sid = Arc::as_ptr(shared) as usize;
                if let Some((csid, crow, vals)) = cache.as_ref() {
                    if *csid == sid && *crow == index {
                        return vals[wid];
                    }
                }
                let vals = shared.eval4_at_row(index, hfrom_cache);
                match cache.as_mut() {
                    Some((csid, crow, cvals)) => {
                        *csid = sid;
                        *crow = index;
                        **cvals = vals;
                    }
                    None => {
                        *cache = Some((sid, index, Box::new(vals)));
                    }
                }
                return vals[wid];
            }
            mle.eval_at_index(index)
        }

        #[cfg(feature = "parallel")]
        let result = (0..domain_half)
            .into_par_iter()
            .fold(scratch, |mut s, b| {
                let idx0 = b << 1;
                let idx1 = (b << 1) | 1;
                for (i, mle) in state.mles.iter().enumerate() {
                    s.vals0[i] = eval_mle_with_cm_cache(
                        mle,
                        idx0,
                        &mut s.cm_cache_even,
                        &mut s.cm_lazy_cache_even,
                        &mut s.hfrom_cache,
                    );
                    s.vals1[i] = eval_mle_with_cm_cache(
                        mle,
                        idx1,
                        &mut s.cm_cache_odd,
                        &mut s.cm_lazy_cache_odd,
                        &mut s.hfrom_cache,
                    );
                }
                s.levals[0] = comb_fn(&s.vals0);
                s.levals[1] = comb_fn(&s.vals1);
                for i in 0..num_polys {
                    s.steps[i] = s.vals1[i] - s.vals0[i];
                    s.vals[i] = s.vals1[i];
                }
                for d in 2..=degree {
                    for i in 0..num_polys {
                        s.vals[i] += s.steps[i];
                    }
                    s.levals[d] = comb_fn(&s.vals);
                }
                for (e, l) in s.evals.iter_mut().zip(s.levals.iter()) {
                    *e += *l;
                }
                s
            })
            .map(|s| s.evals)
            // `reduce_with` avoids repeatedly allocating an identity vector in rayon's reduction tree.
            .reduce_with(|mut acc, evals| {
                for (a, e) in acc.iter_mut().zip(evals) {
                    *a += e;
                }
                acc
            })
            .unwrap_or_else(|| vec![R::ZERO; degree + 1]);

        #[cfg(not(feature = "parallel"))]
        let result = {
            let mut acc = vec![R::ZERO; degree + 1];
            let mut s = scratch();
            for b in 0..domain_half {
                let idx0 = b << 1;
                let idx1 = (b << 1) | 1;
                for (i, mle) in state.mles.iter().enumerate() {
                    s.vals0[i] = eval_mle_with_cm_cache(
                        mle,
                        idx0,
                        &mut s.cm_cache_even,
                        &mut s.cm_lazy_cache_even,
                        &mut s.hfrom_cache,
                    );
                    s.vals1[i] = eval_mle_with_cm_cache(
                        mle,
                        idx1,
                        &mut s.cm_cache_odd,
                        &mut s.cm_lazy_cache_odd,
                        &mut s.hfrom_cache,
                    );
                }
                s.levals[0] = comb_fn(&s.vals0);
                s.levals[1] = comb_fn(&s.vals1);
                for i in 0..num_polys {
                    s.steps[i] = s.vals1[i] - s.vals0[i];
                    s.vals[i] = s.vals1[i];
                }
                for d in 2..=degree {
                    for i in 0..num_polys {
                        s.vals[i] += s.steps[i];
                    }
                    s.levals[d] = comb_fn(&s.vals);
                }
                for (a, l) in acc.iter_mut().zip(s.levals.iter()) {
                    *a += *l;
                }
            }
            acc
        };


        ProverMsg { evaluations: result }
    }

    /// Degree-2 specialized round prover that avoids materializing `vals` for x=2.
    ///
    /// The caller supplies a combiner that can compute `[g(0), g(1), g(2)]` from the paired
    /// MLE evaluations at x=0 and x=1 for the current variable.
    pub fn prove_round_deg2_pairs<R: OverField + PolyRing>(
        state: &mut StreamingSumcheckState<R>,
        v_msg: Option<R::BaseRing>,
        comb_fn2: &(dyn Fn(&[R], &[R]) -> [R; 3] + Sync + Send),
    ) -> ProverMsg<R>
    where
        R::BaseRing: Ring,
    {
        assert_eq!(
            state.max_degree, 2,
            "prove_round_deg2_pairs expects degree=2"
        );

        if let Some(r) = v_msg {
            assert!(state.round > 0);
            state.randomness.push(r);
            if state.round == 1 {
                crate::utils::maybe_print_rss("streaming_sumcheck: fix(start)");
            }
            #[cfg(feature = "parallel")]
            {
                state
                    .mles
                    .par_iter_mut()
                    .for_each(|m| m.fix_variable_in_place_base(r));
            }
            #[cfg(not(feature = "parallel"))]
            {
                for m in state.mles.iter_mut() {
                    m.fix_variable_in_place_base(r);
                }
            }
            if state.round == 1 {
                crate::utils::maybe_print_rss("streaming_sumcheck: fix(done)");
            }
        } else {
            assert!(state.round == 0);
        }

        state.round += 1;
        assert!(state.round <= state.num_vars);

        let nv = state.mles[0].num_vars();
        let domain_half = 1usize << (nv - 1);
        let num_polys = state.mles.len();

        struct Scratch<Rr: OverField + PolyRing>
        where
            Rr::BaseRing: Ring,
        {
            evals: [Rr; 3],
            vals0: Vec<Rr>,
            vals1: Vec<Rr>,
            hfrom_cache: Option<HFromIndexCache<Rr>>,
            // Heap-cached `[R;4]` to avoid large stack frames for big rings (e.g. `FrogRing64`).
            cm_cache_even: Option<(usize, usize, Box<[Rr; 4]>)>,
            cm_cache_odd: Option<(usize, usize, Box<[Rr; 4]>)>,
            cm_lazy_cache_even: Option<(usize, usize, usize, Box<[Rr; 4]>)>,
            cm_lazy_cache_odd: Option<(usize, usize, usize, Box<[Rr; 4]>)>,
        }

        let scratch = || Scratch {
            evals: [R::ZERO; 3],
            vals0: vec![R::ZERO; num_polys],
            vals1: vec![R::ZERO; num_polys],
            hfrom_cache: None,
            cm_cache_even: None,
            cm_cache_odd: None,
            cm_lazy_cache_even: None,
            cm_lazy_cache_odd: None,
        };

        #[inline]
        fn eval_mle_with_cm_cache<R: OverField + PolyRing>(
            mle: &StreamingMleEnum<R>,
            index: usize,
            cache: &mut Option<(usize, usize, Box<[R; 4]>)>,
            lazy_cache: &mut Option<(usize, usize, usize, Box<[R; 4]>)>,
            hfrom_cache: &mut Option<HFromIndexCache<R>>,
        ) -> R
        where
            R::BaseRing: Ring,
        {
            if let StreamingMleEnum::LazyFixed { inner, fixed, weights, .. } = mle {
                if fixed.is_empty() {
                    return eval_mle_with_cm_cache::<R>(
                        inner.as_ref(),
                        index,
                        cache,
                        lazy_cache,
                        hfrom_cache,
                    );
                }
                if let StreamingMleEnum::CmMatVec4Part { shared, which, .. } = inner.as_ref() {
                    let k = fixed.len();
                    let wid = *which as usize;
                    debug_assert!(wid < 4);
                    let sid = Arc::as_ptr(shared) as usize;
                    let base = index << k;

                    if let Some((csid, cbase, ck, vals)) = lazy_cache.as_ref() {
                        if *csid == sid && *cbase == base && *ck == k {
                            return vals[wid];
                        }
                    }

                    let mut acc = [R::ZERO; 4];
                    for (b, &w) in weights.iter().enumerate() {
                        if w == R::BaseRing::ZERO {
                            continue;
                        }
                        let v = shared.eval4_at_row(base | b, hfrom_cache);
                        acc[0] += v[0] * w;
                        acc[1] += v[1] * w;
                        acc[2] += v[2] * w;
                        acc[3] += v[3] * w;
                    }
                    match lazy_cache.as_mut() {
                        Some((csid, cbase, ck, vals)) => {
                            *csid = sid;
                            *cbase = base;
                            *ck = k;
                            **vals = acc;
                        }
                        None => {
                            *lazy_cache = Some((sid, base, k, Box::new(acc)));
                        }
                    }
                    return acc[wid];
                }
            }
            if let StreamingMleEnum::CmMatVec4Part { shared, which, .. } = mle {
                let wid = *which as usize;
                debug_assert!(wid < 4);
                let sid = Arc::as_ptr(shared) as usize;
                if let Some((csid, crow, vals)) = cache.as_ref() {
                    if *csid == sid && *crow == index {
                        return vals[wid];
                    }
                }
                let vals = shared.eval4_at_row(index, hfrom_cache);
                match cache.as_mut() {
                    Some((csid, crow, cvals)) => {
                        *csid = sid;
                        *crow = index;
                        **cvals = vals;
                    }
                    None => {
                        *cache = Some((sid, index, Box::new(vals)));
                    }
                }
                return vals[wid];
            }
            mle.eval_at_index(index)
        }

        #[cfg(feature = "parallel")]
        let result = (0..domain_half)
            .into_par_iter()
            .fold(scratch, |mut s, b| {
                let idx0 = b << 1;
                let idx1 = (b << 1) | 1;
                for (i, mle) in state.mles.iter().enumerate() {
                    s.vals0[i] = eval_mle_with_cm_cache(
                        mle,
                        idx0,
                        &mut s.cm_cache_even,
                        &mut s.cm_lazy_cache_even,
                        &mut s.hfrom_cache,
                    );
                    s.vals1[i] = eval_mle_with_cm_cache(
                        mle,
                        idx1,
                        &mut s.cm_cache_odd,
                        &mut s.cm_lazy_cache_odd,
                        &mut s.hfrom_cache,
                    );
                }
                let le = comb_fn2(&s.vals0, &s.vals1);
                s.evals[0] += le[0];
                s.evals[1] += le[1];
                s.evals[2] += le[2];
                s
            })
            .map(|s| s.evals)
            .reduce_with(|mut acc, evals| {
                acc[0] += evals[0];
                acc[1] += evals[1];
                acc[2] += evals[2];
                acc
            })
            .unwrap_or_else(|| [R::ZERO; 3]);

        #[cfg(not(feature = "parallel"))]
        let result = {
            let mut acc = [R::ZERO; 3];
            let mut s = scratch();
            for b in 0..domain_half {
                let idx0 = b << 1;
                let idx1 = (b << 1) | 1;
                for (i, mle) in state.mles.iter().enumerate() {
                    s.vals0[i] = eval_mle_with_cm_cache(
                        mle,
                        idx0,
                        &mut s.cm_cache_even,
                        &mut s.cm_lazy_cache_even,
                        &mut s.hfrom_cache,
                    );
                    s.vals1[i] = eval_mle_with_cm_cache(
                        mle,
                        idx1,
                        &mut s.cm_cache_odd,
                        &mut s.cm_lazy_cache_odd,
                        &mut s.hfrom_cache,
                    );
                }
                let le = comb_fn2(&s.vals0, &s.vals1);
                acc[0] += le[0];
                acc[1] += le[1];
                acc[2] += le[2];
            }
            acc
        };


        ProverMsg {
            evaluations: vec![result[0], result[1], result[2]],
        }
    }

    /// Base-ring optimized variant of `prove_round` for the case where the entire computation is
    /// constant-coeff (i.e., all MLE values and the combination function live in `R::BaseRing`).
    ///
    /// The returned message is still in `R`, with each evaluation lifted as `R::from(base_eval)`.
    pub fn prove_round_base<R: OverField + PolyRing>(
        state: &mut StreamingSumcheckState<R>,
        v_msg: Option<R::BaseRing>,
        comb_fn0: &(dyn Fn(&[R::BaseRing]) -> R::BaseRing + Sync + Send),
    ) -> ProverMsg<R>
    where
        R::BaseRing: Ring,
    {
        if let Some(r) = v_msg {
            assert!(state.round > 0);
            state.randomness.push(r);
            if state.round == 1 {
                crate::utils::maybe_print_rss("streaming_sumcheck(base): fix(start)");
            }
            #[cfg(feature = "parallel")]
            {
                state
                    .mles
                    .par_iter_mut()
                    .for_each(|m| m.fix_variable_in_place_base(r));
            }
            #[cfg(not(feature = "parallel"))]
            {
                for m in state.mles.iter_mut() {
                    m.fix_variable_in_place_base(r);
                }
            }
            if state.round == 1 {
                crate::utils::maybe_print_rss("streaming_sumcheck(base): fix(done)");
            }
        } else {
            assert!(state.round == 0);
        }

        state.round += 1;
        assert!(state.round <= state.num_vars);

        let nv = state.mles[0].num_vars();
        let degree = state.max_degree;
        let domain_half = 1usize << (nv - 1);
        let num_polys = state.mles.len();

        #[cfg(feature = "parallel")]
        let evals0 = {
            // Avoid per-vertex allocations (critical for performance at 2^27).
            struct Scratch<BR> {
                acc: Vec<BR>,
                vals0: Vec<BR>,
                vals1: Vec<BR>,
                steps: Vec<BR>,
                vals: Vec<BR>,
            }
            let mk_scratch = || Scratch {
                acc: vec![R::BaseRing::ZERO; degree + 1],
                vals0: vec![R::BaseRing::ZERO; num_polys],
                vals1: vec![R::BaseRing::ZERO; num_polys],
                steps: vec![R::BaseRing::ZERO; num_polys],
                vals: vec![R::BaseRing::ZERO; num_polys],
            };
            (0..domain_half)
                .into_par_iter()
                .fold(mk_scratch, |mut s, b| {
                    for (i, mle) in state.mles.iter().enumerate() {
                        s.vals0[i] = mle.eval0_at_index(b << 1);
                        s.vals1[i] = mle.eval0_at_index((b << 1) | 1);
                    }
                    s.acc[0] += comb_fn0(&s.vals0);
                    s.acc[1] += comb_fn0(&s.vals1);
                    for i in 0..num_polys {
                        s.steps[i] = s.vals1[i] - s.vals0[i];
                        s.vals[i] = s.vals1[i];
                    }
                    for d in 2..=degree {
                        for i in 0..num_polys {
                            s.vals[i] += s.steps[i];
                        }
                        s.acc[d] += comb_fn0(&s.vals);
                    }
                    s
                })
                // `reduce_with` avoids repeatedly allocating an identity scratch in rayon's reduction tree.
                .reduce_with(|mut a, b| {
                    for d in 0..=degree {
                        a.acc[d] += b.acc[d];
                    }
                    a
                })
                .unwrap_or_else(mk_scratch)
                .acc
        };

        #[cfg(not(feature = "parallel"))]
        let evals0 = {
            let mut acc = vec![R::BaseRing::ZERO; degree + 1];
            let mut vals0 = vec![R::BaseRing::ZERO; num_polys];
            let mut vals1 = vec![R::BaseRing::ZERO; num_polys];
            let mut steps = vec![R::BaseRing::ZERO; num_polys];
            let mut vals = vec![R::BaseRing::ZERO; num_polys];
            for b in 0..domain_half {
                for (i, mle) in state.mles.iter().enumerate() {
                    vals0[i] = mle.eval0_at_index(b << 1);
                    vals1[i] = mle.eval0_at_index((b << 1) | 1);
                }
                acc[0] += comb_fn0(&vals0);
                acc[1] += comb_fn0(&vals1);
                for i in 0..num_polys {
                    steps[i] = vals1[i] - vals0[i];
                    vals[i] = vals1[i];
                }
                for d in 2..=degree {
                    for i in 0..num_polys {
                        vals[i] += steps[i];
                    }
                    acc[d] += comb_fn0(&vals);
                }
            }
            acc
        };

        ProverMsg {
            evaluations: evals0.into_iter().map(R::from).collect(),
        }
    }

    /// Run streaming sumcheck as a subprotocol (same transcript schedule as LF dense prover).
    ///
    /// Returns:
    /// - `Proof<R>` compatible with `MLSumcheck::verify_as_subprotocol`
    /// - fully sampled verifier randomness vector (length `nvars`)
    /// - final evaluations of all internal MLEs at the sampled point (same ordering as input `mles`)
    pub fn prove_as_subprotocol<R: OverField + PolyRing, T: Transcript<R>>(
        transcript: &mut T,
        mles: Vec<StreamingMleEnum<R>>,
        nvars: usize,
        degree: usize,
        comb_fn: impl Fn(&[R]) -> R + Sync + Send,
    ) -> (Proof<R>, Vec<R::BaseRing>, Vec<R>)
    where
        R::BaseRing: Ring,
    {
        let profile = std::env::var("LF_PLUS_PROFILE").ok().as_deref() == Some("1");
        let t_total = std::time::Instant::now();

        transcript.absorb_field_element(&R::BaseRing::from(nvars as u128));
        transcript.absorb_field_element(&R::BaseRing::from(degree as u128));

        let t_init = std::time::Instant::now();
        let mut state = Self::prover_init(mles, nvars, degree);
        if profile {
            println!(
                "[LF+ streaming_sumcheck] init: {:?} (nvars={}, degree={}, mles={})",
                t_init.elapsed(),
                nvars,
                degree,
                state.mles.len()
            );
        }
        let mut msgs = Vec::with_capacity(nvars);
        let mut v_msg: Option<R::BaseRing> = None;

        let mut t_rounds = std::time::Duration::from_secs(0);
        let mut t_absorb_msgs = std::time::Duration::from_secs(0);
        let mut t_get_chal = std::time::Duration::from_secs(0);
        let mut t_absorb_chal = std::time::Duration::from_secs(0);

        for round in 0..nvars {
            let t_r = std::time::Instant::now();
            let pm = Self::prove_round(&mut state, v_msg, &comb_fn);
            t_rounds += t_r.elapsed();

            let t_a = std::time::Instant::now();
            transcript.absorb_slice(&pm.evaluations);
            t_absorb_msgs += t_a.elapsed();
            msgs.push(pm);

            let t_c = std::time::Instant::now();
            let r = transcript.get_challenge();
            t_get_chal += t_c.elapsed();

            let t_ac = std::time::Instant::now();
            transcript.absorb_field_element(&r);
            t_absorb_chal += t_ac.elapsed();
            v_msg = Some(r);

            if profile && (round == 0 || round + 1 == nvars) {
                println!(
                    "[LF+ streaming_sumcheck] round {}/{} done",
                    round + 1,
                    nvars
                );
            }
        }

        // IMPORTANT: last sampled randomness is not yet applied inside the `nvars` rounds,
        // due to the standard sumcheck schedule (applied at the start of the next round).
        let last_r = v_msg.expect("nvars>0");
        state.randomness.push(last_r);
        let t_fix = std::time::Instant::now();
        state.fix_last_variable(last_r);
        let t_fix_elapsed = t_fix.elapsed();

        let t_final = std::time::Instant::now();
        let final_evals = state.final_evals();
        let t_final_elapsed = t_final.elapsed();

        if profile {
            println!(
                "[LF+ streaming_sumcheck] totals: rounds={:?} absorb_msgs={:?} get_chal={:?} absorb_chal={:?} fix_last={:?} final_evals={:?} total={:?}",
                t_rounds,
                t_absorb_msgs,
                t_get_chal,
                t_absorb_chal,
                t_fix_elapsed,
                t_final_elapsed,
                t_total.elapsed()
            );
        }

        (Proof::new(msgs), state.randomness, final_evals)
    }

    /// Degree-2 specialization of `prove_as_subprotocol` that avoids building `vals(x=2)` vectors.
    ///
    /// The transcript schedule and proof format are identical; only prover-side computation changes.
    pub fn prove_as_subprotocol_deg2_pairs<R: OverField + PolyRing, T: Transcript<R>>(
        transcript: &mut T,
        mles: Vec<StreamingMleEnum<R>>,
        nvars: usize,
        comb_fn2: impl Fn(&[R], &[R]) -> [R; 3] + Sync + Send,
    ) -> (Proof<R>, Vec<R::BaseRing>, Vec<R>)
    where
        R::BaseRing: Ring,
    {
        let profile = std::env::var("LF_PLUS_PROFILE").ok().as_deref() == Some("1");
        let t_total = std::time::Instant::now();
        let degree: usize = 2;

        transcript.absorb_field_element(&R::BaseRing::from(nvars as u128));
        transcript.absorb_field_element(&R::BaseRing::from(degree as u128));

        let t_init = std::time::Instant::now();
        let mut state = Self::prover_init(mles, nvars, degree);
        if profile {
            println!(
                "[LF+ streaming_sumcheck] init: {:?} (nvars={}, degree={}, mles={})",
                t_init.elapsed(),
                nvars,
                degree,
                state.mles.len()
            );
        }
        let mut msgs = Vec::with_capacity(nvars);
        let mut v_msg: Option<R::BaseRing> = None;

        let mut t_rounds = std::time::Duration::from_secs(0);
        let mut t_absorb_msgs = std::time::Duration::from_secs(0);
        let mut t_get_chal = std::time::Duration::from_secs(0);
        let mut t_absorb_chal = std::time::Duration::from_secs(0);

        for round in 0..nvars {
            let t_r = std::time::Instant::now();
            let pm = Self::prove_round_deg2_pairs(&mut state, v_msg, &comb_fn2);
            t_rounds += t_r.elapsed();

            let t_a = std::time::Instant::now();
            transcript.absorb_slice(&pm.evaluations);
            t_absorb_msgs += t_a.elapsed();
            msgs.push(pm);

            let t_c = std::time::Instant::now();
            let r = transcript.get_challenge();
            t_get_chal += t_c.elapsed();

            let t_ac = std::time::Instant::now();
            transcript.absorb_field_element(&r);
            t_absorb_chal += t_ac.elapsed();
            v_msg = Some(r);

            if profile && (round == 0 || round + 1 == nvars) {
                println!(
                    "[LF+ streaming_sumcheck] round {}/{} done",
                    round + 1,
                    nvars
                );
            }
        }

        let last_r = v_msg.expect("nvars>0");
        state.randomness.push(last_r);
        let t_fix = std::time::Instant::now();
        state.fix_last_variable(last_r);
        let t_fix_elapsed = t_fix.elapsed();

        let t_final = std::time::Instant::now();
        let final_evals = state.final_evals();
        let t_final_elapsed = t_final.elapsed();

        if profile {
            println!(
                "[LF+ streaming_sumcheck] totals: rounds={:?} absorb_msgs={:?} get_chal={:?} absorb_chal={:?} fix_last={:?} final_evals={:?} total={:?}",
                t_rounds,
                t_absorb_msgs,
                t_get_chal,
                t_absorb_chal,
                t_fix_elapsed,
                t_final_elapsed,
                t_total.elapsed()
            );
        }

        (Proof::new(msgs), state.randomness, final_evals)
    }

    /// Like `prove_as_subprotocol`, but using the **base-ring optimized** round prover.
    ///
    /// This is valid when:
    /// - every MLE is constant-coeff (so only its base term matters), and
    /// - the combination function lives entirely in `R::BaseRing`.
    ///
    /// The transcript interaction and proof format are identical; only prover-side computation
    /// is faster (avoids lifting to `R` during the hot loop).
    pub fn prove_as_subprotocol_base<R: OverField + PolyRing, T: Transcript<R>>(
        transcript: &mut T,
        mles: Vec<StreamingMleEnum<R>>,
        nvars: usize,
        degree: usize,
        comb_fn0: impl Fn(&[R::BaseRing]) -> R::BaseRing + Sync + Send,
    ) -> (Proof<R>, Vec<R::BaseRing>, Vec<R>)
    where
        R::BaseRing: Ring,
    {
        let profile = std::env::var("LF_PLUS_PROFILE").ok().as_deref() == Some("1");
        let t_total = std::time::Instant::now();

        transcript.absorb_field_element(&R::BaseRing::from(nvars as u128));
        transcript.absorb_field_element(&R::BaseRing::from(degree as u128));
        maybe_print_rss("streaming_sumcheck(base): start");

        let mut state = Self::prover_init(mles, nvars, degree);
        if profile {
            println!(
                "[LF+ streaming_sumcheck] init(base): {:?} (nvars={}, degree={}, mles={})",
                t_total.elapsed(),
                nvars,
                degree,
                state.mles.len()
            );
        }

        let mut msgs = Vec::with_capacity(nvars);
        let mut v_msg: Option<R::BaseRing> = None;

        for _round in 0..nvars {
            let pm = Self::prove_round_base(&mut state, v_msg, &comb_fn0);
            transcript.absorb_slice(&pm.evaluations);
            msgs.push(pm);

            let r = transcript.get_challenge();
            transcript.absorb_field_element(&r);
            v_msg = Some(r);
        }

        // Apply last sampled randomness (standard sumcheck schedule).
        let last_r = v_msg.expect("nvars>0");
        state.randomness.push(last_r);
        state.fix_last_variable(last_r);

        let final_evals = state.final_evals();
        if profile {
            println!(
                "[LF+ streaming_sumcheck] total(base): {:?}",
                t_total.elapsed()
            );
        }
        maybe_print_rss("streaming_sumcheck(base): done");

        (Proof::new(msgs), state.randomness, final_evals)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::RngCore;
    use cyclotomic_rings::rings::FrogRing64 as R;
    use stark_rings::PolyRing;
    use stark_rings_linalg::SparseMatrix;
    use std::sync::Arc;

    #[test]
    fn test_sparse_matvec_const_coeff_matches_ring_matvec() {
        // Small deterministic instance: n = 2^3.
        let n = 8usize;
        let num_vars = 3usize;
        let mut rng = ark_std::test_rng();

        // Random base-scalar witness0.
        let witness0: Vec<<R as PolyRing>::BaseRing> = (0..n)
            .map(|_| <R as PolyRing>::BaseRing::from(rng.next_u64()))
            .collect();
        let witness_ring: Vec<R> = witness0.iter().copied().map(R::from).collect();

        // Random sparse const-coeff matrix (entries are embedded base scalars).
        let mut m = SparseMatrix::<R>::identity(n);
        for row in 0..n {
            m.coeffs[row].clear();
            // ~2 nonzeros per row
            for _ in 0..2 {
                let col = (rng.next_u64() as usize) % n;
                let c0 = <R as PolyRing>::BaseRing::from(rng.next_u64());
                m.coeffs[row].push((R::from(c0), col));
            }
        }
        let matrix = Arc::new(m);

        let mle_ring = StreamingMleEnum::SparseMatVec {
            matrix: matrix.clone(),
            witness: Arc::new(witness_ring),
            num_vars,
        };
        let mle_cc = StreamingMleEnum::SparseMatVecConstCoeff {
            matrix,
            witness0: Arc::new(witness0),
            num_vars,
        };

        for idx in 0..n {
            let a = mle_ring.eval_at_index(idx);
            let b = mle_cc.eval_at_index(idx);
            assert_eq!(a, b, "mismatch at idx={idx}");
        }
    }
}
