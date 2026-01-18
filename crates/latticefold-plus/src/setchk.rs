use ark_std::log2;
use latticefold::{
    transcript::Transcript,
    utils::sumcheck::{
        utils::eq_eval,
        MLSumcheck, Proof, SumCheckError,
    },
};
use stark_rings::{OverField, PolyRing, Ring};
use stark_rings_linalg::{ops::Transpose, Matrix, SparseMatrix};
use thiserror::Error;
use std::sync::Arc;
use std::time::Instant;
use crate::utils::maybe_print_rss;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[inline]
fn is_const_coeff_ring<R: PolyRing>(x: &R) -> bool {
    x.coeffs()
        .iter()
        .skip(1)
        .all(|c| *c == <R as PolyRing>::BaseRing::ZERO)
}

#[inline]
fn is_const_coeff_sparse_matrix<R: PolyRing>(m: &SparseMatrix<R>) -> bool {
    for row in &m.coeffs {
        for (c, _j) in row {
            if !is_const_coeff_ring::<R>(c) {
                return false;
            }
        }
    }
    true
}

// (legacy) build_eq_x_r is no longer used in the streaming prover path

// cM: double commitment, commitment to M
// M: witness matrix of monomials

#[derive(Clone, Debug)]
pub enum MonomialSet<R: PolyRing> {
    /// Legacy sparse representation (kept for unit tests / small cases).
    Matrix(SparseMatrix<R>),
    /// Dense n×d matrix of monomial ring elements (preferred for large instances).
    DenseMatrix(Arc<Matrix<R>>),
    /// Compact monomial matrix represented by small digits:
    /// entry(row,col) = exp(digit_elems[digits[row*ncols + col]]).
    ///
    /// This avoids storing `n×d` full ring elements (which is huge), while keeping prover
    /// transcript / verifier behavior identical.
    DigitsMatrix(Arc<DigitsMatrix<R>>),
    /// Vector set of monomials (Arc-backed to avoid cloning large vectors).
    Vector(Arc<Vec<R>>),
    /// Vector set of monomials stored compactly as indices into an `exp_table`.
    ///
    /// `digits[i]` indexes `exp_table` and the element is `exp_table[digits[i]]`.
    VectorDigits {
        digits: Arc<Vec<u16>>,
        exp_table: Arc<Vec<R>>,
    },
}

/// Compact monomial matrix backed by a digit table.
///
/// - `digits` stores indices into `exp_table` (row-major).
/// - `exp_table[idx]` is the ring monomial element corresponding to that digit.
///
/// NOTE: For large digit bases (e.g. balanced base \(2^{16}\)), `exp_table` may be a fixed-size
/// lookup table (e.g. 65536 entries) and `digits` are pre-mapped into `[0..exp_table.len())`.
#[derive(Clone, Debug)]
pub struct DigitsMatrix<R: PolyRing> {
    pub nrows: usize,
    pub ncols: usize,
    pub digits: DigitsBacking,
    pub exp_table: Arc<Vec<R>>,
}

/// Backing storage for [`DigitsMatrix`].
///
/// For large instances where the underlying witness is known to be **constant-coefficient**
/// (i.e. only `col=0` can be non-trivial), we can store only the digits for `col=0` and
/// treat all `col>0` entries as the digit for zero.
#[derive(Clone, Debug)]
pub enum DigitsBacking {
    /// Full row-major `nrows × ncols` table.
    Full(Arc<Vec<u16>>),
    /// Only `col=0` digits are stored (length `nrows`).
    ConstCol0 {
        col0: Arc<Vec<u16>>,
        zero_idx: u16,
    },
}

impl<R: PolyRing> DigitsMatrix<R> {
    #[inline]
    pub fn digit_idx(&self, row: usize, col: usize) -> usize {
        debug_assert!(row < self.nrows);
        debug_assert!(col < self.ncols);
        match &self.digits {
            DigitsBacking::Full(d) => (d[row * self.ncols + col]) as usize,
            DigitsBacking::ConstCol0 { col0, zero_idx } => {
                if col == 0 {
                    // Allow implicit zero-padding when we only materialize a prefix.
                    col0.get(row).copied().unwrap_or(*zero_idx) as usize
                } else {
                    (*zero_idx) as usize
                }
            }
        }
    }

    #[inline]
    pub fn get(&self, row: usize, col: usize) -> R {
        let idx = self.digit_idx(row, col);
        self.exp_table[idx]
    }
}

#[derive(Clone, Debug)]
pub struct In<R: PolyRing> {
    pub nvars: usize,
    pub sets: Vec<MonomialSet<R>>, // Ms and ms: n x m, or n
}

/// External matrices used by the LF+ set-check (coming from the surrounding protocol).
///
/// In the SP1/const-coeff regime, these matrices are naturally represented over the base ring.
/// We allow passing them either as full ring matrices (`Ring`) or base-ring matrices (`Base`) to
/// avoid catastrophic memory usage when `R::dimension()` is large (e.g. d64).
#[derive(Clone, Copy, Debug)]
pub enum ExternalMats<'a, R: PolyRing> {
    Ring(&'a [Arc<SparseMatrix<R>>]),
    Base(&'a [Arc<SparseMatrix<R::BaseRing>>]),
}

impl<'a, R: PolyRing> ExternalMats<'a, R> {
    #[inline]
    pub fn len(&self) -> usize {
        match self {
            ExternalMats::Ring(m) => m.len(),
            ExternalMats::Base(m) => m.len(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct Out<R: PolyRing> {
    pub nvars: usize,
    pub r: Vec<R::BaseRing>, // log n
    pub sumcheck_proof: Proof<R>,
    pub e: Vec<Vec<Vec<R>>>, // m, matrices outputs
    pub b: Vec<R>,           // vectors outputs
}

#[derive(Debug, Error)]
pub enum SetCheckError<R: Ring> {
    #[error("Sumcheck failed: {0}")]
    Sumcheck(#[from] SumCheckError<R>),
    #[error("Recomputed claim `v` mismatch: expected = {0}, received = {1}")]
    ExpectedEvaluation(R, R),
}

fn ev<R: PolyRing>(r: &R, x: R::BaseRing) -> R::BaseRing {
    r.coeffs()
        .iter()
        .fold(
            (R::BaseRing::ZERO, R::BaseRing::ONE),
            |(mut acc, exp), c| {
                acc += *c * exp;
                (acc, exp * x)
            },
        )
        .0
}

impl<R: OverField + PolyRing> In<R> {
    /// Monomial set check
    ///
    /// Proves sets rings are all unit monomials.
    /// Currently requires k >= 1 monomial matrices sets. TODO support other scenarios.
    /// If k > 1, sumcheck batching is employed.
    pub fn set_check(&self, M: ExternalMats<'_, R>, transcript: &mut impl Transcript<R>) -> Out<R> {
        let profile = std::env::var("LF_PLUS_PROFILE").ok().as_deref() == Some("1");
        let t_total = Instant::now();
        maybe_print_rss("setchk: start");
        let mlen = M.len();

        let Ms_sparse: Vec<&SparseMatrix<R>> = self
            .sets
            .iter()
            .filter_map(|set| match set {
                MonomialSet::Matrix(m) => Some(m),
                _ => None,
            })
            .collect();
        let Ms_dense: Vec<Arc<Matrix<R>>> = self
            .sets
            .iter()
            .filter_map(|set| match set {
                MonomialSet::DenseMatrix(m) => Some(m.clone()),
                _ => None,
            })
            .collect();
        let Ms_digits: Vec<Arc<DigitsMatrix<R>>> = self
            .sets
            .iter()
            .filter_map(|set| match set {
                MonomialSet::DigitsMatrix(m) => Some(m.clone()),
                _ => None,
            })
            .collect();
        enum VecSet<Rr: PolyRing> {
            Dense(Arc<Vec<Rr>>),
            Digits { digits: Arc<Vec<u16>>, exp_table: Arc<Vec<Rr>> },
        }
        let ms: Vec<VecSet<R>> = self
            .sets
            .iter()
            .filter_map(|set| match set {
                MonomialSet::Vector(v) => Some(VecSet::Dense(v.clone())),
                MonomialSet::VectorDigits { digits, exp_table } => {
                    Some(VecSet::Digits { digits: digits.clone(), exp_table: exp_table.clone() })
                }
                _ => None,
            })
            .collect();

        assert!(
            !Ms_sparse.is_empty() || !Ms_dense.is_empty() || !Ms_digits.is_empty(),
            "set_check requires at least one matrix set"
        );
        let (nrows, ncols) = if let Some(m0) = Ms_dense.first() {
            (m0.nrows, m0.ncols)
        } else if let Some(m0) = Ms_digits.first() {
            (m0.nrows, m0.ncols)
        } else {
            (Ms_sparse[0].nrows, Ms_sparse[0].ncols)
        };
        maybe_print_rss("setchk: classified sets");
        let tnvars = log2(nrows.next_power_of_two()) as usize;
        let MTs = Ms_sparse.iter().map(|M| (*M).transpose()).collect::<Vec<_>>();

        // Streaming MLEs (avoid materializing DenseMultilinearExtension tables).
        use crate::streaming_sumcheck::{StreamingMleEnum, StreamingSumcheck};
        let Ms_len = Ms_dense.len() + Ms_digits.len() + Ms_sparse.len();
        let mut mles: Vec<StreamingMleEnum<R>> =
            Vec::with_capacity((Ms_len + ms.len()) * (ncols * 2 + 1));
        let mut alphas = Vec::with_capacity(Ms_len);
        // Track which matrix sets are stored as `DigitsBacking::ConstCol0`.
        // For those, only column 0 varies by row; columns 1..d-1 are fixed to a constant monomial
        // and contribute identically zero to the set-check polynomial. We can omit them from the
        // sumcheck MLE list to save huge memory/time.
        let mut mat_is_constcol0: Vec<bool> = Vec::with_capacity(Ms_len);

        // matrix sets (dense path)
        for (_mi, Md) in Ms_dense.iter().enumerate() {
            let _t_mat = Instant::now();
            // Step 1
            let c0 = transcript.get_challenges(self.nvars);
            let one_minus_c0 = c0.iter().copied().map(|x| R::BaseRing::ONE - x).collect();
            let beta = transcript.get_challenge();

            // Step 2
            // Fast evaluation uses precomputed beta powers (degree = ring dimension).
            let beta_pows = beta_pows::<R>(beta);

            // Avoid materializing full `nrows` tables up front:
            // represent each column as an on-demand MLE that materializes only after the first fix.
            let mat = Md.clone(); // Arc clone (no data copy)
            let beta_pows = Arc::new(beta_pows);
            for col in 0..ncols {
                mles.push(StreamingMleEnum::DenseMatrixColEv {
                    mat: mat.clone(),
                    col,
                    beta_pows: beta_pows.clone(),
                    num_vars: tnvars,
                    square: false,
                });
                mles.push(StreamingMleEnum::DenseMatrixColEv {
                    mat: mat.clone(),
                    col,
                    beta_pows: beta_pows.clone(),
                    num_vars: tnvars,
                    square: true,
                });
            }

            // eq(x,c) as base-ring structured MLE (constant-coeff)
            mles.push(StreamingMleEnum::EqBase {
                scale: R::BaseRing::ONE,
                r: c0,
                one_minus_r: one_minus_c0,
            });

            let alpha = transcript.get_challenge();
            alphas.push(alpha);
            mat_is_constcol0.push(false);
        }

        // matrix sets (digit/oracle path)
        for (_mi, Md) in Ms_digits.iter().enumerate() {
            let _t_mat = Instant::now();
            // Step 1
            let c0 = transcript.get_challenges(self.nvars);
            let one_minus_c0 = c0.iter().copied().map(|x| R::BaseRing::ONE - x).collect();
            let beta = transcript.get_challenge();

            // Step 2
            let beta_pows = beta_pows::<R>(beta);
            let mat = Md.clone();
            let beta_pows = Arc::new(beta_pows);
            let is_constcol0 = matches!(&mat.digits, DigitsBacking::ConstCol0 { .. });
            if is_constcol0 {
                // Only column 0 participates in the polynomial; other columns are constant and contribute 0.
                mles.push(StreamingMleEnum::DigitsMatrixColEv {
                    mat: mat.clone(),
                    col: 0,
                    beta_pows: beta_pows.clone(),
                    num_vars: tnvars,
                    square: false,
                });
                mles.push(StreamingMleEnum::DigitsMatrixColEv {
                    mat: mat.clone(),
                    col: 0,
                    beta_pows: beta_pows.clone(),
                    num_vars: tnvars,
                    square: true,
                });
            } else {
                for col in 0..ncols {
                    mles.push(StreamingMleEnum::DigitsMatrixColEv {
                        mat: mat.clone(),
                        col,
                        beta_pows: beta_pows.clone(),
                        num_vars: tnvars,
                        square: false,
                    });
                    mles.push(StreamingMleEnum::DigitsMatrixColEv {
                        mat: mat.clone(),
                        col,
                        beta_pows: beta_pows.clone(),
                        num_vars: tnvars,
                        square: true,
                    });
                }
            }

            // eq(x,c)
            mles.push(StreamingMleEnum::EqBase {
                scale: R::BaseRing::ONE,
                r: c0,
                one_minus_r: one_minus_c0,
            });

            let alpha = transcript.get_challenge();
            alphas.push(alpha);
            mat_is_constcol0.push(is_constcol0);
        }

        // matrix sets (legacy sparse path)
        for (_mi, M) in Ms_sparse.iter().enumerate() {
            let _t_mat = Instant::now();
            let c0 = transcript.get_challenges(self.nvars);
            let one_minus_c0 = c0.iter().copied().map(|x| R::BaseRing::ONE - x).collect();
            let beta = transcript.get_challenge();

            let MT = (*M).transpose();
            let beta_pows = beta_pows::<R>(beta);

            // This is inherently limited-parallel (ncols is small). Kept for small/unit tests.
            let col_tables: Vec<Arc<Vec<R::BaseRing>>> = (0..ncols)
                .map(|col| {
                    let row = &MT.coeffs[col];
                    let mut v = vec![R::BaseRing::ZERO; M.nrows];
                    for (r_ij, idx) in row.iter() {
                        v[*idx] = ev_fast::<R>(r_ij, &beta_pows);
                    }
                    Arc::new(v)
                })
                .collect();

            for col in 0..ncols {
                let tab = col_tables[col].clone();
                mles.push(StreamingMleEnum::BaseScalarArc { evals: tab.clone(), num_vars: tnvars, square: false });
                mles.push(StreamingMleEnum::BaseScalarArc { evals: tab, num_vars: tnvars, square: true });
            }
            mles.push(StreamingMleEnum::EqBase {
                scale: R::BaseRing::ONE,
                r: c0,
                one_minus_r: one_minus_c0,
            });
            let alpha = transcript.get_challenge();
            alphas.push(alpha);
            mat_is_constcol0.push(false);
        }

        // vector sets
        for (vi, mset) in ms.iter().enumerate() {
            let t_vec = Instant::now();
            // Step 1
            let c0 = transcript.get_challenges(self.nvars);
            let one_minus_c0 = c0.iter().copied().map(|x| R::BaseRing::ONE - x).collect();
            let beta = transcript.get_challenge();

            let beta_pows = beta_pows::<R>(beta);
            let (v0, m_len) = match mset {
                VecSet::Dense(m) => {
                    let mut v0 = vec![R::BaseRing::ZERO; m.len()];
                    for (i, r_i) in m.iter().enumerate() {
                        v0[i] = ev_fast::<R>(r_i, &beta_pows);
                    }
                    (v0, m.len())
                }
                VecSet::Digits { digits, exp_table } => {
                    // Precompute ev(exp_table[d], beta) for the small digit alphabet, then map.
                    let mut ev_tab = vec![R::BaseRing::ZERO; exp_table.len()];
                    for (di, r_di) in exp_table.iter().enumerate() {
                        ev_tab[di] = ev_fast::<R>(r_di, &beta_pows);
                    }
                    let mut v0 = vec![R::BaseRing::ZERO; digits.len()];
                    for (i, &dix) in digits.iter().enumerate() {
                        v0[i] = ev_tab[dix as usize];
                    }
                    (v0, digits.len())
                }
            };
            let tab = Arc::new(v0);
            mles.push(StreamingMleEnum::BaseScalarArc { evals: tab.clone(), num_vars: tnvars, square: false });
            mles.push(StreamingMleEnum::BaseScalarArc { evals: tab, num_vars: tnvars, square: true });
            mles.push(StreamingMleEnum::EqBase {
                scale: R::BaseRing::ONE,
                r: c0,
                one_minus_r: one_minus_c0,
            });

            let alpha = transcript.get_challenge();
            alphas.push(alpha);

            if profile {
                println!(
                    "[LF+ setchk] vector_set[{vi}] build_table: {:?} (len={})",
                    t_vec.elapsed(),
                    m_len
                );
            }
        }
        maybe_print_rss("setchk: after build mles");

        // random linear combinator, for batching
        let rc: Option<R::BaseRing> = (Ms_len > 1).then(|| transcript.get_challenge());

        // Precompute alpha powers for the (rare) non-ConstCol0 matrix sets.
        // This avoids doing `alpha.pow([j])` inside the per-point sumcheck combiner.
        let alpha_pows: Vec<Option<Vec<R::BaseRing>>> = (0..Ms_len)
            .map(|i| {
                if mat_is_constcol0.get(i).copied().unwrap_or(false) {
                    None
                } else {
                    let alpha = alphas[i];
                    let mut p = Vec::with_capacity(ncols);
                    let mut acc = R::BaseRing::ONE; // alpha^0
                    for _ in 0..ncols {
                        p.push(acc);
                        acc *= alpha;
                    }
                    Some(p)
                }
            })
            .collect();

        // Base-ring combiner (all tables here are constant-coeff), so we can use the
        // base-optimized streaming sumcheck prover.
        let comb_fn0 = |vals: &[R::BaseRing]| -> R::BaseRing {
            // When there is only a single term, `rc` is omitted; semantically this means
            // "all terms have weight 1". (We must still include both matrix- and vector-set
            // contributions, otherwise the sumcheck claim won't match Step 3.)
            let mut lc = R::BaseRing::ZERO;

            let mut s = 0usize;
            let mut rc_pow = R::BaseRing::ONE;
            for i in 0..Ms_len {
                let mut res = R::BaseRing::ZERO;
                if mat_is_constcol0.get(i).copied().unwrap_or(false) {
                    // Layout: [m0, m0', eq] (stride 3). alpha^0 = 1.
                    res += vals[s] * vals[s] - vals[s + 1];
                    res *= vals[s + 2]; // eq
                    s += 3;
                } else {
                    // Layout: [m_j, m'_j] for j=0..d-1, then eq (stride 2*d+1).
                    let ap = alpha_pows[i].as_ref().expect("alpha powers must exist for non-ConstCol0");
                    for j in 0..ncols {
                        res += (vals[s + j * 2] * vals[s + j * 2] - vals[s + j * 2 + 1])
                            * ap[j];
                    }
                    res *= vals[s + 2 * ncols]; // eq
                    s += 2 * ncols + 1;
                }
                lc += res * rc_pow;
                if let Some(rc) = rc {
                    rc_pow *= rc;
                }
            }

            for i in 0..ms.len() {
                let s0 = s + i * 3;
                let alpha_idx = Ms_len + i;
                let mut res = R::BaseRing::ZERO;
                res += (vals[s0] * vals[s0] - vals[s0 + 1]) * alphas[alpha_idx];
                res *= vals[s0 + 2]; // eq
                lc += res * rc_pow;
                if let Some(rc) = rc {
                    rc_pow *= rc;
                }
            }
            lc
        };

        let t_sc = Instant::now();
        maybe_print_rss("setchk: before sumcheck");
        let (sumcheck_proof, r, _final_vals) =
            StreamingSumcheck::prove_as_subprotocol_base(transcript, mles, self.nvars, 3, comb_fn0);
        maybe_print_rss("setchk: after sumcheck");
        if profile {
            println!(
                "[LF+ setchk] sumcheck: {:?} (nvars={}, degree=3, ncols={}, Ms={}, ms={})",
                t_sc.elapsed(),
                self.nvars,
                ncols,
                Ms_len,
                ms.len()
            );
        }

        // Step 3
        let t_step3 = Instant::now();
        maybe_print_rss("setchk: step3 start");
        let t_y_mats = Instant::now();
        // Avoid materializing full length-2^n eq table: use (low-table + high-scale).
        let one_minus_r = r.iter().copied().map(|x| R::BaseRing::ONE - x).collect::<Vec<_>>();
        let t_low = choose_t_low(self.nvars);
        let low = build_eq_low_table::<R>(&r[..t_low], &one_minus_r[..t_low]);
        let low_len = 1usize << t_low;
        let low_mask = low_len - 1;
        let scale_high = build_scale_high::<R>(&r, &one_minus_r, t_low);
        let eq_at = |idx: usize| -> R::BaseRing {
            eq_weight_base::<R>(idx, &low, &scale_high, low_mask, t_low)
        };

        // Precompute y_i = M_i^T * eq_r (so eval(M_i * row)(r) = <y_i, row>).
        //
        // CRITICAL (Symphony lesson): if `M` is constant-coeff (SP1 regime), `y_i` is also
        // constant-coeff. Storing it as `Vec<R>` blows up RAM for large `ncols` (e.g. 2^27 with d=64).
        // We therefore store `y_i` in the **base ring** whenever possible.
        enum YMats<Rr: PolyRing> {
            Base(Vec<Vec<Rr::BaseRing>>),
            Ring(Vec<Vec<Rr>>),
        }

        let y_mats: YMats<R> = match M {
            ExternalMats::Base(M0) => {
                #[cfg(feature = "parallel")]
                {
                    YMats::Base(
                        M0.par_iter()
                            .map(|mi| {
                                let mi = mi.as_ref();
                                let mut y0 = vec![R::BaseRing::ZERO; mi.ncols];
                                for (row_idx, row) in mi.coeffs.iter().enumerate() {
                                    let w = eq_at(row_idx);
                                    for (coeff0, col_idx) in row {
                                        y0[*col_idx] += *coeff0 * w;
                                    }
                                }
                                y0
                            })
                            .collect(),
                    )
                }
                #[cfg(not(feature = "parallel"))]
                {
                    YMats::Base(
                        M0.iter()
                            .map(|mi| {
                                let mi = mi.as_ref();
                                let mut y0 = vec![R::BaseRing::ZERO; mi.ncols];
                                for (row_idx, row) in mi.coeffs.iter().enumerate() {
                                    let w = eq_at(row_idx);
                                    for (coeff0, col_idx) in row {
                                        y0[*col_idx] += *coeff0 * w;
                                    }
                                }
                                y0
                            })
                            .collect(),
                    )
                }
            }
            ExternalMats::Ring(Mr) => {
                let mats_const = Mr.iter().all(|m| is_const_coeff_sparse_matrix::<R>(m.as_ref()));
                if mats_const {
                    #[cfg(feature = "parallel")]
                    {
                        YMats::Base(
                            Mr.par_iter()
                                .map(|mi| {
                                    let mi = mi.as_ref();
                                    let mut y0 = vec![R::BaseRing::ZERO; mi.ncols];
                                    for (row_idx, row) in mi.coeffs.iter().enumerate() {
                                        let w = eq_at(row_idx);
                                        for (coeff, col_idx) in row {
                                            // Constant-coeff assumption: use only coeffs()[0].
                                            y0[*col_idx] += coeff.coeffs()[0] * w;
                                        }
                                    }
                                    y0
                                })
                                .collect(),
                        )
                    }
                    #[cfg(not(feature = "parallel"))]
                    {
                        YMats::Base(
                            Mr.iter()
                                .map(|mi| {
                                    let mi = mi.as_ref();
                                    let mut y0 = vec![R::BaseRing::ZERO; mi.ncols];
                                    for (row_idx, row) in mi.coeffs.iter().enumerate() {
                                        let w = eq_at(row_idx);
                                        for (coeff, col_idx) in row {
                                            y0[*col_idx] += coeff.coeffs()[0] * w;
                                        }
                                    }
                                    y0
                                })
                                .collect(),
                        )
                    }
                } else {
                    #[cfg(feature = "parallel")]
                    {
                        YMats::Ring(
                            Mr.par_iter()
                                .map(|mi| {
                                    let mi = mi.as_ref();
                                    let mut y = vec![R::ZERO; mi.ncols];
                                    for (row_idx, row) in mi.coeffs.iter().enumerate() {
                                        let w = R::from(eq_at(row_idx));
                                        for (coeff, col_idx) in row {
                                            y[*col_idx] += *coeff * w;
                                        }
                                    }
                                    y
                                })
                                .collect(),
                        )
                    }
                    #[cfg(not(feature = "parallel"))]
                    {
                        YMats::Ring(
                            Mr.iter()
                                .map(|mi| {
                                    let mi = mi.as_ref();
                                    let mut y = vec![R::ZERO; mi.ncols];
                                    for (row_idx, row) in mi.coeffs.iter().enumerate() {
                                        let w = R::from(eq_at(row_idx));
                                        for (coeff, col_idx) in row {
                                            y[*col_idx] += *coeff * w;
                                        }
                                    }
                                    y
                                })
                                .collect(),
                        )
                    }
                }
            }
        };
        if profile {
            println!("[LF+ setchk] step3(y_mats): {:?}", t_y_mats.elapsed());
        }

        let t_e = Instant::now();
        let e: Vec<Vec<Vec<R>>> = {
            let mut e = Vec::with_capacity(1 + mlen);

            // e0:
            // - for dense matrices: e0[m][col] = Σ_row Md[row][col] * eq(row)
            // - for sparse matrices: keep legacy transpose iteration
            #[cfg(feature = "parallel")]
            let mut e0: Vec<Vec<R>> = Vec::with_capacity(Ms_dense.len() + Ms_digits.len() + MTs.len());
            #[cfg(not(feature = "parallel"))]
            let mut e0: Vec<Vec<R>> = Vec::with_capacity(Ms_dense.len() + Ms_digits.len() + MTs.len());

            // Dense sets
            for Md in &Ms_dense {
                #[cfg(feature = "parallel")]
                let v = (0..nrows)
                    .into_par_iter()
                    .fold(
                        || vec![R::ZERO; ncols],
                        |mut acc, row| {
                            let w = R::from(eq_at(row));
                            let row_vals = &Md.vals[row];
                            for col in 0..ncols {
                                acc[col] += row_vals[col] * w;
                            }
                            acc
                        },
                    )
                    // `reduce_with` avoids repeatedly allocating an identity vector.
                    .reduce_with(|mut a, b| {
                        for col in 0..ncols {
                            a[col] += b[col];
                        }
                        a
                    })
                    .unwrap_or_else(|| vec![R::ZERO; ncols]);
                #[cfg(not(feature = "parallel"))]
                let v = (0..ncols)
                    .map(|col| {
                        let mut acc = R::ZERO;
                        for row in 0..nrows {
                            acc += Md.vals[row][col] * R::from(eq_at(row));
                        }
                        acc
                    })
                    .collect::<Vec<_>>();
                e0.push(v);
            }
            // Digit sets
            if !Ms_digits.is_empty() && Ms_digits
                .iter()
                .all(|md| matches!(md.digits, DigitsBacking::ConstCol0 { .. }))
            {
                // Fast fused path for the common SP1 regime:
                // all digit matrices are `ConstCol0`. Then for each matrix we only need:
                // - acc0 = Σ_row exp(digit[row,0]) * eq(row)
                // - common = exp(0) * Σ_row eq(row)
                //
                // Crucially, `Σ_row eq(row)` is the same across all matrices, so we compute it once.
                let md_count = Ms_digits.len();

                #[cfg(feature = "parallel")]
                let sum_w0: R::BaseRing = (0..nrows)
                    .into_par_iter()
                    .map(|row| eq_at(row))
                    .reduce(|| R::BaseRing::ZERO, |a, b| a + b);
                #[cfg(not(feature = "parallel"))]
                let sum_w0: R::BaseRing = (0..nrows).map(|row| eq_at(row)).fold(R::BaseRing::ZERO, |a, b| a + b);

                #[cfg(feature = "parallel")]
                let acc0s: Vec<R> = {
                    (0..nrows)
                        .into_par_iter()
                        .fold(
                            || vec![R::ZERO; md_count],
                            |mut accs, row| {
                                let w0 = eq_at(row);
                                let rw = R::from(w0);
                                for (k, md) in Ms_digits.iter().enumerate() {
                                    let DigitsBacking::ConstCol0 { col0, zero_idx } = &md.digits else { unreachable!() };
                                    let dix = col0.get(row).copied().unwrap_or(*zero_idx) as usize;
                                    accs[k] += md.exp_table[dix] * rw;
                                }
                                accs
                            },
                        )
                        // `reduce_with` avoids repeatedly allocating an identity vector.
                        .reduce_with(|mut a, b| {
                            for k in 0..md_count {
                                a[k] += b[k];
                            }
                            a
                        })
                        .unwrap_or_else(|| vec![R::ZERO; md_count])
                };
                #[cfg(not(feature = "parallel"))]
                let acc0s: Vec<R> = {
                    let mut accs = vec![R::ZERO; md_count];
                    for row in 0..nrows {
                        let w0 = eq_at(row);
                        let rw = R::from(w0);
                        for (k, md) in Ms_digits.iter().enumerate() {
                            let DigitsBacking::ConstCol0 { col0, zero_idx } = &md.digits else { unreachable!() };
                            let dix = col0.get(row).copied().unwrap_or(*zero_idx) as usize;
                            accs[k] += md.exp_table[dix] * rw;
                        }
                    }
                    accs
                };

                for (k, md) in Ms_digits.iter().enumerate() {
                    let DigitsBacking::ConstCol0 { zero_idx, .. } = &md.digits else { unreachable!() };
                    let exp0 = md.exp_table[*zero_idx as usize];
                    let common = exp0 * R::from(sum_w0);
                    let mut out = vec![common; ncols];
                    out[0] = acc0s[k];
                    e0.push(out);
                }
            } else {
                // General path (mixed digit backings)
                for Md in &Ms_digits {
                    // Fast path: if only col0 digits are stored, then for col>0 the matrix entry is
                    // a constant digit (zero_idx) for all rows, so we can avoid the inner `for col` loop.
                    let v = match &Md.digits {
                        DigitsBacking::ConstCol0 { col0, zero_idx } => {
                            let exp0 = Md.exp_table[*zero_idx as usize];
                            let s0 = exp0; // entry(row,col>0)
                            #[cfg(feature = "parallel")]
                            {
                                let (acc0, sum_w0) = (0..nrows)
                                    .into_par_iter()
                                    .fold(
                                        || (R::ZERO, R::BaseRing::ZERO),
                                        |(mut acc0, mut sum_w0), row| {
                                            let w0 = eq_at(row);
                                            sum_w0 += w0;
                                            let dix = col0.get(row).copied().unwrap_or(*zero_idx) as usize;
                                            acc0 += Md.exp_table[dix] * R::from(w0);
                                            (acc0, sum_w0)
                                        },
                                    )
                                    .reduce(
                                        || (R::ZERO, R::BaseRing::ZERO),
                                        |(a0, aw), (b0, bw)| (a0 + b0, aw + bw),
                                    );
                                let common = s0 * R::from(sum_w0);
                                let mut out = vec![common; ncols];
                                out[0] = acc0;
                                out
                            }
                            #[cfg(not(feature = "parallel"))]
                            {
                                let mut acc0 = R::ZERO;
                                let mut sum_w0 = R::BaseRing::ZERO;
                                for row in 0..nrows {
                                    let w0 = eq_at(row);
                                    sum_w0 += w0;
                                    let dix = col0.get(row).copied().unwrap_or(*zero_idx) as usize;
                                    acc0 += Md.exp_table[dix] * R::from(w0);
                                }
                                let common = s0 * R::from(sum_w0);
                                let mut out = vec![common; ncols];
                                out[0] = acc0;
                                out
                            }
                        }
                        DigitsBacking::Full(_) => {
                            #[cfg(feature = "parallel")]
                            {
                                (0..nrows)
                                    .into_par_iter()
                                    .fold(
                                        || vec![R::ZERO; ncols],
                                        |mut acc, row| {
                                            let w = R::from(eq_at(row));
                                            for col in 0..ncols {
                                                acc[col] += Md.get(row, col) * w;
                                            }
                                            acc
                                        },
                                    )
                                    // `reduce_with` avoids repeatedly allocating an identity vector.
                                    .reduce_with(|mut a, b| {
                                        for col in 0..ncols {
                                            a[col] += b[col];
                                        }
                                        a
                                    })
                                    .unwrap_or_else(|| vec![R::ZERO; ncols])
                            }
                            #[cfg(not(feature = "parallel"))]
                            {
                                (0..ncols)
                                    .map(|col| {
                                        let mut acc = R::ZERO;
                                        for row in 0..nrows {
                                            acc += Md.get(row, col) * R::from(eq_at(row));
                                        }
                                        acc
                                    })
                                    .collect::<Vec<_>>()
                            }
                        }
                    };
                    e0.push(v);
                }
            }
            // Sparse sets
            #[cfg(feature = "parallel")]
            {
                let v = MTs
                    .par_iter()
                    .map(|MT| {
                        MT.coeffs
                            .par_iter()
                            .map(|row| {
                                let mut acc = R::ZERO;
                                for &(rij, idx) in row {
                                    acc += rij * R::from(eq_at(idx));
                                }
                                acc
                            })
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<Vec<R>>>();
                e0.extend(v);
            }
            #[cfg(not(feature = "parallel"))]
            {
                let v = MTs
                    .iter()
                    .map(|MT| {
                        MT.coeffs
                            .iter()
                            .map(|row| {
                                let mut acc = R::ZERO;
                                for &(rij, idx) in row {
                                    acc += rij * R::from(eq_at(idx));
                                }
                                acc
                            })
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<Vec<R>>>();
                e0.extend(v);
            }
            e.push(e0);

            // Mf
            for mi in 0..mlen {
                let mut ei: Vec<Vec<R>> = Vec::with_capacity(Ms_dense.len() + Ms_digits.len() + MTs.len());

                // Dense sets: Σ_row Md[row][col] * y[row]
                for Md in &Ms_dense {
                    #[cfg(feature = "parallel")]
                    let v = (0..nrows)
                        .into_par_iter()
                        .fold(
                            || vec![R::ZERO; ncols],
                            |mut acc, row| {
                                let row_vals = &Md.vals[row];
                                let wy = match &y_mats {
                                    YMats::Base(ys) => R::from(ys[mi][row]),
                                    YMats::Ring(ys) => ys[mi][row],
                                };
                                for col in 0..ncols {
                                    acc[col] += row_vals[col] * wy;
                                }
                                acc
                            },
                        )
                        // `reduce_with` avoids repeatedly allocating an identity vector.
                        .reduce_with(|mut a, b| {
                            for col in 0..ncols {
                                a[col] += b[col];
                            }
                            a
                        })
                        .unwrap_or_else(|| vec![R::ZERO; ncols]);
                    #[cfg(not(feature = "parallel"))]
                    let v = (0..ncols)
                        .map(|col| {
                            let mut acc = R::ZERO;
                            for row in 0..nrows {
                                let wy = match &y_mats {
                                    YMats::Base(ys) => R::from(ys[mi][row]),
                                    YMats::Ring(ys) => ys[mi][row],
                                };
                                acc += Md.vals[row][col] * wy;
                            }
                            acc
                        })
                        .collect::<Vec<_>>();
                    ei.push(v);
                }
                // Digit sets
                for Md in &Ms_digits {
                    let v = match &Md.digits {
                        DigitsBacking::ConstCol0 { col0, zero_idx } => match &y_mats {
                            YMats::Base(ys) => {
                                let exp0 = Md.exp_table[*zero_idx as usize];
                                let s0 = exp0;
                                #[cfg(feature = "parallel")]
                                {
                                    let (acc0, sum_wy0) = (0..nrows)
                                        .into_par_iter()
                                        .fold(
                                            || (R::ZERO, R::BaseRing::ZERO),
                                            |(mut acc0, mut sum_wy0), row| {
                                                let wy0 = ys[mi][row];
                                                sum_wy0 += wy0;
                                                let dix = col0.get(row).copied().unwrap_or(*zero_idx) as usize;
                                                acc0 += Md.exp_table[dix] * R::from(wy0);
                                                (acc0, sum_wy0)
                                            },
                                        )
                                        .reduce(
                                            || (R::ZERO, R::BaseRing::ZERO),
                                            |(a0, aw), (b0, bw)| (a0 + b0, aw + bw),
                                        );
                                    let common = s0 * R::from(sum_wy0);
                                    let mut out = vec![common; ncols];
                                    out[0] = acc0;
                                    out
                                }
                                #[cfg(not(feature = "parallel"))]
                                {
                                    let mut acc0 = R::ZERO;
                                    let mut sum_wy0 = R::BaseRing::ZERO;
                                    for row in 0..nrows {
                                        let wy0 = ys[mi][row];
                                        sum_wy0 += wy0;
                                        let dix = col0.get(row).copied().unwrap_or(*zero_idx) as usize;
                                        acc0 += Md.exp_table[dix] * R::from(wy0);
                                    }
                                    let common = s0 * R::from(sum_wy0);
                                    let mut out = vec![common; ncols];
                                    out[0] = acc0;
                                    out
                                }
                            }
                            YMats::Ring(ys) => {
                                let exp0 = Md.exp_table[*zero_idx as usize];
                                let s0 = exp0;
                                #[cfg(feature = "parallel")]
                                {
                                    let (acc0, sum_wy) = (0..nrows)
                                        .into_par_iter()
                                        .fold(
                                            || (R::ZERO, R::ZERO),
                                            |(mut acc0, mut sum_wy), row| {
                                                let wy = ys[mi][row];
                                                sum_wy += wy;
                                                let dix = col0.get(row).copied().unwrap_or(*zero_idx) as usize;
                                                acc0 += Md.exp_table[dix] * wy;
                                                (acc0, sum_wy)
                                            },
                                        )
                                        .reduce(|| (R::ZERO, R::ZERO), |(a0, aw), (b0, bw)| (a0 + b0, aw + bw));
                                    let common = s0 * sum_wy;
                                    let mut out = vec![common; ncols];
                                    out[0] = acc0;
                                    out
                                }
                                #[cfg(not(feature = "parallel"))]
                                {
                                    let mut acc0 = R::ZERO;
                                    let mut sum_wy = R::ZERO;
                                    for row in 0..nrows {
                                        let wy = ys[mi][row];
                                        sum_wy += wy;
                                        let dix = col0.get(row).copied().unwrap_or(*zero_idx) as usize;
                                        acc0 += Md.exp_table[dix] * wy;
                                    }
                                    let common = s0 * sum_wy;
                                    let mut out = vec![common; ncols];
                                    out[0] = acc0;
                                    out
                                }
                            }
                        },
                        DigitsBacking::Full(_) => {
                            #[cfg(feature = "parallel")]
                            {
                                (0..nrows)
                                    .into_par_iter()
                                    .fold(
                                        || vec![R::ZERO; ncols],
                                        |mut acc, row| {
                                            let wy = match &y_mats {
                                                YMats::Base(ys) => R::from(ys[mi][row]),
                                                YMats::Ring(ys) => ys[mi][row],
                                            };
                                            for col in 0..ncols {
                                                acc[col] += Md.get(row, col) * wy;
                                            }
                                            acc
                                        },
                                    )
                                    // `reduce_with` avoids repeatedly allocating an identity vector.
                                    .reduce_with(|mut a, b| {
                                        for col in 0..ncols {
                                            a[col] += b[col];
                                        }
                                        a
                                    })
                                    .unwrap_or_else(|| vec![R::ZERO; ncols])
                            }
                            #[cfg(not(feature = "parallel"))]
                            {
                                (0..ncols)
                                    .map(|col| {
                                        let mut acc = R::ZERO;
                                        for row in 0..nrows {
                                            let wy = match &y_mats {
                                                YMats::Base(ys) => R::from(ys[mi][row]),
                                                YMats::Ring(ys) => ys[mi][row],
                                            };
                                            acc += Md.get(row, col) * wy;
                                        }
                                        acc
                                    })
                                    .collect::<Vec<_>>()
                            }
                        }
                    };
                    ei.push(v);
                }

                // Sparse sets (legacy)
                #[cfg(feature = "parallel")]
                {
                    let v = MTs
                        .par_iter()
                        .map(|MT| {
                            MT.coeffs
                                .par_iter()
                                .map(|row| {
                                    let mut acc = R::ZERO;
                                    for &(rij, idx) in row {
                                        let wy = match &y_mats {
                                            YMats::Base(ys) => R::from(ys[mi][idx]),
                                            YMats::Ring(ys) => ys[mi][idx],
                                        };
                                        acc += rij * wy;
                                    }
                                    acc
                                })
                                .collect::<Vec<_>>()
                        })
                        .collect::<Vec<Vec<R>>>();
                    ei.extend(v);
                }
                #[cfg(not(feature = "parallel"))]
                {
                    let v = MTs
                        .iter()
                        .map(|MT| {
                            MT.coeffs
                                .iter()
                                .map(|row| {
                                    let mut acc = R::ZERO;
                                    for &(rij, idx) in row {
                                        let wy = match &y_mats {
                                            YMats::Base(ys) => R::from(ys[mi][idx]),
                                            YMats::Ring(ys) => ys[mi][idx],
                                        };
                                        acc += rij * wy;
                                    }
                                    acc
                                })
                                .collect::<Vec<_>>()
                        })
                        .collect::<Vec<Vec<R>>>();
                    ei.extend(v);
                }
                e.push(ei);
            }
            e
        };
        if profile {
            println!("[LF+ setchk] step3(e): {:?}", t_e.elapsed());
        }

        let t_b = Instant::now();
        #[cfg(feature = "parallel")]
        let b: Vec<R> = ms
            .par_iter()
            .map(|mset| match mset {
                // NOTE: `ms.len()` is often 1, so parallelizing over `ms` doesn't help.
                // Parallelize over the vector length (nrows) instead.
                VecSet::Dense(m) => (0..m.len())
                    .into_par_iter()
                    .map(|i| m[i] * R::from(eq_at(i)))
                    .reduce(|| R::ZERO, |a, b| a + b),
                VecSet::Digits { digits, exp_table } => (0..digits.len())
                    .into_par_iter()
                    .map(|i| exp_table[digits[i] as usize] * R::from(eq_at(i)))
                    .reduce(|| R::ZERO, |a, b| a + b),
            })
            .collect();
        #[cfg(not(feature = "parallel"))]
        let b: Vec<R> = ms
            .iter()
            .map(|mset| {
                let mut acc = R::ZERO;
                match mset {
                    VecSet::Dense(m) => {
                        for (i, &mi) in m.iter().enumerate() {
                            acc += mi * R::from(eq_at(i));
                        }
                    }
                    VecSet::Digits { digits, exp_table } => {
                        for (i, &dix) in digits.iter().enumerate() {
                            acc += exp_table[dix as usize] * R::from(eq_at(i));
                        }
                    }
                }
                acc
            })
            .collect();
        if profile {
            println!("[LF+ setchk] step3(b): {:?}", t_b.elapsed());
        }

        let t_absorb = std::time::Instant::now();
        // Prover to Verifier messages
        absorb_evaluations(&e, &b, transcript);
        if profile {
            println!("[LF+ setchk] step3(absorb): {:?}", t_absorb.elapsed());
        }
        maybe_print_rss("setchk: step3 done");

        if profile {
            println!(
                "[LF+ setchk] step3(e,b)+absorb: {:?}  total: {:?}",
                t_step3.elapsed(),
                t_total.elapsed()
            );
        }

        Out {
            nvars: self.nvars,
            e,
            b,
            r,
            sumcheck_proof,
        }
    }
}

#[inline]
fn beta_pows<R: PolyRing>(beta: R::BaseRing) -> Vec<R::BaseRing>
where
    R::BaseRing: Ring,
{
    let d = R::dimension();
    let mut out = Vec::with_capacity(d);
    let mut acc = R::BaseRing::ONE;
    for _ in 0..d {
        out.push(acc);
        acc *= beta;
    }
    out
}

/// Fast `ev(r, beta)`:
/// - if `r` is monomial-like (<=1 nonzero coeff), do O(1) lookup via `beta_pows`
/// - otherwise fall back to full dot product against `beta_pows`
#[inline]
fn ev_fast<R: PolyRing>(r: &R, beta_pows: &[R::BaseRing]) -> R::BaseRing
where
    R::BaseRing: Ring,
{
    let coeffs = r.coeffs();
    debug_assert_eq!(coeffs.len(), beta_pows.len());

    let mut idx: Option<usize> = None;
    let mut c: R::BaseRing = R::BaseRing::ZERO;
    for (i, &ci) in coeffs.iter().enumerate() {
        if ci != R::BaseRing::ZERO {
            if idx.is_some() {
                // not monomial
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

/// Build eq(bits(idx), r) table (little-endian index order), matching latticefold build_eq_x_r_vec.
#[allow(dead_code)]
fn build_eq_table_base<R: PolyRing>(r: &[R::BaseRing]) -> Vec<R::BaseRing>
where
    R::BaseRing: Ring,
{
    let mut buf = vec![R::BaseRing::ONE];
    for &ri in r.iter().rev() {
        let mut res = vec![R::BaseRing::ZERO; buf.len() << 1];
        for (i, out) in res.iter_mut().enumerate() {
            let bi = buf[i >> 1];
            let tmp = ri * bi;
            *out = if (i & 1) == 0 { bi - tmp } else { tmp };
        }
        buf = res;
    }
    buf
}

#[inline]
fn choose_t_low(nvars: usize) -> usize {
    // Keep a tiny 2^t table; stream the remaining high bits.
    nvars.min(12)
}

/// Build eq weights for the low `t` variables (LSB-first), returned as a length-2^t table.
fn build_eq_low_table<R: PolyRing>(
    r_low: &[R::BaseRing],
    one_minus_r_low: &[R::BaseRing],
) -> Vec<R::BaseRing>
where
    R::BaseRing: Ring,
{
    debug_assert_eq!(r_low.len(), one_minus_r_low.len());
    let t = r_low.len();
    let mut buf = vec![R::BaseRing::ONE];
    for i in (0..t).rev() {
        let ri = r_low[i];
        let omi = one_minus_r_low[i];
        let mut res = vec![R::BaseRing::ZERO; buf.len() << 1];
        for (j, out) in res.iter_mut().enumerate() {
            let bi = buf[j >> 1];
            *out = if (j & 1) == 0 { bi * omi } else { bi * ri };
        }
        buf = res;
    }
    buf
}

/// Precompute scale factors for the high bits (t..nvars-1): scale_high[high] = Π_i (bit? r[i] : (1-r[i])).
fn build_scale_high<R: PolyRing>(
    r: &[R::BaseRing],
    one_minus_r: &[R::BaseRing],
    t_low: usize,
) -> Vec<R::BaseRing>
where
    R::BaseRing: Ring,
{
    debug_assert_eq!(r.len(), one_minus_r.len());
    let nvars = r.len();
    let high_bits = nvars - t_low;
    let high_len = 1usize << high_bits;
    let mut out = vec![R::BaseRing::ONE; high_len];
    for h in 0..high_len {
        let mut prod = R::BaseRing::ONE;
        for i in t_low..nvars {
            let bit = ((h >> (i - t_low)) & 1) == 1;
            prod *= if bit { r[i] } else { one_minus_r[i] };
        }
        out[h] = prod;
    }
    out
}

#[inline]
fn eq_weight_base<R: PolyRing>(
    idx: usize,
    low: &[R::BaseRing],
    scale_high: &[R::BaseRing],
    low_mask: usize,
    t_low: usize,
) -> R::BaseRing
where
    R::BaseRing: Ring,
{
    let low_idx = idx & low_mask;
    let high_idx = idx >> t_low;
    scale_high[high_idx] * low[low_idx]
}

impl<R: OverField> Out<R> {
    pub fn verify(&self, transcript: &mut impl Transcript<R>) -> Result<(), SetCheckError<R>> {
        let nclaims = self.e[0].len() + self.b.len();

        let cba: Vec<(Vec<R>, R::BaseRing, R::BaseRing)> = (0..nclaims)
            .map(|_| {
                let c: Vec<R> = transcript
                    .get_challenges(self.nvars)
                    .into_iter()
                    .map(|x| x.into())
                    .collect();
                let beta = transcript.get_challenge();
                let alpha = transcript.get_challenge();
                (c, beta, alpha)
            })
            .collect();

        let rc: Option<R::BaseRing> = (self.e[0].len() > 1).then(|| transcript.get_challenge());

        let subclaim = MLSumcheck::verify_as_subprotocol(
            transcript,
            self.nvars,
            3,
            R::zero(),
            &self.sumcheck_proof,
        )?;

        let r: Vec<R> = subclaim.point.into_iter().map(|x| x.into()).collect();

        let v = subclaim.expected_evaluation;

        // Prover to Verifier messages
        absorb_evaluations(&self.e, &self.b, transcript);

        use ark_std::One;
        let mut ver = R::zero();
        // Avoid pow() in verification as well: use precomputed alpha powers and iterative rc powers.
        let mut rc_pow: R::BaseRing = R::BaseRing::one();
        for (i, e) in self.e[0].iter().enumerate() {
            let c = &cba[i].0;
            let beta = &cba[i].1;
            let alpha = &cba[i].2;
            let eq = eq_eval(c, &r).unwrap();
            let mut alpha_pows: Vec<R::BaseRing> = Vec::with_capacity(e.len());
            {
                let mut ap = R::BaseRing::one(); // alpha^0
                for _ in 0..e.len() {
                    alpha_pows.push(ap);
                    ap *= *alpha;
                }
            }
            let e_sum = e
                .iter()
                .enumerate()
                .map(|(j, e_j)| {
                    let ev1 = R::from(ev(e_j, *beta));
                    let ev2 = R::from(ev(e_j, *beta * beta));
                    (ev1 * ev1 - ev2) * alpha_pows[j]
                })
                .sum::<R>();
            ver += eq * e_sum * R::from(rc_pow);
            if let Some(rc) = rc {
                rc_pow *= rc;
            }
        }
        for (i, b) in self.b.iter().enumerate() {
            let offset = self.e[0].len();
            let c = &cba[i + offset].0;
            let beta = &cba[i + offset].1;
            let alpha = &cba[i + offset].2;
            let eq = eq_eval(c, &r).unwrap();
            let b_claim = {
                let ev1 = R::from(ev(b, *beta));
                let ev2 = R::from(ev(b, *beta * *beta));
                ev1 * ev1 - ev2
            };
            ver += eq * *alpha * b_claim * R::from(rc_pow);
            if let Some(rc) = rc {
                rc_pow *= rc;
            }
        }

        (ver == v)
            .then_some(())
            .ok_or(SetCheckError::ExpectedEvaluation(ver, v))?;

        Ok(())
    }
}

fn absorb_evaluations<R: OverField>(
    e: &[Vec<Vec<R>>],
    b: &[R],
    transcript: &mut impl Transcript<R>,
) {
    for ek in e {
        for ej in ek {
            transcript.absorb_slice(ej);
        }
    }
    transcript.absorb_slice(b);
}

#[cfg(test)]
mod tests {
    use ark_std::One;
    use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
    use stark_rings::{cyclotomic_ring::models::frog_ring::RqPoly as R, unit_monomial};
    use stark_rings_linalg::SparseMatrix;
    use std::sync::Arc;

    use super::*;
    use crate::transcript::PoseidonTranscript;

    #[test]
    fn test_set_check() {
        let n = 4;
        let M = SparseMatrix::<R>::identity(n);

        let scin = In {
            sets: vec![MonomialSet::Matrix(M)],
            nvars: log2(n) as usize,
        };

        let mut ts = PoseidonTranscript::empty::<PC>();
        let empty_m: Vec<Arc<SparseMatrix<R>>> = Vec::new();
        let out = scin.set_check(ExternalMats::Ring(&empty_m), &mut ts);

        let mut ts = PoseidonTranscript::empty::<PC>();
        out.verify(&mut ts).unwrap();
    }

    #[test]
    fn test_set_check_bad() {
        let n = 4;
        let mut M = SparseMatrix::<R>::identity(n);
        // 1 + X, not a monomial
        let mut onepx = R::one();
        onepx.coeffs_mut()[1] = 1u128.into();
        M.coeffs[0][0].0 = onepx;

        let scin = In {
            sets: vec![MonomialSet::Matrix(M)],
            nvars: log2(n) as usize,
        };

        let mut ts = PoseidonTranscript::empty::<PC>();
        let empty_m: Vec<Arc<SparseMatrix<R>>> = Vec::new();
        let out = scin.set_check(ExternalMats::Ring(&empty_m), &mut ts);

        let mut ts = PoseidonTranscript::empty::<PC>();
        assert!(out.verify(&mut ts).is_err());
    }

    #[test]
    fn test_set_check_batched() {
        let n = 4;
        let M0 = SparseMatrix::<R>::identity(n);
        let M1 = SparseMatrix::<R>::identity(n);

        let scin = In {
            sets: vec![MonomialSet::Matrix(M0), MonomialSet::Matrix(M1)],
            nvars: log2(n) as usize,
        };

        let mut ts = PoseidonTranscript::empty::<PC>();
        let empty_m: Vec<Arc<SparseMatrix<R>>> = Vec::new();
        let out = scin.set_check(ExternalMats::Ring(&empty_m), &mut ts);

        let mut ts = PoseidonTranscript::empty::<PC>();
        out.verify(&mut ts).unwrap();
    }

    #[test]
    fn test_set_check_batched_bad() {
        let n = 4;
        let M0 = SparseMatrix::<R>::identity(n);
        let mut M1 = SparseMatrix::<R>::identity(n);
        // 1 + X, not a monomial
        let mut onepx = R::one();
        onepx.coeffs_mut()[1] = 1u128.into();
        M1.coeffs[0][0].0 = onepx;

        let scin = In {
            sets: vec![MonomialSet::Matrix(M0), MonomialSet::Matrix(M1)],
            nvars: log2(n) as usize,
        };

        let mut ts = PoseidonTranscript::empty::<PC>();
        let empty_m: Vec<Arc<SparseMatrix<R>>> = Vec::new();
        let out = scin.set_check(ExternalMats::Ring(&empty_m), &mut ts);

        let mut ts = PoseidonTranscript::empty::<PC>();
        assert!(out.verify(&mut ts).is_err());
    }

    #[test]
    fn test_set_check_mix() {
        let n = 4;
        let M0 = SparseMatrix::<R>::identity(n);
        let M1 = SparseMatrix::<R>::identity(n);
        let m0 = Arc::new(vec![R::one(); n]);
        let m1 = Arc::new(vec![unit_monomial(2); n]);

        let scin = In {
            sets: vec![
                MonomialSet::Matrix(M0),
                MonomialSet::Matrix(M1),
                MonomialSet::Vector(m0),
                MonomialSet::Vector(m1),
            ],
            nvars: log2(n) as usize,
        };

        let mut ts = PoseidonTranscript::empty::<PC>();
        let empty_m: Vec<Arc<SparseMatrix<R>>> = Vec::new();
        let out = scin.set_check(ExternalMats::Ring(&empty_m), &mut ts);

        let mut ts = PoseidonTranscript::empty::<PC>();
        out.verify(&mut ts).unwrap();
    }

    #[test]
    fn test_set_check_mix_bad() {
        let n = 4;
        let M0 = SparseMatrix::<R>::identity(n);
        let M1 = SparseMatrix::<R>::identity(n);
        let mut m0 = vec![R::one(); n];
        let mut onepx = R::one();
        onepx.coeffs_mut()[1] = 1u128.into();
        m0[0] = onepx;
        let m0 = Arc::new(m0);

        let scin = In {
            sets: vec![
                MonomialSet::Matrix(M0),
                MonomialSet::Matrix(M1),
                MonomialSet::Vector(m0),
            ],
            nvars: log2(n) as usize,
        };

        let mut ts = PoseidonTranscript::empty::<PC>();
        let empty_m: Vec<Arc<SparseMatrix<R>>> = Vec::new();
        let out = scin.set_check(ExternalMats::Ring(&empty_m), &mut ts);

        let mut ts = PoseidonTranscript::empty::<PC>();
        assert!(out.verify(&mut ts).is_err());
    }

    #[test]
    fn test_set_check_digits_bad() {
        // Ensure the digit-backed (oracle) path rejects non-monomial entries.
        //
        // This is important for the SP1 regime where we heavily rely on DigitsMatrix/VectorDigits
        // representations for performance.
        let n = 8usize;
        let nvars = ark_std::log2(n) as usize;
        let d = R::dimension();

        // exp_table[0] = 1 (monomial), exp_table[1] = 1 + X (NOT a monomial).
        let exp0 = R::one();
        let mut exp1 = R::one();
        exp1.coeffs_mut()[1] = 1u128.into();
        let exp_table: Arc<Vec<R>> = Arc::new(vec![exp0, exp1]);

        // Put the bad entry in col=0 at one row; all other cols are treated as exp(0).
        let mut col0 = vec![0u16; n];
        col0[3] = 1u16;

        let dm = DigitsMatrix::<R> {
            nrows: n,
            ncols: d,
            digits: DigitsBacking::ConstCol0 { col0: Arc::new(col0), zero_idx: 0u16 },
            exp_table,
        };

        let scin = In {
            sets: vec![MonomialSet::DigitsMatrix(Arc::new(dm))],
            nvars,
        };

        let mut ts = PoseidonTranscript::empty::<PC>();
        let empty_m: Vec<Arc<SparseMatrix<R>>> = Vec::new();
        let out = scin.set_check(ExternalMats::Ring(&empty_m), &mut ts);

        let mut ts = PoseidonTranscript::empty::<PC>();
        assert!(out.verify(&mut ts).is_err());
    }
}
