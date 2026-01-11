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

#[cfg(feature = "parallel")]
use rayon::prelude::*;

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
            StreamingMleEnum::DenseMatrixColEv { num_vars, .. } => *num_vars,
            StreamingMleEnum::DigitsMatrixColEv { num_vars, .. } => *num_vars,
            StreamingMleEnum::EqBase { r, .. } => r.len(),
            StreamingMleEnum::SparseMatVec { num_vars, .. } => *num_vars,
            StreamingMleEnum::SparseMatVecConstCoeff { num_vars, .. } => *num_vars,
            StreamingMleEnum::Tensor4Padded { num_vars, .. } => *num_vars,
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
        }
    }

    #[inline]
    pub fn fix_variable_in_place_base(&mut self, r0: R::BaseRing) {
        let nv = self.num_vars();
        assert!(nv > 0);
        let half = 1usize << (nv - 1);
        let one_minus0 = R::BaseRing::ONE - r0;
        let r_ring = R::from(r0);
        match self {
            StreamingMleEnum::DenseOwned { evals, num_vars } => {
                let one_minus = R::ONE - r_ring;
                for i in 0..half {
                    let a = evals[i << 1];
                    let b = evals[(i << 1) | 1];
                    evals[i] = one_minus * a + r_ring * b;
                }
                evals.truncate(half);
                *num_vars -= 1;
            }
            StreamingMleEnum::DenseArc { .. } => {
                let next = self.fix_variable(r_ring);
                *self = next;
            }
            StreamingMleEnum::BaseScalarOwned { evals, num_vars } => {
                for i in 0..half {
                    let a = evals[i << 1];
                    let b = evals[(i << 1) | 1];
                    evals[i] = one_minus0 * a + r0 * b;
                }
                evals.truncate(half);
                *num_vars -= 1;
            }
            StreamingMleEnum::BaseScalarArc {
                evals,
                num_vars,
                square,
            } => {
                // Take ownership of the Arc if possible; otherwise clone.
                let arc = std::mem::take(evals);
                let mut owned = match Arc::try_unwrap(arc) {
                    Ok(v) => v,
                    Err(a) => (*a).clone(),
                };
                if *square {
                    // Vertex-wise squares: square BEFORE combining.
                    for i in 0..half {
                        let mut a = owned[i << 1];
                        let mut b = owned[(i << 1) | 1];
                        a *= a;
                        b *= b;
                        owned[i] = one_minus0 * a + r0 * b;
                    }
                } else {
                    for i in 0..half {
                        let a = owned[i << 1];
                        let b = owned[(i << 1) | 1];
                        owned[i] = one_minus0 * a + r0 * b;
                    }
                }
                owned.truncate(half);
                // After fixing, the table is now the correct MLE values (square semantics consumed).
                *self = StreamingMleEnum::BaseScalarOwned {
                    evals: owned,
                    num_vars: *num_vars - 1,
                };
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
            StreamingMleEnum::Tensor4Padded { .. } => {
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
            StreamingMleEnum::Tensor4Padded { .. } => {
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
        } else {
            assert!(state.round == 0);
        }

        state.round += 1;
        assert!(state.round <= state.num_vars);

        let nv = state.mles[0].num_vars();
        let degree = state.max_degree;
        let domain_half = 1usize << (nv - 1);
        let num_polys = state.mles.len();

        struct Scratch<Rr> {
            evals: Vec<Rr>,
            steps: Vec<Rr>,
            vals0: Vec<Rr>,
            vals1: Vec<Rr>,
            vals: Vec<Rr>,
            levals: Vec<Rr>,
        }

        let scratch = || Scratch {
            evals: vec![R::ZERO; degree + 1],
            steps: vec![R::ZERO; num_polys],
            vals0: vec![R::ZERO; num_polys],
            vals1: vec![R::ZERO; num_polys],
            vals: vec![R::ZERO; num_polys],
            levals: vec![R::ZERO; degree + 1],
        };

        #[cfg(feature = "parallel")]
        let result = (0..domain_half)
            .into_par_iter()
            .fold(scratch, |mut s, b| {
                for (i, mle) in state.mles.iter().enumerate() {
                    s.vals0[i] = mle.eval_at_index(b << 1);
                    s.vals1[i] = mle.eval_at_index((b << 1) | 1);
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
            .reduce(|| vec![R::ZERO; degree + 1], |mut acc, evals| {
                for (a, e) in acc.iter_mut().zip(evals) {
                    *a += e;
                }
                acc
            });

        #[cfg(not(feature = "parallel"))]
        let result = {
            let mut acc = vec![R::ZERO; degree + 1];
            let mut s = scratch();
            for b in 0..domain_half {
                for (i, mle) in state.mles.iter().enumerate() {
                    s.vals0[i] = mle.eval_at_index(b << 1);
                    s.vals1[i] = mle.eval_at_index((b << 1) | 1);
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
        let evals0 = (0..domain_half)
            .into_par_iter()
            .map(|b| {
                let mut vals0 = vec![R::BaseRing::ZERO; num_polys];
                let mut vals1 = vec![R::BaseRing::ZERO; num_polys];
                let mut steps = vec![R::BaseRing::ZERO; num_polys];
                let mut vals = vec![R::BaseRing::ZERO; num_polys];
                let mut out = vec![R::BaseRing::ZERO; degree + 1];

                for (i, mle) in state.mles.iter().enumerate() {
                    vals0[i] = mle.eval0_at_index(b << 1);
                    vals1[i] = mle.eval0_at_index((b << 1) | 1);
                }
                out[0] = comb_fn0(&vals0);
                out[1] = comb_fn0(&vals1);
                for i in 0..num_polys {
                    steps[i] = vals1[i] - vals0[i];
                    vals[i] = vals1[i];
                }
                for d in 2..=degree {
                    for i in 0..num_polys {
                        vals[i] += steps[i];
                    }
                    out[d] = comb_fn0(&vals);
                }
                out
            })
            .reduce(|| vec![R::BaseRing::ZERO; degree + 1], |mut acc, v| {
                for (a, e) in acc.iter_mut().zip(v) {
                    *a += e;
                }
                acc
            });

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
        let profile_detail = std::env::var("LF_PLUS_PROFILE_DETAIL").ok().as_deref() == Some("1");
        let t_total = std::time::Instant::now();

        transcript.absorb_field_element(&R::BaseRing::from(nvars as u128));
        transcript.absorb_field_element(&R::BaseRing::from(degree as u128));

        let t_init = std::time::Instant::now();
        let mut state = Self::prover_init(mles, nvars, degree);
        if profile && profile_detail {
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

            if profile && profile_detail && (round == 0 || round + 1 == nvars) {
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

        if profile && profile_detail {
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
}

