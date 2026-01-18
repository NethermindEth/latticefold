//! Streaming Sumcheck Prover for Symphony.
//!
//! This implements the sumcheck protocol without materializing O(2^n) evaluation tables.
//! Evaluations are computed on-demand via index-based lookup functions.
//!
//! Per HackMD (<https://hackmd.io/@zkpunk/SkleumdxZl>), the key insight is that
//! structured MLEs can be evaluated at hypercube indices without storing all values.

use ark_std::vec::Vec;
use rayon::prelude::*;
use std::sync::Arc;
use stark_rings::{OverField, PolyRing, Ring};
use stark_rings_linalg::SparseMatrix;

use latticefold::transcript::Transcript;

/// Concrete streaming MLE variants (no dynamic dispatch in hot loops).
#[derive(Clone)]
pub enum StreamingMleEnum<R: OverField + PolyRing>
where
    R::BaseRing: Ring,
{
    Dense(DenseStreamingMle<R>),
    /// Arc-wrapped dense table - avoids cloning when sharing across multiple consumers.
    DenseArc {
        evals: Arc<Vec<R>>,
        num_vars: usize,
    },
    /// y[row] = (M * w)[row], computed from sparse rows (only used before the first fix).
    SparseMatVec {
        matrix: Arc<SparseMatrix<R>>,
        witness: Arc<Vec<R>>,
        num_vars: usize,
    },
    /// Constant-coefficient specialization of `SparseMatVec`.
    ///
    /// When the matrix and witness are embedded base-field scalars inside the ring (i.e., all
    /// non-constant coefficients are zero), the product `M*w` stays constant-coefficient too.
    /// We can then evaluate this MLE purely in the base ring and lift to `R` at the end.
    ///
    /// This is expected for many applications (including SP1-style R1CS), and matches the
    /// “base-field embedded in the ring” setting used in the Symphony paper streaming notes.
    SparseMatVecConstCoeff {
        matrix: Arc<SparseMatrix<R>>,
        witness0: Arc<Vec<R::BaseRing>>,
        num_vars: usize,
    },
    /// A base-ring scalar table (stored as BaseRing), lifted to R via R::from.
    BaseScalarVec {
        evals: Arc<Vec<R::BaseRing>>,
        num_vars: usize,
    },
    /// Owned base-scalar table used after the first fix so we can fix variables in-place
    /// (avoids allocating a fresh `Vec` every round).
    BaseScalarVecOwned {
        evals: Vec<R::BaseRing>,
        num_vars: usize,
    },
    /// A compact base-scalar table that is periodic in the "row" dimension.
///
    /// This is tailored for Symphony's monomial side, where `m_j` values repeat across `rep`
    /// blocks, so we can store only `m_j * d` base scalars instead of `m * d`.
    ///
    /// Storage layout is column-major: `evals[col * m_j + row]`.
    PeriodicBaseScalarVec {
        evals: Arc<Vec<R::BaseRing>>,
        num_vars: usize,
        /// Current row count for the full domain (power-of-two).
        m: usize,
        /// Current projected row count (power-of-two, <= m).
        m_j: usize,
        /// Current number of columns (power-of-two).
        d: usize,
        /// If true, this MLE represents **vertex-wise squares**:
        /// it evaluates \(\mathrm{MLE}(m^{\circ 2})(r)\), where \(m^{\circ 2}\) denotes squaring
        /// each hypercube vertex value before multilinear interpolation.
        ///
        /// This is intentionally **not** the same as \(\mathrm{MLE}(m)(r)^2\) at non-Boolean points,
        /// and the non-commutation is exactly what monomial/set-check style claims rely on.
        ///
        /// Implementation detail: we square at lookup for unfixed variables, and when fixing variables
        /// we square the paired vertex values *before* combining (see `fix_variable`).
        square: bool,
    },
    /// eq(bits(index), r) evaluated in the base ring, then lifted to R.
    EqBase {
        scale: R::BaseRing,
        r: Vec<R::BaseRing>,
        one_minus_r: Vec<R::BaseRing>,
    },
    /// Seeded Ajtai column: a_ρ[j] = ⟨ρ, A[:,j]⟩ computed on-demand via PRG.
    /// 
    /// This avoids materializing the O(n) table upfront. Instead, each `eval_at_index`
    /// call computes the value via κ PRG evaluations.
    /// 
    /// **Const-coeff optimization**: `rho_ring` stores precomputed `R::from(rho0[i])` to
    /// avoid repeated lifting in the hot loop. The multiplication is still O(d) because
    /// `rho_ring[i]` is a constant-coefficient ring element.
    /// 
    /// Uses Arc to avoid per-digit cloning (κ is small, but Arc is cleaner).
    /// 
    /// After the first `fix_variable`, this transitions to `Dense` (since we need to
    /// actually combine paired values, which requires materializing).
    AjtaiRho {
        domain: Vec<u8>,
        seed: [u8; 32],
        kappa: usize,
        /// Precomputed ring elements: rho_ring[i] = R::from(rho0[i])
        rho_ring: Arc<Vec<R>>,
        num_vars: usize,
    },
}

impl<R: OverField + PolyRing> StreamingMleEnum<R>
where
    R::BaseRing: Ring,
{
    pub fn dense(evals: Vec<R>) -> Self {
        Self::Dense(DenseStreamingMle::new(evals))
    }

    /// Create a dense MLE from an Arc-wrapped vector (avoids cloning for shared data).
    pub fn dense_arc(evals: Arc<Vec<R>>) -> Self {
        let len = evals.len();
        assert!(len.is_power_of_two(), "Length must be power of 2");
        let num_vars = len.trailing_zeros() as usize;
        Self::DenseArc { evals, num_vars }
    }

    pub fn sparse_mat_vec(matrix: Arc<SparseMatrix<R>>, witness: Arc<Vec<R>>) -> Self {
        let nrows = matrix.nrows;
        assert!(nrows.is_power_of_two(), "nrows must be power-of-two");
        let num_vars = nrows.trailing_zeros() as usize;
        Self::SparseMatVec {
            matrix,
            witness,
            num_vars,
        }
    }

    pub fn sparse_mat_vec_const_coeff(matrix: Arc<SparseMatrix<R>>, witness0: Arc<Vec<R::BaseRing>>) -> Self {
        let nrows = matrix.nrows;
        assert!(nrows.is_power_of_two(), "nrows must be power-of-two");
        let num_vars = nrows.trailing_zeros() as usize;
        Self::SparseMatVecConstCoeff {
            matrix,
            witness0,
            num_vars,
        }
    }

    pub fn base_scalar_vec(num_vars: usize, evals: Arc<Vec<R::BaseRing>>) -> Self {
        Self::BaseScalarVec { evals, num_vars }
    }

    #[inline]
    pub fn base_scalar_vec_owned(num_vars: usize, evals: Vec<R::BaseRing>) -> Self {
        Self::BaseScalarVecOwned { evals, num_vars }
    }

    pub fn periodic_base_scalar_vec(
        num_vars: usize,
        m: usize,
        m_j: usize,
        d: usize,
        evals: Arc<Vec<R::BaseRing>>,
        square: bool,
    ) -> Self {
        debug_assert!(m.is_power_of_two());
        debug_assert!(m_j.is_power_of_two());
        debug_assert!(d.is_power_of_two());
        debug_assert!(m_j <= m);
        debug_assert_eq!(evals.len(), m_j * d);
        Self::PeriodicBaseScalarVec {
            evals,
            num_vars,
            m,
            m_j,
            d,
            square,
        }
    }

    pub fn eq_base(r: Vec<R::BaseRing>) -> Self {
        let one_minus_r = r.iter().copied().map(|x| R::BaseRing::ONE - x).collect();
        Self::EqBase {
            scale: R::BaseRing::ONE,
            r,
            one_minus_r,
        }
    }

    /// Create a streaming MLE for `a_ρ[j] = ⟨ρ, A[:,j]⟩` where A is a seeded Ajtai matrix.
    /// 
    /// This computes values on-demand via PRG, avoiding O(n) memory allocation upfront.
    /// After the first variable fix, it transitions to dense (required for combining pairs).
    /// 
    /// **Const-coeff optimization**: takes `rho0: Arc<Vec<R::BaseRing>>` and precomputes
    /// `rho_ring = R::from(rho0[i])` once to avoid repeated lifting in `eval_at_index`.
    /// Arc avoids per-digit cloning.
    pub fn ajtai_rho(domain: Vec<u8>, seed: [u8; 32], kappa: usize, rho0: Arc<Vec<R::BaseRing>>, num_vars: usize) -> Self {
        debug_assert_eq!(rho0.len(), kappa);
        // Precompute rho_ring once here (κ is small, ~8 elements)
        let rho_ring: Vec<R> = rho0.iter().map(|&r| R::from(r)).collect();
        Self::AjtaiRho {
            domain,
            seed,
            kappa,
            rho_ring: Arc::new(rho_ring),
            num_vars,
        }
    }

    #[inline]
    pub fn num_vars(&self) -> usize {
        match self {
            StreamingMleEnum::Dense(m) => m.num_vars,
            StreamingMleEnum::DenseArc { num_vars, .. } => *num_vars,
            StreamingMleEnum::SparseMatVec { num_vars, .. } => *num_vars,
            StreamingMleEnum::SparseMatVecConstCoeff { num_vars, .. } => *num_vars,
            StreamingMleEnum::BaseScalarVec { num_vars, .. } => *num_vars,
            StreamingMleEnum::BaseScalarVecOwned { num_vars, .. } => *num_vars,
            StreamingMleEnum::PeriodicBaseScalarVec { num_vars, .. } => *num_vars,
            StreamingMleEnum::EqBase { r, .. } => r.len(),
            StreamingMleEnum::AjtaiRho { num_vars, .. } => *num_vars,
        }
    }

    #[inline]
    pub fn eval_at_index(&self, index: usize) -> R {
        match self {
            StreamingMleEnum::Dense(m) => m.evals[index],
            StreamingMleEnum::DenseArc { evals, .. } => evals[index],
            StreamingMleEnum::SparseMatVec {
                matrix, witness, ..
            } => {
                let mut sum = R::ZERO;
                for (coeff, col_idx) in &matrix.coeffs[index] {
                    if *col_idx < witness.len() {
                        sum += *coeff * witness[*col_idx];
                    }
                }
                sum
            }
            StreamingMleEnum::SparseMatVecConstCoeff {
                matrix, witness0, ..
            } => {
                let mut sum = R::BaseRing::ZERO;
                for (coeff, col_idx) in &matrix.coeffs[index] {
                    if *col_idx < witness0.len() {
                        // Safety: caller guarantees constant-coefficient embedding.
                        sum += coeff.coeffs()[0] * witness0[*col_idx];
                    }
                }
                R::from(sum)
            }
            StreamingMleEnum::BaseScalarVec { evals, .. } => R::from(evals[index]),
            StreamingMleEnum::BaseScalarVecOwned { evals, .. } => R::from(evals[index]),
            StreamingMleEnum::PeriodicBaseScalarVec {
                evals,
                m,
                m_j,
                d,
                square,
                ..
            } => {
                debug_assert_eq!(*m * *d, 1usize << self.num_vars());
                let col = index / *m;
                let r = index % *m;
                let row = if *m_j == 0 { 0 } else { r % *m_j };
                let mut v = evals[col * *m_j + row];
                if *square {
                    v *= v;
                }
                R::from(v)
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
            StreamingMleEnum::AjtaiRho {
                domain,
                seed,
                kappa,
                rho_ring,
                ..
            } => {
                // Compute a_ρ[index] = ⟨ρ, A[:,index]⟩ on-demand via PRG
                // rho_ring[row] is const-coeff (= R::from(scalar)), so this multiplication
                // is inherently O(d) not O(d²) in coefficient complexity
                use rand::SeedableRng;
                use rand_chacha::ChaCha20Rng;
                use latticefold::commitment::AjtaiCommitmentScheme;
                
                let col_seed = AjtaiCommitmentScheme::<R>::derive_col_seed(domain, seed, index as u64);
                let mut rng = ChaCha20Rng::from_seed(col_seed);
                
                let mut sum = R::ZERO;
                for row in 0..*kappa {
                    let aij = R::rand(&mut rng);
                    // const-coeff × ring: even if ring impl doesn't optimize,
                    // the convolution naturally skips non-contributing terms
                    sum += rho_ring[row] * aij;
                }
                sum
            }
        }
    }

    /// Evaluate this MLE at `index`, returning only the constant coefficient (base ring element).
    ///
    /// Used by the constant-coefficient Π_had fast path to avoid constructing full ring elements
    /// in the hot sumcheck loop.
    #[inline]
    pub fn eval0_at_index(&self, index: usize) -> R::BaseRing {
        match self {
            StreamingMleEnum::Dense(m) => m.evals[index].coeffs()[0],
            StreamingMleEnum::DenseArc { evals, .. } => evals[index].coeffs()[0],
            StreamingMleEnum::SparseMatVec { .. } => self.eval_at_index(index).coeffs()[0],
            StreamingMleEnum::SparseMatVecConstCoeff {
                matrix, witness0, ..
            } => {
                let mut sum = R::BaseRing::ZERO;
                for (coeff, col_idx) in &matrix.coeffs[index] {
                    if *col_idx < witness0.len() {
                        sum += coeff.coeffs()[0] * witness0[*col_idx];
                    }
                }
                sum
            }
            StreamingMleEnum::BaseScalarVec { evals, .. } => evals[index],
            StreamingMleEnum::BaseScalarVecOwned { evals, .. } => evals[index],
            StreamingMleEnum::PeriodicBaseScalarVec {
                evals,
                m,
                m_j,
                d: _d,
                square,
                ..
            } => {
                let col = index / *m;
                let r = index % *m;
                let row = if *m_j == 0 { 0 } else { r % *m_j };
                let mut v = evals[col * *m_j + row];
                if *square {
                    v *= v;
                }
                v
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
            StreamingMleEnum::AjtaiRho { .. } => {
                // AjtaiRho produces full ring elements, take constant coefficient
                self.eval_at_index(index).coeffs()[0]
            }
        }
    }

    /// Fix the next variable (LSB-first) and return a new MLE over one fewer variable.
    pub fn fix_variable(&self, r: R) -> StreamingMleEnum<R> {
        let nv = self.num_vars();
        assert!(nv > 0);
        let half = 1 << (nv - 1);
        match self {
            // Dense: materialize smaller dense.
            StreamingMleEnum::Dense(m) => {
                let new_evals: Vec<R> = (0..half)
                    .map(|i| (R::ONE - r) * m.evals[i << 1] + r * m.evals[(i << 1) | 1])
                    .collect();
                StreamingMleEnum::dense(new_evals)
            }
            // DenseArc: same as Dense, but read from Arc
            StreamingMleEnum::DenseArc { evals, .. } => {
                let new_evals: Vec<R> = (0..half)
                    .map(|i| (R::ONE - r) * evals[i << 1] + r * evals[(i << 1) | 1])
                    .collect();
                StreamingMleEnum::dense(new_evals)
            }
            // Sparse mat-vec: materialize to dense once fixed (avoids exponential recursion).
            StreamingMleEnum::SparseMatVec { .. } => {
                let new_evals: Vec<R> = (0..half)
                    .map(|i| {
                        let f0 = self.eval_at_index(i << 1);
                        let f1 = self.eval_at_index((i << 1) | 1);
                        (R::ONE - r) * f0 + r * f1
                    })
                    .collect();
                StreamingMleEnum::dense(new_evals)
            }
            // Constant-coeff sparse mat-vec: keep it base-scalar after fixing.
            StreamingMleEnum::SparseMatVecConstCoeff { num_vars, .. } => {
                let r0 = r.coeffs()[0];
                let one_minus = R::BaseRing::ONE - r0;
                let new_evals: Vec<R::BaseRing> = (0..half)
                    .map(|i| {
                        // Avoid constructing full ring elements; evaluate directly in base ring.
                        let f0 = self.eval0_at_index(i << 1);
                        let f1 = self.eval0_at_index((i << 1) | 1);
                        one_minus * f0 + r0 * f1
                    })
                    .collect();
                StreamingMleEnum::BaseScalarVecOwned { evals: new_evals, num_vars: num_vars - 1 }
            }
            // Base-scalar table: keep it base-scalar after fixing (critical for memory).
            //
            // Materializing this into `Vec<R>` blows up memory by ~dim(R)× because `R` stores `d`
            // base-ring coefficients. We only need constant-coefficient values here.
            StreamingMleEnum::BaseScalarVec { evals, num_vars } => {
                let r0 = r.coeffs()[0];
                let one_minus = R::BaseRing::ONE - r0;
                let new_evals: Vec<R::BaseRing> = (0..half)
                    .map(|i| one_minus * evals[i << 1] + r0 * evals[(i << 1) | 1])
                    .collect();
                StreamingMleEnum::BaseScalarVecOwned { evals: new_evals, num_vars: num_vars - 1 }
            }
            StreamingMleEnum::BaseScalarVecOwned { evals, num_vars } => {
                let r0 = r.coeffs()[0];
                let one_minus = R::BaseRing::ONE - r0;
                let mut evals = evals.clone();
                for i in 0..half {
                    let a = evals[i << 1];
                    let b = evals[(i << 1) | 1];
                    evals[i] = one_minus * a + r0 * b;
                }
                evals.truncate(half);
                StreamingMleEnum::BaseScalarVecOwned {
                    evals,
                    num_vars: num_vars - 1,
                }
            }
            StreamingMleEnum::PeriodicBaseScalarVec {
                evals,
                num_vars,
                m,
                m_j,
                d,
                square,
            } => {
                let r0 = r.coeffs()[0];
                let one_minus = R::BaseRing::ONE - r0;

                // We are fixing variables LSB-first. As long as m > 1, the LSB is a row bit.
                // Once m == 1, remaining bits are column bits.
                if *m > 1 {
                    let new_m = *m / 2;
                    // If m_j > 1, the periodicity depends on this bit and we must compress m_j too.
                    // If m_j == 1, the function is already independent of row bits.
                    if *m_j > 1 {
                        let new_mj = *m_j / 2;
                        let mut new_evals = vec![R::BaseRing::ZERO; new_mj * *d];
                        for col in 0..*d {
                            let base = col * *m_j;
                            let out_base = col * new_mj;
                            for i in 0..new_mj {
                                let mut a = evals[base + (i << 1)];
                                let mut b = evals[base + (i << 1) | 1];
                                // IMPORTANT: if this MLE represents vertex-wise squares, we must
                                // fix variables on the squared values (not square after fixing).
                                if *square {
                                    a *= a;
                                    b *= b;
                                }
                                new_evals[out_base + i] = one_minus * a + r0 * b;
                            }
                        }
                        StreamingMleEnum::PeriodicBaseScalarVec {
                            evals: Arc::new(new_evals),
                            num_vars: num_vars - 1,
                            m: new_m,
                            m_j: new_mj,
                            d: *d,
                            square: false,
                        }
                    } else {
                        // Independent of this row bit: combining equal values yields same table.
                        StreamingMleEnum::PeriodicBaseScalarVec {
                            evals: if *square {
                                // Materialize squared values once so future fixes are correct.
                                Arc::new(evals.iter().map(|x| (*x) * (*x)).collect())
                            } else {
                                evals.clone()
                            },
                            num_vars: num_vars - 1,
                            m: new_m,
                            m_j: *m_j,
                            d: *d,
                            square: false,
                        }
                    }
                } else {
                    // Column bits.
                    if *d == 1 {
                        // Fully fixed after this; treat as Dense of length 1.
                        let mut v = evals[0];
                        if *square { v *= v; }
                        StreamingMleEnum::base_scalar_vec(0, Arc::new(vec![v]))
                    } else {
                        let new_d = *d / 2;
                        let mut new_evals = vec![R::BaseRing::ZERO; *m_j * new_d];
                        // m_j should be 1 once m==1 in our use-case, but keep it generic.
                        for col in 0..new_d {
                            for row in 0..*m_j {
                                let mut a = evals[(col << 1) * *m_j + row];
                                let mut b = evals[((col << 1) | 1) * *m_j + row];
                                if *square {
                                    a *= a;
                                    b *= b;
                                }
                                new_evals[col * *m_j + row] = one_minus * a + r0 * b;
                            }
                        }
                        StreamingMleEnum::PeriodicBaseScalarVec {
                            evals: Arc::new(new_evals),
                            num_vars: num_vars - 1,
                            m: *m,
                            m_j: *m_j,
                            d: new_d,
                            square: false,
                        }
                    }
                }
            }
            // EqBase: keep structural (base-field) with updated scale and shortened r vector.
            StreamingMleEnum::EqBase {
                scale,
                r: rr,
                one_minus_r,
            } => {
                let r0 = r.coeffs()[0];
                let eq_factor = (R::BaseRing::ONE - r0) * one_minus_r[0] + r0 * rr[0];
                let new_scale = *scale * eq_factor;
                let new_r = rr[1..].to_vec();
                let new_om = one_minus_r[1..].to_vec();
                StreamingMleEnum::EqBase {
                    scale: new_scale,
                    r: new_r,
                    one_minus_r: new_om,
                }
            }
            // AjtaiRho: materialize to dense after first fix (no nice structural property)
            StreamingMleEnum::AjtaiRho { .. } => {
                // Materialize all values, then combine pairs
                let new_evals: Vec<R> = (0..half)
                    .map(|i| {
                        let v0 = self.eval_at_index(i << 1);
                        let v1 = self.eval_at_index((i << 1) | 1);
                        (R::ONE - r) * v0 + r * v1
                    })
                    .collect();
                StreamingMleEnum::dense(new_evals)
            }
        }
    }

    /// Fix the next variable (LSB-first) in-place for variants that can reuse buffers.
    pub fn fix_variable_in_place_base(&mut self, r0: R::BaseRing) {
        let nv = self.num_vars();
        assert!(nv > 0);
        let half = 1 << (nv - 1);
        let one_minus0 = R::BaseRing::ONE - r0;
        let r_ring = R::from(r0);
        match self {
            StreamingMleEnum::Dense(m) => {
                let one_minus = R::ONE - r_ring;
                for i in 0..half {
                    let a = m.evals[i << 1];
                    let b = m.evals[(i << 1) | 1];
                    m.evals[i] = one_minus * a + r_ring * b;
                }
                m.evals.truncate(half);
                m.num_vars -= 1;
            }
            StreamingMleEnum::DenseArc { .. } => {
                // Can't modify Arc in-place; fall back to fix_variable
                let new = self.fix_variable(r_ring);
                *self = new;
            }
            StreamingMleEnum::SparseMatVec { .. } => {
                let new = self.fix_variable(r_ring);
                *self = new;
            }
            StreamingMleEnum::SparseMatVecConstCoeff { .. } => {
                // Fall back to `fix_variable` (which already does the right thing); we overwrite
                // `self` afterward so future fixes can be in-place.
                let next = self.fix_variable(r_ring);
                *self = next;
            }
            StreamingMleEnum::BaseScalarVec { evals, num_vars } => {
                // Take ownership of the Arc, then try to unwrap; fall back to cloning if shared.
                let arc = std::mem::take(evals);
                let nv = *num_vars;
                let mut owned = match Arc::try_unwrap(arc) {
                    Ok(v) => v,
                    Err(a) => (*a).clone(),
                };
                for i in 0..half {
                    let a = owned[i << 1];
                    let b = owned[(i << 1) | 1];
                    owned[i] = one_minus0 * a + r0 * b;
                }
                owned.truncate(half);
                *self = StreamingMleEnum::BaseScalarVecOwned {
                    evals: owned,
                    num_vars: nv - 1,
                };
            }
            StreamingMleEnum::BaseScalarVecOwned { evals, num_vars } => {
                for i in 0..half {
                    let a = evals[i << 1];
                    let b = evals[(i << 1) | 1];
                    evals[i] = one_minus0 * a + r0 * b;
                }
                evals.truncate(half);
                *num_vars -= 1;
            }
            StreamingMleEnum::PeriodicBaseScalarVec { .. } => {
                let new = self.fix_variable(r_ring);
                *self = new;
            }
            StreamingMleEnum::EqBase {
                scale,
                r,
                one_minus_r,
            } => {
                let eq_factor = one_minus0 * one_minus_r[0] + r0 * r[0];
                *scale *= eq_factor;
                // Keep order; lengths are small (<= ~24), so shifting cost is negligible.
                r.remove(0);
                one_minus_r.remove(0);
            }
            StreamingMleEnum::AjtaiRho { .. } => {
                // Materialize to dense after first fix
                let new = self.fix_variable(r_ring);
                *self = new;
            }
        }
    }
}

// ============================================================================
// Concrete Streaming MLE Implementations
// ============================================================================

/// Dense MLE (fallback, stores all evaluations).
#[derive(Clone)]
pub struct DenseStreamingMle<R: OverField> {
    evals: Vec<R>,
    num_vars: usize,
}

impl<R: OverField> DenseStreamingMle<R> {
    pub fn new(evals: Vec<R>) -> Self {
        let len = evals.len();
        assert!(len.is_power_of_two(), "Length must be power of 2");
        let num_vars = len.trailing_zeros() as usize;
        Self { evals, num_vars }
}

    pub fn num_vars(&self) -> usize {
        self.num_vars
    }

    pub fn eval_at_index(&self, index: usize) -> R {
        self.evals[index]
    }
}

// ============================================================================
// Streaming Sumcheck Prover
// ============================================================================

/// Prover message for streaming sumcheck.
#[derive(Clone, Debug)]
pub struct StreamingProverMsg<R: OverField> {
    /// Evaluations at 0, 1, 2, ..., degree
    pub evaluations: Vec<R>,
}

/// Streaming sumcheck proof.
#[derive(Clone, Debug)]
pub struct StreamingProof<R: OverField>(pub Vec<StreamingProverMsg<R>>);

/// Streaming sumcheck prover state.
pub struct StreamingSumcheckState<R: OverField + PolyRing>
where
    R::BaseRing: Ring,
{
    mles: Vec<StreamingMleEnum<R>>,
    pub randomness: Vec<R::BaseRing>,
    num_vars: usize,
    max_degree: usize,
    round: usize,
}

impl<R: OverField + PolyRing> StreamingSumcheckState<R>
where
    R::BaseRing: Ring,
{
    /// Number of variables remaining in the internal MLEs.
    #[inline]
    pub fn remaining_vars(&self) -> usize {
        self.mles[0].num_vars()
    }

    /// Apply one additional verifier challenge by fixing the next variable in all internal MLEs.
    ///
    /// This is useful when the caller samples the final challenge but does not call `prove_round`
    /// again (e.g., shared-schedule drivers that end immediately after the last transcript squeeze).
    pub fn fix_last_variable(&mut self, r: R::BaseRing) {
        let nv = self.remaining_vars();
        assert!(
            nv == 1,
            "fix_last_variable expects exactly 1 remaining var, got {nv}"
        );
        for m in self.mles.iter_mut() {
            m.fix_variable_in_place_base(r);
        }
    }

    /// Evaluate all internal MLEs at the fully fixed point (after all variables are fixed).
    pub fn final_evals(&self) -> Vec<R> {
        let nv = self.remaining_vars();
        assert!(
            nv == 0,
            "final_evals expects 0 remaining vars (fully fixed), got {nv}"
        );
        self.mles.iter().map(|m| m.eval_at_index(0)).collect()
    }
}

/// Streaming sumcheck prover.
pub struct StreamingSumcheck;

impl StreamingSumcheck {
    /// Initialize prover with streaming MLEs.
    pub fn prover_init<R: OverField + PolyRing>(
        mles: Vec<StreamingMleEnum<R>>,
        nvars: usize,
        degree: usize,
    ) -> StreamingSumcheckState<R> {
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

    /// Prove one round of sumcheck.
    pub fn prove_round<R: OverField + PolyRing>(
        state: &mut StreamingSumcheckState<R>,
        v_msg: Option<R::BaseRing>,
        comb_fn: &(dyn Fn(&[R]) -> R + Sync + Send),
    ) -> StreamingProverMsg<R> {
        // Handle previous round's randomness
        if let Some(r) = v_msg {
            assert!(state.round > 0);
            state.randomness.push(r);

            // Fix variable in all MLEs
            for m in state.mles.iter_mut() {
                m.fix_variable_in_place_base(r);
            }
        } else {
            assert!(state.round == 0);
        }

        state.round += 1;
        assert!(state.round <= state.num_vars);

        let nv = state.mles[0].num_vars();
        let degree = state.max_degree;
        let domain_half = 1 << (nv - 1);
        let num_polys = state.mles.len();

        // Compute round polynomial evaluations at 0, 1, ..., degree
        // Using the Jolt-style optimization: for each b in {0,1}^{nv-1},
        // evaluate combiner at (0, b) and (1, b), then interpolate for higher degrees.

        struct Scratch<R> {
            evals: Vec<R>,
            steps: Vec<R>,
            vals0: Vec<R>,
            vals1: Vec<R>,
            vals: Vec<R>,
            levals: Vec<R>,
        }

        let scratch = || Scratch {
            evals: vec![R::ZERO; degree + 1],
            steps: vec![R::ZERO; num_polys],
            vals0: vec![R::ZERO; num_polys],
            vals1: vec![R::ZERO; num_polys],
            vals: vec![R::ZERO; num_polys],
            levals: vec![R::ZERO; degree + 1],
        };

        let result = (0..domain_half)
            .into_par_iter()
            .fold(scratch, |mut s, b| {
                // Evaluate all MLEs at indices 2*b (for x_i=0) and 2*b+1 (for x_i=1)
                for (i, mle) in state.mles.iter().enumerate() {
                    s.vals0[i] = mle.eval_at_index(b << 1);
                    s.vals1[i] = mle.eval_at_index((b << 1) | 1);
                }

                s.levals[0] = comb_fn(&s.vals0);
                s.levals[1] = comb_fn(&s.vals1);

                // Compute steps for interpolation: step[i] = vals1[i] - vals0[i]
                for i in 0..num_polys {
                    s.steps[i] = s.vals1[i] - s.vals0[i];
                    s.vals[i] = s.vals1[i];
                }

                // Extrapolate to higher degrees
                for d in 2..=degree {
                    for i in 0..num_polys {
                        s.vals[i] += s.steps[i];
                    }
                    s.levals[d] = comb_fn(&s.vals);
                }

                // Accumulate
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

        StreamingProverMsg { evaluations: result }
    }

    /// Prove one round of sumcheck, specialized for the constant-coefficient (base ring) case.
    ///
    /// The combiner operates over base-ring scalars and the produced univariate evaluations are
    /// lifted back into the ring as constant-coefficient elements.
    pub fn prove_round_base<R: OverField + PolyRing>(
        state: &mut StreamingSumcheckState<R>,
        v_msg: Option<R::BaseRing>,
        comb_fn0: &(dyn Fn(&[R::BaseRing]) -> R::BaseRing + Sync + Send),
    ) -> StreamingProverMsg<R>
    where
        R::BaseRing: Ring,
    {
        if let Some(r) = v_msg {
            assert!(state.round > 0);
            state.randomness.push(r);
            for m in state.mles.iter_mut() {
                m.fix_variable_in_place_base(r);
            }
        } else {
            assert!(state.round == 0);
        }

        state.round += 1;
        assert!(state.round <= state.num_vars);

        let nv = state.mles[0].num_vars();
        let degree = state.max_degree;
        let domain_half = 1 << (nv - 1);
        let num_polys = state.mles.len();

        struct Scratch<F> {
            evals: Vec<F>,
            steps: Vec<F>,
            vals0: Vec<F>,
            vals1: Vec<F>,
            vals: Vec<F>,
            levals: Vec<F>,
        }
        let scratch = || Scratch {
            evals: vec![R::BaseRing::ZERO; degree + 1],
            steps: vec![R::BaseRing::ZERO; num_polys],
            vals0: vec![R::BaseRing::ZERO; num_polys],
            vals1: vec![R::BaseRing::ZERO; num_polys],
            vals: vec![R::BaseRing::ZERO; num_polys],
            levals: vec![R::BaseRing::ZERO; degree + 1],
        };

        let result0 = (0..domain_half)
            .into_par_iter()
            .fold(scratch, |mut s, b| {
                for (i, mle) in state.mles.iter().enumerate() {
                    s.vals0[i] = mle.eval0_at_index(b << 1);
                    s.vals1[i] = mle.eval0_at_index((b << 1) | 1);
                }

                s.levals[0] = comb_fn0(&s.vals0);
                s.levals[1] = comb_fn0(&s.vals1);

                for i in 0..num_polys {
                    s.steps[i] = s.vals1[i] - s.vals0[i];
                    s.vals[i] = s.vals1[i];
                }

                for d in 2..=degree {
                    for i in 0..num_polys {
                        s.vals[i] += s.steps[i];
                    }
                    s.levals[d] = comb_fn0(&s.vals);
                }

                for (e, l) in s.evals.iter_mut().zip(s.levals.iter()) {
                    *e += *l;
                }
                s
            })
            .map(|s| s.evals)
            .reduce(|| vec![R::BaseRing::ZERO; degree + 1], |mut acc, evals| {
                for (a, e) in acc.iter_mut().zip(evals) {
                    *a += e;
                }
                acc
            });

        StreamingProverMsg {
            evaluations: result0.into_iter().map(R::from).collect(),
        }
    }

    /// Run streaming sumcheck as subprotocol.
    pub fn prove_as_subprotocol<R: OverField + PolyRing, T: Transcript<R>>(
        transcript: &mut T,
        mles: Vec<StreamingMleEnum<R>>,
        nvars: usize,
        degree: usize,
        comb_fn: impl Fn(&[R]) -> R + Sync + Send,
    ) -> (StreamingProof<R>, Vec<R::BaseRing>) {
        // Use absorb_field_element for scalars (saves ~(d-1) base-field absorbs each).
        transcript.absorb_field_element(&R::BaseRing::from(nvars as u128));
        transcript.absorb_field_element(&R::BaseRing::from(degree as u128));

        let mut state = Self::prover_init(mles, nvars, degree);
        let mut msgs = Vec::with_capacity(nvars);
        let mut v_msg: Option<R::BaseRing> = None;

        for _ in 0..nvars {
            let pm = Self::prove_round(&mut state, v_msg, &comb_fn);
            transcript.absorb_slice(&pm.evaluations);
            msgs.push(pm);

            let r = transcript.get_challenge();
            transcript.absorb_field_element(&r);  // scalar absorb
            v_msg = Some(r);
        }

        state.randomness.push(v_msg.unwrap());
        (StreamingProof(msgs), state.randomness)
    }

    /// Run streaming sumcheck with hook after `hook_round` challenges.
    pub fn prove_as_subprotocol_with_hook<R: OverField + PolyRing, T: Transcript<R>>(
        transcript: &mut T,
        mles: Vec<StreamingMleEnum<R>>,
        nvars: usize,
        degree: usize,
        comb_fn: impl Fn(&[R]) -> R + Sync + Send,
        hook_round: usize,
        mut hook: impl FnMut(&mut T, &[R::BaseRing]),
    ) -> (StreamingProof<R>, Vec<R::BaseRing>) {
        assert!(hook_round <= nvars);

        // Use absorb_field_element for scalars (saves ~(d-1) base-field absorbs each).
        transcript.absorb_field_element(&R::BaseRing::from(nvars as u128));
        transcript.absorb_field_element(&R::BaseRing::from(degree as u128));

        let mut state = Self::prover_init(mles, nvars, degree);
        let mut msgs = Vec::with_capacity(nvars);
        let mut v_msg: Option<R::BaseRing> = None;
        let mut sampled: Vec<R::BaseRing> = Vec::with_capacity(nvars);

        for _ in 0..nvars {
            let pm = Self::prove_round(&mut state, v_msg, &comb_fn);
            transcript.absorb_slice(&pm.evaluations);
            msgs.push(pm);

            let r = transcript.get_challenge();
            transcript.absorb_field_element(&r);  // scalar absorb
            sampled.push(r);

            if hook_round != 0 && sampled.len() == hook_round {
                hook(transcript, &sampled);
            }

            v_msg = Some(r);
        }

        state.randomness.push(v_msg.unwrap());
        (StreamingProof(msgs), state.randomness)
    }

    /// Run TWO streaming sumchecks with shared randomness (Symphony schedule).
    pub fn prove_two_shared_with_hook<R: OverField + PolyRing, T: Transcript<R>>(
        transcript: &mut T,
        // A
        mles_a: Vec<StreamingMleEnum<R>>,
        nvars_a: usize,
        degree_a: usize,
        comb_fn_a: impl Fn(&[R]) -> R + Sync + Send,
        // B
        mles_b: Vec<StreamingMleEnum<R>>,
        nvars_b: usize,
        degree_b: usize,
        comb_fn_b: impl Fn(&[R]) -> R + Sync + Send,
        // hook
        hook_round: usize,
        mut hook: impl FnMut(&mut T, &[R::BaseRing]),
    ) -> (
        (StreamingProof<R>, Vec<R::BaseRing>),
        (StreamingProof<R>, Vec<R::BaseRing>),
    ) {
        assert!(nvars_a > 0 && nvars_b > 0);
        assert!(hook_round <= nvars_b);

        // Use absorb_field_element for scalars (saves ~(d-1) base-field absorbs each).
        transcript.absorb_field_element(&R::BaseRing::from(nvars_a as u128));
        transcript.absorb_field_element(&R::BaseRing::from(degree_a as u128));
        transcript.absorb_field_element(&R::BaseRing::from(nvars_b as u128));
        transcript.absorb_field_element(&R::BaseRing::from(degree_b as u128));

        let mut state_a = Self::prover_init(mles_a, nvars_a, degree_a);
        let mut state_b = Self::prover_init(mles_b, nvars_b, degree_b);

        let mut msgs_a = Vec::with_capacity(nvars_a);
        let mut msgs_b = Vec::with_capacity(nvars_b);
        let mut sampled: Vec<R::BaseRing> = Vec::with_capacity(nvars_a.max(nvars_b));
        let mut v_msg: Option<R::BaseRing> = None;

        let rounds = nvars_a.max(nvars_b);
        for round_idx in 0..rounds {
            // Prover messages
            if round_idx < nvars_a {
                let pm_a = Self::prove_round(&mut state_a, v_msg, &comb_fn_a);
                transcript.absorb_slice(&pm_a.evaluations);
                msgs_a.push(pm_a);
            }
            if round_idx < nvars_b {
                let pm_b = Self::prove_round(&mut state_b, v_msg, &comb_fn_b);
                transcript.absorb_slice(&pm_b.evaluations);
                msgs_b.push(pm_b);
            }

            // Sample shared randomness
            let r = transcript.get_challenge();
            transcript.absorb_field_element(&r);  // scalar absorb
            sampled.push(r);

            if hook_round != 0 && sampled.len() == hook_round {
                hook(transcript, &sampled);
            }

            v_msg = Some(r);
        }

        // Push final randomness
        let last_a = sampled[nvars_a - 1];
        let last_b = sampled[nvars_b - 1];
        state_a.randomness.push(last_a);
        state_b.randomness.push(last_b);

        (
            (StreamingProof(msgs_a), state_a.randomness),
            (StreamingProof(msgs_b), state_b.randomness),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_streaming_sumcheck_matches_dense_on_linear_combiner() {
        use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as R;
        use stark_rings::PolyRing;
        use stark_rings::Ring;
        use latticefold::utils::sumcheck::prover::ProverMsg;
        use latticefold::utils::sumcheck::IPForMLSumcheck;

        // 3 variables => 8 evals
        let nvars = 3usize;
        let len = 1usize << nvars;
        let f1 = (0..len).map(|i| R::from((i + 1) as u128)).collect::<Vec<_>>();
        let f2 = (0..len)
            .map(|i| R::from(((i * 7 + 3) % 97) as u128))
            .collect::<Vec<_>>();

        let m1 = StreamingMleEnum::dense(f1.clone());
        let m2 = StreamingMleEnum::dense(f2.clone());

        let comb = |vals: &[R]| -> R { vals[0] * R::from(13u128) + vals[1] * R::from(29u128) };

        // Streaming: run the protocol directly with fixed challenges (BaseRing elements).
        type K = <R as PolyRing>::BaseRing;
        let challenges = vec![K::from(11u128), K::from(22u128), K::from(33u128)];

        // Asserted sum over {0,1}^n
        let mut st = StreamingSumcheck::prover_init(vec![m1, m2], nvars, 1);
        let mut v_msg = None;
        let mut msgs = Vec::new();
        for r in &challenges {
            let pm = StreamingSumcheck::prove_round(&mut st, v_msg, &comb);
            msgs.push(pm);
            v_msg = Some(*r);
        }

        // Verify the sumcheck identity on the produced univariates:
        // at each round, g_i(0)+g_i(1) must equal the previous claim.
        let mut claim = {
            let mut acc = R::ZERO;
            for idx in 0..len {
                acc += comb(&[f1[idx], f2[idx]]);
            }
            acc
        };
        for (round, pm) in msgs.iter().enumerate() {
            assert_eq!(
                pm.evaluations[0] + pm.evaluations[1],
                claim,
                "round {round} failed"
            );
            // next claim is g_i(r_i) (degree=1 so interpolate linearly)
            let r = R::from(challenges[round]);
            claim = (R::ONE - r) * pm.evaluations[0] + r * pm.evaluations[1];
        }

        // Stronger check: feed the same prover messages into latticefold's *real* verifier state
        // with the chosen randomness, and ensure the produced subclaim matches evaluation at that point.
        let mut vs = IPForMLSumcheck::<R, crate::transcript::PoseidonTranscript<R>>::verifier_init(nvars, 1);
        for (pm, r) in msgs.iter().zip(challenges.iter().copied()) {
            IPForMLSumcheck::<R, crate::transcript::PoseidonTranscript<R>>::verify_round_with_randomness(
                ProverMsg::new(pm.evaluations.clone()),
                &mut vs,
                r,
            );
        }
        let sub = IPForMLSumcheck::<R, crate::transcript::PoseidonTranscript<R>>::check_and_generate_subclaim(vs, {
            let mut acc = R::ZERO;
            for idx in 0..len {
                acc += comb(&[f1[idx], f2[idx]]);
            }
            acc
        })
        .expect("verifier should accept correct proof");

        // Compute expected evaluation by fixing the same point in the prover state.
        //
        // IMPORTANT: In the sumcheck schedule, the verifier samples `r_i` *after* the prover sends
        // the i-th univariate. The prover only applies `r_i` at the *start of the next round*.
        // So after `nvars` rounds, `st.mles` has only applied r_0..r_{nvars-2}; we must still
        // apply the last challenge r_{nvars-1} to reach the full evaluation point.
        let last_r = challenges[nvars - 1];
        let fixed_last = st
            .mles
            .iter()
            .map(|m| m.fix_variable(R::from(last_r)))
            .collect::<Vec<_>>();
        let vals_at_point: Vec<R> = fixed_last.iter().map(|m| m.eval_at_index(0)).collect();
        let expected = comb(&vals_at_point);
        assert_eq!(sub.expected_evaluation, expected, "subclaim mismatch");
    }

    #[test]
    fn test_dense_streaming_mle() {
        use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as R;
        let evals = vec![R::from(1u128), R::from(2u128), R::from(3u128), R::from(4u128)];
        let mle = DenseStreamingMle::new(evals.clone());
        assert_eq!(mle.num_vars(), 2);
        for i in 0..4 {
            assert_eq!(mle.eval_at_index(i), evals[i]);
        }
    }

    #[test]
    fn test_periodic_base_scalar_vec_square_is_vertex_square_not_eval_square() {
        // Pins the intended meaning of the `square` flag:
        //   square=false => evaluates MLE(m)(r)
        //   square=true  => evaluates MLE(m^{∘2})(r), i.e. squares *vertex values* before interpolation.
        //
        // For a non-Boolean point r, we generally have:
        //   MLE(m)(r)^2 != MLE(m^{∘2})(r)
        // and this non-commutation is exactly what monomial/set-check tests rely on.
        use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as R;
        use stark_rings::PolyRing;

        type F = <R as PolyRing>::BaseRing;

        // One variable: m=2 rows, d=1 column => g_len = m*d = 2 => nvars=1.
        let m = 2usize;
        let m_j = 2usize;
        let d = 1usize;
        let nvars = 1usize;

        // Vertex values: [a, b]
        let a = F::from(1u128);
        let b = F::from(2u128);
        let evals = Arc::new(vec![a, b]); // length m_j*d = 2

        let mle = StreamingMleEnum::<R>::periodic_base_scalar_vec(nvars, m, m_j, d, evals.clone(), false);
        let mle_sq = StreamingMleEnum::<R>::periodic_base_scalar_vec(nvars, m, m_j, d, evals, true);

        // Non-Boolean point.
        let r = F::from(2u128);

        // Evaluate by fixing variables (LSB-first).
        let eval0 = |mut m: StreamingMleEnum<R>| -> F {
            m.fix_variable_in_place_base(r);
            assert_eq!(m.num_vars(), 0);
            m.eval0_at_index(0)
        };

        let v = eval0(mle);
        let v_sq = eval0(mle_sq);

        // Expected:
        // MLE([a,b])(r) = (1-r)*a + r*b
        // MLE([a^2,b^2])(r) = (1-r)*a^2 + r*b^2
        let one_minus = F::ONE - r;
        let expected = one_minus * a + r * b;
        let expected_sq = one_minus * (a * a) + r * (b * b);

        assert_eq!(v, expected);
        assert_eq!(v_sq, expected_sq);

        // And in general, these differ at non-Boolean r.
        assert_ne!(v_sq, v * v);
    }
}
