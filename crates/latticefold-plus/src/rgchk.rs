use ark_std::iter::once;
use latticefold::transcript::Transcript;
use latticefold::commitment::AjtaiCommitmentScheme;
use stark_rings::{
    balanced_decomposition::Decompose,
    exp, psi, CoeffRing, OverField, PolyRing, Ring, Zq,
};
use stark_rings_linalg::{Matrix, SparseMatrix};
use std::sync::Arc;
use thiserror::Error;

use crate::{
    setchk::{DigitsMatrix, In, MonomialSet, Out},
    utils::split,
};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[inline]
fn digit_lookup_or_panic<BR>(
    digit_elems: &[BR],
    dig: BR,
    digit_abs_max: i128,
    b: u128,
    k: usize,
    ctx: &'static str,
) -> u16
where
    BR: Zq + Ring + Copy + PartialEq,
{
    if let Some(pos) = digit_elems.iter().position(|&x| x == dig) {
        return pos as u16;
    }

    // Diagnostic: try to interpret `dig` as a signed small integer under Zq centering.
    let canon = dig.to_u64().ok();
    let mag = dig.center().to_u64().ok();
    let is_neg = dig.sign() == -BR::ONE;
    let signed_mag: Option<i128> = mag.map(|m| if is_neg { -(m as i128) } else { m as i128 });
    panic!(
        "digit not in alphabet: \
         ctx={} b={} k={} expected alphabet=[-D,D] => [-{},{}]; \
         dig: canon_u64={:?} centered_mag_u64={:?} sign_is_neg={} signed_mag={:?} \
         (check |signed_mag| <= D={} )",
        ctx,
        b,
        k,
        digit_abs_max,
        digit_abs_max,
        canon,
        mag,
        is_neg,
        signed_mag,
        digit_abs_max
    );
}

#[inline]
fn div_floor_i128(a: i128, b: i128) -> i128 {
    debug_assert!(b > 0);
    let q = a / b;
    let r = a % b;
    if r != 0 && a < 0 { q - 1 } else { q }
}

#[inline]
fn div_ceil_i128(a: i128, b: i128) -> i128 {
    debug_assert!(b > 0);
    let q = a / b;
    let r = a % b;
    if r != 0 && a > 0 { q + 1 } else { q }
}

#[inline]
fn div_round_i128(a: i128, b: i128) -> i128 {
    debug_assert!(b > 0);
    let q = div_floor_i128(a, b);
    let r = a - q * b; // 0..b-1
    if r * 2 >= b { q + 1 } else { q }
}

#[inline]
fn br_from_i128<BR: Ring + From<u128>>(x: i128) -> BR {
    if x >= 0 {
        BR::from(x as u128)
    } else {
        -BR::from((-x) as u128)
    }
}

/// Deterministically decompose a *small integer* representative into `k` digits in base `b`,
/// with each digit in `[-digit_abs_max, digit_abs_max]`.
///
/// This is the decomposition model that matches what the verifier can enforce today via
/// unit-monomial exponents: for Frog(d=64), `digit_abs_max = d/2 - 1 = 31`.
#[inline]
fn bounded_decompose_to_digits<BR>(
    x: BR,
    b: i128,
    digit_abs_max: i128,
    pow_b: &[i128],
    rem_bound: &[i128],
    out: &mut [BR],
    row_idx: usize,
    ctx: &'static str,
) where
    BR: Zq + Ring + Copy + From<u128>,
{
    debug_assert!(b >= 2);
    debug_assert_eq!(pow_b.len(), out.len());
    debug_assert_eq!(rem_bound.len(), out.len());

    // Interpret as signed integer via Zq centering (must fit in i128 for our parameter regimes).
    let mag_u64 = x
        .center()
        .to_u64()
        .expect("centered magnitude should fit in u64") as i128;
    let is_neg = x.sign() == -BR::ONE;
    let mut cur: i128 = if is_neg { -mag_u64 } else { mag_u64 };

    // Choose digits from high-to-low to maintain the invariant that the remaining tail is representable.
    for i in (0..out.len()).rev() {
        let bi = pow_b[i];
        let r = rem_bound[i]; // representable bound using lower i digits

        // Need: |cur - d*bi| <= r  =>  (cur-r)/bi <= d <= (cur+r)/bi
        let mut lo = div_ceil_i128(cur - r, bi);
        let mut hi = div_floor_i128(cur + r, bi);
        lo = lo.max(-digit_abs_max);
        hi = hi.min(digit_abs_max);
        if lo > hi {
            panic!(
                "bounded decomposition failed: ctx={} row_idx={} cur={} b={} digit_abs_max={} (lo={}, hi={})",
                ctx, row_idx, cur, b, digit_abs_max, lo, hi
            );
        }

        let mut d = div_round_i128(cur, bi);
        if d < lo {
            d = lo;
        } else if d > hi {
            d = hi;
        }

        cur -= d * bi;
        out[i] = br_from_i128::<BR>(d);
    }

    if cur != 0 {
        panic!(
            "bounded decomposition remainder nonzero: ctx={} row_idx={} rem={} (b={}, digit_abs_max={}, k={})",
            ctx,
            row_idx,
            cur,
            b,
            digit_abs_max,
            out.len()
        );
    }
}

// D_f: decomposed cf(f), Z n x dk
// M_f: EXP(D_f)

#[derive(Clone, Debug)]
pub struct DecompParameters {
    pub b: u128,
    pub k: usize,
    pub l: usize,
}

#[derive(Clone, Debug)]
pub struct FComs<R> {
    pub cm_f: Vec<R>,
    pub C_Mf: Vec<R>,
    pub cm_mtau: Vec<R>,
}

#[derive(Clone, Debug)]
pub struct Rg<R: PolyRing> {
    pub nvars: usize,
    pub instances: Vec<RgInstance<R>>, // L instances
    pub dparams: DecompParameters,
}

#[derive(Clone, Debug)]
pub struct RgInstance<R: PolyRing> {
    /// Monomial matrices in compact digit form (k matrices, each n×d).
    pub M_f: Vec<Arc<DigitsMatrix<R>>>,
    pub tau: Arc<Vec<R::BaseRing>>, // n
    pub m_tau: MonomialVec<R>,      // n, monomials (compact)
    pub f: WitnessVec<R>,           // n
    pub comM_f: Vec<Matrix<R>>,
    pub fcoms: FComs<R>,
}

/// Compact representation for a monomial witness vector `m_tau`.
///
/// In LF+/SP1 regimes, `m_tau[i]` is a unit monomial, so storing a full `Vec<R>` is extremely
/// memory-inefficient (O(n*d) coefficients). We store per-entry indices into a small `exp_table`
/// instead (O(n) u16), and reconstruct `R` on demand where needed.
#[derive(Clone, Debug)]
pub enum MonomialVec<R: PolyRing> {
    Dense(Arc<Vec<R>>),
    Digits {
        digits: Arc<Vec<u16>>,
        exp_table: Arc<Vec<R>>,
    },
}

impl<R: PolyRing> MonomialVec<R> {
    #[inline]
    pub fn len(&self) -> usize {
        match self {
            MonomialVec::Dense(v) => v.len(),
            MonomialVec::Digits { digits, .. } => digits.len(),
        }
    }

    #[inline]
    pub fn get(&self, idx: usize) -> R {
        match self {
            MonomialVec::Dense(v) => v.get(idx).copied().unwrap_or(R::ZERO),
            MonomialVec::Digits { digits, exp_table } => {
                let di = digits.get(idx).copied().unwrap_or(0) as usize;
                exp_table[di]
            }
        }
    }

    #[inline]
    pub fn as_dense_arc(&self) -> Option<Arc<Vec<R>>> {
        match self {
            MonomialVec::Dense(v) => Some(v.clone()),
            MonomialVec::Digits { .. } => None,
        }
    }
}

/// Witness vector representation used by the prover.
///
/// This is a **prover-only** representation choice; verifier behavior/proofs are unchanged.
#[derive(Clone, Debug)]
pub enum WitnessVec<R: PolyRing> {
    /// Fully materialized ring vector.
    Ring(Arc<Vec<R>>),
    /// Constant-coefficient embedding stored as base scalars (avoids allocating `Vec<R>`),
    /// with an explicit *domain length* (typically the Ajtai width / padded `ncols`).
    ///
    /// The `values` are interpreted as a prefix of the witness; indices `j >= values.len()`
    /// are treated as zero. This is CRITICAL for SP1 where `ncols` can be huge but the
    /// nontrivial witness prefix is much smaller.
    ConstCoeffBase {
        values: Arc<Vec<R::BaseRing>>,
        domain_len: usize,
    },
}

impl<R: PolyRing> WitnessVec<R> {
    #[inline]
    pub fn len(&self) -> usize {
        match self {
            WitnessVec::Ring(v) => v.len(),
            WitnessVec::ConstCoeffBase { domain_len, .. } => *domain_len,
        }
    }

    #[inline]
    pub fn values_len(&self) -> usize {
        match self {
            WitnessVec::Ring(v) => v.len(),
            WitnessVec::ConstCoeffBase { values, .. } => values.len(),
        }
    }

    #[inline]
    pub fn as_ring_arc(&self) -> Option<Arc<Vec<R>>> {
        match self {
            WitnessVec::Ring(v) => Some(v.clone()),
            WitnessVec::ConstCoeffBase { .. } => None,
        }
    }

    #[inline]
    pub fn as_const_coeff_base_arc(&self) -> Option<Arc<Vec<R::BaseRing>>> {
        match self {
            WitnessVec::ConstCoeffBase { values, .. } => Some(values.clone()),
            WitnessVec::Ring(_) => None,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Dcom<R: PolyRing> {
    pub evals: Vec<DcomEvals<R>>, // L evals
    pub fcoms: Vec<FComs<R>>,     // L commitments
    pub out: Out<R>,              // set checks
    pub dparams: DecompParameters,
}

#[derive(Clone, Debug)]
pub struct DcomEvals<R: PolyRing> {
    pub v: Vec<R::BaseRing>, // eval over M_f
    pub a: Vec<R::BaseRing>, // eval over tau
    pub b: Vec<R>,           // eval over m_tau
    pub c: Vec<R>,           // eval over f
}

#[derive(Debug, Error)]
pub enum RangeCheckError<R: PolyRing> {
    #[error("Psi check failed: a = {0}, b = {1}")]
    PsiCheckAB(R::BaseRing, R),
    #[error("Psi check failed: v = {0}, u-comb = {1}")]
    PsiCheckVU(Vec<R::BaseRing>, Vec<R>),
}

impl<R: CoeffRing> Rg<R>
where
    R::BaseRing: Zq,
{
    /// Range checks
    ///
    /// Support for `L` [`RgInstance`]s mapped to the corresponding [`DcomEvals`].
    pub fn range_check(
        &self,
        M: &[Arc<SparseMatrix<R>>],
        transcript: &mut impl Transcript<R>,
    ) -> Dcom<R> {
        let profile = std::env::var("LF_PLUS_PROFILE").ok().as_deref() == Some("1");
        let t_total = std::time::Instant::now();

        let mut sets =
            Vec::with_capacity(self.instances.len() * (self.instances[0].M_f.len() + 1));
        for inst in &self.instances {
            inst.M_f.iter().for_each(|m| {
                sets.push(MonomialSet::DigitsMatrix(m.clone()));
            });
        }
        for inst in &self.instances {
            sets.push(match &inst.m_tau {
                MonomialVec::Dense(v) => MonomialSet::Vector(v.clone()),
                MonomialVec::Digits { digits, exp_table } => {
                    MonomialSet::VectorDigits { digits: digits.clone(), exp_table: exp_table.clone() }
                }
            });
        }

        let in_rel = In {
            sets,
            nvars: self.nvars,
        };
        let out_rel = in_rel.set_check(crate::setchk::ExternalMats::Ring(M), transcript);

        // Avoid allocating a full eq-table of size 2^nvars.
        // We instead stream eq-weights in small blocks in the evaluation routines below.
        let one_minus_r = out_rel
            .r
            .iter()
            .copied()
            .map(|x| R::BaseRing::ONE - x)
            .collect::<Vec<_>>();
        if profile {
            println!(
                "[LF+ Rg::range_check] set_check: {:?} (nvars={})",
                t_total.elapsed(),
                self.nvars
            );
        }

        let evals = self
            .instances
            .iter()
            .enumerate()
            .map(|(l, inst)| {
                let mut a = Vec::with_capacity(1 + M.len());
                let mut b = Vec::with_capacity(1 + M.len());
                // Let `c` be the evaluation of `f` over r
                let mut c = Vec::with_capacity(1 + M.len());

                // v: coefficient-wise evaluation of f at out_rel.r
                let v =
                    eval_vec_coeffs_at_point_streaming_witness::<R>(&inst.f, &out_rel.r, &one_minus_r);

                a.push(dot_base_streaming::<R>(
                    inst.tau.as_ref(),
                    &out_rel.r,
                    &one_minus_r,
                ));
                b.push(out_rel.b[l]);
                c.push(dot_ring_streaming_witness::<R>(&inst.f, &out_rel.r, &one_minus_r));

                // Evaluate M * tau / m_tau / f at out_rel.r *without materializing length-n vectors*.
                //
                // For each matrix M:
                //   eval(M*w)(r) = Σ_row eq[row] * (Σ_{(coeff,col) in row} coeff * w[col])
                for m in M {
                    a.push(sparse_mat_vec_eval_ct_streaming::<R>(
                        m,
                        inst.tau.as_ref(),
                        &out_rel.r,
                        &one_minus_r,
                    ));
                    b.push(match &inst.m_tau {
                        MonomialVec::Dense(v) => sparse_mat_vec_eval_ring_streaming::<R>(
                            m,
                            v.as_ref(),
                            &out_rel.r,
                            &one_minus_r,
                        ),
                        MonomialVec::Digits { digits, exp_table } => {
                            sparse_mat_vec_eval_ring_streaming_monomial_digits::<R>(
                                m,
                                digits.as_ref(),
                                exp_table.as_ref(),
                                &out_rel.r,
                                &one_minus_r,
                            )
                        }
                    });
                    c.push(sparse_mat_vec_eval_ring_streaming_witness::<R>(
                        m,
                        &inst.f,
                        &out_rel.r,
                        &one_minus_r,
                    ));
                }
                DcomEvals { v, a, b, c }
            })
            .collect::<Vec<_>>();

        absorb_evaluations(&evals, transcript);

        if profile {
            println!("[LF+ Rg::range_check] evals+absorb: {:?}", t_total.elapsed());
        }

        Dcom {
            evals,
            fcoms: self
                .instances
                .iter()
                .map(|inst| inst.fcoms.clone())
                .collect(),
            out: out_rel,
            dparams: self.dparams.clone(),
        }
    }

    /// Range checks, but with external matrices represented over the **base ring**.
    ///
    /// This is the natural representation for SP1/R1LF chunks (const-coeff by construction) and
    /// avoids materializing `SparseMatrix<R>` which is catastrophic at large `R::dimension()`.
    pub fn range_check_base(
        &self,
        M0: &[Arc<SparseMatrix<R::BaseRing>>],
        transcript: &mut impl Transcript<R>,
    ) -> Dcom<R> {
        let profile = std::env::var("LF_PLUS_PROFILE").ok().as_deref() == Some("1");
        let t_total = std::time::Instant::now();

        let mut sets =
            Vec::with_capacity(self.instances.len() * (self.instances[0].M_f.len() + 1));
        for inst in &self.instances {
            inst.M_f.iter().for_each(|m| {
                sets.push(MonomialSet::DigitsMatrix(m.clone()));
            });
        }
        for inst in &self.instances {
            sets.push(match &inst.m_tau {
                MonomialVec::Dense(v) => MonomialSet::Vector(v.clone()),
                MonomialVec::Digits { digits, exp_table } => {
                    MonomialSet::VectorDigits { digits: digits.clone(), exp_table: exp_table.clone() }
                }
            });
        }

        let in_rel = In {
            sets,
            nvars: self.nvars,
        };
        let out_rel = in_rel.set_check(crate::setchk::ExternalMats::Base(M0), transcript);

        let one_minus_r = out_rel
            .r
            .iter()
            .copied()
            .map(|x| R::BaseRing::ONE - x)
            .collect::<Vec<_>>();
        if profile {
            println!(
                "[LF+ Rg::range_check] set_check: {:?} (nvars={})",
                t_total.elapsed(),
                self.nvars
            );
        }

        let evals = self
            .instances
            .iter()
            .enumerate()
            .map(|(l, inst)| {
                let mut a = Vec::with_capacity(1 + M0.len());
                let mut b = Vec::with_capacity(1 + M0.len());
                let mut c = Vec::with_capacity(1 + M0.len());

                let v =
                    eval_vec_coeffs_at_point_streaming_witness::<R>(&inst.f, &out_rel.r, &one_minus_r);

                a.push(dot_base_streaming::<R>(
                    inst.tau.as_ref(),
                    &out_rel.r,
                    &one_minus_r,
                ));
                b.push(out_rel.b[l]);
                c.push(dot_ring_streaming_witness::<R>(&inst.f, &out_rel.r, &one_minus_r));

                for m0 in M0 {
                    a.push(sparse_mat0_vec_eval_ct_streaming::<R>(
                        m0,
                        inst.tau.as_ref(),
                        &out_rel.r,
                        &one_minus_r,
                    ));
                    b.push(match &inst.m_tau {
                        MonomialVec::Dense(v) => sparse_mat0_vec_eval_ring_streaming::<R>(
                            m0,
                            v.as_ref(),
                            &out_rel.r,
                            &one_minus_r,
                        ),
                        MonomialVec::Digits { digits, exp_table } => {
                            sparse_mat0_vec_eval_ring_streaming_monomial_digits::<R>(
                                m0,
                                digits.as_ref(),
                                exp_table.as_ref(),
                                &out_rel.r,
                                &one_minus_r,
                            )
                        }
                    });
                    c.push(sparse_mat0_vec_eval_ring_streaming_witness::<R>(
                        m0,
                        &inst.f,
                        &out_rel.r,
                        &one_minus_r,
                    ));
                }
                DcomEvals { v, a, b, c }
            })
            .collect::<Vec<_>>();

        absorb_evaluations(&evals, transcript);

        if profile {
            println!("[LF+ Rg::range_check] evals+absorb: {:?}", t_total.elapsed());
        }

        Dcom {
            evals,
            fcoms: self
                .instances
                .iter()
                .map(|inst| inst.fcoms.clone())
                .collect(),
            out: out_rel,
            dparams: self.dparams.clone(),
        }
    }
}

impl<R: CoeffRing> Dcom<R>
where
    R::BaseRing: Zq,
{
    pub fn verify(&self, transcript: &mut impl Transcript<R>) -> Result<(), RangeCheckError<R>> {
        self.out.verify(transcript).unwrap(); //.map_err(|_| ())?;

        absorb_evaluations(&self.evals, transcript);

        for (l, eval) in self.evals.iter().enumerate() {
            // ct(psi b) =? a
            for (&a_i, b_i) in eval.a.iter().zip(eval.b.iter()) {
                ((psi::<R>() * b_i).ct() == a_i)
                    .then_some(())
                    .ok_or(RangeCheckError::PsiCheckAB(a_i, *b_i))?;
            }

            let d = R::dimension();
            let base = self.dparams.b;
            for (ni, _) in self.out.e.iter().enumerate() {
                let base_br = R::BaseRing::from(base);
                let mut d_ppow = R::BaseRing::ONE;
                let u_comb = self.out.e[ni]
                    .iter()
                    .skip(self.dparams.k * l)
                    .take(self.dparams.k)
                    .fold(vec![R::zero(); d], |mut acc, u_i| {
                        // Avoid pow() in verification: maintain `base^i` incrementally.
                        u_i.iter()
                            .zip(acc.iter_mut())
                            .for_each(|(u_ij, a_j)| *a_j += *u_ij * d_ppow);
                        d_ppow *= base_br;
                        acc
                    });

                // ct(psi (sum d^i u_i)) =? v
                let v_rec = u_comb
                    .iter()
                    .map(|&uc| (psi::<R>() * uc).ct())
                    .collect::<Vec<_>>();

                if ni == 0 {
                    (eval.v == v_rec)
                        .then_some(())
                        .ok_or(RangeCheckError::PsiCheckVU(v_rec, u_comb))?;
                } else {
                    (eval.c[ni].coeffs() == v_rec)
                        .then_some(())
                        .ok_or(RangeCheckError::PsiCheckVU(v_rec, u_comb))?;
                }
            }
        }

        Ok(())
    }
}

impl<R: PolyRing> RgInstance<R> {
    /// Construct monomial sets from `M_f` and `m_tau`
    pub fn sets(&self) -> Vec<MonomialSet<R>> {
        self.M_f
            .iter()
            .map(|m| MonomialSet::DigitsMatrix(m.clone()))
            .chain(once(match &self.m_tau {
                MonomialVec::Dense(v) => MonomialSet::Vector(v.clone()),
                MonomialVec::Digits { digits, exp_table } => {
                    MonomialSet::VectorDigits { digits: digits.clone(), exp_table: exp_table.clone() }
                }
            }))
            .collect()
    }
}

impl<R: CoeffRing> RgInstance<R>
where
    R::BaseRing: Zq + From<u128>,
    R: Decompose,
{
    pub fn from_f(f: Vec<R>, A: &Matrix<R>, decomp: &DecompParameters) -> Self {
        let profile = std::env::var("LF_PLUS_PROFILE").ok().as_deref() == Some("1");
        let t_total = std::time::Instant::now();

        let n = f.len();

        // Build compact digit matrices for cf(f) decomposition.
        //
        // Previous code used `row.decompose_to_vec(...)` which allocates a `Vec` of length `k`
        // for every coefficient => ~16M allocations at n=1M,d=16,k=4. This dominated runtime.
        let d = R::dimension();
        let k = decomp.k;
        let t = std::time::Instant::now();
        // Digit encoding (monomial path):
        //
        // The verifier enforces that each digit is a unit monomial exponent (via psi/exp),
        // which implies a conservative per-digit bound of |digit| <= D where D = d/2 - 1.
        // We therefore build the digit alphabet as [-D, D] and decompose coefficients into that set.
        let digit_abs_max: i128 = (R::dimension() as i128) / 2 - 1;
        assert!(digit_abs_max >= 1, "ring dimension too small for monomial digits");
        let b_i128: i128 = decomp.b as i128;
        assert!(b_i128 >= 2, "decomposition base must be >= 2");
        let digit_elems: Vec<R::BaseRing> = (-digit_abs_max..=digit_abs_max)
            .map(|x| br_from_i128::<R::BaseRing>(x))
            .collect();
        assert!(
            digit_elems.len() <= (u16::MAX as usize),
            "digit alphabet too large for u16 indices (len={})",
            digit_elems.len()
        );
        let exp_table: Arc<Vec<R>> = Arc::new(
            digit_elems
                .iter()
                .map(|&x| exp::<R>(x).unwrap())
                .collect::<Vec<_>>(),
        );
        let digit_elems = Arc::new(digit_elems);
        let b_u128 = decomp.b;
        let ctx: &'static str = "RgInstance::from_f";
        let map_digit_to_idx: Box<dyn Fn(R::BaseRing) -> u16 + Send + Sync> =
            Box::new(move |dig: R::BaseRing| -> u16 {
                digit_lookup_or_panic(&digit_elems, dig, digit_abs_max, b_u128, k, ctx)
            });

        // Precompute base powers and tail-bounds for bounded decomposition.
        let mut pow_b: Vec<i128> = vec![1; k];
        for i in 1..k {
            pow_b[i] = pow_b[i - 1] * b_i128;
        }
        let mut rem_bound: Vec<i128> = vec![0; k];
        // rem_bound[i] = D * (b^i - 1)/(b - 1)
        let denom = b_i128 - 1;
        for i in 0..k {
            rem_bound[i] = digit_abs_max * (pow_b[i].saturating_sub(1)) / denom;
        }

        // If `f` is constant-coefficient (only col=0 can be nonzero), we can store only the
        // `col=0` digit table and treat all other coeff columns as digit=0.
        let is_const_coeff = f.iter().all(|fi| fi.coeffs().iter().skip(1).all(|c| *c == R::BaseRing::ZERO));
        let zero_idx: u16 = (map_digit_to_idx)(R::BaseRing::ZERO);

        // Allocate digit tables:
        // - const-coeff: `k × n` (only col0)
        // - general:     `k × (n*d)` (row-major full)
        let mut digits_tables: Vec<Vec<u16>> = if is_const_coeff {
            (0..k).map(|_| vec![zero_idx; n]).collect()
        } else {
            (0..k).map(|_| vec![0u16; n * d]).collect()
        };
        #[cfg(feature = "parallel")]
        {
            use rayon::prelude::*;
            digits_tables
                .par_iter_mut()
                .enumerate()
                .for_each(|(k_i, table)| {
                    // Each thread gets its own scratch to avoid races.
                    let mut tmp_local = vec![R::BaseRing::ZERO; k];
                    for (row_idx, fi) in f.iter().enumerate() {
                        let coeffs = fi.coeffs();
                        debug_assert_eq!(coeffs.len(), d);
                        if is_const_coeff {
                            // Only constant coefficient matters; other columns are fixed to digit=0.
                            bounded_decompose_to_digits(
                                coeffs[0],
                                b_i128,
                                digit_abs_max,
                                &pow_b,
                                &rem_bound,
                                &mut tmp_local,
                                row_idx,
                                ctx,
                            );
                            table[row_idx] = (map_digit_to_idx)(tmp_local[k_i]);
                        } else {
                            for (col_idx, &c) in coeffs.iter().enumerate() {
                                bounded_decompose_to_digits(
                                    c,
                                    b_i128,
                                    digit_abs_max,
                                    &pow_b,
                                    &rem_bound,
                                    &mut tmp_local,
                                    row_idx,
                                    ctx,
                                );
                                table[row_idx * d + col_idx] = (map_digit_to_idx)(tmp_local[k_i]);
                            }
                        }
                    }
                });
        }
        #[cfg(not(feature = "parallel"))]
        {
            let mut tmp = vec![R::BaseRing::ZERO; k];
            for (row_idx, fi) in f.iter().enumerate() {
                let coeffs = fi.coeffs();
                debug_assert_eq!(coeffs.len(), d);
                if is_const_coeff {
                    bounded_decompose_to_digits(
                        coeffs[0],
                        b_i128,
                        digit_abs_max,
                        &pow_b,
                        &rem_bound,
                        &mut tmp,
                        row_idx,
                        ctx,
                    );
                    for k_i in 0..k {
                        digits_tables[k_i][row_idx] = (map_digit_to_idx)(tmp[k_i]);
                    }
                } else {
                    for (col_idx, &c) in coeffs.iter().enumerate() {
                        // Writes into tmp[0..k] in-place.
                        bounded_decompose_to_digits(
                            c,
                            b_i128,
                            digit_abs_max,
                            &pow_b,
                            &rem_bound,
                            &mut tmp,
                            row_idx,
                            ctx,
                        );
                        for k_i in 0..k {
                            digits_tables[k_i][row_idx * d + col_idx] = (map_digit_to_idx)(tmp[k_i]);
                        }
                    }
                }
            }
        }
        if profile {
            println!(
                "[LF+ RgInstance::from_f] decompose_to (no-alloc): {:?} (n={}, d={}, k={})",
                t.elapsed(),
                n,
                d,
                k
            );
        }

        let t = std::time::Instant::now();
        // Commit monomial matrices: comM_f[k_i] = A * M_f[k_i] without materializing full `M_f`.
        //
        // `A.try_mul_mat` appears to under-utilize CPU (kappa small, not parallelized),
        // so we explicitly parallelize over columns + rows (rayon reduction).
        fn commit_digits_matrix<Rr>(a: &Matrix<Rr>, m: &DigitsMatrix<Rr>) -> Matrix<Rr>
        where
            Rr: CoeffRing,
            Rr::BaseRing: Zq,
        {
            let kappa = a.nrows;
            let n = a.ncols;
            debug_assert_eq!(m.nrows, n);
            let dcols = m.ncols;
            #[cfg(feature = "parallel")]
            {
                // Compute one commitment vector per column in parallel, then assemble.
                let cols: Vec<Vec<Rr>> = (0..dcols)
                    .into_par_iter()
                    .map(|col| {
                        (0..n)
                            .into_par_iter()
                            .fold(
                                || vec![Rr::ZERO; kappa],
                                |mut acc, i| {
                                    let mi = m.get(i, col);
                                    for r in 0..kappa {
                                        acc[r] += a.vals[r][i] * mi;
                                    }
                                    acc
                                },
                            )
                            .reduce(
                                || vec![Rr::ZERO; kappa],
                                |mut a0, b0| {
                                    for r in 0..kappa {
                                        a0[r] += b0[r];
                                    }
                                    a0
                                },
                            )
                    })
                    .collect();
                let mut out = Matrix::zero(kappa, dcols);
                for col in 0..dcols {
                    for r in 0..kappa {
                        out.vals[r][col] = cols[col][r];
                    }
                }
                out
            }
            #[cfg(not(feature = "parallel"))]
            {
                let mut out = Matrix::zero(kappa, dcols);
                for col in 0..dcols {
                    for r in 0..kappa {
                        let mut acc = Rr::ZERO;
                        for i in 0..n {
                            acc += a.vals[r][i] * m.get(i, col);
                        }
                        out.vals[r][col] = acc;
                    }
                }
                out
            }
        }

        let M_f: Vec<Arc<DigitsMatrix<R>>> = digits_tables
            .into_iter()
            .map(|digits| {
                Arc::new(DigitsMatrix {
                    nrows: n,
                    ncols: d,
                    digits: if is_const_coeff {
                        crate::setchk::DigitsBacking::ConstCol0 {
                            col0: Arc::new(digits),
                            zero_idx,
                        }
                    } else {
                        crate::setchk::DigitsBacking::Full(Arc::new(digits))
                    },
                    exp_table: exp_table.clone(),
                })
            })
            .collect();

        let comM_f = M_f
            .iter()
            .map(|M| commit_digits_matrix(A, M.as_ref()))
            .collect::<Vec<_>>();
        let com = Matrix::hconcat(&comM_f).unwrap();
        if profile {
            println!(
                "[LF+ RgInstance::from_f] commit monomial mats (A*M_f): {:?} (kappa×(k*d) = {}×{})",
                t.elapsed(),
                A.nrows,
                decomp.k * d
            );
        }

        let t = std::time::Instant::now();
        let tau = split(&com, n, (R::dimension() / 2) as u128, decomp.l);
        if profile {
            println!("[LF+ RgInstance::from_f] split tau: {:?}", t.elapsed());
        }

        let t = std::time::Instant::now();
        let m_tau = tau
            .iter()
            .map(|c| exp::<R>(*c).unwrap())
            .collect::<Vec<_>>();
        if profile {
            println!("[LF+ RgInstance::from_f] build m_tau via exp: {:?}", t.elapsed());
        }

        let t = std::time::Instant::now();
        let cm_f = A.try_mul_vec(&f).unwrap();
        let C_Mf = A
            .try_mul_vec(&tau.iter().map(|z| R::from(*z)).collect::<Vec<R>>())
            .unwrap();
        let cm_mtau = A.try_mul_vec(&m_tau).unwrap();
        if profile {
            println!("[LF+ RgInstance::from_f] commit f/tau/m_tau: {:?}", t.elapsed());
            println!("[LF+ RgInstance::from_f] total: {:?}", t_total.elapsed());
        }
        let fcoms = FComs {
            cm_f,
            C_Mf,
            cm_mtau,
        };

        Self {
            M_f,
            tau: Arc::new(tau),
            m_tau: MonomialVec::Dense(Arc::new(m_tau)),
            f: WitnessVec::Ring(Arc::new(f)),
            comM_f,
            fcoms,
        }
    }

    /// Construct an [`RgInstance`] from witness `f`, using a **seeded implicit Ajtai matrix**.
    ///
    /// This avoids materializing a `kappa × n` dense matrix in memory. Outputs are identical
    /// in distribution (up to the Ajtai matrix generation procedure), and verifier behavior is unchanged.
    pub fn from_f_seeded(f: Vec<R>, scheme: &AjtaiCommitmentScheme<R>, decomp: &DecompParameters) -> Self {
        let profile = std::env::var("LF_PLUS_PROFILE").ok().as_deref() == Some("1");
        let t_total = std::time::Instant::now();

        let n = f.len();
        let kappa = scheme.kappa();
        let d = R::dimension();
        let k = decomp.k;

        // Reuse the same digit-table logic as `from_f`, including the const-coeff optimization.
        // (This keeps transcript behavior and setchk wiring identical.)
        let digit_abs_max: i128 = (R::dimension() as i128) / 2 - 1;
        assert!(digit_abs_max >= 1, "ring dimension too small for monomial digits");
        let b_i128: i128 = decomp.b as i128;
        assert!(b_i128 >= 2, "decomposition base must be >= 2");
        let digit_elems: Vec<R::BaseRing> = (-digit_abs_max..=digit_abs_max)
            .map(|x| br_from_i128::<R::BaseRing>(x))
            .collect();
        let exp_table: Arc<Vec<R>> = Arc::new(
            digit_elems
                .iter()
                .map(|&x| exp::<R>(x).unwrap())
                .collect::<Vec<_>>(),
        );
        let digit_elems = Arc::new(digit_elems);
        let b_u128 = decomp.b;
        let ctx: &'static str = "RgInstance::from_f_seeded";
        let map_digit_to_idx: Box<dyn Fn(R::BaseRing) -> u16 + Send + Sync> =
            Box::new(move |dig: R::BaseRing| -> u16 {
                digit_lookup_or_panic(&digit_elems, dig, digit_abs_max, b_u128, k, ctx)
            });

        let is_const_coeff = f.iter().all(|fi| fi.coeffs().iter().skip(1).all(|c| *c == R::BaseRing::ZERO));
        let zero_idx: u16 = (map_digit_to_idx)(R::BaseRing::ZERO);

        let mut pow_b: Vec<i128> = vec![1; k];
        for i in 1..k {
            pow_b[i] = pow_b[i - 1] * b_i128;
        }
        let mut rem_bound: Vec<i128> = vec![0; k];
        let denom = b_i128 - 1;
        for i in 0..k {
            rem_bound[i] = digit_abs_max * (pow_b[i].saturating_sub(1)) / denom;
        }
        let mut digits_tables: Vec<Vec<u16>> = if is_const_coeff {
            (0..k).map(|_| vec![zero_idx; n]).collect()
        } else {
            (0..k).map(|_| vec![0u16; n * d]).collect()
        };
        #[cfg(feature = "parallel")]
        {
            use rayon::prelude::*;
            digits_tables
                .par_iter_mut()
                .enumerate()
                .for_each(|(k_i, table)| {
                    let mut tmp_local = vec![R::BaseRing::ZERO; k];
                    for (row_idx, fi) in f.iter().enumerate() {
                        let coeffs = fi.coeffs();
                        debug_assert_eq!(coeffs.len(), d);
                        if is_const_coeff {
                            bounded_decompose_to_digits(
                                coeffs[0],
                                b_i128,
                                digit_abs_max,
                                &pow_b,
                                &rem_bound,
                                &mut tmp_local,
                                row_idx,
                                ctx,
                            );
                            table[row_idx] = (map_digit_to_idx)(tmp_local[k_i]);
                        } else {
                            for (col_idx, &c) in coeffs.iter().enumerate() {
                                bounded_decompose_to_digits(
                                    c,
                                    b_i128,
                                    digit_abs_max,
                                    &pow_b,
                                    &rem_bound,
                                    &mut tmp_local,
                                    row_idx,
                                    ctx,
                                );
                                table[row_idx * d + col_idx] = (map_digit_to_idx)(tmp_local[k_i]);
                            }
                        }
                    }
                });
        }
        #[cfg(not(feature = "parallel"))]
        {
            let mut tmp = vec![R::BaseRing::ZERO; k];
            for (row_idx, fi) in f.iter().enumerate() {
                let coeffs = fi.coeffs();
                debug_assert_eq!(coeffs.len(), d);
                if is_const_coeff {
                    bounded_decompose_to_digits(
                        coeffs[0],
                        b_i128,
                        digit_abs_max,
                        &pow_b,
                        &rem_bound,
                        &mut tmp,
                        row_idx,
                        ctx,
                    );
                    for k_i in 0..k {
                        digits_tables[k_i][row_idx] = (map_digit_to_idx)(tmp[k_i]);
                    }
                } else {
                    for (col_idx, &c) in coeffs.iter().enumerate() {
                        bounded_decompose_to_digits(
                            c,
                            b_i128,
                            digit_abs_max,
                            &pow_b,
                            &rem_bound,
                            &mut tmp,
                            row_idx,
                            ctx,
                        );
                        for k_i in 0..k {
                            digits_tables[k_i][row_idx * d + col_idx] = (map_digit_to_idx)(tmp[k_i]);
                        }
                    }
                }
            }
        }

        let M_f: Vec<Arc<DigitsMatrix<R>>> = digits_tables
            .into_iter()
            .map(|digits| {
                Arc::new(DigitsMatrix {
                    nrows: n,
                    ncols: d,
                    digits: if is_const_coeff {
                        crate::setchk::DigitsBacking::ConstCol0 {
                            col0: Arc::new(digits),
                            zero_idx,
                        }
                    } else {
                        crate::setchk::DigitsBacking::Full(Arc::new(digits))
                    },
                    exp_table: exp_table.clone(),
                })
            })
            .collect();

        // Commit monomial matrices: comM_f[k_i] is kappa×d, column-wise Ajtai commitments under the same scheme.
        //
        // Optimization (critical for SP1 const-coeff regime):
        // If `DigitsBacking::ConstCol0` is used, then only column 0 varies by row; all columns 1..d-1
        // are the constant `exp(0)` column. We can commit that constant column once and reuse it,
        // turning a ~d× scan of the Ajtai matrix into ~2× (col0 + const column).
        let t = std::time::Instant::now();
        let const_exp0_commit: Option<Vec<R>> = if is_const_coeff && d > 1 {
            let exp0 = M_f[0].exp_table[zero_idx as usize];
            let c = scheme
                .commit_many_with(n, 1, move |_row, out| {
                    out[0] = exp0;
                })
                .expect("commit_many_with const exp(0) col");
            Some(c[0].as_ref().to_vec())
        } else {
            None
        };

        let comM_f = M_f
            .iter()
            .map(|dm| {
                let mut mat = Matrix::zero(kappa, d);

                // Commit column 0 (the only non-constant column in const-coeff mode).
                let c0 = scheme
                    .commit_many_with(n, 1, {
                        let dm = dm.clone();
                        move |row, out| {
                            out[0] = dm.get(row, 0);
                        }
                    })
                    .expect("commit_many_with M_f col0");
                let c0 = c0[0].as_ref();
                for r in 0..kappa {
                    mat.vals[r][0] = c0[r];
                }

                if let Some(cc) = const_exp0_commit.as_ref() {
                    // Reuse constant exp(0) commitment for columns 1..d-1.
                    for col in 1..d {
                        for r in 0..kappa {
                            mat.vals[r][col] = cc[r];
                        }
                    }
                } else if d > 1 {
                    // General path: commit remaining columns in one batch.
                    let cs = scheme
                        .commit_many_with(n, d - 1, {
                            let dm = dm.clone();
                            move |row, out| {
                                for col in 1..d {
                                    out[col - 1] = dm.get(row, col);
                                }
                            }
                        })
                        .expect("commit_many_with M_f cols1..d-1");
                    for col in 1..d {
                        let ccol = cs[col - 1].as_ref();
                        for r in 0..kappa {
                            mat.vals[r][col] = ccol[r];
                        }
                    }
                }
                mat
            })
            .collect::<Vec<_>>();
        if profile {
            println!(
                "[LF+ RgInstance::from_f_seeded] commit monomial mats (Ajtai seeded): {:?} (kappa×(k*d) = {}×{})",
                t.elapsed(),
                kappa,
                decomp.k * d
            );
        }

        let com = Matrix::hconcat(&comM_f).unwrap();

        let t = std::time::Instant::now();
        let tau = split(&com, n, (R::dimension() / 2) as u128, decomp.l);
        if profile {
            println!("[LF+ RgInstance::from_f_seeded] split tau: {:?}", t.elapsed());
        }

        let t = std::time::Instant::now();
        let m_tau = tau.iter().map(|c| exp::<R>(*c).unwrap()).collect::<Vec<_>>();
        if profile {
            println!("[LF+ RgInstance::from_f_seeded] build m_tau via exp: {:?}", t.elapsed());
        }

        let t = std::time::Instant::now();
        // Batch the two constant-coeff commitments (f and tau) in a single pass over columns.
        //
        // This matches Symphony's `commit_many_with` pattern and avoids an extra full scan / RNG stream.
        let f0: Arc<Vec<R::BaseRing>> = Arc::new(f.iter().map(|x| x.coeffs()[0]).collect());
        // `tau` is huge (length n); avoid cloning it. Wrap once and share the Arc.
        let tau0: Arc<Vec<R::BaseRing>> = Arc::new(tau);
        let cm_pair = scheme
            .commit_many_const_coeff_base_fast(n, 2, {
                let f0 = f0.clone();
                let tau0 = tau0.clone();
                move |j, out| {
                    // out.len()==2
                    out[0] = f0[j];
                    out[1] = tau0[j];
                }
            })
            .expect("commit_many_const_coeff_base_fast (f,tau)");
        let cm_f = cm_pair[0].as_ref().to_vec();
        let C_Mf = cm_pair[1].as_ref().to_vec();
        let cm_mtau = scheme.commit(&m_tau).expect("commit m_tau").as_ref().to_vec();
        if profile {
            println!("[LF+ RgInstance::from_f_seeded] commit f/tau/m_tau: {:?}", t.elapsed());
            println!("[LF+ RgInstance::from_f_seeded] total: {:?}", t_total.elapsed());
        }
        let fcoms = FComs { cm_f, C_Mf, cm_mtau };

        Self {
            M_f,
            tau: tau0.clone(),
            m_tau: MonomialVec::Dense(Arc::new(m_tau)),
            f: WitnessVec::Ring(Arc::new(f)),
            comM_f,
            fcoms,
        }
    }

    /// Construct an [`RgInstance`] from a **constant-coefficient** witness `f0` (base scalars),
    /// using a seeded implicit Ajtai matrix.
    ///
    /// This avoids allocating the huge `Vec<R>` when the witness is embedded as constant-coeff
    /// ring elements (the SP1 production regime).
    pub fn from_f0_seeded(
        f0: Arc<Vec<R::BaseRing>>,
        scheme: &AjtaiCommitmentScheme<R>,
        decomp: &DecompParameters,
    ) -> Self {
        let profile = std::env::var("LF_PLUS_PROFILE").ok().as_deref() == Some("1");
        let t_total = std::time::Instant::now();

        // Domain length must match the Ajtai width (padded `ncols`), even if the witness is only a prefix.
        let n = scheme.width();
        let prefix_len = f0.len();
        let kappa = scheme.kappa();
        let d = R::dimension();
        let k = decomp.k;
        if profile {
            println!(
                "[LF+ RgInstance::from_f0_seeded] start: n(domain)={} prefix_len={} kappa={} d={} k={}",
                n, prefix_len, kappa, d, k
            );
        }

        let digit_abs_max: i128 = (R::dimension() as i128) / 2 - 1;
        assert!(digit_abs_max >= 1, "ring dimension too small for monomial digits");
        let b_i128: i128 = decomp.b as i128;
        assert!(b_i128 >= 2, "decomposition base must be >= 2");
        let digit_elems: Vec<R::BaseRing> = (-digit_abs_max..=digit_abs_max)
            .map(|x| br_from_i128::<R::BaseRing>(x))
            .collect();
        let exp_table: Arc<Vec<R>> = Arc::new(
            digit_elems
                .iter()
                .map(|&x| exp::<R>(x).unwrap())
                .collect::<Vec<_>>(),
        );
        let digit_elems = Arc::new(digit_elems);
        let digit_elems_for_idx = digit_elems.clone();
        let b_u128 = decomp.b;
        let ctx: &'static str = "RgInstance::from_f0_seeded";
        let map_digit_to_idx: Box<dyn Fn(R::BaseRing) -> u16 + Send + Sync> =
            Box::new(move |dig: R::BaseRing| -> u16 {
                digit_lookup_or_panic(&digit_elems_for_idx, dig, digit_abs_max, b_u128, k, ctx)
            });

        // Tau digits from `split` are within ±(d/4), which is always within our witness digit alphabet [-D,D].
        let tau_b: u128 = (R::dimension() / 2) as u128;
        assert!(tau_b >= 2 && (tau_b % 2 == 0), "tau decomposition base must be even");
        // Map tau digits using the same alphabet/exp-table.
        let ctx_tau: &'static str = "RgInstance::from_f0_seeded::m_tau";
        let map_tau_digit_to_idx: Box<dyn Fn(R::BaseRing) -> u16 + Send + Sync> = Box::new({
            let digit_elems = digit_elems.clone();
            move |dig: R::BaseRing| -> u16 { digit_lookup_or_panic(&digit_elems, dig, digit_abs_max, tau_b, 0, ctx_tau) }
        });

        // Const-coeff witness: store only col0 digit table.
        let zero_idx: u16 = (map_digit_to_idx)(R::BaseRing::ZERO);
        // Only materialize digits for the prefix; rows beyond `prefix_len` are implicitly zero digits.
        let mut digits_tables: Vec<Vec<u16>> = (0..k).map(|_| vec![zero_idx; prefix_len]).collect();
        let t_digits = std::time::Instant::now();
        let mut pow_b: Vec<i128> = vec![1; k];
        for i in 1..k {
            pow_b[i] = pow_b[i - 1] * b_i128;
        }
        let mut rem_bound: Vec<i128> = vec![0; k];
        let denom = b_i128 - 1;
        for i in 0..k {
            rem_bound[i] = digit_abs_max * (pow_b[i].saturating_sub(1)) / denom;
        }
        #[cfg(feature = "parallel")]
        {
            use rayon::prelude::*;
            // Don't parallelize over `k` (k=11) because it only uses ~k threads and also recomputes the
            // full decomposition k times per row (O(k^2) work). Instead decompose once per row and
            // write all k digits, parallelizing over rows.
            //
            // Safety: each parallel task writes a distinct `row_idx` across all tables. The `Vec<u16>`
            // allocations are fixed-size and won't reallocate during the fill.
            struct TablePtrs {
                // Store as `usize` to avoid `Sync` capture issues with raw pointers.
                ptrs_addr: usize,
                len: usize,
            }
            // We only ever write to disjoint indices (row-wise), so sharing these pointers is safe.
            unsafe impl Sync for TablePtrs {}
            unsafe impl Send for TablePtrs {}

            let ptrs: Vec<*mut u16> = digits_tables.iter_mut().map(|t| t.as_mut_ptr()).collect();
            let tbl = TablePtrs {
                ptrs_addr: ptrs.as_ptr() as usize,
                len: ptrs.len(),
            };
            f0.par_iter()
                .enumerate()
                .for_each_init(|| vec![R::BaseRing::ZERO; k], |tmp, (row_idx, &c0)| {
                    bounded_decompose_to_digits(
                        c0,
                        b_i128,
                        digit_abs_max,
                        &pow_b,
                        &rem_bound,
                        tmp,
                        row_idx,
                        ctx,
                    );
                    for k_i in 0..k {
                        let dig = (map_digit_to_idx)(tmp[k_i]);
                        unsafe {
                            debug_assert!(k_i < tbl.len);
                            let base = tbl.ptrs_addr as *const *mut u16;
                            *(*base.add(k_i)).add(row_idx) = dig;
                        }
                    }
                });
        }
        #[cfg(not(feature = "parallel"))]
        {
            let mut tmp = vec![R::BaseRing::ZERO; k];
            for (row_idx, &c0) in f0.iter().enumerate() {
                bounded_decompose_to_digits(
                    c0,
                    b_i128,
                    digit_abs_max,
                    &pow_b,
                    &rem_bound,
                    &mut tmp,
                    row_idx,
                    ctx,
                );
                for k_i in 0..k {
                    digits_tables[k_i][row_idx] = (map_digit_to_idx)(tmp[k_i]);
                }
            }
        }
        if profile {
            println!(
                "[LF+ RgInstance::from_f0_seeded] build digit tables (prefix only): {:?}",
                t_digits.elapsed()
            );
        }

        let M_f: Vec<Arc<DigitsMatrix<R>>> = digits_tables
            .into_iter()
            .map(|digits| {
                Arc::new(DigitsMatrix {
                    nrows: n,
                    ncols: d,
                    digits: crate::setchk::DigitsBacking::ConstCol0 {
                        col0: Arc::new(digits),
                        zero_idx,
                    },
                    exp_table: exp_table.clone(),
                })
            })
            .collect();

        // Commit monomial matrices (Ajtai seeded).
        let t = std::time::Instant::now();
        // Major speed win: batch ALL k col0 vectors + the constant exp(0) vector in one Ajtai pass,
        // using the **monomial-digit** specialization (rotation instead of full ring mul).
        let mfs = M_f.clone();
        let exp_table = exp_table.clone();
        let cs = scheme
            .commit_many_with_monomial_digits(n, k + 1, exp_table.clone(), move |row, out| {
                for ki in 0..k {
                    let dm = mfs[ki].as_ref();
                    match &dm.digits {
                        crate::setchk::DigitsBacking::ConstCol0 { col0, zero_idx } => {
                            out[ki] = col0.get(row).copied().unwrap_or(*zero_idx);
                        }
                        crate::setchk::DigitsBacking::Full(_) => unreachable!("from_f0_seeded uses ConstCol0"),
                    }
                }
                // constant exp(0) digit
                out[k] = zero_idx;
            })
            .expect("commit_many_with_monomial_digits (M_f col0 batch + const exp0)");

        let const_exp0_commit = cs[k].as_ref().to_vec();
        let comM_f = (0..k)
            .map(|ki| {
                let mut mat = Matrix::zero(kappa, d);
                let c0 = cs[ki].as_ref();
                for r in 0..kappa {
                    mat.vals[r][0] = c0[r];
                }
                for col in 1..d {
                    for r in 0..kappa {
                        mat.vals[r][col] = const_exp0_commit[r];
                    }
                }
                mat
            })
            .collect::<Vec<_>>();
        if profile {
            println!(
                "[LF+ RgInstance::from_f0_seeded] commit monomial mats (Ajtai seeded): {:?} (kappa×(k*d) = {}×{})",
                t.elapsed(),
                kappa,
                decomp.k * d
            );
        }

        let com = Matrix::hconcat(&comM_f).unwrap();

        let t = std::time::Instant::now();
        let tau = split(&com, n, tau_b, decomp.l);
        if profile {
            println!("[LF+ RgInstance::from_f0_seeded] split tau: {:?}", t.elapsed());
        }

        let t = std::time::Instant::now();
        // Build compact m_tau digits instead of materializing Vec<R> (saves ~70GiB at n=2^27,d=64).
        let m_tau_digits: Arc<Vec<u16>> = {
            #[cfg(feature = "parallel")]
            {
                use rayon::prelude::*;
                Arc::new(
                    tau.par_iter()
                        .map(|&c| (map_tau_digit_to_idx)(c))
                        .collect::<Vec<u16>>(),
                )
            }
            #[cfg(not(feature = "parallel"))]
            {
                Arc::new(tau.iter().copied().map(|c| (map_tau_digit_to_idx)(c)).collect::<Vec<u16>>())
            }
        };
        if profile {
            println!(
                "[LF+ RgInstance::from_f0_seeded] build m_tau digits: {:?}",
                t.elapsed()
            );
        }

        let t = std::time::Instant::now();
        // `tau` is huge (length n); avoid cloning it. Wrap once and share the Arc.
        let tau0: Arc<Vec<R::BaseRing>> = Arc::new(tau);
        let cm_pair = scheme
            .commit_many_const_coeff_base_fast(n, 2, {
                let f0 = f0.clone();
                let tau0 = tau0.clone();
                move |j, out| {
                    // f0 is a prefix; missing entries are implicit zeros.
                    out[0] = f0.get(j).copied().unwrap_or(R::BaseRing::ZERO);
                    out[1] = tau0[j];
                }
            })
            .expect("commit_many_const_coeff_base_fast (f0,tau)");
        let cm_f = cm_pair[0].as_ref().to_vec();
        let C_Mf = cm_pair[1].as_ref().to_vec();
        let cm_mtau = scheme
            .commit_many_with_monomial_digits(n, 1, exp_table.clone(), {
                let digits = m_tau_digits.clone();
                move |j, out| {
                    out[0] = digits[j];
                }
            })
            .expect("commit m_tau (digits)")[0]
            .as_ref()
            .to_vec();
        if profile {
            println!("[LF+ RgInstance::from_f0_seeded] commit f/tau/m_tau: {:?}", t.elapsed());
            println!("[LF+ RgInstance::from_f0_seeded] total: {:?}", t_total.elapsed());
        }
        let fcoms = FComs { cm_f, C_Mf, cm_mtau };

        Self {
            M_f,
            tau: tau0.clone(),
            m_tau: MonomialVec::Digits { digits: m_tau_digits, exp_table: exp_table.clone() },
            f: WitnessVec::ConstCoeffBase {
                values: f0,
                domain_len: n,
            },
            comM_f,
            fcoms,
        }
    }
}

fn absorb_evaluations<R: OverField>(evals: &[DcomEvals<R>], transcript: &mut impl Transcript<R>) {
    evals.iter().for_each(|eval| {
        transcript.absorb_slice(&eval.a.iter().map(|z| R::from(*z)).collect::<Vec<R>>());
        transcript.absorb_slice(&eval.c);
    });
}

/// Precompute eq weights for the first `t` (low) variables (LSB-first).
fn build_eq_low_table<R: PolyRing>(r_low: &[R::BaseRing], one_minus_r_low: &[R::BaseRing]) -> Vec<R::BaseRing>
where
    R::BaseRing: Ring,
{
    debug_assert_eq!(r_low.len(), one_minus_r_low.len());
    let t = r_low.len();
    let mut buf = vec![R::BaseRing::ONE];
    // Expand in the same LSB-first convention used elsewhere.
    // For low bits, we can fold from high-to-low within this slice.
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

#[inline]
fn eq_scale_for_high_bits<R: PolyRing>(
    high: usize,
    r: &[R::BaseRing],
    one_minus_r: &[R::BaseRing],
    t_low: usize,
) -> R::BaseRing
where
    R::BaseRing: Ring,
{
    let mut prod = R::BaseRing::ONE;
    for i in t_low..r.len() {
        let bit = ((high >> (i - t_low)) & 1) == 1;
        prod *= if bit { r[i] } else { one_minus_r[i] };
    }
    prod
}

#[inline]
fn choose_t_low(nvars: usize) -> usize {
    // Keep a tiny table (<= 2^12 = 4096) to avoid big allocations across many chunks.
    nvars.min(12)
}

fn eval_vec_coeffs_at_point_streaming_witness<R: PolyRing>(
    v: &WitnessVec<R>,
    r: &[R::BaseRing],
    one_minus_r: &[R::BaseRing],
) -> Vec<R::BaseRing>
where
    R::BaseRing: Ring,
{
    match v {
        WitnessVec::Ring(vr) => eval_vec_coeffs_at_point_streaming::<R>(vr.as_ref(), r, one_minus_r),
        WitnessVec::ConstCoeffBase { values: v0, .. } => {
            let mut out = vec![R::BaseRing::ZERO; R::dimension()];
            out[0] = dot_base_streaming::<R>(v0.as_ref(), r, one_minus_r);
            out
        }
    }
}

fn dot_ring_streaming_witness<R>(
    v: &WitnessVec<R>,
    r: &[R::BaseRing],
    one_minus_r: &[R::BaseRing],
) -> R
where
    R: PolyRing + From<R::BaseRing>,
    R::BaseRing: Ring,
{
    match v {
        WitnessVec::Ring(vr) => dot_ring_streaming::<R>(vr.as_ref(), r, one_minus_r),
        WitnessVec::ConstCoeffBase { values: v0, .. } => {
            R::from(dot_base_streaming::<R>(v0.as_ref(), r, one_minus_r))
        }
    }
}

fn sparse_mat_vec_eval_ring_streaming_witness<R>(
    m: &SparseMatrix<R>,
    witness: &WitnessVec<R>,
    r: &[R::BaseRing],
    one_minus_r: &[R::BaseRing],
) -> R
where
    R: PolyRing + From<R::BaseRing>,
    R::BaseRing: Ring,
{
    match witness {
        WitnessVec::Ring(vr) => sparse_mat_vec_eval_ring_streaming::<R>(m, vr.as_ref(), r, one_minus_r),
        WitnessVec::ConstCoeffBase { values: v0, .. } => {
            R::from(sparse_mat_vec_eval_ct_streaming::<R>(m, v0.as_ref(), r, one_minus_r))
        }
    }
}

fn dot_base_streaming<R: PolyRing>(
    v: &[R::BaseRing],
    r: &[R::BaseRing],
    one_minus_r: &[R::BaseRing],
) -> R::BaseRing
where
    R::BaseRing: Ring,
{
    debug_assert_eq!(r.len(), one_minus_r.len());
    let nvars = r.len();
    let n = v.len();
    let t = choose_t_low(nvars);
    let low = build_eq_low_table::<R>(&r[..t], &one_minus_r[..t]);
    let low_len = 1usize << t;
    let high_bits = nvars - t;
    // Cap at the actual table length; anything beyond `n` is implicit zero padding.
    let high_len = ((n + low_len - 1) / low_len).min(1usize << high_bits);
    #[cfg(feature = "parallel")]
    {
        (0..high_len)
            .into_par_iter()
            .map(|h| {
                let scale = eq_scale_for_high_bits::<R>(h, r, one_minus_r, t);
                let base = h * low_len;
                let mut acc = R::BaseRing::ZERO;
                for i in 0..low_len {
                    let idx = base + i;
                    if idx >= n {
                        break;
                    }
                    acc += v[idx] * (scale * low[i]);
                }
                acc
            })
            .reduce(|| R::BaseRing::ZERO, |a, b| a + b)
    }
    #[cfg(not(feature = "parallel"))]
    {
        let mut acc = R::BaseRing::ZERO;
        for h in 0..high_len {
            let scale = eq_scale_for_high_bits::<R>(h, r, one_minus_r, t);
            let base = h * low_len;
            for i in 0..low_len {
                let idx = base + i;
                if idx >= n {
                    break;
                }
                acc += v[idx] * (scale * low[i]);
            }
        }
        acc
    }
}

fn dot_ring_streaming<R>(v: &[R], r: &[R::BaseRing], one_minus_r: &[R::BaseRing]) -> R
where
    R: PolyRing + From<R::BaseRing>,
    R::BaseRing: Ring,
{
    debug_assert_eq!(r.len(), one_minus_r.len());
    let nvars = r.len();
    let n = v.len();
    let t = choose_t_low(nvars);
    let low = build_eq_low_table::<R>(&r[..t], &one_minus_r[..t]);
    let low_len = 1usize << t;
    let high_bits = nvars - t;
    let high_len = ((n + low_len - 1) / low_len).min(1usize << high_bits);
    #[cfg(feature = "parallel")]
    {
        (0..high_len)
            .into_par_iter()
            .map(|h| {
                let scale = eq_scale_for_high_bits::<R>(h, r, one_minus_r, t);
                let base = h * low_len;
                let mut acc = R::ZERO;
                for i in 0..low_len {
                    let idx = base + i;
                    if idx >= n {
                        break;
                    }
                    let w = scale * low[i];
                    acc += v[idx] * R::from(w);
                }
                acc
            })
            .reduce(|| R::ZERO, |a, b| a + b)
    }
    #[cfg(not(feature = "parallel"))]
    {
        let mut acc = R::ZERO;
        for h in 0..high_len {
            let scale = eq_scale_for_high_bits::<R>(h, r, one_minus_r, t);
            let base = h * low_len;
            for i in 0..low_len {
                let idx = base + i;
                if idx >= n {
                    break;
                }
                let w = scale * low[i];
                acc += v[idx] * R::from(w);
            }
        }
        acc
    }
}

fn eval_vec_coeffs_at_point_streaming<R: PolyRing>(
    v: &[R],
    r: &[R::BaseRing],
    one_minus_r: &[R::BaseRing],
) -> Vec<R::BaseRing>
where
    R::BaseRing: Ring,
{
    debug_assert_eq!(r.len(), one_minus_r.len());
    let nvars = r.len();
    let n = v.len();
    let t = choose_t_low(nvars);
    let low = build_eq_low_table::<R>(&r[..t], &one_minus_r[..t]);
    let low_len = 1usize << t;
    let high_bits = nvars - t;
    let high_len = ((n + low_len - 1) / low_len).min(1usize << high_bits);
    let d = R::dimension();
    #[cfg(feature = "parallel")]
    {
        (0..d)
            .into_par_iter()
            .map(|j| {
                (0..high_len)
                    .into_par_iter()
                    .map(|h| {
                        let scale = eq_scale_for_high_bits::<R>(h, r, one_minus_r, t);
                        let base = h * low_len;
                        let mut acc = R::BaseRing::ZERO;
                        for i in 0..low_len {
                            let w = scale * low[i];
                            let idx = base + i;
                            if idx >= n {
                                break;
                            }
                            acc += v[idx].coeffs()[j] * w;
                        }
                        acc
                    })
                    .reduce(|| R::BaseRing::ZERO, |a, b| a + b)
            })
            .collect()
    }
    #[cfg(not(feature = "parallel"))]
    {
        let mut out = vec![R::BaseRing::ZERO; d];
        for h in 0..high_len {
            let scale = eq_scale_for_high_bits::<R>(h, r, one_minus_r, t);
            let base = h * low_len;
            for i in 0..low_len {
                let w = scale * low[i];
                let idx = base + i;
                if idx >= n {
                    break;
                }
                let x = &v[idx];
                for j in 0..d {
                    out[j] += x.coeffs()[j] * w;
                }
            }
        }
        out
    }
}

fn sparse_mat_vec_eval_ct_streaming<R: PolyRing>(
    m: &SparseMatrix<R>,
    witness0: &[R::BaseRing],
    r: &[R::BaseRing],
    one_minus_r: &[R::BaseRing],
) -> R::BaseRing
where
    R::BaseRing: Ring,
{
    debug_assert_eq!(r.len(), one_minus_r.len());
    let nvars = r.len();
    let n = m.nrows;
    let t = choose_t_low(nvars);
    let low = build_eq_low_table::<R>(&r[..t], &one_minus_r[..t]);
    let low_len = 1usize << t;
    let high_bits = nvars - t;
    let high_len = ((n + low_len - 1) / low_len).min(1usize << high_bits);
    #[cfg(feature = "parallel")]
    {
        (0..high_len)
            .into_par_iter()
            .map(|h| {
                let scale = eq_scale_for_high_bits::<R>(h, r, one_minus_r, t);
                let base = h * low_len;
                let mut acc = R::BaseRing::ZERO;
                for i in 0..low_len {
                    let row_idx = base + i;
                    if row_idx >= n {
                        break;
                    }
                    let w_row = scale * low[i];
                    let row = &m.coeffs[row_idx];
                    let mut sum0 = R::BaseRing::ZERO;
                    for (coeff, col_idx) in row {
                        if *col_idx < witness0.len() {
                            sum0 += coeff.coeffs()[0] * witness0[*col_idx];
                        }
                    }
                    acc += sum0 * w_row;
                }
                acc
            })
            .reduce(|| R::BaseRing::ZERO, |a, b| a + b)
    }
    #[cfg(not(feature = "parallel"))]
    {
        let mut acc = R::BaseRing::ZERO;
        for h in 0..high_len {
            let scale = eq_scale_for_high_bits::<R>(h, r, one_minus_r, t);
            let base = h * low_len;
            for i in 0..low_len {
                let row_idx = base + i;
                if row_idx >= n {
                    break;
                }
                let w_row = scale * low[i];
                let row = &m.coeffs[row_idx];
                let mut sum0 = R::BaseRing::ZERO;
                for (coeff, col_idx) in row {
                    if *col_idx < witness0.len() {
                        sum0 += coeff.coeffs()[0] * witness0[*col_idx];
                    }
                }
                acc += sum0 * w_row;
            }
        }
        acc
    }
}

fn sparse_mat0_vec_eval_ct_streaming<R: PolyRing>(
    m0: &SparseMatrix<R::BaseRing>,
    witness0: &[R::BaseRing],
    r: &[R::BaseRing],
    one_minus_r: &[R::BaseRing],
) -> R::BaseRing
where
    R::BaseRing: Ring,
{
    debug_assert_eq!(r.len(), one_minus_r.len());
    let nvars = r.len();
    let n = m0.nrows;
    let t = choose_t_low(nvars);
    let low = build_eq_low_table::<R>(&r[..t], &one_minus_r[..t]);
    let low_len = 1usize << t;
    let high_bits = nvars - t;
    let high_len = ((n + low_len - 1) / low_len).min(1usize << high_bits);
    #[cfg(feature = "parallel")]
    {
        (0..high_len)
            .into_par_iter()
            .map(|h| {
                let scale = eq_scale_for_high_bits::<R>(h, r, one_minus_r, t);
                let base = h * low_len;
                let mut acc = R::BaseRing::ZERO;
                for i in 0..low_len {
                    let row_idx = base + i;
                    if row_idx >= n {
                        break;
                    }
                    let w_row = scale * low[i];
                    let row = &m0.coeffs[row_idx];
                    let mut sum0 = R::BaseRing::ZERO;
                    for (coeff0, col_idx) in row {
                        if *col_idx < witness0.len() {
                            sum0 += *coeff0 * witness0[*col_idx];
                        }
                    }
                    acc += sum0 * w_row;
                }
                acc
            })
            .reduce(|| R::BaseRing::ZERO, |a, b| a + b)
    }
    #[cfg(not(feature = "parallel"))]
    {
        let mut acc = R::BaseRing::ZERO;
        for h in 0..high_len {
            let scale = eq_scale_for_high_bits::<R>(h, r, one_minus_r, t);
            let base = h * low_len;
            for i in 0..low_len {
                let row_idx = base + i;
                if row_idx >= n {
                    break;
                }
                let w_row = scale * low[i];
                let row = &m0.coeffs[row_idx];
                let mut sum0 = R::BaseRing::ZERO;
                for (coeff0, col_idx) in row {
                    if *col_idx < witness0.len() {
                        sum0 += *coeff0 * witness0[*col_idx];
                    }
                }
                acc += sum0 * w_row;
            }
        }
        acc
    }
}

fn sparse_mat_vec_eval_ring_streaming<R>(
    m: &SparseMatrix<R>,
    witness: &[R],
    r: &[R::BaseRing],
    one_minus_r: &[R::BaseRing],
) -> R
where
    R: PolyRing + From<R::BaseRing>,
    R::BaseRing: Ring,
{
    debug_assert_eq!(r.len(), one_minus_r.len());
    let nvars = r.len();
    let n = m.nrows;
    let t = choose_t_low(nvars);
    let low = build_eq_low_table::<R>(&r[..t], &one_minus_r[..t]);
    let low_len = 1usize << t;
    let high_bits = nvars - t;
    let high_len = ((n + low_len - 1) / low_len).min(1usize << high_bits);
    #[cfg(feature = "parallel")]
    {
        (0..high_len)
            .into_par_iter()
            .map(|h| {
                let scale = eq_scale_for_high_bits::<R>(h, r, one_minus_r, t);
                let base = h * low_len;
                let mut acc = R::ZERO;
                for i in 0..low_len {
                    let row_idx = base + i;
                    if row_idx >= n {
                        break;
                    }
                    let w_row = scale * low[i];
                    let row = &m.coeffs[row_idx];
                    let mut row_dot = R::ZERO;
                    for (coeff, col_idx) in row {
                        if *col_idx < witness.len() {
                            row_dot += *coeff * witness[*col_idx];
                        }
                    }
                    acc += row_dot * R::from(w_row);
                }
                acc
            })
            .reduce(|| R::ZERO, |a, b| a + b)
    }
    #[cfg(not(feature = "parallel"))]
    {
        let mut acc = R::ZERO;
        for h in 0..high_len {
            let scale = eq_scale_for_high_bits::<R>(h, r, one_minus_r, t);
            let base = h * low_len;
            for i in 0..low_len {
                let row_idx = base + i;
                if row_idx >= n {
                    break;
                }
                let w_row = R::from(scale * low[i]);
                let row = &m.coeffs[row_idx];
                let mut row_dot = R::ZERO;
                for (coeff, col_idx) in row {
                    if *col_idx < witness.len() {
                        row_dot += *coeff * witness[*col_idx];
                    }
                }
                acc += row_dot * w_row;
            }
        }
        acc
    }
}

fn sparse_mat0_vec_eval_ring_streaming<R>(
    m0: &SparseMatrix<R::BaseRing>,
    witness: &[R],
    r: &[R::BaseRing],
    one_minus_r: &[R::BaseRing],
) -> R
where
    R: PolyRing + From<R::BaseRing>,
    R::BaseRing: Ring,
{
    debug_assert_eq!(r.len(), one_minus_r.len());
    let nvars = r.len();
    let n = m0.nrows;
    let t = choose_t_low(nvars);
    let low = build_eq_low_table::<R>(&r[..t], &one_minus_r[..t]);
    let low_len = 1usize << t;
    let high_bits = nvars - t;
    let high_len = ((n + low_len - 1) / low_len).min(1usize << high_bits);
    #[cfg(feature = "parallel")]
    {
        (0..high_len)
            .into_par_iter()
            .map(|h| {
                let scale = eq_scale_for_high_bits::<R>(h, r, one_minus_r, t);
                let base = h * low_len;
                let mut acc = R::ZERO;
                for i in 0..low_len {
                    let row_idx = base + i;
                    if row_idx >= n {
                        break;
                    }
                    let w_row0 = scale * low[i];
                    let row = &m0.coeffs[row_idx];
                    let mut row_dot = R::ZERO;
                    for (coeff0, col_idx) in row {
                        if *col_idx < witness.len() {
                            row_dot += witness[*col_idx] * R::from(*coeff0);
                        }
                    }
                    acc += row_dot * R::from(w_row0);
                }
                acc
            })
            .reduce(|| R::ZERO, |a, b| a + b)
    }
    #[cfg(not(feature = "parallel"))]
    {
        let mut acc = R::ZERO;
        for h in 0..high_len {
            let scale = eq_scale_for_high_bits::<R>(h, r, one_minus_r, t);
            let base = h * low_len;
            for i in 0..low_len {
                let row_idx = base + i;
                if row_idx >= n {
                    break;
                }
                let w_row = R::from(scale * low[i]);
                let row = &m0.coeffs[row_idx];
                let mut row_dot = R::ZERO;
                for (coeff0, col_idx) in row {
                    if *col_idx < witness.len() {
                        row_dot += witness[*col_idx] * *coeff0;
                    }
                }
                acc += row_dot * w_row;
            }
        }
        acc
    }
}

fn sparse_mat_vec_eval_ring_streaming_monomial_digits<R>(
    m: &SparseMatrix<R>,
    digits: &[u16],
    exp_table: &[R],
    r: &[R::BaseRing],
    one_minus_r: &[R::BaseRing],
) -> R
where
    R: PolyRing + From<R::BaseRing>,
    R::BaseRing: Ring,
{
    debug_assert_eq!(r.len(), one_minus_r.len());
    let nvars = r.len();
    let n = m.nrows;
    let t = choose_t_low(nvars);
    let low = build_eq_low_table::<R>(&r[..t], &one_minus_r[..t]);
    let low_len = 1usize << t;
    let high_bits = nvars - t;
    let high_len = ((n + low_len - 1) / low_len).min(1usize << high_bits);
    #[cfg(feature = "parallel")]
    {
        (0..high_len)
            .into_par_iter()
            .map(|h| {
                let scale = eq_scale_for_high_bits::<R>(h, r, one_minus_r, t);
                let base = h * low_len;
                let mut acc = R::ZERO;
                for i in 0..low_len {
                    let row_idx = base + i;
                    if row_idx >= n {
                        break;
                    }
                    let w_row = scale * low[i];
                    let row = &m.coeffs[row_idx];
                    let mut row_dot = R::ZERO;
                    for (coeff, col_idx) in row {
                        let cj = *col_idx;
                        if cj < digits.len() {
                            row_dot += *coeff * exp_table[digits[cj] as usize];
                        }
                    }
                    acc += row_dot * R::from(w_row);
                }
                acc
            })
            .reduce(|| R::ZERO, |a, b| a + b)
    }
    #[cfg(not(feature = "parallel"))]
    {
        let mut acc = R::ZERO;
        for h in 0..high_len {
            let scale = eq_scale_for_high_bits::<R>(h, r, one_minus_r, t);
            let base = h * low_len;
            for i in 0..low_len {
                let row_idx = base + i;
                if row_idx >= n {
                    break;
                }
                let w_row = R::from(scale * low[i]);
                let row = &m.coeffs[row_idx];
                let mut row_dot = R::ZERO;
                for (coeff, col_idx) in row {
                    let cj = *col_idx;
                    if cj < digits.len() {
                        row_dot += *coeff * exp_table[digits[cj] as usize];
                    }
                }
                acc += row_dot * w_row;
            }
        }
        acc
    }
}

fn sparse_mat0_vec_eval_ring_streaming_monomial_digits<R>(
    m0: &SparseMatrix<R::BaseRing>,
    digits: &[u16],
    exp_table: &[R],
    r: &[R::BaseRing],
    one_minus_r: &[R::BaseRing],
) -> R
where
    R: PolyRing + From<R::BaseRing>,
    R::BaseRing: Ring,
{
    debug_assert_eq!(r.len(), one_minus_r.len());
    let nvars = r.len();
    let n = m0.nrows;
    let t = choose_t_low(nvars);
    let low = build_eq_low_table::<R>(&r[..t], &one_minus_r[..t]);
    let low_len = 1usize << t;
    let high_bits = nvars - t;
    let high_len = ((n + low_len - 1) / low_len).min(1usize << high_bits);
    #[cfg(feature = "parallel")]
    {
        (0..high_len)
            .into_par_iter()
            .map(|h| {
                let scale = eq_scale_for_high_bits::<R>(h, r, one_minus_r, t);
                let base = h * low_len;
                let mut acc = R::ZERO;
                for i in 0..low_len {
                    let row_idx = base + i;
                    if row_idx >= n {
                        break;
                    }
                    let w_row0 = scale * low[i];
                    let row = &m0.coeffs[row_idx];
                    let mut row_dot = R::ZERO;
                    for (coeff0, col_idx) in row {
                        let cj = *col_idx;
                        if cj < digits.len() {
                            row_dot += exp_table[digits[cj] as usize] * R::from(*coeff0);
                        }
                    }
                    acc += row_dot * R::from(w_row0);
                }
                acc
            })
            .reduce(|| R::ZERO, |a, b| a + b)
    }
    #[cfg(not(feature = "parallel"))]
    {
        let mut acc = R::ZERO;
        for h in 0..high_len {
            let scale = eq_scale_for_high_bits::<R>(h, r, one_minus_r, t);
            let base = h * low_len;
            for i in 0..low_len {
                let row_idx = base + i;
                if row_idx >= n {
                    break;
                }
                let w_row = R::from(scale * low[i]);
                let row = &m0.coeffs[row_idx];
                let mut row_dot = R::ZERO;
                for (coeff0, col_idx) in row {
                    let cj = *col_idx;
                    if cj < digits.len() {
                        row_dot += exp_table[digits[cj] as usize] * *coeff0;
                    }
                }
                acc += row_dot * w_row;
            }
        }
        acc
    }
}

fn sparse_mat0_vec_eval_ring_streaming_witness<R>(
    m0: &SparseMatrix<R::BaseRing>,
    witness: &WitnessVec<R>,
    r: &[R::BaseRing],
    one_minus_r: &[R::BaseRing],
) -> R
where
    R: PolyRing + From<R::BaseRing>,
    R::BaseRing: Ring,
{
    match witness {
        WitnessVec::Ring(vr) => sparse_mat0_vec_eval_ring_streaming::<R>(m0, vr.as_ref(), r, one_minus_r),
        WitnessVec::ConstCoeffBase { values: v0, .. } => {
            R::from(sparse_mat0_vec_eval_ct_streaming::<R>(m0, v0.as_ref(), r, one_minus_r))
        }
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::PrimeField;
    use ark_std::{log2, Zero};
    use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
    use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as R;
    use std::sync::Arc;

    use super::*;
    use crate::transcript::PoseidonTranscript;

    #[test]
    fn test_range_check() {
        // f: [
        // 2 + 5X
        // 4 + X^2
        // ]
        let mut f = vec![R::zero(); 1 << 15];
        f[0].coeffs_mut()[0] = 2u128.into();
        f[0].coeffs_mut()[1] = 5u128.into();
        f[1].coeffs_mut()[0] = 4u128.into();
        f[1].coeffs_mut()[2] = 1u128.into();

        let n = f.len();
        let kappa = 1;
        let b = (R::dimension() / 2) as u128;
        let k = 2;
        // log_d' (q)
        let l = ((<<R as PolyRing>::BaseRing>::MODULUS.0[0] as f64).ln()
            / ((R::dimension() / 2) as f64).ln())
        .ceil() as usize;

        let A = Matrix::<R>::rand(&mut ark_std::test_rng(), kappa, n);

        let dparams = DecompParameters { b, k, l };
        let instance = RgInstance::from_f(f.clone(), &A, &dparams);

        let rg = Rg {
            nvars: log2(n) as usize,
            instances: vec![instance],
            dparams,
        };

        let mut ts = PoseidonTranscript::empty::<PC>();
        let dcom = rg.range_check(&Vec::<Arc<SparseMatrix<R>>>::new(), &mut ts);

        let mut ts = PoseidonTranscript::empty::<PC>();
        dcom.verify(&mut ts).unwrap();
    }

    #[test]
    fn test_range_check_mm() {
        // f: [
        // 2 + 5X
        // 4 + X^2
        // ]
        let n = 1 << 15;
        let mut f = vec![R::zero(); n];
        f[0].coeffs_mut()[0] = 2u128.into();
        f[0].coeffs_mut()[1] = 5u128.into();
        f[1].coeffs_mut()[0] = 4u128.into();
        f[1].coeffs_mut()[2] = 1u128.into();

        let mut m = SparseMatrix::identity(n);
        m.coeffs[0][0].0 = 2u128.into();
        let M = vec![m];

        let kappa = 1;
        let b = (R::dimension() / 2) as u128;
        let k = 2;
        // log_d' (q)
        let l = ((<<R as PolyRing>::BaseRing>::MODULUS.0[0] as f64).ln()
            / ((R::dimension() / 2) as f64).ln())
        .ceil() as usize;

        let A = Matrix::<R>::rand(&mut ark_std::test_rng(), kappa, n);

        let dparams = DecompParameters { b, k, l };
        let instance = RgInstance::from_f(f.clone(), &A, &dparams);

        let rg = Rg {
            nvars: log2(n) as usize,
            instances: vec![instance],
            dparams,
        };

        let mut ts = PoseidonTranscript::empty::<PC>();
        let M: Vec<Arc<SparseMatrix<R>>> = M.into_iter().map(Arc::new).collect();
        let dcom = rg.range_check(&M, &mut ts);

        let mut ts = PoseidonTranscript::empty::<PC>();
        dcom.verify(&mut ts).unwrap();
    }
}
