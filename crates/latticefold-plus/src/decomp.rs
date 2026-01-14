use ark_std::log2;
use ark_ff::{BigInteger, PrimeField};
use latticefold::commitment::AjtaiCommitmentScheme;
use stark_rings::{
    balanced_decomposition::{recompose, Decompose},
    OverField, PolyRing, Ring, Zq,
};
use stark_rings_linalg::{Matrix, SparseMatrix};
use std::time::Instant;
use std::sync::Arc;

use crate::lin::{LinB, LinBX};
use crate::rgchk::WitnessVec;
use crate::utils::maybe_print_rss;

pub type RxR<R> = (R, R);

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// Packed representation for a length-`n` ring vector whose coefficients are small.
///
/// Layout: `coeffs[j * d + k]` stores the `k`-th base coefficient for entry `j`.
/// We use this to avoid allocating `Vec<R>` for digit vectors during decomposition.
#[derive(Clone, Debug)]
enum PackedDigitVec<BR: Ring> {
    I16 { coeffs: Vec<i16>, d: usize, n: usize, _br: core::marker::PhantomData<BR> },
    I32 { coeffs: Vec<i32>, d: usize, n: usize, _br: core::marker::PhantomData<BR> },
}

impl<BR: Ring> PackedDigitVec<BR> {
    #[inline]
    fn new_i16(n: usize, d: usize) -> Self {
        PackedDigitVec::I16 { coeffs: vec![0i16; n * d], d, n, _br: core::marker::PhantomData }
    }

    #[inline]
    fn upgrade_to_i32(&mut self) {
        if let PackedDigitVec::I16 { coeffs, d, n, .. } = self {
            let mut out = Vec::<i32>::with_capacity(coeffs.len());
            out.extend(coeffs.iter().map(|&x| x as i32));
            *self = PackedDigitVec::I32 { coeffs: out, d: *d, n: *n, _br: core::marker::PhantomData };
        }
    }

    #[inline]
    fn dims(&self) -> (usize, usize) {
        match self {
            PackedDigitVec::I16 { d, n, .. } => (*d, *n),
            PackedDigitVec::I32 { d, n, .. } => (*d, *n),
        }
    }

    #[inline]
    fn get_i64(&self, j: usize, k: usize) -> i64 {
        match self {
            PackedDigitVec::I16 { coeffs, d, .. } => coeffs[j * *d + k] as i64,
            PackedDigitVec::I32 { coeffs, d, .. } => coeffs[j * *d + k] as i64,
        }
    }

    #[inline]
    fn fill_ring_at<Rr>(&self, j: usize, out: &mut Rr)
    where
        Rr: PolyRing<BaseRing = BR>,
        BR: Ring + Copy,
    {
        let oc = out.coeffs_mut();
        let d = oc.len();
        let (_, n) = self.dims();
        debug_assert!(j < n);
        debug_assert_eq!(d, self.dims().0);
        for k in 0..d {
            let v = self.get_i64(j, k);
            if v >= 0 {
                oc[k] = BR::from(v as u128);
            } else {
                oc[k] = -BR::from((-v) as u128);
            }
        }
    }

    /// dst += scale * F1[j], where `F1` is represented in packed coefficients.
    #[inline]
    fn mul_add_into_ring<Rr>(&self, dst: &mut Rr, j: usize, scale: BR)
    where
        Rr: PolyRing<BaseRing = BR>,
        BR: Ring + Copy,
    {
        if scale == BR::ZERO {
            return;
        }
        let dc = dst.coeffs_mut();
        let d = dc.len();
        debug_assert_eq!(d, self.dims().0);
        for k in 0..d {
            let v = self.get_i64(j, k);
            if v >= 0 {
                dc[k] += BR::from(v as u128) * scale;
            } else {
                dc[k] -= BR::from((-v) as u128) * scale;
            }
        }
    }

    #[inline]
    fn set_coeffs_from_ring<Rr>(&mut self, j: usize, x: &Rr)
    where
        Rr: PolyRing<BaseRing = BR>,
        BR: PrimeField,
        <BR as PrimeField>::BigInt: BigInteger,
    {
        #[inline]
        fn br_to_i64_balanced<BRr: PrimeField>(x: BRr) -> i64
        where
            <BRr as PrimeField>::BigInt: BigInteger,
        {
            let rep = x.into_bigint();
            let limbs = rep.as_ref();
            // Convert canonical rep to u128 (we only expect small moduli in this codebase).
            let mut v: u128 = 0;
            let take = core::cmp::min(limbs.len(), 2);
            for i in 0..take {
                v |= (limbs[i] as u128) << (64 * i);
            }
            for i in take..limbs.len() {
                debug_assert_eq!(limbs[i], 0, "unexpected large modulus limb in coefficient rep");
            }

            let modulus = BRr::MODULUS;
            let mod_limbs = modulus.as_ref();
            let mut q: u128 = 0;
            let take_q = core::cmp::min(mod_limbs.len(), 2);
            for i in 0..take_q {
                q |= (mod_limbs[i] as u128) << (64 * i);
            }
            for i in take_q..mod_limbs.len() {
                debug_assert_eq!(mod_limbs[i], 0, "unexpected large modulus limb");
            }

            // Map to balanced representative in (-q/2, q/2].
            let half = q >> 1;
            if v <= half {
                v as i64
            } else {
                (v as i128 - q as i128) as i64
            }
        }

        let xc = x.coeffs();
        loop {
            match self {
                PackedDigitVec::I16 { coeffs, d, .. } => {
                    debug_assert_eq!(xc.len(), *d);
                    let base = j * *d;
                    let mut ok = true;
                    for k in 0..*d {
                        let v = br_to_i64_balanced::<BR>(xc[k]);
                        if v < i16::MIN as i64 || v > i16::MAX as i64 {
                            ok = false;
                            break;
                        }
                        coeffs[base + k] = v as i16;
                    }
                    if ok {
                        return;
                    }
                    self.upgrade_to_i32();
                }
                PackedDigitVec::I32 { coeffs, d, .. } => {
                    debug_assert_eq!(xc.len(), *d);
                    let base = j * *d;
                    for k in 0..*d {
                        let v = br_to_i64_balanced::<BR>(xc[k]);
                        if v < i32::MIN as i64 || v > i32::MAX as i64 {
                            panic!("packed digit overflow: coefficient doesn't fit i32");
                        }
                        coeffs[base + k] = v as i32;
                    }
                    return;
                }
            }
        }
    }
}

#[inline]
fn choose_t_low(nvars: usize) -> usize {
    // Keep small tables; high part is 2^(nvars-t_low).
    nvars.min(12)
}

/// Build eq weights for the low `t` variables (LSB-first), returned as a length-2^t table.
fn build_eq_low_table_base<BR: Ring>(r_low: &[BR], one_minus_r_low: &[BR]) -> Vec<BR> {
    debug_assert_eq!(r_low.len(), one_minus_r_low.len());
    let t = r_low.len();
    let mut buf = vec![BR::ONE];
    for i in (0..t).rev() {
        let ri = r_low[i];
        let omi = one_minus_r_low[i];
        let mut res = vec![BR::ZERO; buf.len() << 1];
        for (j, out) in res.iter_mut().enumerate() {
            let bi = buf[j >> 1];
            *out = if (j & 1) == 0 { bi * omi } else { bi * ri };
        }
        buf = res;
    }
    buf
}

/// Precompute scale factors for the high bits (t..nvars-1): scale_high[high] = Π_i (bit? r[i] : (1-r[i])).
fn build_scale_high_base<BR: Ring>(r: &[BR], one_minus_r: &[BR], t_low: usize) -> Vec<BR> {
    debug_assert_eq!(r.len(), one_minus_r.len());
    let nvars = r.len();
    let high_bits = nvars - t_low;
    let high_len = 1usize << high_bits;
    let mut out = vec![BR::ONE; high_len];
    for h in 0..high_len {
        let mut prod = BR::ONE;
        for i in t_low..nvars {
            let bit = ((h >> (i - t_low)) & 1) == 1;
            prod *= if bit { r[i] } else { one_minus_r[i] };
        }
        out[h] = prod;
    }
    out
}

#[inline]
fn eq_at_base<BR: Ring>(idx: usize, low: &[BR], scale_high: &[BR], low_mask: usize, t_low: usize) -> BR {
    let low_idx = idx & low_mask;
    let high_idx = idx >> t_low;
    scale_high[high_idx] * low[low_idx]
}

#[derive(Debug)]
pub struct Decomp<'a, R> {
    pub f: Vec<R>,
    pub r: Vec<(R, R)>,
    pub M: &'a [Arc<SparseMatrix<R>>],
}

/// Decomposition prover for the const-coeff/SP1 regime where the external matrices are stored over
/// the base ring.
#[derive(Debug)]
pub struct DecompBase<'a, R: PolyRing> {
    pub f: Vec<R>,
    pub r: Vec<(R, R)>,
    pub M0: &'a [Arc<SparseMatrix<R::BaseRing>>],
}

#[derive(Clone, Debug)]
pub struct DecompProof<R> {
    /// C = com(F)
    pub C: (Vec<R>, Vec<R>), // kappa x 2
    pub v: (Vec<RxR<R>>, Vec<RxR<R>>), // (v(0), v(1))
}

impl<R: PolyRing> Decomp<'_, R>
where
    R: Decompose + OverField,
    R::BaseRing: Zq,
{
    pub fn decompose(self, A: &Matrix<R>, B: u128) -> ((LinB<R>, LinB<R>), DecompProof<R>) {
        let profile = std::env::var("LF_PLUS_PROFILE").ok().as_deref() == Some("1");
        let t_total = Instant::now();

        let nvars = log2(A.ncols) as usize;
        // In-place 2-digit decomposition:
        // reuse the original `f` buffer for F0 to avoid holding `f + F0 + F1` at once.
        let mut F0 = self.f;
        let n = F0.len();
        let mut F1 = vec![R::ZERO; n];
        #[cfg(feature = "parallel")]
        {
            const CHUNK: usize = 1 << 14;
            F0.par_chunks_mut(CHUNK)
                .zip(F1.par_chunks_mut(CHUNK))
                .for_each_init(|| vec![R::ZERO; 2], |tmp, (c0, c1)| {
                    for i in 0..c0.len() {
                        let orig = std::mem::replace(&mut c0[i], R::ZERO);
                        orig.decompose_to(B, tmp);
                        c0[i] = tmp[0];
                        c1[i] = tmp[1];
                    }
                });
        }
        #[cfg(not(feature = "parallel"))]
        {
            let mut tmp = vec![R::ZERO; 2];
            for i in 0..n {
                let orig = std::mem::replace(&mut F0[i], R::ZERO);
                orig.decompose_to(B, &mut tmp);
                F0[i] = tmp[0];
                F1[i] = tmp[1];
            }
        }

        let r_a = self.r.iter().map(|rr| rr.0).collect::<Vec<_>>();
        let r_b = self.r.iter().map(|rr| rr.1).collect::<Vec<_>>();

        #[inline]
        fn is_identity_matrix<Rr: PolyRing>(m: &SparseMatrix<Rr>) -> bool {
            if m.nrows != m.ncols {
                return false;
            }
            // Fast reject: identity must have exactly one entry per row.
            if m.coeffs.len() != m.nrows {
                return false;
            }
            for (i, row) in m.coeffs.iter().enumerate() {
                if row.len() != 1 {
                    return false;
                }
                let (c, j) = row[0];
                if j != i {
                    return false;
                }
                if c != Rr::ONE {
                    return false;
                }
            }
            true
        }

        // Build the multilinear “equality” weights for evaluating a vector of length 2^n at point r.
        // For a multilinear extension with evaluations `f[x]` (x in {0,1}^n), we have:
        //   f(r) = Σ_x f[x] * eq_r[x],
        // where eq_r[x] = Π_j (x_j ? r_j : (1 - r_j)).
        #[inline]
        fn eq_weights<Rr: PolyRing>(r: &[Rr]) -> Vec<Rr> {
            // Match `DenseMultilinearExtension::evaluate` variable ordering: it folds evaluations
            // by combining consecutive pairs per coordinate, which corresponds to iterating the
            // point coordinates from last to first (LSB-first indexing in the evaluation vector).
            let nvars = r.len();
            let n = 1usize << nvars;

            // Double-buffer to preserve the *interleaved* layout:
            // after each variable, weights are [w0*(1-r), w0*r, w1*(1-r), w1*r, ...],
            // matching the evaluation indexing used by DenseMultilinearExtension.
            let mut cur = vec![Rr::ZERO; n];
            let mut next = vec![Rr::ZERO; n];
            cur[0] = Rr::ONE;

            let mut len = 1usize;
            let mut cur_is_cur = true;
            for &rj in r.iter().rev() {
                let om = Rr::ONE - rj;
                let (src, dst) = if cur_is_cur {
                    (&cur[..len], &mut next[..(2 * len)])
                } else {
                    (&next[..len], &mut cur[..(2 * len)])
                };

                #[cfg(feature = "parallel")]
                {
                    dst.par_chunks_mut(2)
                        .zip(src.par_iter())
                        .for_each(|(pair, &wi)| {
                            pair[0] = wi * om;
                            pair[1] = wi * rj;
                        });
                }
                #[cfg(not(feature = "parallel"))]
                {
                    for (i, &wi) in src.iter().enumerate() {
                        dst[2 * i] = wi * om;
                        dst[2 * i + 1] = wi * rj;
                    }
                }

                len <<= 1;
                cur_is_cur = !cur_is_cur;
            }

            if cur_is_cur { cur } else { next }
        }

        #[inline]
        fn dot_with_eq<Rr: PolyRing>(f: &[Rr], eq: &[Rr]) -> Rr {
            debug_assert_eq!(f.len(), eq.len());
            #[cfg(feature = "parallel")]
            {
                f.par_iter()
                    .zip(eq.par_iter())
                    .map(|(&fx, &wx)| fx * wx)
                    .reduce(|| Rr::ZERO, |a, b| a + b)
            }
            #[cfg(not(feature = "parallel"))]
            {
                f.iter()
                    .zip(eq.iter())
                    .fold(Rr::ZERO, |acc, (&fx, &wx)| acc + fx * wx)
            }
        }


        // Precompute eq-weights once per point; shared across both Fi branches.
        let t_eq = Instant::now();
        let eq_a = eq_weights::<R>(&r_a);
        let eq_b = eq_weights::<R>(&r_b);
        if profile {
            println!("[LF+ Decomp::decompose] eq_weights: {:?}", t_eq.elapsed());
        }

        #[inline]
        fn eval_sparse_mat_two_vecs_at_two_points<Rr: PolyRing>(
            m: &SparseMatrix<Rr>,
            f0: &[Rr],
            f1: &[Rr],
            eq_a: &[Rr],
            eq_b: &[Rr],
        ) -> ((Rr, Rr), (Rr, Rr)) {
            debug_assert_eq!(m.ncols, f0.len());
            debug_assert_eq!(m.ncols, f1.len());
            debug_assert_eq!(m.nrows, eq_a.len());
            debug_assert_eq!(m.nrows, eq_b.len());

            #[cfg(feature = "parallel")]
            {
                m.coeffs
                    .par_iter()
                    .enumerate()
                    .map(|(row_idx, row)| {
                        let mut row_dot0 = Rr::ZERO;
                        let mut row_dot1 = Rr::ZERO;
                        for (coeff, col_idx) in row {
                            if *col_idx < f0.len() {
                                let c = *coeff;
                                let j = *col_idx;
                                row_dot0 += c * f0[j];
                                row_dot1 += c * f1[j];
                            }
                        }
                        let wa = eq_a[row_idx];
                        let wb = eq_b[row_idx];
                        ((wa * row_dot0, wb * row_dot0), (wa * row_dot1, wb * row_dot1))
                    })
                    .reduce(
                        || ((Rr::ZERO, Rr::ZERO), (Rr::ZERO, Rr::ZERO)),
                        |((a00, b00), (a10, b10)), ((a01, b01), (a11, b11))| {
                            ((a00 + a01, b00 + b01), (a10 + a11, b10 + b11))
                        },
                    )
            }
            #[cfg(not(feature = "parallel"))]
            {
                m.coeffs
                    .iter()
                    .enumerate()
                    .fold(
                        ((Rr::ZERO, Rr::ZERO), (Rr::ZERO, Rr::ZERO)),
                        |((a00, b00), (a10, b10)), (row_idx, row)| {
                            let mut row_dot0 = Rr::ZERO;
                            let mut row_dot1 = Rr::ZERO;
                            for (coeff, col_idx) in row {
                                if *col_idx < f0.len() {
                                    let c = *coeff;
                                    let j = *col_idx;
                                    row_dot0 += c * f0[j];
                                    row_dot1 += c * f1[j];
                                }
                            }
                            let wa = eq_a[row_idx];
                            let wb = eq_b[row_idx];
                            (
                                (a00 + wa * row_dot0, b00 + wb * row_dot0),
                                (a10 + wa * row_dot1, b10 + wb * row_dot1),
                            )
                        },
                    )
            }
        }

        // Variant that computes both v0 and v1 in one pass over matrices (better cache reuse).
        let vi_calc_pair = || -> (Vec<(R, R)>, Vec<(R, R)>) {
            let t_fv = Instant::now();
            let fv0 = (dot_with_eq::<R>(&F0, &eq_a), dot_with_eq::<R>(&F0, &eq_b));
            let fv1 = (dot_with_eq::<R>(&F1, &eq_a), dot_with_eq::<R>(&F1, &eq_b));
            if profile {
                println!("[LF+ Decomp::decompose] fv(dot_with_eq) both: {:?}", t_fv.elapsed());
            }

            let mut v0 = Vec::with_capacity(1 + self.M.len());
            let mut v1 = Vec::with_capacity(1 + self.M.len());
            v0.push(fv0);
            v1.push(fv1);

            let t_mats = Instant::now();
            for M_i in self.M {
                let M_i = M_i.as_ref();
                if is_identity_matrix::<R>(M_i) {
                    v0.push(fv0);
                    v1.push(fv1);
                } else {
                    let (m0, m1) = eval_sparse_mat_two_vecs_at_two_points::<R>(
                        M_i, &F0, &F1, &eq_a, &eq_b,
                    );
                    v0.push(m0);
                    v1.push(m1);
                }
            }
            if profile {
                println!(
                    "[LF+ Decomp::decompose] mats(eval_sparse_mat_two_vecs_at_two_points): {:?} (Mlen={})",
                    t_mats.elapsed(),
                    self.M.len()
                );
            }
            (v0, v1)
        };

        if profile {
            println!(
                "[LF+ Decomp::decompose] setup+split: {:?} (nvars={}, Mlen={})",
                t_total.elapsed(),
                nvars,
                self.M.len()
            );
        }

        let t = Instant::now();
        let (v0, v1) = vi_calc_pair();
        if profile {
            println!("[LF+ Decomp::decompose] compute v0/v1: {:?}", t.elapsed());
        }

        let t = Instant::now();
        let (C0, C1) = {
            #[cfg(feature = "parallel")]
            {
                rayon::join(|| A.try_mul_vec(&F0).unwrap(), || A.try_mul_vec(&F1).unwrap())
            }
            #[cfg(not(feature = "parallel"))]
            {
                (A.try_mul_vec(&F0).unwrap(), A.try_mul_vec(&F1).unwrap())
            }
        };
        if profile {
            println!("[LF+ Decomp::decompose] commitments C0/C1: {:?}", t.elapsed());
            println!("[LF+ Decomp::decompose] total: {:?}", t_total.elapsed());
        }

        let linb0 = LinB {
            x: LinBX {
                cm_f: C0.clone(),
                r: self.r.clone(),
                v: v0.clone(),
            },
            f: WitnessVec::Ring(Arc::new(F0)),
        };
        let linb1 = LinB {
            x: LinBX {
                cm_f: C1.clone(),
                r: self.r.clone(),
                v: v1.clone(),
            },
            f: WitnessVec::Ring(Arc::new(F1)),
        };
        let proof = DecompProof {
            C: (C0, C1),
            v: (v0, v1),
        };

        ((linb0, linb1), proof)
    }

    /// Same as [`Decomp::decompose`], but commits using a seeded implicit Ajtai matrix.
    ///
    /// This avoids materializing a `kappa × n` dense matrix. The verifier-side checks are unchanged.
    pub fn decompose_seeded(
        self,
        scheme: &AjtaiCommitmentScheme<R>,
        B: u128,
    ) -> ((LinB<R>, LinB<R>), DecompProof<R>) {
        let profile = std::env::var("LF_PLUS_PROFILE").ok().as_deref() == Some("1");
        let t_total = Instant::now();
        maybe_print_rss("decomp_seeded: start");

        let nvars = log2(scheme.width()) as usize;
        // In-place 2-digit decomposition:
        // reuse the original `f` buffer for F0 to avoid holding `f + F0 + F1` at once.
        let mut F0 = self.f;
        let n = F0.len();
        let mut F1 = vec![R::ZERO; n];
        #[cfg(feature = "parallel")]
        {
            const CHUNK: usize = 1 << 14;
            F0.par_chunks_mut(CHUNK)
                .zip(F1.par_chunks_mut(CHUNK))
                .for_each_init(|| vec![R::ZERO; 2], |tmp, (c0, c1)| {
                    for i in 0..c0.len() {
                        let orig = std::mem::replace(&mut c0[i], R::ZERO);
                        orig.decompose_to(B, tmp);
                        c0[i] = tmp[0];
                        c1[i] = tmp[1];
                    }
                });
        }
        #[cfg(not(feature = "parallel"))]
        {
            let mut tmp = vec![R::ZERO; 2];
            for i in 0..n {
                let orig = std::mem::replace(&mut F0[i], R::ZERO);
                orig.decompose_to(B, &mut tmp);
                F0[i] = tmp[0];
                F1[i] = tmp[1];
            }
        }
        maybe_print_rss("decomp_seeded: after decompose_to_vec");

        let r_a = self.r.iter().map(|rr| rr.0).collect::<Vec<_>>();
        let r_b = self.r.iter().map(|rr| rr.1).collect::<Vec<_>>();

        #[inline]
        fn is_identity_matrix<Rr: PolyRing>(m: &SparseMatrix<Rr>) -> bool {
            if m.nrows != m.ncols {
                return false;
            }
            if m.coeffs.len() != m.nrows {
                return false;
            }
            for (i, row) in m.coeffs.iter().enumerate() {
                if row.len() != 1 {
                    return false;
                }
                let (c, j) = row[0];
                if j != i {
                    return false;
                }
                if c != Rr::ONE {
                    return false;
                }
            }
            true
        }

        #[inline]
        fn eq_weights<Rr: PolyRing>(r: &[Rr]) -> Vec<Rr> {
            let nvars = r.len();
            let n = 1usize << nvars;
            let mut cur = vec![Rr::ZERO; n];
            let mut next = vec![Rr::ZERO; n];
            cur[0] = Rr::ONE;

            let mut len = 1usize;
            let mut cur_is_cur = true;
            for &rj in r.iter().rev() {
                let om = Rr::ONE - rj;
                let (src, dst) = if cur_is_cur {
                    (&cur[..len], &mut next[..(2 * len)])
                } else {
                    (&next[..len], &mut cur[..(2 * len)])
                };

                #[cfg(feature = "parallel")]
                {
                    dst.par_chunks_mut(2)
                        .zip(src.par_iter())
                        .for_each(|(pair, &wi)| {
                            pair[0] = wi * om;
                            pair[1] = wi * rj;
                        });
                }
                #[cfg(not(feature = "parallel"))]
                {
                    for (i, &wi) in src.iter().enumerate() {
                        dst[2 * i] = wi * om;
                        dst[2 * i + 1] = wi * rj;
                    }
                }

                len <<= 1;
                cur_is_cur = !cur_is_cur;
            }
            if cur_is_cur { cur } else { next }
        }

        #[inline]
        fn is_const_coeff_ring<Rr: PolyRing>(x: &Rr) -> bool {
            x.coeffs()
                .iter()
                .skip(1)
                .all(|c| *c == <Rr as PolyRing>::BaseRing::ZERO)
        }

        #[inline]
        fn eval_sparse_mat_two_vecs_at_two_points<Rr: PolyRing>(
            m: &SparseMatrix<Rr>,
            f0: &[Rr],
            f1: &[Rr],
            eq_a: &[Rr],
            eq_b: &[Rr],
        ) -> (RxR<Rr>, RxR<Rr>) {
            debug_assert_eq!(m.ncols, f0.len());
            debug_assert_eq!(m.ncols, f1.len());
            debug_assert_eq!(m.nrows, eq_a.len());
            debug_assert_eq!(m.nrows, eq_b.len());

            let mut out0 = (Rr::ZERO, Rr::ZERO);
            let mut out1 = (Rr::ZERO, Rr::ZERO);
            for (row, terms) in m.coeffs.iter().enumerate() {
                let wa = eq_a[row];
                let wb = eq_b[row];
                if wa == Rr::ZERO && wb == Rr::ZERO {
                    continue;
                }
                let mut s0 = Rr::ZERO;
                let mut s1 = Rr::ZERO;
                for (c, j) in terms {
                    s0 += *c * f0[*j];
                    s1 += *c * f1[*j];
                }
                out0.0 += wa * s0;
                out0.1 += wb * s0;
                out1.0 += wa * s1;
                out1.1 += wb * s1;
            }
            (out0, out1)
        }

        let vi_calc_pair = || {
   
            let t_eq = Instant::now();
            let (eq_a_ring, eq_b_ring, eq_plan_a, eq_plan_b) =
                if r_a.iter().all(is_const_coeff_ring::<R>) && r_b.iter().all(is_const_coeff_ring::<R>) {
                    let r_a0 = r_a.iter().map(|x| x.coeffs()[0]).collect::<Vec<_>>();
                    let r_b0 = r_b.iter().map(|x| x.coeffs()[0]).collect::<Vec<_>>();

                    let one_minus_a0 = r_a0.iter().copied().map(|x| R::BaseRing::ONE - x).collect::<Vec<_>>();
                    let one_minus_b0 = r_b0.iter().copied().map(|x| R::BaseRing::ONE - x).collect::<Vec<_>>();
                    let t_low = choose_t_low(r_a0.len());

                    let low_a = build_eq_low_table_base::<R::BaseRing>(&r_a0[..t_low], &one_minus_a0[..t_low]);
                    let low_b = build_eq_low_table_base::<R::BaseRing>(&r_b0[..t_low], &one_minus_b0[..t_low]);
                    let low_mask = (1usize << t_low) - 1;
                    let scale_a = build_scale_high_base::<R::BaseRing>(&r_a0, &one_minus_a0, t_low);
                    let scale_b = build_scale_high_base::<R::BaseRing>(&r_b0, &one_minus_b0, t_low);

                    (None, None, Some((t_low, low_mask, low_a, scale_a)), Some((t_low, low_mask, low_b, scale_b)))
                } else {
                    (Some(eq_weights::<R>(&r_a)), Some(eq_weights::<R>(&r_b)), None, None)
                };
            maybe_print_rss("decomp_seeded: after eq_weights");
            if profile {
                println!(
                    "[LF+ Decomp::decompose_seeded] eq_weights: {:?} (nvars={})",
                    t_eq.elapsed(),
                    nvars
                );
            }

            #[inline]
            fn dot_with_eq<Rr: PolyRing>(f: &[Rr], eq: &[Rr]) -> Rr {
                debug_assert_eq!(f.len(), eq.len());
                #[cfg(feature = "parallel")]
                {
                    f.par_iter()
                        .zip(eq.par_iter())
                        .map(|(&fx, &wx)| fx * wx)
                        .reduce(|| Rr::ZERO, |a, b| a + b)
                }
                #[cfg(not(feature = "parallel"))]
                {
                    f.iter()
                        .zip(eq.iter())
                        .fold(Rr::ZERO, |acc, (&fx, &wx)| acc + fx * wx)
                }
            }

            let t_fv = Instant::now();
            // Base term corresponds to the "no-matrix" entry in `vo`: evaluation of g itself.
            // We need both evaluation points, so we compute both dot-products for both Fi.
            let (fv0, fv1) = if let (Some((t_low, low_mask, low_a, scale_a)), Some((_, _, low_b, scale_b))) =
                (eq_plan_a.as_ref(), eq_plan_b.as_ref())
            {
                // Stream eq weights on the fly: no length-2^n table allocation.
                // IMPORTANT: fuse all four dot-products into a single pass over i to avoid
                // traversing 2^n entries four times.
                let n = F0.len();
                let eval_at = |idx: usize, low: &[R::BaseRing], scale: &[R::BaseRing]| -> R::BaseRing {
                    eq_at_base::<R::BaseRing>(idx, low, scale, *low_mask, *t_low)
                };
                #[cfg(feature = "parallel")]
                let (f0_a, f0_b, f1_a, f1_b) = (0..n)
                    .into_par_iter()
                    .fold(
                        || (R::ZERO, R::ZERO, R::ZERO, R::ZERO),
                        |(mut a0, mut b0, mut a1, mut b1), i| {
                            let wa = eval_at(i, low_a, scale_a);
                            let wb = eval_at(i, low_b, scale_b);
                            a0 += F0[i] * wa;
                            b0 += F0[i] * wb;
                            a1 += F1[i] * wa;
                            b1 += F1[i] * wb;
                            (a0, b0, a1, b1)
                        },
                    )
                    .reduce(
                        || (R::ZERO, R::ZERO, R::ZERO, R::ZERO),
                        |(a0, b0, a1, b1), (c0, d0, c1, d1)| (a0 + c0, b0 + d0, a1 + c1, b1 + d1),
                    );
                #[cfg(not(feature = "parallel"))]
                let (f0_a, f0_b, f1_a, f1_b) = {
                    let mut a0 = R::ZERO;
                    let mut b0 = R::ZERO;
                    let mut a1 = R::ZERO;
                    let mut b1 = R::ZERO;
                    for i in 0..n {
                        let wa = eval_at(i, low_a, scale_a);
                        let wb = eval_at(i, low_b, scale_b);
                        a0 += F0[i] * wa;
                        b0 += F0[i] * wb;
                        a1 += F1[i] * wa;
                        b1 += F1[i] * wb;
                    }
                    (a0, b0, a1, b1)
                };

                ((f0_a, f0_b), (f1_a, f1_b))
            } else {
                let ea = eq_a_ring.as_ref().unwrap();
                let eb = eq_b_ring.as_ref().unwrap();
                (
                    (dot_with_eq::<R>(&F0, ea), dot_with_eq::<R>(&F0, eb)),
                    (dot_with_eq::<R>(&F1, ea), dot_with_eq::<R>(&F1, eb)),
                )
            };
            if profile {
                println!(
                    "[LF+ Decomp::decompose_seeded] fv(dot_with_eq) both: {:?}",
                    t_fv.elapsed()
                );
            }

            let t_mats = Instant::now();
            let mut v0 = Vec::with_capacity(1 + self.M.len());
            let mut v1 = Vec::with_capacity(1 + self.M.len());
            v0.push(fv0);
            v1.push(fv1);
            for M_i in self.M.iter().map(|m| m.as_ref()) {
                if is_identity_matrix::<R>(M_i) {
                    v0.push(fv0);
                    v1.push(fv1);
                } else {
                    let (m0, m1) = if let (Some((t_low, low_mask, low_a, scale_a)), Some((_, _, low_b, scale_b))) =
                        (eq_plan_a.as_ref(), eq_plan_b.as_ref())
                    {
                        let eval_at = |idx: usize, low: &[R::BaseRing], scale: &[R::BaseRing]| -> R::BaseRing {
                            eq_at_base::<R::BaseRing>(idx, low, scale, *low_mask, *t_low)
                        };
                        #[cfg(feature = "parallel")]
                        let (out0, out1) = M_i
                            .coeffs
                            .par_iter()
                            .enumerate()
                            .fold(
                                || ((R::ZERO, R::ZERO), (R::ZERO, R::ZERO)),
                                |(mut o0, mut o1), (row, terms)| {
                                    let wa0 = eval_at(row, low_a, scale_a);
                                    let wb0 = eval_at(row, low_b, scale_b);
                                    if wa0 == R::BaseRing::ZERO && wb0 == R::BaseRing::ZERO {
                                        return (o0, o1);
                                    }
                                    let mut s0 = R::ZERO;
                                    let mut s1 = R::ZERO;
                                    for (c, j) in terms {
                                        s0 += *c * F0[*j];
                                        s1 += *c * F1[*j];
                                    }
                                    o0.0 += s0 * wa0;
                                    o0.1 += s0 * wb0;
                                    o1.0 += s1 * wa0;
                                    o1.1 += s1 * wb0;
                                    (o0, o1)
                                },
                            )
                            .reduce(
                                || ((R::ZERO, R::ZERO), (R::ZERO, R::ZERO)),
                                |(a0, a1), (b0, b1)| {
                                    ((a0.0 + b0.0, a0.1 + b0.1), (a1.0 + b1.0, a1.1 + b1.1))
                                },
                            );
                        #[cfg(not(feature = "parallel"))]
                        let (out0, out1) = {
                            let mut out0 = (R::ZERO, R::ZERO);
                            let mut out1 = (R::ZERO, R::ZERO);
                            for (row, terms) in M_i.coeffs.iter().enumerate() {
                                let wa0 = eval_at(row, low_a, scale_a);
                                let wb0 = eval_at(row, low_b, scale_b);
                                if wa0 == R::BaseRing::ZERO && wb0 == R::BaseRing::ZERO {
                                    continue;
                                }
                                let mut s0 = R::ZERO;
                                let mut s1 = R::ZERO;
                                for (c, j) in terms {
                                    s0 += *c * F0[*j];
                                    s1 += *c * F1[*j];
                                }
                                out0.0 += s0 * wa0;
                                out0.1 += s0 * wb0;
                                out1.0 += s1 * wa0;
                                out1.1 += s1 * wb0;
                            }
                            (out0, out1)
                        };
                        (out0, out1)
                    } else {
                        let ea = eq_a_ring.as_ref().unwrap();
                        let eb = eq_b_ring.as_ref().unwrap();
                        eval_sparse_mat_two_vecs_at_two_points::<R>(M_i, &F0, &F1, ea, eb)
                    };
                    v0.push(m0);
                    v1.push(m1);
                }
            }
            if profile {
                println!(
                    "[LF+ Decomp::decompose_seeded] mats(eval_sparse_mat_two_vecs_at_two_points): {:?} (Mlen={})",
                    t_mats.elapsed(),
                    self.M.len()
                );
            }
            maybe_print_rss("decomp_seeded: after v0/v1 mats");
            (v0, v1)
        };

        if profile {
            println!(
                "[LF+ Decomp::decompose_seeded] setup+split: {:?} (nvars={}, Mlen={})",
                t_total.elapsed(),
                nvars,
                self.M.len()
            );
        }

        let t = Instant::now();
        let (v0, v1) = vi_calc_pair();
        if profile {
            println!("[LF+ Decomp::decompose_seeded] compute v0/v1: {:?}", t.elapsed());
        }
        maybe_print_rss("decomp_seeded: after compute v0/v1");

        let t = Instant::now();
        let (C0, C1) = {
            #[cfg(feature = "parallel")]
            {
                rayon::join(
                    || scheme.commit(&F0).unwrap().as_ref().to_vec(),
                    || scheme.commit(&F1).unwrap().as_ref().to_vec(),
                )
            }
            #[cfg(not(feature = "parallel"))]
            {
                (
                    scheme.commit(&F0).unwrap().as_ref().to_vec(),
                    scheme.commit(&F1).unwrap().as_ref().to_vec(),
                )
            }
        };
        if profile {
            println!("[LF+ Decomp::decompose_seeded] commitments C0/C1: {:?}", t.elapsed());
            println!("[LF+ Decomp::decompose_seeded] total: {:?}", t_total.elapsed());
        }
        maybe_print_rss("decomp_seeded: done");

        let linb0 = LinB {
            x: LinBX {
                cm_f: C0.clone(),
                r: self.r.clone(),
                v: v0.clone(),
            },
            f: WitnessVec::Ring(Arc::new(F0)),
        };
        let linb1 = LinB {
            x: LinBX {
                cm_f: C1.clone(),
                r: self.r.clone(),
                v: v1.clone(),
            },
            f: WitnessVec::Ring(Arc::new(F1)),
        };
        let proof = DecompProof { C: (C0, C1), v: (v0, v1) };

        ((linb0, linb1), proof)
    }

}

impl<R: PolyRing> DecompBase<'_, R>
where
    R: Decompose + OverField,
    R::BaseRing: Zq,
{
    /// One-shot decomposition proof for the SP1/base regime.
    ///
    /// This is equivalent to `decompose_seeded_base`, but avoids materializing the second digit
    /// vector `F1` as `Vec<R>`. Instead, we store `F1` in a packed coefficient form and decode
    /// only when needed. This reduces peak RSS substantially and does not change the proof or
    /// verifier transcript.
    pub fn decompose_seeded_base_one_shot(
        self,
        scheme: &AjtaiCommitmentScheme<R>,
        B: u128,
    ) -> DecompProof<R>
    where
        R::BaseRing: PrimeField,
        <R::BaseRing as PrimeField>::BigInt: BigInteger,
    {
        let profile = std::env::var("LF_PLUS_PROFILE").ok().as_deref() == Some("1");
        let t_total = Instant::now();
        maybe_print_rss("decomp_seeded(one_shot): start");

        let nvars = log2(scheme.width()) as usize;
        let mut F0 = self.f;
        let n = F0.len();
        let d = R::dimension();

        let use_i32 = B > (i16::MAX as u128);
        let mut F1_packed: PackedDigitVec<R::BaseRing> = if use_i32 {
            PackedDigitVec::I32 {
                coeffs: vec![0i32; n * d],
                d,
                n,
                _br: core::marker::PhantomData,
            }
        } else {
            PackedDigitVec::new_i16(n, d)
        };

        #[cfg(feature = "parallel")]
        {
            const CHUNK: usize = 1 << 14;
            match &mut F1_packed {
                PackedDigitVec::I16 { coeffs, d: dd, .. } => {
                    debug_assert_eq!(*dd, d);
                    F0.par_chunks_mut(CHUNK)
                        .zip(coeffs.par_chunks_mut(CHUNK * d))
                        .for_each_init(|| vec![R::ZERO; 2], |tmp, (c0, c1_flat)| {
                            for i in 0..c0.len() {
                                let orig = std::mem::replace(&mut c0[i], R::ZERO);
                                orig.decompose_to(B, tmp);
                                c0[i] = tmp[0];
                                let c1c = tmp[1].coeffs();
                                let base = i * d;
                                for k in 0..d {
                                    let v = {
                                        // Interpret coefficient as signed balanced integer.
                                        let rep = c1c[k].into_bigint();
                                        let limbs = rep.as_ref();
                                        let mut vv: u128 = 0;
                                        let take = core::cmp::min(limbs.len(), 2);
                                        for t in 0..take {
                                            vv |= (limbs[t] as u128) << (64 * t);
                                        }
                                        let modulus = <R::BaseRing as PrimeField>::MODULUS;
                                        let mod_limbs = modulus.as_ref();
                                        let mut q: u128 = 0;
                                        let take_q = core::cmp::min(mod_limbs.len(), 2);
                                        for t in 0..take_q {
                                            q |= (mod_limbs[t] as u128) << (64 * t);
                                        }
                                        let half = q >> 1;
                                        if vv <= half {
                                            vv as i64
                                        } else {
                                            (vv as i128 - q as i128) as i64
                                        }
                                    };
                                    debug_assert!(
                                        v >= i16::MIN as i64 && v <= i16::MAX as i64,
                                        "packed F1 coefficient out of i16 range"
                                    );
                                    c1_flat[base + k] = v as i16;
                                }
                            }
                        });
                }
                PackedDigitVec::I32 { .. } => {
                    // Rare case (very large B): fall back to the existing materialized implementation.
                    // This preserves correctness/transcript but doesn't reduce memory.
                    drop(F1_packed);
                    let ((_linb0, _linb1), proof) = DecompBase { f: F0, r: self.r, M0: self.M0 }
                        .decompose_seeded_base(scheme, B);
                    return proof;
                }
            }
        }
        #[cfg(not(feature = "parallel"))]
        {
            let mut tmp = vec![R::ZERO; 2];
            for i in 0..n {
                let orig = std::mem::replace(&mut F0[i], R::ZERO);
                orig.decompose_to(B, &mut tmp);
                F0[i] = tmp[0];
                F1_packed.set_coeffs_from_ring(i, &tmp[1]);
            }
        }
        maybe_print_rss("decomp_seeded(one_shot): after decompose_to_packed");

        let r_a = self.r.iter().map(|rr| rr.0).collect::<Vec<_>>();
        let r_b = self.r.iter().map(|rr| rr.1).collect::<Vec<_>>();

        #[inline]
        fn is_identity_matrix_base<BR: Ring>(m: &SparseMatrix<BR>) -> bool {
            if m.nrows != m.ncols {
                return false;
            }
            if m.coeffs.len() != m.nrows {
                return false;
            }
            for (i, row) in m.coeffs.iter().enumerate() {
                if row.len() != 1 {
                    return false;
                }
                let (c, j) = row[0];
                if j != i {
                    return false;
                }
                if c != BR::ONE {
                    return false;
                }
            }
            true
        }

        #[inline]
        fn eq_weights<Rr: PolyRing>(r: &[Rr]) -> Vec<Rr> {
            let nvars = r.len();
            let n = 1usize << nvars;
            let mut cur = vec![Rr::ZERO; n];
            let mut next = vec![Rr::ZERO; n];
            cur[0] = Rr::ONE;
            let mut len = 1usize;
            let mut cur_is_cur = true;
            for &rj in r.iter().rev() {
                let om = Rr::ONE - rj;
                let (src, dst) = if cur_is_cur {
                    (&cur[..len], &mut next[..(2 * len)])
                } else {
                    (&next[..len], &mut cur[..(2 * len)])
                };
                for (i, &wi) in src.iter().enumerate() {
                    dst[2 * i] = wi * om;
                    dst[2 * i + 1] = wi * rj;
                }
                len <<= 1;
                cur_is_cur = !cur_is_cur;
            }
            if cur_is_cur { cur } else { next }
        }

        #[inline]
        fn is_const_coeff_ring<Rr: PolyRing>(x: &Rr) -> bool {
            x.coeffs()
                .iter()
                .skip(1)
                .all(|c| *c == <Rr as PolyRing>::BaseRing::ZERO)
        }

        #[inline]
        fn dot_with_eq_ring<Rr: PolyRing>(f: &[Rr], eq: &[Rr]) -> Rr {
            debug_assert_eq!(f.len(), eq.len());
            #[cfg(feature = "parallel")]
            {
                f.par_iter()
                    .zip(eq.par_iter())
                    .map(|(&fx, &wx)| fx * wx)
                    .reduce(|| Rr::ZERO, |a, b| a + b)
            }
            #[cfg(not(feature = "parallel"))]
            {
                f.iter()
                    .zip(eq.iter())
                    .fold(Rr::ZERO, |acc, (&fx, &wx)| acc + fx * wx)
            }
        }

        #[inline]
        fn dot_with_eq_packed<Rr: PolyRing>(
            f0: &[Rr],
            f1_packed: &PackedDigitVec<Rr::BaseRing>,
            eq: &[Rr],
        ) -> Rr
        where
            Rr::BaseRing: Ring + Copy,
        {
            debug_assert_eq!(f0.len(), eq.len());
            let mut acc = Rr::ZERO;
            let mut tmp = Rr::ZERO;
            for i in 0..f0.len() {
                // f0 part not used here, but keep signature similar if we need both later.
                let _ = &f0[i];
                f1_packed.fill_ring_at(i, &mut tmp);
                acc += tmp * eq[i];
            }
            acc
        }

        let vi_calc_pair = || {
            let t_eq = Instant::now();
            let (eq_a_ring, eq_b_ring, eq_plan_a, eq_plan_b) = if r_a
                .iter()
                .all(is_const_coeff_ring::<R>)
                && r_b.iter().all(is_const_coeff_ring::<R>)
            {
                let r_a0 = r_a.iter().map(|x| x.coeffs()[0]).collect::<Vec<_>>();
                let r_b0 = r_b.iter().map(|x| x.coeffs()[0]).collect::<Vec<_>>();

                let one_minus_a0 =
                    r_a0.iter().copied().map(|x| R::BaseRing::ONE - x).collect::<Vec<_>>();
                let one_minus_b0 =
                    r_b0.iter().copied().map(|x| R::BaseRing::ONE - x).collect::<Vec<_>>();
                let t_low = choose_t_low(r_a0.len());

                let low_a =
                    build_eq_low_table_base::<R::BaseRing>(&r_a0[..t_low], &one_minus_a0[..t_low]);
                let low_b =
                    build_eq_low_table_base::<R::BaseRing>(&r_b0[..t_low], &one_minus_b0[..t_low]);
                let low_mask = (1usize << t_low) - 1;
                let scale_a = build_scale_high_base::<R::BaseRing>(&r_a0, &one_minus_a0, t_low);
                let scale_b = build_scale_high_base::<R::BaseRing>(&r_b0, &one_minus_b0, t_low);

                (
                    None,
                    None,
                    Some((t_low, low_mask, low_a, scale_a)),
                    Some((t_low, low_mask, low_b, scale_b)),
                )
            } else {
                (Some(eq_weights::<R>(&r_a)), Some(eq_weights::<R>(&r_b)), None, None)
            };
            maybe_print_rss("decomp_seeded(one_shot): after eq_weights");
            if profile {
                println!(
                    "[LF+ Decomp::decompose_seeded_base_one_shot] eq_weights: {:?} (nvars={})",
                    t_eq.elapsed(),
                    nvars
                );
            }

            let t_fv = Instant::now();
            let (fv0, fv1) = if let (
                Some((t_low, low_mask, low_a, scale_a)),
                Some((_, _, low_b, scale_b)),
            ) = (eq_plan_a.as_ref(), eq_plan_b.as_ref())
            {
                let n = F0.len();
                let eval_at = |idx: usize, low: &[R::BaseRing], scale: &[R::BaseRing]| -> R::BaseRing {
                    eq_at_base::<R::BaseRing>(idx, low, scale, *low_mask, *t_low)
                };
                #[cfg(feature = "parallel")]
                let (f0_a, f0_b, f1_a, f1_b) = (0..n)
                    .into_par_iter()
                    .fold(
                        || (R::ZERO, R::ZERO, R::ZERO, R::ZERO),
                        |(mut a0, mut b0, mut a1, mut b1), i| {
                            let wa = eval_at(i, low_a, scale_a);
                            let wb = eval_at(i, low_b, scale_b);
                            a0 += F0[i] * wa;
                            b0 += F0[i] * wb;
                            F1_packed.mul_add_into_ring(&mut a1, i, wa);
                            F1_packed.mul_add_into_ring(&mut b1, i, wb);
                            (a0, b0, a1, b1)
                        },
                    )
                    .reduce(
                        || (R::ZERO, R::ZERO, R::ZERO, R::ZERO),
                        |(a0, b0, a1, b1), (c0, d0, c1, d1)| (a0 + c0, b0 + d0, a1 + c1, b1 + d1),
                    );
                #[cfg(not(feature = "parallel"))]
                let (f0_a, f0_b, f1_a, f1_b) = {
                    let mut a0 = R::ZERO;
                    let mut b0 = R::ZERO;
                    let mut a1 = R::ZERO;
                    let mut b1 = R::ZERO;
                    for i in 0..n {
                        let wa = eval_at(i, low_a, scale_a);
                        let wb = eval_at(i, low_b, scale_b);
                        a0 += F0[i] * wa;
                        b0 += F0[i] * wb;
                        F1_packed.mul_add_into_ring(&mut a1, i, wa);
                        F1_packed.mul_add_into_ring(&mut b1, i, wb);
                    }
                    (a0, b0, a1, b1)
                };
                ((f0_a, f0_b), (f1_a, f1_b))
            } else {
                let ea = eq_a_ring.as_ref().unwrap();
                let eb = eq_b_ring.as_ref().unwrap();
                (
                    (dot_with_eq_ring(&F0, ea), dot_with_eq_ring(&F0, eb)),
                    (dot_with_eq_packed(&F0, &F1_packed, ea), dot_with_eq_packed(&F0, &F1_packed, eb)),
                )
            };
            if profile {
                println!(
                    "[LF+ Decomp::decompose_seeded_base_one_shot] fv both: {:?}",
                    t_fv.elapsed()
                );
            }

            let t_mats = Instant::now();
            let mut v0 = Vec::with_capacity(1 + self.M0.len());
            let mut v1 = Vec::with_capacity(1 + self.M0.len());
            v0.push(fv0);
            v1.push(fv1);

            for M_i in self.M0.iter().map(|x| x.as_ref()) {
                if is_identity_matrix_base::<R::BaseRing>(M_i) {
                    v0.push(fv0);
                    v1.push(fv1);
                    continue;
                }
                let (m0, m1) = if let (
                    Some((t_low, low_mask, low_a, scale_a)),
                    Some((_, _, low_b, scale_b)),
                ) = (eq_plan_a.as_ref(), eq_plan_b.as_ref())
                {
                    let eval_at = |idx: usize, low: &[R::BaseRing], scale: &[R::BaseRing]| -> R::BaseRing {
                        eq_at_base::<R::BaseRing>(idx, low, scale, *low_mask, *t_low)
                    };
                    #[cfg(feature = "parallel")]
                    let (out0, out1) = M_i
                        .coeffs
                        .par_iter()
                        .enumerate()
                        .fold(
                            || ((R::ZERO, R::ZERO), (R::ZERO, R::ZERO)),
                            |(mut o0, mut o1), (row, terms)| {
                                let wa0 = eval_at(row, low_a, scale_a);
                                let wb0 = eval_at(row, low_b, scale_b);
                                if wa0 == R::BaseRing::ZERO && wb0 == R::BaseRing::ZERO {
                                    return (o0, o1);
                                }
                                let mut s0 = R::ZERO;
                                let mut s1 = R::ZERO;
                                for (c0, j) in terms {
                                    s0 += F0[*j] * *c0;
                                    F1_packed.mul_add_into_ring(&mut s1, *j, *c0);
                                }
                                o0.0 += s0 * wa0;
                                o0.1 += s0 * wb0;
                                o1.0 += s1 * wa0;
                                o1.1 += s1 * wb0;
                                (o0, o1)
                            },
                        )
                        .reduce(
                            || ((R::ZERO, R::ZERO), (R::ZERO, R::ZERO)),
                            |(a0, a1), (b0, b1)| {
                                ((a0.0 + b0.0, a0.1 + b0.1), (a1.0 + b1.0, a1.1 + b1.1))
                            },
                        );
                    #[cfg(not(feature = "parallel"))]
                    let (out0, out1) = {
                        let mut out0 = (R::ZERO, R::ZERO);
                        let mut out1 = (R::ZERO, R::ZERO);
                        for (row, terms) in M_i.coeffs.iter().enumerate() {
                            let wa0 = eval_at(row, low_a, scale_a);
                            let wb0 = eval_at(row, low_b, scale_b);
                            if wa0 == R::BaseRing::ZERO && wb0 == R::BaseRing::ZERO {
                                continue;
                            }
                            let mut s0 = R::ZERO;
                            let mut s1 = R::ZERO;
                            for (c0, j) in terms {
                                s0 += F0[*j] * *c0;
                                F1_packed.mul_add_into_ring(&mut s1, *j, *c0);
                            }
                            out0.0 += s0 * wa0;
                            out0.1 += s0 * wb0;
                            out1.0 += s1 * wa0;
                            out1.1 += s1 * wb0;
                        }
                        (out0, out1)
                    };
                    (out0, out1)
                } else {
                    let ea = eq_a_ring.as_ref().unwrap();
                    let eb = eq_b_ring.as_ref().unwrap();
                    // Same as `eval_sparse_mat0_two_vecs_at_two_points`, but with `F1` packed.
                    let mut out0 = (R::ZERO, R::ZERO);
                    let mut out1 = (R::ZERO, R::ZERO);
                    for (row, terms) in M_i.coeffs.iter().enumerate() {
                        let wa = ea[row];
                        let wb = eb[row];
                        if wa == R::ZERO && wb == R::ZERO {
                            continue;
                        }
                        let mut s0 = R::ZERO;
                        let mut s1 = R::ZERO;
                        for (c0, j) in terms {
                            s0 += F0[*j] * *c0;
                            F1_packed.mul_add_into_ring(&mut s1, *j, *c0);
                        }
                        out0.0 += wa * s0;
                        out0.1 += wb * s0;
                        out1.0 += wa * s1;
                        out1.1 += wb * s1;
                    }
                    (out0, out1)
                };
                v0.push(m0);
                v1.push(m1);
            }
            if profile {
                println!(
                    "[LF+ Decomp::decompose_seeded_base_one_shot] mats: {:?} (Mlen={})",
                    t_mats.elapsed(),
                    self.M0.len()
                );
            }
            maybe_print_rss("decomp_seeded(one_shot): after v0/v1 mats");
            (v0, v1)
        };

        if profile {
            println!(
                "[LF+ Decomp::decompose_seeded_base_one_shot] setup+split: {:?} (nvars={}, Mlen={})",
                t_total.elapsed(),
                nvars,
                self.M0.len()
            );
        }

        let t = Instant::now();
        let (v0, v1) = vi_calc_pair();
        if profile {
            println!(
                "[LF+ Decomp::decompose_seeded_base_one_shot] compute v0/v1: {:?}",
                t.elapsed()
            );
        }
        maybe_print_rss("decomp_seeded(one_shot): after compute v0/v1");

        let t = Instant::now();
        let C0 = scheme.commit(&F0).unwrap().as_ref().to_vec();
        let C1 = {
            let mut out = scheme
                .commit_many_with(n, 1, |j, vals| {
                    let mut tmp = R::ZERO;
                    F1_packed.fill_ring_at(j, &mut tmp);
                    vals[0] = tmp;
                })
                .unwrap();
            out.remove(0).as_ref().to_vec()
        };
        if profile {
            println!(
                "[LF+ Decomp::decompose_seeded_base_one_shot] commitments C0/C1: {:?}",
                t.elapsed()
            );
            println!(
                "[LF+ Decomp::decompose_seeded_base_one_shot] total: {:?}",
                t_total.elapsed()
            );
        }
        maybe_print_rss("decomp_seeded(one_shot): done");

        DecompProof { C: (C0, C1), v: (v0, v1) }
    }

    /// Same as [`Decomp::decompose_seeded`], but with external matrices represented over the base ring.
    pub fn decompose_seeded_base(
        self,
        scheme: &AjtaiCommitmentScheme<R>,
        B: u128,
    ) -> ((LinB<R>, LinB<R>), DecompProof<R>) {
        let profile = std::env::var("LF_PLUS_PROFILE").ok().as_deref() == Some("1");
        let t_total = Instant::now();
        maybe_print_rss("decomp_seeded: start");

        let nvars = log2(scheme.width()) as usize;
        let mut F0 = self.f;
        let n = F0.len();
        let mut F1 = vec![R::ZERO; n];
        #[cfg(feature = "parallel")]
        {
            const CHUNK: usize = 1 << 14;
            F0.par_chunks_mut(CHUNK)
                .zip(F1.par_chunks_mut(CHUNK))
                .for_each_init(|| vec![R::ZERO; 2], |tmp, (c0, c1)| {
                    for i in 0..c0.len() {
                        let orig = std::mem::replace(&mut c0[i], R::ZERO);
                        orig.decompose_to(B, tmp);
                        c0[i] = tmp[0];
                        c1[i] = tmp[1];
                    }
                });
        }
        #[cfg(not(feature = "parallel"))]
        {
            let mut tmp = vec![R::ZERO; 2];
            for i in 0..n {
                let orig = std::mem::replace(&mut F0[i], R::ZERO);
                orig.decompose_to(B, &mut tmp);
                F0[i] = tmp[0];
                F1[i] = tmp[1];
            }
        }
        maybe_print_rss("decomp_seeded: after decompose_to_vec");

        let r_a = self.r.iter().map(|rr| rr.0).collect::<Vec<_>>();
        let r_b = self.r.iter().map(|rr| rr.1).collect::<Vec<_>>();

        #[inline]
        fn is_identity_matrix_base<BR: Ring>(m: &SparseMatrix<BR>) -> bool {
            if m.nrows != m.ncols {
                return false;
            }
            if m.coeffs.len() != m.nrows {
                return false;
            }
            for (i, row) in m.coeffs.iter().enumerate() {
                if row.len() != 1 {
                    return false;
                }
                let (c, j) = row[0];
                if j != i {
                    return false;
                }
                if c != BR::ONE {
                    return false;
                }
            }
            true
        }

        #[inline]
        fn eq_weights<Rr: PolyRing>(r: &[Rr]) -> Vec<Rr> {
            let nvars = r.len();
            let n = 1usize << nvars;
            let mut cur = vec![Rr::ZERO; n];
            let mut next = vec![Rr::ZERO; n];
            cur[0] = Rr::ONE;

            let mut len = 1usize;
            let mut cur_is_cur = true;
            for &rj in r.iter().rev() {
                let om = Rr::ONE - rj;
                let (src, dst) = if cur_is_cur {
                    (&cur[..len], &mut next[..(2 * len)])
                } else {
                    (&next[..len], &mut cur[..(2 * len)])
                };

                #[cfg(feature = "parallel")]
                {
                    dst.par_chunks_mut(2)
                        .zip(src.par_iter())
                        .for_each(|(pair, &wi)| {
                            pair[0] = wi * om;
                            pair[1] = wi * rj;
                        });
                }
                #[cfg(not(feature = "parallel"))]
                {
                    for (i, &wi) in src.iter().enumerate() {
                        dst[2 * i] = wi * om;
                        dst[2 * i + 1] = wi * rj;
                    }
                }

                len <<= 1;
                cur_is_cur = !cur_is_cur;
            }
            if cur_is_cur { cur } else { next }
        }

        #[inline]
        fn is_const_coeff_ring<Rr: PolyRing>(x: &Rr) -> bool {
            x.coeffs()
                .iter()
                .skip(1)
                .all(|c| *c == <Rr as PolyRing>::BaseRing::ZERO)
        }

        #[inline]
        fn eval_sparse_mat0_two_vecs_at_two_points<Rr: PolyRing>(
            m0: &SparseMatrix<Rr::BaseRing>,
            f0: &[Rr],
            f1: &[Rr],
            eq_a: &[Rr],
            eq_b: &[Rr],
        ) -> (RxR<Rr>, RxR<Rr>)
        where
            Rr::BaseRing: Ring,
            Rr: core::ops::Mul<Rr::BaseRing, Output = Rr>,
        {
            debug_assert_eq!(m0.ncols, f0.len());
            debug_assert_eq!(m0.ncols, f1.len());
            debug_assert_eq!(m0.nrows, eq_a.len());
            debug_assert_eq!(m0.nrows, eq_b.len());

            let mut out0 = (Rr::ZERO, Rr::ZERO);
            let mut out1 = (Rr::ZERO, Rr::ZERO);
            for (row, terms) in m0.coeffs.iter().enumerate() {
                let wa = eq_a[row];
                let wb = eq_b[row];
                if wa == Rr::ZERO && wb == Rr::ZERO {
                    continue;
                }
                let mut s0 = Rr::ZERO;
                let mut s1 = Rr::ZERO;
                for (c0, j) in terms {
                    s0 += f0[*j] * *c0;
                    s1 += f1[*j] * *c0;
                }
                out0.0 += wa * s0;
                out0.1 += wb * s0;
                out1.0 += wa * s1;
                out1.1 += wb * s1;
            }
            (out0, out1)
        }

        #[inline]
        fn dot_with_eq<Rr: PolyRing>(f: &[Rr], eq: &[Rr]) -> Rr {
            debug_assert_eq!(f.len(), eq.len());
            #[cfg(feature = "parallel")]
            {
                f.par_iter()
                    .zip(eq.par_iter())
                    .map(|(&fx, &wx)| fx * wx)
                    .reduce(|| Rr::ZERO, |a, b| a + b)
            }
            #[cfg(not(feature = "parallel"))]
            {
                f.iter()
                    .zip(eq.iter())
                    .fold(Rr::ZERO, |acc, (&fx, &wx)| acc + fx * wx)
            }
        }

        let vi_calc_pair = || {
            let t_eq = Instant::now();
            let (eq_a_ring, eq_b_ring, eq_plan_a, eq_plan_b) = if r_a
                .iter()
                .all(is_const_coeff_ring::<R>)
                && r_b.iter().all(is_const_coeff_ring::<R>)
            {
                let r_a0 = r_a.iter().map(|x| x.coeffs()[0]).collect::<Vec<_>>();
                let r_b0 = r_b.iter().map(|x| x.coeffs()[0]).collect::<Vec<_>>();

                let one_minus_a0 =
                    r_a0.iter().copied().map(|x| R::BaseRing::ONE - x).collect::<Vec<_>>();
                let one_minus_b0 =
                    r_b0.iter().copied().map(|x| R::BaseRing::ONE - x).collect::<Vec<_>>();
                let t_low = choose_t_low(r_a0.len());

                let low_a =
                    build_eq_low_table_base::<R::BaseRing>(&r_a0[..t_low], &one_minus_a0[..t_low]);
                let low_b =
                    build_eq_low_table_base::<R::BaseRing>(&r_b0[..t_low], &one_minus_b0[..t_low]);
                let low_mask = (1usize << t_low) - 1;
                let scale_a = build_scale_high_base::<R::BaseRing>(&r_a0, &one_minus_a0, t_low);
                let scale_b = build_scale_high_base::<R::BaseRing>(&r_b0, &one_minus_b0, t_low);

                (None, None, Some((t_low, low_mask, low_a, scale_a)), Some((t_low, low_mask, low_b, scale_b)))
            } else {
                (Some(eq_weights::<R>(&r_a)), Some(eq_weights::<R>(&r_b)), None, None)
            };
            maybe_print_rss("decomp_seeded: after eq_weights");
            if profile {
                println!(
                    "[LF+ Decomp::decompose_seeded_base] eq_weights: {:?} (nvars={})",
                    t_eq.elapsed(),
                    nvars
                );
            }

            let t_fv = Instant::now();
            let (fv0, fv1) = if let (Some((t_low, low_mask, low_a, scale_a)), Some((_, _, low_b, scale_b))) =
                (eq_plan_a.as_ref(), eq_plan_b.as_ref())
            {
                let n = F0.len();
                let eval_at = |idx: usize, low: &[R::BaseRing], scale: &[R::BaseRing]| -> R::BaseRing {
                    eq_at_base::<R::BaseRing>(idx, low, scale, *low_mask, *t_low)
                };
                #[cfg(feature = "parallel")]
                let (f0_a, f0_b, f1_a, f1_b) = (0..n)
                    .into_par_iter()
                    .fold(
                        || (R::ZERO, R::ZERO, R::ZERO, R::ZERO),
                        |(mut a0, mut b0, mut a1, mut b1), i| {
                            let wa = eval_at(i, low_a, scale_a);
                            let wb = eval_at(i, low_b, scale_b);
                            a0 += F0[i] * wa;
                            b0 += F0[i] * wb;
                            a1 += F1[i] * wa;
                            b1 += F1[i] * wb;
                            (a0, b0, a1, b1)
                        },
                    )
                    .reduce(
                        || (R::ZERO, R::ZERO, R::ZERO, R::ZERO),
                        |(a0, b0, a1, b1), (c0, d0, c1, d1)| (a0 + c0, b0 + d0, a1 + c1, b1 + d1),
                    );
                #[cfg(not(feature = "parallel"))]
                let (f0_a, f0_b, f1_a, f1_b) = {
                    let mut a0 = R::ZERO;
                    let mut b0 = R::ZERO;
                    let mut a1 = R::ZERO;
                    let mut b1 = R::ZERO;
                    for i in 0..n {
                        let wa = eval_at(i, low_a, scale_a);
                        let wb = eval_at(i, low_b, scale_b);
                        a0 += F0[i] * wa;
                        b0 += F0[i] * wb;
                        a1 += F1[i] * wa;
                        b1 += F1[i] * wb;
                    }
                    (a0, b0, a1, b1)
                };
                ((f0_a, f0_b), (f1_a, f1_b))
            } else {
                let ea = eq_a_ring.as_ref().unwrap();
                let eb = eq_b_ring.as_ref().unwrap();
                (
                    (dot_with_eq(&F0, ea), dot_with_eq(&F0, eb)),
                    (dot_with_eq(&F1, ea), dot_with_eq(&F1, eb)),
                )
            };
            if profile {
                println!(
                    "[LF+ Decomp::decompose_seeded_base] fv(dot_with_eq) both: {:?}",
                    t_fv.elapsed()
                );
            }

            let t_mats = Instant::now();
            let mut v0 = Vec::with_capacity(1 + self.M0.len());
            let mut v1 = Vec::with_capacity(1 + self.M0.len());
            v0.push(fv0);
            v1.push(fv1);

            for M_i in self.M0.iter().map(|x| x.as_ref()) {
                if is_identity_matrix_base::<R::BaseRing>(M_i) {
                    v0.push(fv0);
                    v1.push(fv1);
                    continue;
                }
                let (m0, m1) = if let (Some((t_low, low_mask, low_a, scale_a)), Some((_, _, low_b, scale_b))) =
                    (eq_plan_a.as_ref(), eq_plan_b.as_ref())
                {
                    let eval_at = |idx: usize, low: &[R::BaseRing], scale: &[R::BaseRing]| -> R::BaseRing {
                        eq_at_base::<R::BaseRing>(idx, low, scale, *low_mask, *t_low)
                    };
                    #[cfg(feature = "parallel")]
                    let (out0, out1) = M_i
                        .coeffs
                        .par_iter()
                        .enumerate()
                        .fold(
                            || ((R::ZERO, R::ZERO), (R::ZERO, R::ZERO)),
                            |(mut o0, mut o1), (row, terms)| {
                                let wa0 = eval_at(row, low_a, scale_a);
                                let wb0 = eval_at(row, low_b, scale_b);
                                if wa0 == R::BaseRing::ZERO && wb0 == R::BaseRing::ZERO {
                                    return (o0, o1);
                                }
                                let mut s0 = R::ZERO;
                                let mut s1 = R::ZERO;
                                for (c0, j) in terms {
                                    s0 += F0[*j] * *c0;
                                    s1 += F1[*j] * *c0;
                                }
                                o0.0 += s0 * wa0;
                                o0.1 += s0 * wb0;
                                o1.0 += s1 * wa0;
                                o1.1 += s1 * wb0;
                                (o0, o1)
                            },
                        )
                        .reduce(
                            || ((R::ZERO, R::ZERO), (R::ZERO, R::ZERO)),
                            |(a0, a1), (b0, b1)| {
                                ((a0.0 + b0.0, a0.1 + b0.1), (a1.0 + b1.0, a1.1 + b1.1))
                            },
                        );
                    #[cfg(not(feature = "parallel"))]
                    let (out0, out1) = {
                        let mut out0 = (R::ZERO, R::ZERO);
                        let mut out1 = (R::ZERO, R::ZERO);
                        for (row, terms) in M_i.coeffs.iter().enumerate() {
                            let wa0 = eval_at(row, low_a, scale_a);
                            let wb0 = eval_at(row, low_b, scale_b);
                            if wa0 == R::BaseRing::ZERO && wb0 == R::BaseRing::ZERO {
                                continue;
                            }
                            let mut s0 = R::ZERO;
                            let mut s1 = R::ZERO;
                            for (c0, j) in terms {
                                s0 += F0[*j] * *c0;
                                s1 += F1[*j] * *c0;
                            }
                            out0.0 += s0 * wa0;
                            out0.1 += s0 * wb0;
                            out1.0 += s1 * wa0;
                            out1.1 += s1 * wb0;
                        }
                        (out0, out1)
                    };
                    (out0, out1)
                } else {
                    let ea = eq_a_ring.as_ref().unwrap();
                    let eb = eq_b_ring.as_ref().unwrap();
                    eval_sparse_mat0_two_vecs_at_two_points::<R>(M_i, &F0, &F1, ea, eb)
                };
                v0.push(m0);
                v1.push(m1);
            }
            if profile {
                println!(
                    "[LF+ Decomp::decompose_seeded_base] mats(eval_sparse_mat_two_vecs_at_two_points): {:?} (Mlen={})",
                    t_mats.elapsed(),
                    self.M0.len()
                );
            }
            maybe_print_rss("decomp_seeded: after v0/v1 mats");
            (v0, v1)
        };

        if profile {
            println!(
                "[LF+ Decomp::decompose_seeded_base] setup+split: {:?} (nvars={}, Mlen={})",
                t_total.elapsed(),
                nvars,
                self.M0.len()
            );
        }

        let t = Instant::now();
        let (v0, v1) = vi_calc_pair();
        if profile {
            println!("[LF+ Decomp::decompose_seeded_base] compute v0/v1: {:?}", t.elapsed());
        }
        maybe_print_rss("decomp_seeded: after compute v0/v1");

        let t = Instant::now();
        let (C0, C1) = {
            #[cfg(feature = "parallel")]
            {
                rayon::join(
                    || scheme.commit(&F0).unwrap().as_ref().to_vec(),
                    || scheme.commit(&F1).unwrap().as_ref().to_vec(),
                )
            }
            #[cfg(not(feature = "parallel"))]
            {
                (
                    scheme.commit(&F0).unwrap().as_ref().to_vec(),
                    scheme.commit(&F1).unwrap().as_ref().to_vec(),
                )
            }
        };
        if profile {
            println!("[LF+ Decomp::decompose_seeded_base] commitments C0/C1: {:?}", t.elapsed());
            println!("[LF+ Decomp::decompose_seeded_base] total: {:?}", t_total.elapsed());
        }
        maybe_print_rss("decomp_seeded: done");

        let linb0 = LinB {
            x: LinBX {
                cm_f: C0.clone(),
                r: self.r.clone(),
                v: v0.clone(),
            },
            f: WitnessVec::Ring(Arc::new(F0)),
        };
        let linb1 = LinB {
            x: LinBX {
                cm_f: C1.clone(),
                r: self.r.clone(),
                v: v1.clone(),
            },
            f: WitnessVec::Ring(Arc::new(F1)),
        };
        let proof = DecompProof { C: (C0, C1), v: (v0, v1) };

        ((linb0, linb1), proof)
    }
}

impl<R: PolyRing> DecompProof<R> {
    pub fn verify(&self, cm_f: &[R], v: &[(R, R)], B: u128) {
        let Br = R::from(B);
        let rec_cm = self
            .C
            .0
            .iter()
            .zip(self.C.1.iter())
            .map(|(&r0, &r1)| recompose(&[r0, r1], Br))
            .collect::<Vec<R>>();

        let rec_v = self
            .v
            .0
            .iter()
            .zip(self.v.1.iter())
            .map(|(v0, v1)| (recompose(&[v0.0, v1.0], Br), recompose(&[v0.1, v1.1], Br)))
            .collect::<Vec<(R, R)>>();

        assert_eq!(rec_cm, cm_f);
        assert_eq!(rec_v, v);
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::PrimeField;
    use ark_std::One;
    use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
    use latticefold::arith::r1cs::R1CS;
    use stark_rings::{
        balanced_decomposition::GadgetDecompose, cyclotomic_ring::models::frog_ring::RqPoly as R,
    };
    use stark_rings_linalg::SparseMatrix;

    use super::*;
    use crate::{
        lin::{LinParameters, Linearize, LinearizedVerify},
        mlin::Mlin,
        r1cs::{r1cs_decomposed_square, ComR1CS},
        rgchk::DecompParameters,
        transcript::PoseidonTranscript,
    };

    fn identity_cs(n: usize) -> (R1CS<R>, Vec<R>) {
        let r1cs = R1CS::<R> {
            l: 1,
            A: SparseMatrix::identity(n),
            B: SparseMatrix::identity(n),
            C: SparseMatrix::identity(n),
        };
        let z = vec![R::one(); n];
        (r1cs, z)
    }

    #[test]
    fn test_decomp_r1cs() {
        let B = 50u128;
        let kappa = 2;
        let n = 1 << 15;
        let k = 4;

        let (mut r1cs, z) = identity_cs(n / k);
        r1cs.A.coeffs[0][0].0 = 2u128.into();
        r1cs.C.coeffs[0][0].0 = 2u128.into();
        let r1cs = r1cs_decomposed_square(r1cs, n, 2, k);

        let A = Matrix::<R>::rand(&mut ark_std::test_rng(), kappa, n);

        let cr1cs = ComR1CS::new(r1cs, z, 1, 2, k, &A);

        let M = cr1cs.x.matrices_arc();

        let mut ts = PoseidonTranscript::empty::<PC>();
        let (linb, lproof) = cr1cs.linearize(&mut ts);

        let mut ts = PoseidonTranscript::empty::<PC>();
        lproof.verify(&mut ts);

        let r = lproof.r.iter().map(|&r| (r, r)).collect::<Vec<_>>();

        let decomp = Decomp {
            // Decomp currently expects an owned witness vector.
            // This test is small; cloning is fine here.
            f: cr1cs
                .f
                .as_ring_arc()
                .expect("test uses ring witness")
                .as_ref()
                .clone(),
            r,
            M: &M,
        };

        let ((_linb0, _linb1), proof) = decomp.decompose(&A, B);

        proof.verify(&cr1cs.x.cm_f, &linb.x.v, B);
    }

    #[test]
    fn test_decomp_g() {
        let B = (<<R as PolyRing>::BaseRing>::MODULUS.0[0] as f64)
            .sqrt()
            .ceil() as u128
            + 1;
        let n = 1 << 15;
        let k = 2;
        let kappa = 2;
        let b = (R::dimension() / 2) as u128;
        // log_d' (q)
        let l = ((<<R as PolyRing>::BaseRing>::MODULUS.0[0] as f64).ln()
            / ((R::dimension() / 2) as f64).ln())
        .ceil() as usize;

        let params = LinParameters {
            kappa,
            decomp: DecompParameters { b, k, l },
        };

        let z0 = vec![R::one(); n / k];
        let mut z1 = vec![R::one(); n / k];
        z1[0] = R::from(0u128);

        let mut r1cs = R1CS::<R> {
            l: 1,
            A: SparseMatrix::identity(n / k),
            B: SparseMatrix::identity(n / k),
            C: SparseMatrix::identity(n / k),
        };

        r1cs.A.coeffs[0][0].0 = 2u128.into();
        r1cs.C.coeffs[0][0].0 = 2u128.into();

        r1cs.A = r1cs.A.gadget_decompose(2, k);
        r1cs.B = r1cs.B.gadget_decompose(2, k);
        r1cs.C = r1cs.C.gadget_decompose(2, k);
        r1cs.A.pad_rows(n);
        r1cs.B.pad_rows(n);
        r1cs.C.pad_rows(n);

        let f0 = z0.gadget_decompose(2, k);
        let f1 = z1.gadget_decompose(2, k);
        r1cs.check_relation(&f0).unwrap();
        r1cs.check_relation(&f1).unwrap();

        let A = Matrix::<R>::rand(&mut ark_std::test_rng(), params.kappa, n);

        let cr1cs0 = ComR1CS::new(r1cs.clone(), z0, 1, B, k, &A);
        let cr1cs1 = ComR1CS::new(r1cs, z1, 1, B, k, &A);

        let mut ts = PoseidonTranscript::empty::<PC>();
        let (linb0, lproof0) = cr1cs0.linearize(&mut ts);
        let (linb1, lproof1) = cr1cs1.linearize(&mut ts);

        let M = cr1cs0.x.matrices_arc();

        let mlin = Mlin {
            lins: vec![linb0, linb1],
            params,
        };

        let (linb2, cmproof) = mlin.mlin(&A, &M, &mut ts);

        let mut ts = PoseidonTranscript::empty::<PC>();
        lproof0.verify(&mut ts);
        lproof1.verify(&mut ts);
        cmproof.verify(&M, &mut ts).unwrap();

        let decomp = Decomp {
            f: linb2.g,
            r: linb2.x.ro,
            M: &M,
        };

        let (_linb, proof) = decomp.decompose(&A, B);

        proof.verify(&linb2.x.cm_g, &linb2.x.vo, B);
    }

    #[test]
    fn test_decomp_base_one_shot_matches_materialized() {
        use ark_std::UniformRand;
        let mut rng = ark_std::test_rng();
        let kappa = 2;
        let nvars = 12usize;
        let n = 1usize << nvars;
        let B = (R::dimension() / 2) as u128;

        // Explicit Ajtai matrix (small test).
        let A = Matrix::<R>::rand(&mut rng, kappa, n);
        let scheme = AjtaiCommitmentScheme::new(A);

        // Base-ring matrices (identity is enough to exercise the code).
        let M0: Vec<Arc<SparseMatrix<<R as PolyRing>::BaseRing>>> =
            vec![Arc::new(SparseMatrix::<<R as PolyRing>::BaseRing>::identity(n))];

        // Random witness and transcript-derived r (const-coeff path).
        let f = (0..n).map(|_| R::rand(&mut rng)).collect::<Vec<_>>();
        let r = (0..nvars)
            .map(|i| {
                let x = R::from((i as u128) + 7);
                (x, x)
            })
            .collect::<Vec<_>>();

        let decomp_old = DecompBase {
            f: f.clone(),
            r: r.clone(),
            M0: &M0,
        };
        let ((_linb0, _linb1), proof_old) = decomp_old.decompose_seeded_base(&scheme, B);

        let decomp_new = DecompBase { f, r, M0: &M0 };
        let proof_new = decomp_new.decompose_seeded_base_one_shot(&scheme, B);

        assert_eq!(proof_old.C.0, proof_new.C.0);
        assert_eq!(proof_old.C.1, proof_new.C.1);
        assert_eq!(proof_old.v.0, proof_new.v.0);
        assert_eq!(proof_old.v.1, proof_new.v.1);
    }

    #[test]
    fn test_decomp_base_one_shot_matches_materialized_non_const_r() {
        use ark_std::UniformRand;
        let mut rng = ark_std::test_rng();
        let kappa = 2;
        let nvars = 10usize;
        let n = 1usize << nvars;
        let B = (R::dimension() / 2) as u128;

        let A = Matrix::<R>::rand(&mut rng, kappa, n);
        let scheme = AjtaiCommitmentScheme::new(A);
        let M0: Vec<Arc<SparseMatrix<<R as PolyRing>::BaseRing>>> =
            vec![Arc::new(SparseMatrix::<<R as PolyRing>::BaseRing>::identity(n))];

        let f = (0..n).map(|_| R::rand(&mut rng)).collect::<Vec<_>>();
        // Non-const ring challenges to force the ring-eq branch.
        let r = (0..nvars)
            .map(|_| {
                let x = R::rand(&mut rng);
                (x, x)
            })
            .collect::<Vec<_>>();

        let decomp_old = DecompBase {
            f: f.clone(),
            r: r.clone(),
            M0: &M0,
        };
        let ((_linb0, _linb1), proof_old) = decomp_old.decompose_seeded_base(&scheme, B);

        let decomp_new = DecompBase { f, r, M0: &M0 };
        let proof_new = decomp_new.decompose_seeded_base_one_shot(&scheme, B);

        assert_eq!(proof_old.C.0, proof_new.C.0);
        assert_eq!(proof_old.C.1, proof_new.C.1);
        assert_eq!(proof_old.v.0, proof_new.v.0);
        assert_eq!(proof_old.v.1, proof_new.v.1);
    }
}
