use ark_std::log2;
use latticefold::commitment::AjtaiCommitmentScheme;
use stark_rings::{
    balanced_decomposition::{recompose, Decompose, DecomposeToVec},
    PolyRing, Zq,
};
use stark_rings_linalg::{ops::Transpose, Matrix, SparseMatrix};
use std::time::Instant;
use std::sync::Arc;

use crate::lin::{LinB, LinBX};
use crate::rgchk::WitnessVec;
use crate::utils::maybe_print_rss;

pub type RxR<R> = (R, R);

#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[derive(Debug)]
pub struct Decomp<'a, R> {
    pub f: Vec<R>,
    pub r: Vec<(R, R)>,
    pub M: &'a [Arc<SparseMatrix<R>>],
}

#[derive(Clone, Debug)]
pub struct DecompProof<R> {
    /// C = com(F)
    pub C: (Vec<R>, Vec<R>), // kappa x 2
    pub v: (Vec<RxR<R>>, Vec<RxR<R>>), // (v(0), v(1))
}

impl<R: PolyRing> Decomp<'_, R>
where
    R: Decompose,
    R::BaseRing: Zq,
{
    pub fn decompose(&self, A: &Matrix<R>, B: u128) -> ((LinB<R>, LinB<R>), DecompProof<R>) {
        let profile = std::env::var("LF_PLUS_PROFILE").ok().as_deref() == Some("1");
        let t_total = Instant::now();

        let nvars = log2(A.ncols) as usize;
        let mut F = self.f.decompose_to_vec(B, 2).transpose().into_iter();
        let F0 = F.next().unwrap();
        let F1 = F.next().unwrap();

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

        let detail = std::env::var("LF_PLUS_PROFILE_DETAIL").ok().as_deref() == Some("1");

        // Precompute eq-weights once per point; shared across both Fi branches.
        let t_eq = Instant::now();
        let eq_a = eq_weights::<R>(&r_a);
        let eq_b = eq_weights::<R>(&r_b);
        if profile && detail {
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
            if profile && detail {
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
            if profile && detail {
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
        &self,
        scheme: &AjtaiCommitmentScheme<R>,
        B: u128,
    ) -> ((LinB<R>, LinB<R>), DecompProof<R>) {
        let profile = std::env::var("LF_PLUS_PROFILE").ok().as_deref() == Some("1");
        let t_total = Instant::now();
        maybe_print_rss("decomp_seeded: start");

        let nvars = log2(scheme.width()) as usize;
        let mut F = self.f.decompose_to_vec(B, 2).transpose().into_iter();
        let F0 = F.next().unwrap();
        let F1 = F.next().unwrap();
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
            let detail = std::env::var("LF_PLUS_PROFILE_DETAIL").ok().as_deref() == Some("1");

            let t_eq = Instant::now();
            let eq_a = eq_weights::<R>(&r_a);
            let eq_b = eq_weights::<R>(&r_b);
            maybe_print_rss("decomp_seeded: after eq_weights");
            if profile && detail {
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
            let fv0 = (dot_with_eq::<R>(&F0, &eq_a), dot_with_eq::<R>(&F0, &eq_b));
            let fv1 = (dot_with_eq::<R>(&F1, &eq_a), dot_with_eq::<R>(&F1, &eq_b));
            if profile && detail {
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
                    let (m0, m1) =
                        eval_sparse_mat_two_vecs_at_two_points::<R>(M_i, &F0, &F1, &eq_a, &eq_b);
                    v0.push(m0);
                    v1.push(m1);
                }
            }
            if profile && detail {
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
}
