use ark_std::log2;
use stark_rings::{
    balanced_decomposition::{recompose, Decompose, DecomposeToVec},
    PolyRing, Zq,
};
use stark_rings_linalg::{ops::Transpose, Matrix, SparseMatrix};
use std::time::Instant;
use std::sync::Arc;

use crate::lin::{LinB, LinBX};

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
            let mut w = vec![Rr::ZERO; n];
            w[0] = Rr::ONE;
            let mut len = 1usize;
            for &rj in r.iter().rev() {
                let om = Rr::ONE - rj;
                // Split `[0..2*len)` into low/high halves and update in-place:
                //   high[i] = old_low[i] * rj
                //   low[i]  = old_low[i] * (1-rj)
                let (low, high) = w[..(2 * len)].split_at_mut(len);
                #[cfg(feature = "parallel")]
                {
                    low.par_iter_mut()
                        .zip(high.par_iter_mut())
                        .for_each(|(lo, hi)| {
                            let old = *lo;
                            *hi = old * rj;
                            *lo = old * om;
                        });
                }
                #[cfg(not(feature = "parallel"))]
                {
                    for i in 0..len {
                        let old = low[i];
                        high[i] = old * rj;
                        low[i] = old * om;
                    }
                }
                len <<= 1;
            }
            w
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

        // Evaluate (M * f)(r_a) and (M * f)(r_b) without materializing M*f or building per-matrix MLEs:
        //   (M*f)(r) = Σ_row eq_r[row] * (Σ_{(c,col)∈row} c * f[col]).
        #[inline]
        fn eval_sparse_mat_times_vec_at_two_points<Rr: PolyRing>(
            m: &SparseMatrix<Rr>,
            f: &[Rr],
            eq_a: &[Rr],
            eq_b: &[Rr],
        ) -> (Rr, Rr) {
            debug_assert_eq!(m.ncols, f.len());
            debug_assert_eq!(m.nrows, eq_a.len());
            debug_assert_eq!(m.nrows, eq_b.len());
            #[cfg(feature = "parallel")]
            {
                m.coeffs
                    .par_iter()
                    .enumerate()
                    .map(|(row_idx, row)| {
                        let mut row_dot = Rr::ZERO;
                        for (coeff, col_idx) in row {
                            if *col_idx < f.len() {
                                row_dot += *coeff * f[*col_idx];
                            }
                        }
                        (eq_a[row_idx] * row_dot, eq_b[row_idx] * row_dot)
                    })
                    .reduce(|| (Rr::ZERO, Rr::ZERO), |(a0, b0), (a1, b1)| (a0 + a1, b0 + b1))
            }
            #[cfg(not(feature = "parallel"))]
            {
                m.coeffs
                    .iter()
                    .enumerate()
                    .fold((Rr::ZERO, Rr::ZERO), |(aa, bb), (row_idx, row)| {
                        let mut row_dot = Rr::ZERO;
                        for (coeff, col_idx) in row {
                            if *col_idx < f.len() {
                                row_dot += *coeff * f[*col_idx];
                            }
                        }
                        (aa + eq_a[row_idx] * row_dot, bb + eq_b[row_idx] * row_dot)
                    })
            }
        }

        let detail = std::env::var("LF_PLUS_PROFILE_DETAIL").ok().as_deref() == Some("1");

        let vi_calc = |Fi: &[R]| -> Vec<(R, R)> {
            // Precompute eq-weights once per point; reuse for fv and all M_i*Fi evaluations.
            let t_eq = Instant::now();
            let eq_a = eq_weights::<R>(&r_a);
            let eq_b = eq_weights::<R>(&r_b);
            if profile && detail {
                println!("[LF+ Decomp::decompose] eq_weights: {:?}", t_eq.elapsed());
            }

            let t_fv = Instant::now();
            let fv = (dot_with_eq::<R>(Fi, &eq_a), dot_with_eq::<R>(Fi, &eq_b));
            if profile && detail {
                println!("[LF+ Decomp::decompose] fv(dot_with_eq): {:?}", t_fv.elapsed());
            }

            let mut vi = vec![fv];
            let t_mats = Instant::now();
            for M_i in self.M {
                let M_i = M_i.as_ref();
                // Hot-path for the common identity/permutation test harness case: M_i == I implies M_i * Fi == Fi.
                if is_identity_matrix::<R>(M_i) {
                    vi.push(fv);
                } else {
                    vi.push(eval_sparse_mat_times_vec_at_two_points::<R>(
                        M_i, Fi, &eq_a, &eq_b,
                    ));
                }
            }
            if profile && detail {
                println!(
                    "[LF+ Decomp::decompose] mats(eval_sparse_mat_times_vec_at_two_points): {:?} (Mlen={})",
                    t_mats.elapsed(),
                    self.M.len()
                );
            }
            vi
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
        let v0 = vi_calc(&F0);
        let v1 = vi_calc(&F1);
        if profile {
            println!("[LF+ Decomp::decompose] compute v0/v1: {:?}", t.elapsed());
        }

        let t = Instant::now();
        let C0 = A.try_mul_vec(&F0).unwrap();
        let C1 = A.try_mul_vec(&F1).unwrap();
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
            f: F0,
        };
        let linb1 = LinB {
            x: LinBX {
                cm_f: C1.clone(),
                r: self.r.clone(),
                v: v1.clone(),
            },
            f: F1,
        };
        let proof = DecompProof {
            C: (C0, C1),
            v: (v0, v1),
        };

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
            f: cr1cs.f,
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
