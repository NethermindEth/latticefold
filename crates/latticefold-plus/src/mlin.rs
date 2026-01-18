use ark_std::log2;
use latticefold::commitment::AjtaiCommitmentScheme;
use latticefold::transcript::Transcript;
use stark_rings::{
    balanced_decomposition::{convertible_ring::ConvertibleRing, Decompose},
    CoeffRing, PolyRing, Zq,
};
use stark_rings_linalg::{Matrix, SparseMatrix};
use std::sync::Arc;
use std::time::Instant;

use crate::{
    cm::{Cm, CmProof},
    lin::{LinB, LinParameters},
    rgchk::{Rg, RgInstance},
};

#[derive(Clone, Debug)]
pub struct Mlin<R: PolyRing> {
    pub lins: Vec<LinB<R>>,
    pub params: LinParameters,
}

#[derive(Clone, Debug)]
pub struct LinB2X<R: PolyRing> {
    pub cm_g: Vec<R>,
    pub ro: Vec<(R, R)>,
    pub vo: Vec<(R, R)>,
}

#[derive(Clone, Debug)]
pub struct LinB2<R: PolyRing> {
    pub g: Vec<R>,
    pub x: LinB2X<R>,
}

impl<R: CoeffRing + PolyRing> Mlin<R>
where
    R::BaseRing: ConvertibleRing + Decompose + Zq,
    R: Decompose,
{
    /// Πmlin protocol
    ///
    /// Folds L `LinB` instances.
    pub fn mlin(
        &self,
        A: &Matrix<R>,
        M: &[Arc<SparseMatrix<R>>],
        transcript: &mut impl Transcript<R>,
    ) -> (LinB2<R>, CmProof<R>) {
        let profile = std::env::var("LF_PLUS_PROFILE").ok().as_deref() == Some("1");
        let t_total = Instant::now();
        let n = self.lins[0].f.len();

        let t = Instant::now();
        let instances = self
            .lins
            .iter()
            .map(|lin| match &lin.f {
                crate::rgchk::WitnessVec::Ring(vr) => {
                    // Non-seeded path is used only for small/unit tests; cloning here is OK.
                    RgInstance::from_f(vr.as_ref().clone(), A, &self.params.decomp)
                }
                crate::rgchk::WitnessVec::ConstCoeffBase { values: v0, domain_len } => {
                    // Fallback for non-seeded path: materialize ring elements (zero-padded).
                    let mut f_ring = vec![R::ZERO; *domain_len];
                    for (i, &x) in v0.iter().enumerate() {
                        f_ring[i] = R::from(x);
                    }
                    RgInstance::from_f(f_ring, A, &self.params.decomp)
                }
            })
            .collect::<Vec<_>>();
        if profile {
            println!(
                "[LF+ Mlin::mlin] build instances: {:?} (L={}, n={}, kappa={})",
                t.elapsed(),
                self.lins.len(),
                n,
                self.params.kappa
            );
        }

        let rg = Rg {
            nvars: log2(n) as usize,
            instances,
            dparams: self.params.decomp.clone(),
        };

        let cm = Cm { rg };

        let t = Instant::now();
        let (com, proof) = cm.prove(M, transcript);
        if profile {
            println!("[LF+ Mlin::mlin] Cm::prove: {:?}", t.elapsed());
        }

        // IMPORTANT (peak RSS):
        // Avoid allocating a fresh `vec![0; n]` accumulator for large vectors (notably `g`).
        // For L=1 this would double peak memory. Instead, reuse the first vector as the accumulator
        // and add remaining vectors into it.
        let crate::cm::Com { g: g_vecs, x: com_x } = com;

        let mut cmg_iter = com_x.cm_g.into_iter();
        let mut cm_g = cmg_iter
            .next()
            .unwrap_or_else(|| vec![R::zero(); self.params.kappa]);
        for cm in cmg_iter {
            for (acc_r, cm_r) in cm_g.iter_mut().zip(cm.into_iter()) {
                *acc_r += cm_r;
            }
        }

        let mut vo_iter = com_x.vo.into_iter();
        let mut vo = vo_iter
            .next()
            .unwrap_or_else(|| Vec::<(R, R)>::new());
        for v in vo_iter {
            debug_assert_eq!(v.len(), vo.len());
            for (acc, vv) in vo.iter_mut().zip(v.into_iter()) {
                acc.0 += vv.0;
                acc.1 += vv.1;
            }
        }

        let x = LinB2X {
            cm_g,
            ro: com_x.ro,
            vo,
        };

        let mut g_iter = g_vecs.into_iter();
        let mut g = g_iter.next().unwrap_or_else(|| vec![R::zero(); n]);
        for gi in g_iter {
            debug_assert_eq!(gi.len(), g.len());
            for (acc_r, gi_r) in g.iter_mut().zip(gi.into_iter()) {
                *acc_r += gi_r;
            }
        }
        let linb2 = LinB2 { g, x };

        if profile {
            println!("[LF+ Mlin::mlin] total: {:?}", t_total.elapsed());
        }

        (linb2, proof)
    }

    /// Πmlin protocol, but using an **implicitly-defined** Ajtai matrix (seeded) instead of an explicit dense `Matrix<R>`.
    ///
    /// This is intended for very large `n` where materializing a `kappa × n` Ajtai matrix is not viable.
    /// The verifier-side protocol and proof objects are unchanged; only prover-side commitment computation differs.
    pub fn mlin_seeded(
        &self,
        scheme: &AjtaiCommitmentScheme<R>,
        M: &[Arc<SparseMatrix<R>>],
        transcript: &mut impl Transcript<R>,
    ) -> (LinB2<R>, CmProof<R>) {
        let profile = std::env::var("LF_PLUS_PROFILE").ok().as_deref() == Some("1");
        let t_total = Instant::now();
        let n = self.lins[0].f.len();

        let t = Instant::now();
        let mut instances = Vec::with_capacity(self.lins.len());
        let mut n_f0 = 0usize;
        let mut n_ring = 0usize;
        for (i, lin) in self.lins.iter().enumerate() {
            match &lin.f {
                crate::rgchk::WitnessVec::ConstCoeffBase { values: v0, .. } => {
                    n_f0 += 1;
                    if profile {
                        println!(
                            "[LF+ Mlin::mlin_seeded] instance[{i}] witness=ConstCoeffBase(len={}) -> RgInstance::from_f0_seeded",
                            v0.len()
                        );
                    }
                    instances.push(RgInstance::from_f0_seeded(
                        v0.clone(),
                        scheme,
                        &self.params.decomp,
                    ));
                }
                crate::rgchk::WitnessVec::Ring(vr) => {
                    n_ring += 1;
                    if profile {
                        println!(
                            "[LF+ Mlin::mlin_seeded] instance[{i}] witness=Ring(len={}) -> RgInstance::from_f_seeded",
                            vr.len()
                        );
                    }
                    // Fallback for small cases: keep existing ring-vector seeded constructor.
                    instances.push(RgInstance::from_f_seeded(
                        vr.as_ref().clone(),
                        scheme,
                        &self.params.decomp,
                    ));
                }
            }
        }
        if profile {
            println!(
                "[LF+ Mlin::mlin_seeded] build instances: {:?} (L={}, n={}, kappa={}, f0_instances={}, ring_instances={})",
                t.elapsed(),
                self.lins.len(),
                n,
                self.params.kappa,
                n_f0,
                n_ring
            );
        }

        let rg = Rg {
            nvars: log2(n) as usize,
            instances,
            dparams: self.params.decomp.clone(),
        };
        let cm = Cm { rg };

        let t = Instant::now();
        let (com, proof) = cm.prove(M, transcript);
        if profile {
            println!("[LF+ Mlin::mlin_seeded] Cm::prove: {:?}", t.elapsed());
        }

        let crate::cm::Com { g: g_vecs, x: com_x } = com;

        let mut cmg_iter = com_x.cm_g.into_iter();
        let mut cm_g = cmg_iter
            .next()
            .unwrap_or_else(|| vec![R::zero(); self.params.kappa]);
        for cm in cmg_iter {
            for (acc_r, cm_r) in cm_g.iter_mut().zip(cm.into_iter()) {
                *acc_r += cm_r;
            }
        }

        let mut vo_iter = com_x.vo.into_iter();
        let mut vo = vo_iter
            .next()
            .unwrap_or_else(|| Vec::<(R, R)>::new());
        for v in vo_iter {
            debug_assert_eq!(v.len(), vo.len());
            for (acc, vv) in vo.iter_mut().zip(v.into_iter()) {
                acc.0 += vv.0;
                acc.1 += vv.1;
            }
        }

        let x = LinB2X {
            cm_g,
            ro: com_x.ro,
            vo,
        };

        let mut g_iter = g_vecs.into_iter();
        let mut g = g_iter.next().unwrap_or_else(|| vec![R::zero(); n]);
        for gi in g_iter {
            debug_assert_eq!(gi.len(), g.len());
            for (acc_r, gi_r) in g.iter_mut().zip(gi.into_iter()) {
                *acc_r += gi_r;
            }
        }
        let linb2 = LinB2 { g, x };

        if profile {
            println!("[LF+ Mlin::mlin_seeded] total: {:?}", t_total.elapsed());
        }
        (linb2, proof)
    }

    /// Πmlin protocol, but with external matrices represented over the **base ring**.
    ///
    /// This avoids materializing `SparseMatrix<R>` when the matrices are constant-coefficient by
    /// construction (SP1/R1LF regime), which is catastrophic at large ring dimension (e.g. d64).
    pub fn mlin_seeded_base(
        &self,
        scheme: &AjtaiCommitmentScheme<R>,
        M0: &[Arc<SparseMatrix<R::BaseRing>>],
        transcript: &mut impl Transcript<R>,
    ) -> (LinB2<R>, CmProof<R>) {
        let profile = std::env::var("LF_PLUS_PROFILE").ok().as_deref() == Some("1");
        let t_total = Instant::now();
        let n = self.lins[0].f.len();

        let t = Instant::now();
        let mut instances = Vec::with_capacity(self.lins.len());
        let mut n_f0 = 0usize;
        let mut n_ring = 0usize;
        for (i, lin) in self.lins.iter().enumerate() {
            match &lin.f {
                crate::rgchk::WitnessVec::ConstCoeffBase { values: v0, .. } => {
                    n_f0 += 1;
                    if profile {
                        println!(
                            "[LF+ Mlin::mlin_seeded_base] instance[{i}] witness=ConstCoeffBase(len={}) -> RgInstance::from_f0_seeded",
                            v0.len()
                        );
                    }
                    instances.push(RgInstance::from_f0_seeded(
                        v0.clone(),
                        scheme,
                        &self.params.decomp,
                    ));
                }
                crate::rgchk::WitnessVec::Ring(vr) => {
                    n_ring += 1;
                    if profile {
                        println!(
                            "[LF+ Mlin::mlin_seeded_base] instance[{i}] witness=Ring(len={}) -> RgInstance::from_f_seeded",
                            vr.len()
                        );
                    }
                    instances.push(RgInstance::from_f_seeded(
                        vr.as_ref().clone(),
                        scheme,
                        &self.params.decomp,
                    ));
                }
            }
        }
        if profile {
            println!(
                "[LF+ Mlin::mlin_seeded_base] build instances: {:?} (L={}, n={}, kappa={}, f0_instances={}, ring_instances={})",
                t.elapsed(),
                self.lins.len(),
                n,
                self.params.kappa,
                n_f0,
                n_ring
            );
        }

        let rg = Rg {
            nvars: log2(n) as usize,
            instances,
            dparams: self.params.decomp.clone(),
        };
        let cm = Cm { rg };

        let t = Instant::now();
        let (com, proof) = cm.prove_base(M0, transcript);
        if profile {
            println!("[LF+ Mlin::mlin_seeded_base] Cm::prove_base: {:?}", t.elapsed());
        }

        let crate::cm::Com { g: g_vecs, x: com_x } = com;

        let mut cmg_iter = com_x.cm_g.into_iter();
        let mut cm_g = cmg_iter
            .next()
            .unwrap_or_else(|| vec![R::zero(); self.params.kappa]);
        for cm in cmg_iter {
            for (acc_r, cm_r) in cm_g.iter_mut().zip(cm.into_iter()) {
                *acc_r += cm_r;
            }
        }

        let mut vo_iter = com_x.vo.into_iter();
        let mut vo = vo_iter
            .next()
            .unwrap_or_else(|| Vec::<(R, R)>::new());
        for v in vo_iter {
            debug_assert_eq!(v.len(), vo.len());
            for (acc, vv) in vo.iter_mut().zip(v.into_iter()) {
                acc.0 += vv.0;
                acc.1 += vv.1;
            }
        }

        let x = LinB2X {
            cm_g,
            ro: com_x.ro,
            vo,
        };

        let mut g_iter = g_vecs.into_iter();
        let mut g = g_iter.next().unwrap_or_else(|| vec![R::zero(); n]);
        for gi in g_iter {
            debug_assert_eq!(gi.len(), g.len());
            for (acc_r, gi_r) in g.iter_mut().zip(gi.into_iter()) {
                *acc_r += gi_r;
            }
        }
        let linb2 = LinB2 { g, x };

        if profile {
            println!("[LF+ Mlin::mlin_seeded_base] total: {:?}", t_total.elapsed());
        }

        (linb2, proof)
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::PrimeField;
    use ark_std::One;
    use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
    use latticefold::arith::r1cs::R1CS;
    use stark_rings::{cyclotomic_ring::models::frog_ring::RqPoly as R, PolyRing};
    use stark_rings_linalg::SparseMatrix;

    use super::*;
    use crate::{
        lin::{Linearize, LinearizedVerify},
        r1cs::{r1cs_decomposed_square, ComR1CS},
        rgchk::DecompParameters,
        transcript::PoseidonTranscript,
    };

    #[test]
    fn test_mlin() {
        let n = 1 << 15;
        let k = 2;
        let B = 2;
        let b = (R::dimension() / 2) as u128;
        let kappa = 2;
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

        let mut r1cs = r1cs_decomposed_square(
            R1CS::<R> {
                l: 1,
                A: SparseMatrix::identity(n / k),
                B: SparseMatrix::identity(n / k),
                C: SparseMatrix::identity(n / k),
            },
            n,
            B,
            k,
        );

        r1cs.A.coeffs[0][0].0 = 2u128.into();
        r1cs.C.coeffs[0][0].0 = 2u128.into();

        let A = Matrix::<R>::rand(&mut ark_std::test_rng(), params.kappa, n);

        let cr1cs0 = ComR1CS::new(r1cs.clone(), z0, 1, B, k, &A);
        let cr1cs1 = ComR1CS::new(r1cs, z1, 1, B, k, &A);

        let mut ts = PoseidonTranscript::empty::<PC>();
        let (linb0, lproof0) = cr1cs0.linearize(&mut ts);
        let (linb1, lproof1) = cr1cs1.linearize(&mut ts);

        let M = cr1cs0.x.matrices_arc();

        let A = Matrix::<R>::rand(&mut ark_std::test_rng(), params.kappa, n);

        let mlin = Mlin {
            lins: vec![linb0, linb1],
            params,
        };

        let (_linb2, cmproof) = mlin.mlin(&A, &M, &mut ts);

        let mut ts = PoseidonTranscript::empty::<PC>();
        lproof0.verify(&mut ts);
        lproof1.verify(&mut ts);
        cmproof.verify(&M, &mut ts).unwrap();
    }
}
