use ark_std::{
    log2,
    ops::{Mul, Sub},
    One,
};
use latticefold::{
    transcript::Transcript,
    utils::sumcheck::{
        utils::eq_eval,
        MLSumcheck, Proof, SumCheckError,
    },
};
use stark_rings::{unit_monomial, CoeffRing, OverField, PolyRing, Ring, Zq};
use stark_rings_linalg::SparseMatrix;
use std::sync::Arc;
use std::time::Instant;

use crate::{
    rgchk::{Dcom, Rg},
    streaming_sumcheck::{HCol0Precomp, StreamingMleEnum, StreamingSumcheck},
    utils::{maybe_print_rss, short_challenge, tensor, tensor_product},
};

use crate::rgchk::WitnessVec;

#[inline]
fn is_const_coeff_ring<R: PolyRing>(x: &R) -> bool {
    let coeffs = x.coeffs();
    // constant term can be anything; all higher coeffs must be zero
    coeffs.iter().skip(1).all(|c| *c == <R as PolyRing>::BaseRing::ZERO)
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

fn try_as_base_scalars<R: PolyRing>(v: &[R]) -> Option<Vec<R::BaseRing>> {
    let mut out = Vec::with_capacity(v.len());
    for x in v {
        if !is_const_coeff_ring::<R>(x) {
            return None;
        }
        out.push(x.coeffs()[0]);
    }
    Some(out)
}

#[derive(Clone, Debug)]
pub struct Cm<R: PolyRing> {
    pub rg: Rg<R>,
}

// eval over r_o of [tau (a), m_tau (b), f (c), h (u)] over 1 + n_lin
#[derive(Clone, Debug)]
pub struct InstanceEvals<R>(Vec<[R; 4]>);

impl<R> InstanceEvals<R> {
    pub(crate) fn rows(&self) -> &[[R; 4]] {
        &self.0
    }
}

#[derive(Clone, Debug)]
pub struct CmProof<R: PolyRing> {
    pub dcom: Dcom<R>,
    pub comh: Vec<Vec<R>>,
    pub sumcheck_proofs: (Proof<R>, Proof<R>),
    pub evals: (Vec<InstanceEvals<R>>, Vec<InstanceEvals<R>>),
}

#[derive(Clone, Debug)]
pub struct Com<R> {
    pub g: Vec<Vec<R>>,
    pub x: ComX<R>,
}

#[derive(Clone, Debug)]
pub struct ComX<R> {
    pub cm_g: Vec<Vec<R>>,
    pub ro: Vec<(R, R)>,
    pub vo: Vec<Vec<(R, R)>>,
}

impl<R: CoeffRing> Cm<R>
where
    R::BaseRing: Zq,
{
    pub fn prove(
        &self,
        M: &[Arc<SparseMatrix<R>>],
        transcript: &mut impl Transcript<R>,
    ) -> (Com<R>, CmProof<R>) {
        let profile = std::env::var("LF_PLUS_PROFILE").ok().as_deref() == Some("1");
        let t_total = Instant::now();
        maybe_print_rss("cm: prove start");

        let k = self.rg.dparams.k;
        let d = R::dimension();
        let dp = R::dimension() / 2;
        let l = self.rg.dparams.l;
        let n = self.rg.instances[0].tau.len();

        if profile {
            #[cfg(feature = "parallel")]
            println!(
                "[LF+ Cm::prove] start: n={} nvars={} Mlen={} rayon_threads={}",
                n,
                self.rg.nvars,
                M.len(),
                rayon::current_num_threads()
            );
            #[cfg(not(feature = "parallel"))]
            println!(
                "[LF+ Cm::prove] start: n={} nvars={} Mlen={} rayon_threads=DISABLED(feature=parallel)",
                n,
                self.rg.nvars,
                M.len(),
            );
        }

        let t = Instant::now();
        let dcom = self.rg.range_check(M, transcript);
        if profile {
            println!("[LF+ Cm::prove] range_check: {:?}", t.elapsed());
        }
        maybe_print_rss("cm: after range_check");

        let s = (0..3)
            .map(|_| short_challenge(128, transcript))
            .collect::<Vec<R>>();

        let s_prime = (0..k)
            .map(|_| {
                (0..d)
                    .map(|_| short_challenge(128, transcript))
                    .collect::<Vec<R>>()
            })
            .collect::<Vec<_>>();
        let s_prime_flat = s_prime.clone().into_iter().flatten().collect::<Vec<R>>();

        // Optional: avoid materializing `h` as a full length-2^n ring vector.
        // This reduces peak RSS substantially (especially for d64), at the cost of extra compute
        // when evaluating `h`/`M*h` during CM sumcheck.
        let stream_h = std::env::var("LF_PLUS_CM_STREAM_H").ok().as_deref() == Some("1");
        let t = Instant::now();
        maybe_print_rss("cm: build_h start");

        // Build an *immutable* on-demand `h` MLE for each instance (used for g and for `M*h` MLEs).
        // If any `M_f` is not `ConstCol0`, we fall back to materializing `h`.
        let mut h_mles_full: Option<Vec<Arc<StreamingMleEnum<R>>>> = if stream_h {
            Some(Vec::with_capacity(self.rg.instances.len()))
        } else {
            None
        };

        // Materialized `h` (legacy) if `stream_h` is off or if we can't stream due to `Full` backing.
        let mut h_vecs: Option<Vec<Vec<R>>> = if stream_h { None } else { Some(Vec::with_capacity(self.rg.instances.len())) };

        #[inline]
        fn mon_info<Rr: PolyRing>(mono: &Rr) -> Option<(usize, Rr::BaseRing)>
        where
            Rr::BaseRing: Ring,
        {
            let coeffs = mono.coeffs();
            let mut found: Option<(usize, Rr::BaseRing)> = None;
            for (i, &ci) in coeffs.iter().enumerate() {
                if ci != Rr::BaseRing::ZERO {
                    if found.is_some() {
                        return None;
                    }
                    found = Some((i, ci));
                }
            }
            found
        }
        #[inline]
        fn mul_negacyclic_by_monomial<Rr>(a: &Rr, shift: usize, scale: Rr::BaseRing) -> Rr
        where
            Rr: PolyRing,
            Rr::BaseRing: Ring + Copy,
        {
            if scale == Rr::BaseRing::ZERO {
                return Rr::ZERO;
            }
            let ac = a.coeffs();
            let d = ac.len();
            if shift == 0 && scale == Rr::BaseRing::ONE {
                return *a;
            }
            let mut out = Rr::ZERO;
            let outc = out.coeffs_mut();
            for i in 0..d {
                let v = ac[i] * scale;
                if v == Rr::BaseRing::ZERO {
                    continue;
                }
                let j = i + shift;
                if j < d {
                    outc[j] += v;
                } else {
                    outc[j - d] -= v;
                }
            }
            out
        }

        for inst in self.rg.instances.iter() {
            let n = 1 << self.rg.nvars;
            maybe_print_rss("cm: build_h one inst start");

            if let Some(hs) = h_mles_full.as_mut() {
                // Streaming path requires all M_f to be ConstCol0.
                let mut precomps: Vec<HCol0Precomp<R>> = Vec::with_capacity(inst.M_f.len());
                let mut ok = true;
                for (M, s_i) in inst.M_f.iter().zip(s_prime.iter()) {
                    match &M.digits {
                        crate::setchk::DigitsBacking::ConstCol0 { col0, zero_idx } => {
                            let mut mi_tab = Vec::with_capacity(M.exp_table.len());
                            for r in M.exp_table.iter() {
                                mi_tab.push(mon_info::<R>(r).expect("exp_table entry must be monomial"));
                            }
                            let s0 = s_i[0];
                            let rest_sum = s_i.iter().skip(1).copied().sum::<R>();
                            let term0_tab: Arc<Vec<R>> = Arc::new(
                                mi_tab
                                    .iter()
                                    .map(|(shift, scale)| mul_negacyclic_by_monomial::<R>(&s0, *shift, *scale))
                                    .collect::<Vec<_>>(),
                            );
                            let (shift0, scale0) = mi_tab[*zero_idx as usize];
                            let term_rest = if shift0 == 0 && scale0 == R::BaseRing::ONE {
                                rest_sum
                            } else {
                                mul_negacyclic_by_monomial::<R>(&rest_sum, shift0, scale0)
                            };
                            precomps.push(HCol0Precomp {
                                col0: col0.clone(),
                                zero_idx: *zero_idx,
                                term0_tab,
                                term_rest,
                            });
                        }
                        crate::setchk::DigitsBacking::Full(_) => {
                            ok = false;
                            break;
                        }
                    }
                }
                if ok {
                    let mle = StreamingMleEnum::HFromMfDigitsConstCol0 {
                        precomps: Arc::new(precomps),
                        num_vars: self.rg.nvars,
                    };
                    hs.push(Arc::new(mle));
                    maybe_print_rss("cm: build_h one inst done");
                    continue;
                }
                // If streaming is requested but unsupported, fall back to materializing for all instances.
                h_mles_full = None;
                h_vecs = Some(Vec::with_capacity(self.rg.instances.len()));
            }

            // Materialize `h` (legacy)
            let mut h = vec![R::ZERO; n];
            #[cfg(feature = "parallel")]
            {
                use rayon::prelude::*;
                h.par_iter_mut().enumerate().for_each(|(row, out)| {
                    let mut acc = R::ZERO;
                    for (i, M) in inst.M_f.iter().enumerate() {
                        let s_i = &s_prime[i];
                        for col in 0..M.ncols {
                            acc += M.get(row, col) * s_i[col];
                        }
                    }
                    *out = acc;
                });
            }
            #[cfg(not(feature = "parallel"))]
            {
                for row in 0..n {
                    let mut acc = R::ZERO;
                    for (i, M) in inst.M_f.iter().enumerate() {
                        let s_i = &s_prime[i];
                        for col in 0..M.ncols {
                            acc += M.get(row, col) * s_i[col];
                        }
                    }
                    h[row] = acc;
                }
            }
            h_vecs.as_mut().unwrap().push(h);
            maybe_print_rss("cm: build_h one inst done");
        }

        if profile {
            println!("[LF+ Cm::prove] build h: {:?}", t.elapsed());
        }
        maybe_print_rss("cm: build_h done");

        let t = Instant::now();
        let comh: Vec<Vec<R>> = self
            .rg
            .instances
            .iter()
            .map(|inst| {
                let comh_vectors = inst
                    .comM_f
                    .iter()
                    .zip(s_prime.iter())
                    .map(|(comM_f_i, s_i)| comM_f_i.try_mul_vec(s_i).unwrap())
                    .collect::<Vec<_>>();

                let mut comh = vec![R::zero(); inst.comM_f[0].nrows];
                for v in comh_vectors {
                    for (i, val) in v.iter().enumerate() {
                        comh[i] += *val;
                    }
                }
                comh
            })
            .collect();
        if profile {
            println!("[LF+ Cm::prove] build comh: {:?}", t.elapsed());
        }

        absorb_comh(&comh, transcript);

        let kappa = comh[0].len();
        let log_kappa = log2(kappa) as usize;

        let c = (0..2)
            .map(|_| {
                transcript
                    .get_challenges(log_kappa)
                    .into_iter()
                    .map(|x| x.into())
                    .collect::<Vec<R>>()
            })
            .collect::<Vec<_>>();
        // `h` input for sumcheck:
        // - either as materialized vectors, or
        // - as on-demand MLEs (streaming path).
        let h: Vec<Arc<Vec<R>>> = match h_vecs {
            Some(vs) => vs.into_iter().map(Arc::new).collect(),
            None => Vec::new(),
        };

        let dpp = (0..l)
            .map(|i| R::from(R::BaseRing::from(dp as u128).pow([i as u64])))
            .collect::<Vec<_>>();
        let xp = (0..d).map(|i| unit_monomial::<R>(i)).collect::<Vec<_>>();

        // Build *structured* tensor tables without materializing O(n) vectors.
        let t = Instant::now();
        let tensor_c0 = crate::utils::tensor(&c[0]);
        let tensor_c1 = crate::utils::tensor(&c[1]);
        let tensor_len = tensor_c0.len() * s_prime_flat.len() * dpp.len() * xp.len();
        assert_eq!(tensor_c0.len(), tensor_c1.len());
        if tensor_len > n {
            panic!("t(z) tensor_len {} > n {}", tensor_len, n);
        }
        let t0_mle = StreamingMleEnum::Tensor4Padded {
            t1: Arc::new(tensor_c0),
            t2: Arc::new(s_prime_flat.clone()),
            t3: Arc::new(dpp.clone()),
            t4: Arc::new(xp.clone()),
            tensor_len,
            num_vars: self.rg.nvars,
        };
        let t1_mle = StreamingMleEnum::Tensor4Padded {
            t1: Arc::new(tensor_c1),
            t2: Arc::new(s_prime_flat.clone()),
            t3: Arc::new(dpp.clone()),
            t4: Arc::new(xp.clone()),
            tensor_len,
            num_vars: self.rg.nvars,
        };
        if profile {
            println!(
                "[LF+ Cm::prove] build t(z) streaming: {:?} (tensor_len={}, padded_to_n={})",
                t.elapsed(),
                tensor_len,
                n
            );
        }

        // Share `M` matrices across both sumchecks (avoid cloning them twice).
        let t_m_arcs = Instant::now();
        // NOTE: `M` is already Arc-wrapped by the caller, so this is cheap (Arc refcount clones only).
        let m_arcs: Vec<Arc<SparseMatrix<R>>> = M.to_vec();
        if profile {
            println!(
                "[LF+ Cm::prove] build shared m_arcs: {:?} (Mlen={})",
                t_m_arcs.elapsed(),
                M.len()
            );
        }
        let mats_const = M.iter().all(|m| is_const_coeff_sparse_matrix::<R>(m.as_ref()));

        let (proof_a, evals_a, ro_a, proof_b, evals_b, ro_b, g) = if let Some(hm) = h_mles_full {
            // Streaming-h path: build g using on-demand h.
            let g = self
                .rg
                .instances
                .iter()
                .enumerate()
                .map(|(i, inst)| {
                    let n = inst.tau.len();
                    let h_mle = hm[i].clone();
                    #[cfg(feature = "parallel")]
                    {
                        use rayon::prelude::*;
                        (0..n)
                            .into_par_iter()
                            .map(|j| {
                                let r_tau = inst.tau[j];
                                let r_mtau = inst.m_tau.get(j);
                                let r_f = match &inst.f {
                                    WitnessVec::Ring(vr) => vr[j],
                                    WitnessVec::ConstCoeffBase { values: v0, .. } => {
                                        R::from(v0.get(j).copied().unwrap_or(R::BaseRing::ZERO))
                                    }
                                };
                                let r_h = h_mle.eval_at_index(j);
                                (s[0] * R::from(r_tau)) + (s[1] * r_mtau) + (s[2] * r_f) + r_h
                            })
                            .collect::<Vec<R>>()
                    }
                    #[cfg(not(feature = "parallel"))]
                    {
                        (0..n)
                            .map(|j| {
                                let r_tau = inst.tau[j];
                                let r_mtau = inst.m_tau.get(j);
                                let r_f = match &inst.f {
                                    WitnessVec::Ring(vr) => vr[j],
                                    WitnessVec::ConstCoeffBase { values: v0, .. } => {
                                        R::from(v0.get(j).copied().unwrap_or(R::BaseRing::ZERO))
                                    }
                                };
                                let r_h = h_mle.eval_at_index(j);
                                (s[0] * R::from(r_tau)) + (s[1] * r_mtau) + (s[2] * r_f) + r_h
                            })
                            .collect::<Vec<R>>()
                    }
                })
                .collect::<Vec<_>>();

            let (p0, e0, r0) = self.sumchecker_streaming_hmle(
                &dcom,
                &hm,
                &t0_mle,
                &t1_mle,
                &m_arcs,
                mats_const,
                transcript,
                profile,
            );
            let (p1, e1, r1) = self.sumchecker_streaming_hmle(
                &dcom,
                &hm,
                &t0_mle,
                &t1_mle,
                &m_arcs,
                mats_const,
                transcript,
                profile,
            );
            (p0, e0, r0, p1, e1, r1, g)
        } else {
            let (p0, e0, r0) = self.sumchecker_streaming(
                &dcom,
                &h,
                &t0_mle,
                &t1_mle,
                &m_arcs,
                mats_const,
                transcript,
                profile,
            );
            let (p1, e1, r1) = self.sumchecker_streaming(
                &dcom,
                &h,
                &t0_mle,
                &t1_mle,
                &m_arcs,
                mats_const,
                transcript,
                profile,
            );
            // Step 7 (legacy uses materialized h)
            let g = self
                .rg
                .instances
                .iter()
                .enumerate()
                .map(|(i, inst)| {
                    let n = inst.tau.len();
                    #[cfg(feature = "parallel")]
                    {
                        use rayon::prelude::*;
                        (0..n)
                            .into_par_iter()
                            .map(|j| {
                                let r_tau = inst.tau[j];
                                let r_mtau = inst.m_tau.get(j);
                                let r_f = match &inst.f {
                                    WitnessVec::Ring(vr) => vr[j],
                                    WitnessVec::ConstCoeffBase { values: v0, .. } => {
                                        R::from(v0.get(j).copied().unwrap_or(R::BaseRing::ZERO))
                                    }
                                };
                                let r_h = h[i][j];
                                (s[0] * R::from(r_tau)) + (s[1] * r_mtau) + (s[2] * r_f) + r_h
                            })
                            .collect::<Vec<R>>()
                    }
                    #[cfg(not(feature = "parallel"))]
                    {
                        (0..n)
                            .map(|j| {
                                let r_tau = inst.tau[j];
                                let r_mtau = inst.m_tau.get(j);
                                let r_f = match &inst.f {
                                    WitnessVec::Ring(vr) => vr[j],
                                    WitnessVec::ConstCoeffBase { values: v0, .. } => {
                                        R::from(v0.get(j).copied().unwrap_or(R::BaseRing::ZERO))
                                    }
                                };
                                let r_h = h[i][j];
                                (s[0] * R::from(r_tau)) + (s[1] * r_mtau) + (s[2] * r_f) + r_h
                            })
                            .collect::<Vec<R>>()
                    }
                })
                .collect::<Vec<_>>();
            (p0, e0, r0, p1, e1, r1, g)
        };

        let proof = CmProof {
            dcom,
            comh,
            sumcheck_proofs: (proof_a, proof_b),
            evals: (evals_a, evals_b),
        };

        let ro = ro_a.into_iter().zip(ro_b).collect::<Vec<_>>();

        let x = proof.x(&s, ro);

        let com = Com { g, x };

        if profile {
            println!("[LF+ Cm::prove] total: {:?}", t_total.elapsed());
        }

        (com, proof)
    }

    /// Prove, but with external matrices represented over the **base ring**.
    ///
    /// This is the natural representation for SP1/R1LF chunks (const-coeff by construction) and
    /// avoids materializing `SparseMatrix<R>` which is catastrophic at large `R::dimension()`.
    pub fn prove_base(
        &self,
        M0: &[Arc<SparseMatrix<R::BaseRing>>],
        transcript: &mut impl Transcript<R>,
    ) -> (Com<R>, CmProof<R>) {
        let profile = std::env::var("LF_PLUS_PROFILE").ok().as_deref() == Some("1");
        let t_total = Instant::now();
        maybe_print_rss("cm: prove start");

        let k = self.rg.dparams.k;
        let d = R::dimension();
        let dp = R::dimension() / 2;
        let l = self.rg.dparams.l;
        let n = self.rg.instances[0].tau.len();

        if profile {
            #[cfg(feature = "parallel")]
            println!(
                "[LF+ Cm::prove] start: n={} nvars={} Mlen={} rayon_threads={}",
                n,
                self.rg.nvars,
                M0.len(),
                rayon::current_num_threads()
            );
            #[cfg(not(feature = "parallel"))]
            println!(
                "[LF+ Cm::prove] start: n={} nvars={} Mlen={} rayon_threads=DISABLED(feature=parallel)",
                n,
                self.rg.nvars,
                M0.len(),
            );
        }

        let t = Instant::now();
        let dcom = self.rg.range_check_base(M0, transcript);
        if profile {
            println!("[LF+ Cm::prove] range_check: {:?}", t.elapsed());
        }
        maybe_print_rss("cm: after range_check");

        let s = (0..3)
            .map(|_| short_challenge(128, transcript))
            .collect::<Vec<R>>();

        let s_prime = (0..k)
            .map(|_| {
                (0..d)
                    .map(|_| short_challenge(128, transcript))
                    .collect::<Vec<R>>()
            })
            .collect::<Vec<_>>();
        let s_prime_flat = s_prime.clone().into_iter().flatten().collect::<Vec<R>>();

        let t = Instant::now();
        maybe_print_rss("cm: build_h start");
        let h_vecs: Vec<Vec<R>> = self
            .rg
            .instances
            .iter()
            .map(|inst| {
                let n = 1 << self.rg.nvars;
                maybe_print_rss("cm: build_h one inst start");
                // Same build-h as `prove` (does not depend on external M representation).
                #[inline]
                fn mon_info<Rr: PolyRing>(mono: &Rr) -> Option<(usize, Rr::BaseRing)>
                where
                    Rr::BaseRing: Ring,
                {
                    let coeffs = mono.coeffs();
                    let mut found: Option<(usize, Rr::BaseRing)> = None;
                    for (i, &ci) in coeffs.iter().enumerate() {
                        if ci != Rr::BaseRing::ZERO {
                            if found.is_some() {
                                return None;
                            }
                            found = Some((i, ci));
                        }
                    }
                    found
                }
                #[inline]
                fn mul_negacyclic_by_monomial<Rr>(a: &Rr, shift: usize, scale: Rr::BaseRing) -> Rr
                where
                    Rr: PolyRing,
                    Rr::BaseRing: Ring + Copy,
                {
                    if scale == Rr::BaseRing::ZERO {
                        return Rr::ZERO;
                    }
                    let ac = a.coeffs();
                    let d = ac.len();
                    if shift == 0 && scale == Rr::BaseRing::ONE {
                        return *a;
                    }
                    let mut out = Rr::ZERO;
                    let outc = out.coeffs_mut();
                    for i in 0..d {
                        let v = ac[i] * scale;
                        if v == Rr::BaseRing::ZERO {
                            continue;
                        }
                        let j = i + shift;
                        if j < d {
                            outc[j] += v;
                        } else {
                            outc[j - d] -= v;
                        }
                    }
                    out
                }

                struct Col0Precomp<Rr: PolyRing> {
                    col0: Arc<Vec<u16>>,
                    zero_idx: u16,
                    term0_tab: Vec<Rr>,
                    term_rest: Rr,
                }

                let mut pre = Vec::<Option<Col0Precomp<R>>>::with_capacity(inst.M_f.len());
                for (M, s_i) in inst.M_f.iter().zip(s_prime.iter()) {
                    match &M.digits {
                        crate::setchk::DigitsBacking::ConstCol0 { col0, zero_idx } => {
                            let mut mi_tab = Vec::with_capacity(M.exp_table.len());
                            for r in M.exp_table.iter() {
                                mi_tab.push(mon_info::<R>(r).expect("exp_table entry must be monomial"));
                            }
                            let s0 = s_i[0];
                            let rest_sum = s_i.iter().skip(1).copied().sum::<R>();
                            let term0_tab = mi_tab
                                .iter()
                                .map(|(shift, scale)| mul_negacyclic_by_monomial::<R>(&s0, *shift, *scale))
                                .collect::<Vec<_>>();
                            let (shift0, scale0) = mi_tab[*zero_idx as usize];
                            let term_rest = if shift0 == 0 && scale0 == R::BaseRing::ONE {
                                rest_sum
                            } else {
                                mul_negacyclic_by_monomial::<R>(&rest_sum, shift0, scale0)
                            };
                            pre.push(Some(Col0Precomp {
                                col0: col0.clone(),
                                zero_idx: *zero_idx,
                                term0_tab,
                                term_rest,
                            }));
                        }
                        crate::setchk::DigitsBacking::Full(_) => pre.push(None),
                    }
                }

                #[cfg(feature = "parallel")]
                let h = {
                    use rayon::prelude::*;
                    (0..n)
                        .into_par_iter()
                        .map(|row| {
                            let mut acc = R::ZERO;
                            for (i, M) in inst.M_f.iter().enumerate() {
                                if let Some(p) = &pre[i] {
                                    let dix = p.col0.get(row).copied().unwrap_or(p.zero_idx) as usize;
                                    acc += p.term0_tab[dix] + p.term_rest;
                                } else {
                                    let s_i = &s_prime[i];
                                    for col in 0..M.ncols {
                                        acc += M.get(row, col) * s_i[col];
                                    }
                                }
                            }
                            acc
                        })
                        .collect::<Vec<_>>()
                };
                #[cfg(not(feature = "parallel"))]
                let h = {
                    let mut h = vec![R::ZERO; n];
                    for row in 0..n {
                        let mut acc = R::ZERO;
                        for (i, M) in inst.M_f.iter().enumerate() {
                            if let Some(p) = &pre[i] {
                                let dix = p.col0.get(row).copied().unwrap_or(p.zero_idx) as usize;
                                acc += p.term0_tab[dix] + p.term_rest;
                            } else {
                                let s_i = &s_prime[i];
                                for col in 0..M.ncols {
                                    acc += M.get(row, col) * s_i[col];
                                }
                            }
                        }
                        h[row] = acc;
                    }
                    h
                };
                maybe_print_rss("cm: build_h one inst done");
                h
            })
            .collect();
        if profile {
            println!("[LF+ Cm::prove] build h: {:?}", t.elapsed());
        }
        maybe_print_rss("cm: build_h done");

        let t = Instant::now();
        let comh: Vec<Vec<R>> = self
            .rg
            .instances
            .iter()
            .map(|inst| {
                let comh_vectors = inst
                    .comM_f
                    .iter()
                    .zip(s_prime.iter())
                    .map(|(comM_f_i, s_i)| comM_f_i.try_mul_vec(s_i).unwrap())
                    .collect::<Vec<_>>();

                let mut comh = vec![R::zero(); inst.comM_f[0].nrows];
                for v in comh_vectors {
                    for (i, val) in v.iter().enumerate() {
                        comh[i] += *val;
                    }
                }
                comh
            })
            .collect();
        if profile {
            println!("[LF+ Cm::prove] build comh: {:?}", t.elapsed());
        }

        absorb_comh(&comh, transcript);

        let kappa = comh[0].len();
        let log_kappa = log2(kappa) as usize;

        let c = (0..2)
            .map(|_| {
                transcript
                    .get_challenges(log_kappa)
                    .into_iter()
                    .map(|x| x.into())
                    .collect::<Vec<R>>()
            })
            .collect::<Vec<_>>();
        let h: Vec<Arc<Vec<R>>> = h_vecs.into_iter().map(Arc::new).collect();

        let dpp = (0..l)
            .map(|i| R::from(R::BaseRing::from(dp as u128).pow([i as u64])))
            .collect::<Vec<_>>();
        let xp = (0..d).map(|i| unit_monomial::<R>(i)).collect::<Vec<_>>();

        let t = Instant::now();
        let tensor_c0 = crate::utils::tensor(&c[0]);
        let tensor_c1 = crate::utils::tensor(&c[1]);
        let tensor_len = tensor_c0.len() * s_prime_flat.len() * dpp.len() * xp.len();
        assert_eq!(tensor_c0.len(), tensor_c1.len());
        if tensor_len > n {
            panic!("t(z) tensor_len {} > n {}", tensor_len, n);
        }
        let t0_mle = StreamingMleEnum::Tensor4Padded {
            t1: Arc::new(tensor_c0),
            t2: Arc::new(s_prime_flat.clone()),
            t3: Arc::new(dpp.clone()),
            t4: Arc::new(xp.clone()),
            tensor_len,
            num_vars: self.rg.nvars,
        };
        let t1_mle = StreamingMleEnum::Tensor4Padded {
            t1: Arc::new(tensor_c1),
            t2: Arc::new(s_prime_flat.clone()),
            t3: Arc::new(dpp.clone()),
            t4: Arc::new(xp.clone()),
            tensor_len,
            num_vars: self.rg.nvars,
        };
        if profile {
            println!(
                "[LF+ Cm::prove] build t(z) streaming: {:?} (tensor_len={}, padded_to_n={})",
                t.elapsed(),
                tensor_len,
                n
            );
        }

        let t_m_arcs = Instant::now();
        let m_arcs0: Vec<Arc<SparseMatrix<R::BaseRing>>> = M0.to_vec();
        if profile {
            println!(
                "[LF+ Cm::prove] build shared m_arcs: {:?} (Mlen={})",
                t_m_arcs.elapsed(),
                M0.len()
            );
        }
        let mats_const = true;

        let (proof_a, evals_a, ro_a) = self.sumchecker_streaming_base(
            &dcom,
            &h,
            &t0_mle,
            &t1_mle,
            &m_arcs0,
            mats_const,
            transcript,
            profile,
        );
        let (proof_b, evals_b, ro_b) = self.sumchecker_streaming_base(
            &dcom,
            &h,
            &t0_mle,
            &t1_mle,
            &m_arcs0,
            mats_const,
            transcript,
            profile,
        );

        // Step 7 (unchanged)
        let g = self
            .rg
            .instances
            .iter()
            .enumerate()
            .map(|(i, inst)| {
                let n = inst.tau.len();
                debug_assert_eq!(inst.m_tau.len(), n);
                debug_assert_eq!(h[i].len(), n);
                debug_assert_eq!(inst.f.len(), n);

                #[cfg(feature = "parallel")]
                {
                    use rayon::prelude::*;
                    (0..n)
                        .into_par_iter()
                        .map(|j| {
                            let r_tau = inst.tau[j];
                            let r_mtau = inst.m_tau.get(j);
                            let r_f = match &inst.f {
                                WitnessVec::Ring(vr) => vr[j],
                                WitnessVec::ConstCoeffBase { values: v0, .. } => {
                                    R::from(v0.get(j).copied().unwrap_or(R::BaseRing::ZERO))
                                }
                            };
                            let r_h = h[i][j];
                            (s[0] * R::from(r_tau)) + (s[1] * r_mtau) + (s[2] * r_f) + r_h
                        })
                        .collect::<Vec<R>>()
                }
                #[cfg(not(feature = "parallel"))]
                {
                    (0..n)
                        .map(|j| {
                            let r_tau = inst.tau[j];
                            let r_mtau = inst.m_tau.get(j);
                            let r_f = match &inst.f {
                                WitnessVec::Ring(vr) => vr[j],
                                WitnessVec::ConstCoeffBase { values: v0, .. } => {
                                    R::from(v0.get(j).copied().unwrap_or(R::BaseRing::ZERO))
                                }
                            };
                            let r_h = h[i][j];
                            (s[0] * R::from(r_tau)) + (s[1] * r_mtau) + (s[2] * r_f) + r_h
                        })
                        .collect::<Vec<R>>()
                }
            })
            .collect::<Vec<_>>();

        let proof = CmProof {
            dcom,
            comh,
            sumcheck_proofs: (proof_a, proof_b),
            evals: (evals_a, evals_b),
        };

        let ro = ro_a.into_iter().zip(ro_b).collect::<Vec<_>>();
        let x = proof.x(&s, ro);
        let com = Com { g, x };

        if profile {
            println!("[LF+ Cm::prove] total: {:?}", t_total.elapsed());
        }
        (com, proof)
    }

    fn sumchecker_streaming(
        &self,
        dcom: &Dcom<R>,
        h: &[Arc<Vec<R>>],
        t0_mle: &StreamingMleEnum<R>,
        t1_mle: &StreamingMleEnum<R>,
        m_arcs: &[Arc<SparseMatrix<R>>],
        mats_const: bool,
        transcript: &mut impl Transcript<R>,
        profile: bool,
    ) -> (Proof<R>, Vec<InstanceEvals<R>>, Vec<R>) {
        let t_sumcheck = Instant::now();
        let nvars = self.rg.nvars;

        let rc = transcript.get_challenge();

        let L = self.rg.instances.len();

        let mut mles = Vec::with_capacity(
            1 // eq
            + L * (
                4  // [tau, m_tau, f, h]
                + 4 * m_arcs.len() // M * [tau, ...]
            )
            + 2, // t(z)
        );

        // eq table as structured base-ring MLE.
        let r0 = dcom.out.r.clone();
        let one_minus_r0 = r0.iter().copied().map(|x| R::BaseRing::ONE - x).collect();
        mles.push(StreamingMleEnum::EqBase {
            scale: R::BaseRing::ONE,
            r: r0,
            one_minus_r: one_minus_r0,
        });

        // Symphony-style conditional fast path: if the matrix coefficients AND the relevant witness
        // vectors are constant-coeff, use base-scalar mat-vec MLEs (cheaper eval_at_index and avoids
        // materializing tau as a full ring vector).
        //
        // IMPORTANT: this must be a sound detector; if false-positives occur it breaks correctness.
        // `mats_const` is computed once in `prove` and threaded through to cover both sumchecks.

        let t_build_mles = Instant::now();
        let cm_lazy: usize = std::env::var("LF_PLUS_CM_LAZY_FIX")
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(4);
        for (i, inst) in self.rg.instances.iter().enumerate() {
            // Build the base-scalar tables once and share them across:
            // - the direct MLEs for (tau, m_tau, f, h), and
            // - the const-coeff sparse mat-vec MLEs (M * vec).
            //
            // This is the only path we care about for SP1 (const-coeff matrices), and it is
            // algebraically identical to using ring tables: BaseScalarArc evaluates to `R::from(scalar)`.
            let tau0_arc: Arc<Vec<R::BaseRing>> = inst.tau.clone();
            let mtau0_arc: Option<Arc<Vec<R::BaseRing>>> = if mats_const {
                inst.m_tau
                    .as_dense_arc()
                    .and_then(|v| try_as_base_scalars::<R>(v.as_ref()).map(Arc::new))
            } else {
                None
            };
            let f0_arc: Option<Arc<Vec<R::BaseRing>>> = if mats_const {
                match &inst.f {
                    WitnessVec::ConstCoeffBase { values: v0, .. } => Some(v0.clone()),
                    WitnessVec::Ring(vr) => try_as_base_scalars::<R>(vr.as_ref()).map(Arc::new),
                }
            } else {
                None
            };
            let h0_arc: Option<Arc<Vec<R::BaseRing>>> =
                if mats_const { try_as_base_scalars::<R>(h[i].as_ref()).map(Arc::new) } else { None };

            // We apply the const-coeff optimization **per vector**, not all-or-nothing.
            //
            // - `tau` is always base-scalars by construction.
            // - `f` is const-coeff for SP1 (witness embedded as constant-coeff ring elements).
            // - `m_tau` and `h` are typically **not** const-coeff (monomials / mixed ring challenges),
            //   so insisting on them would disable the optimization in the real production regime.
            let tau_cc = mats_const; // matrix must be const-coeff to use SparseMatVecConstCoeff
            let mtau_cc = mats_const && mtau0_arc.is_some();
            let f_cc = mats_const && f0_arc.is_some();
            let h_cc = mats_const && h0_arc.is_some();

            // Direct tables (tau, m_tau, f, h):
            // Use BaseScalarArc whenever available; otherwise use DenseArc.
            mles.push(StreamingMleEnum::BaseScalarArc {
                evals: tau0_arc.clone(),
                num_vars: nvars,
                square: false,
            });

            let f_arc_ring: Option<Arc<Vec<R>>> = inst.f.as_ring_arc();
            let h_arc_ring: Arc<Vec<R>> = h[i].clone();

            if mtau_cc {
                mles.push(StreamingMleEnum::BaseScalarArc {
                    evals: mtau0_arc.as_ref().unwrap().clone(),
                    num_vars: nvars,
                    square: false,
                });
            } else {
                match &inst.m_tau {
                    crate::rgchk::MonomialVec::Dense(v) => {
                        let mle = StreamingMleEnum::DenseArc {
                            evals: v.clone(),
                            num_vars: nvars,
                        };
                        if cm_lazy > 0 {
                            mles.push(StreamingMleEnum::LazyFixed {
                                inner: Box::new(mle),
                                num_vars: nvars,
                                fixed: Vec::new(),
                                weights: vec![R::BaseRing::ONE],
                                max_lazy: cm_lazy,
                            });
                        } else {
                            mles.push(mle);
                        }
                    }
                    crate::rgchk::MonomialVec::Digits { digits, exp_table } => {
                        let mle = StreamingMleEnum::MonomialDigitsArc {
                            digits: digits.clone(),
                            exp_table: exp_table.clone(),
                            num_vars: nvars,
                        };
                        if cm_lazy > 0 {
                            mles.push(StreamingMleEnum::LazyFixed {
                                inner: Box::new(mle),
                                num_vars: nvars,
                                fixed: Vec::new(),
                                weights: vec![R::BaseRing::ONE],
                                max_lazy: cm_lazy,
                            });
                        } else {
                            mles.push(mle);
                        }
                    }
                }
            }
            if f_cc {
                mles.push(StreamingMleEnum::BaseScalarArc {
                    evals: f0_arc.as_ref().unwrap().clone(),
                    num_vars: nvars,
                    square: false,
                });
            } else {
                mles.push(StreamingMleEnum::DenseArc {
                    evals: f_arc_ring
                        .as_ref()
                        .expect("Ring witness required when f_cc is false")
                        .clone(),
                    num_vars: nvars,
                });
            }
            if h_cc {
                mles.push(StreamingMleEnum::BaseScalarArc {
                    evals: h0_arc.as_ref().unwrap().clone(),
                    num_vars: nvars,
                    square: false,
                });
            } else {
                let mle = StreamingMleEnum::DenseArc {
                    evals: h_arc_ring.clone(),
                    num_vars: nvars,
                };
                if cm_lazy > 0 {
                    mles.push(StreamingMleEnum::LazyFixed {
                        inner: Box::new(mle),
                        num_vars: nvars,
                        fixed: Vec::new(),
                        weights: vec![R::BaseRing::ONE],
                        max_lazy: cm_lazy,
                    });
                } else {
                    mles.push(mle);
                }
            }

            if profile {
                println!(
                    "[LF+ Cm::sumchecker_streaming] const-coeff mat-vec flags (L_idx={}): mats_const={} tau_cc={} mtau_cc={} f_cc={} h_cc={}",
                    i, mats_const, tau_cc, mtau_cc, f_cc, h_cc
                );
            }

            // Only materialize `tau` as a ring vector if we cannot use base-scalar mat-vec for it.
            let tau_ring: Option<Arc<Vec<R>>> = if tau_cc {
                None
            } else {
                // Materialize tau as ring only once for sparse mat-vec evaluation.
                // This is O(n) and can dominate wall time for large n; parallelize the conversion.
                #[cfg(feature = "parallel")]
                let v: Vec<R> = {
                    use rayon::prelude::*;
                    inst.tau.par_iter().copied().map(R::from).collect()
                };
                #[cfg(not(feature = "parallel"))]
                let v: Vec<R> = inst.tau.iter().copied().map(R::from).collect();
                Some(Arc::new(v))
            };

            for m in m_arcs {
                if tau_cc {
                    mles.push(StreamingMleEnum::SparseMatVecConstCoeff {
                        matrix: m.clone(),
                        witness0: tau0_arc.clone(),
                        num_vars: nvars,
                    });
                } else {
                    let tau_ring = tau_ring
                        .as_ref()
                        .expect("tau_ring must exist when tau_cc is false");
                    mles.push(StreamingMleEnum::SparseMatVec {
                        matrix: m.clone(),
                        witness: tau_ring.clone(),
                        num_vars: nvars,
                    });
                }

                if mtau_cc {
                    mles.push(StreamingMleEnum::SparseMatVecConstCoeff {
                        matrix: m.clone(),
                        witness0: mtau0_arc.as_ref().unwrap().clone(),
                        num_vars: nvars,
                    });
                } else {
                    match &inst.m_tau {
                        crate::rgchk::MonomialVec::Dense(v) => {
                            let mle = StreamingMleEnum::SparseMatVec {
                                matrix: m.clone(),
                                witness: v.clone(),
                                num_vars: nvars,
                            };
                            if cm_lazy > 0 {
                                mles.push(StreamingMleEnum::LazyFixed {
                                    inner: Box::new(mle),
                                    num_vars: nvars,
                                    fixed: Vec::new(),
                                    weights: vec![R::BaseRing::ONE],
                                    max_lazy: cm_lazy,
                                });
                            } else {
                                mles.push(mle);
                            }
                        }
                        crate::rgchk::MonomialVec::Digits { digits, exp_table } => {
                            let mle = StreamingMleEnum::SparseMatVecMonomialDigits {
                                matrix: m.clone(),
                                digits: digits.clone(),
                                exp_table: exp_table.clone(),
                                num_vars: nvars,
                            };
                            if cm_lazy > 0 {
                                mles.push(StreamingMleEnum::LazyFixed {
                                    inner: Box::new(mle),
                                    num_vars: nvars,
                                    fixed: Vec::new(),
                                    weights: vec![R::BaseRing::ONE],
                                    max_lazy: cm_lazy,
                                });
                            } else {
                                mles.push(mle);
                            }
                        }
                    }
                }

                if f_cc {
                    mles.push(StreamingMleEnum::SparseMatVecConstCoeff {
                        matrix: m.clone(),
                        witness0: f0_arc.as_ref().unwrap().clone(),
                        num_vars: nvars,
                    });
                } else {
                    mles.push(StreamingMleEnum::SparseMatVec {
                        matrix: m.clone(),
                        witness: f_arc_ring
                            .as_ref()
                            .expect("Ring witness required when f_cc is false")
                            .clone(),
                        num_vars: nvars,
                    });
                }

                if h_cc {
                    mles.push(StreamingMleEnum::SparseMatVecConstCoeff {
                        matrix: m.clone(),
                        witness0: h0_arc.as_ref().unwrap().clone(),
                        num_vars: nvars,
                    });
                } else {
                    let mle = StreamingMleEnum::SparseMatVec {
                        matrix: m.clone(),
                        witness: h_arc_ring.clone(),
                        num_vars: nvars,
                    };
                    if cm_lazy > 0 {
                        mles.push(StreamingMleEnum::LazyFixed {
                            inner: Box::new(mle),
                            num_vars: nvars,
                            fixed: Vec::new(),
                            weights: vec![R::BaseRing::ONE],
                            max_lazy: cm_lazy,
                        });
                    } else {
                        mles.push(mle);
                    }
                }
            }
        }
        if profile {
            println!(
                "[LF+ Cm::sumchecker_streaming] build mles: {:?} (mles={})",
                t_build_mles.elapsed(),
                mles.len()
            );
        }

        mles.push(t0_mle.clone());
        mles.push(t1_mle.clone());

        let Mlen = m_arcs.len();

        // Pre-compute random-combinator powers
        let t_rcps = Instant::now();
        let mut rcps = vec![];
        let mut rcp = R::BaseRing::ONE;
        for _ in 0..L {
            // [tau, m_tau, f, h]
            for _ in 0..4 {
                rcps.push(rcp);
                rcp *= rc;
            }
            for _ in 0..Mlen {
                // M_i * [tau, m_tau, f, h]
                for _ in 0..4 {
                    rcps.push(rcp);
                    rcp *= rc;
                }
            }
        }
        rcps.push(rcp); // t(0)
        rcp *= rc;
        rcps.push(rcp); // t(1)
        if profile {
            println!(
                "[LF+ Cm::sumchecker_streaming] build rc powers: {:?} (len={})",
                t_rcps.elapsed(),
                rcps.len()
            );
        }

        let comb_fn = |vals: &[R]| -> R {
            (0..L)
                .map(|l| {
                    let l_idx = 1 + l * (4 + 4 * Mlen);
                    vals[0] * ( // eq
                    vals[l_idx] * rcps[l_idx - 1]  // tau
                    + vals[l_idx + 1] * rcps[l_idx] // m_tau
                    + vals[l_idx + 2] * rcps[l_idx + 1] // f
                    + vals[l_idx + 3] * rcps[l_idx + 2] // h
                    + (0..Mlen).map(|i| {
                        let idx = l_idx + 4 + i * 4;
                        vals[idx] * rcps[idx - 1] // M_i * tau
                        + vals[idx + 1] * rcps[idx] // M_i * m_tau
                        + vals[idx + 2] * rcps[idx + 1] // M_i * f
                        + vals[idx + 3] * rcps[idx + 2] // M_i * h
                     }).sum::<R>()
                )
            + (vals[l_idx] * vals[vals.len()-2]) * rcps[vals.len() - 3] // t(0)
            + (vals[l_idx] * vals[vals.len()-1]) * rcps[vals.len() - 2] // t(1)
                })
                .sum::<R>()
        };

        let t_sc = Instant::now();
        let (sumcheck_proof, randomness, final_vals) =
            StreamingSumcheck::prove_as_subprotocol(transcript, mles, nvars, 2, comb_fn);
        if profile {
            println!(
                "[LF+ Cm::sumchecker_streaming] streaming sumcheck: {:?}",
                t_sc.elapsed()
            );
        }

        let ro = randomness.into_iter().map(|x| x.into()).collect::<Vec<R>>();

        let t_evals = Instant::now();
        let evals = (0..L)
            .map(|l| {
                let mut e = Vec::with_capacity(1 + Mlen);
                let l_idx = 1 + l * (4 + 4 * Mlen);
                e.push([
                    final_vals[l_idx],
                    final_vals[l_idx + 1],
                    final_vals[l_idx + 2],
                    final_vals[l_idx + 3],
                ]);
                for i in 0..Mlen {
                    let idx = l_idx + 4 + i * 4;
                    e.push([
                        final_vals[idx],
                        final_vals[idx + 1],
                        final_vals[idx + 2],
                        final_vals[idx + 3],
                    ]);
                }
                InstanceEvals(e)
            })
            .collect::<Vec<_>>();
        if profile {
            println!(
                "[LF+ Cm::sumchecker_streaming] build evals structs: {:?}",
                t_evals.elapsed()
            );
        }

        let t_absorb = Instant::now();
        absorb_evaluations(&evals, transcript);
        if profile {
            println!(
                "[LF+ Cm::sumchecker_streaming] absorb evals: {:?}",
                t_absorb.elapsed()
            );
        }

        if profile {
            println!(
                "[LF+ Cm::sumchecker_streaming] sumcheck+evals: {:?} (mles={}, L={}, Mlen={})",
                t_sumcheck.elapsed(),
                final_vals.len(),
                L,
                Mlen
            );
        }

        (sumcheck_proof, evals, ro)
    }

    fn sumchecker_streaming_hmle(
        &self,
        dcom: &Dcom<R>,
        h_mles_full: &[Arc<StreamingMleEnum<R>>],
        t0_mle: &StreamingMleEnum<R>,
        t1_mle: &StreamingMleEnum<R>,
        m_arcs: &[Arc<SparseMatrix<R>>],
        mats_const: bool,
        transcript: &mut impl Transcript<R>,
        profile: bool,
    ) -> (Proof<R>, Vec<InstanceEvals<R>>, Vec<R>) {
        let t_sumcheck = Instant::now();
        let nvars = self.rg.nvars;

        let rc = transcript.get_challenge();
        let L = self.rg.instances.len();
        let Mlen = m_arcs.len();

        let mut mles = Vec::with_capacity(1 + L * (4 + 4 * Mlen) + 2);

        let r0 = dcom.out.r.clone();
        let one_minus_r0 = r0.iter().copied().map(|x| R::BaseRing::ONE - x).collect();
        mles.push(StreamingMleEnum::EqBase {
            scale: R::BaseRing::ONE,
            r: r0,
            one_minus_r: one_minus_r0,
        });

        let t_build_mles = Instant::now();
        let cm_lazy: usize = std::env::var("LF_PLUS_CM_LAZY_FIX")
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(4);

        for (i, inst) in self.rg.instances.iter().enumerate() {
            let tau0_arc: Arc<Vec<R::BaseRing>> = inst.tau.clone();
            let mtau0_arc: Option<Arc<Vec<R::BaseRing>>> = if mats_const {
                inst.m_tau
                    .as_dense_arc()
                    .and_then(|v| try_as_base_scalars::<R>(v.as_ref()).map(Arc::new))
            } else {
                None
            };
            let f0_arc: Option<Arc<Vec<R::BaseRing>>> = if mats_const {
                match &inst.f {
                    WitnessVec::ConstCoeffBase { values: v0, .. } => Some(v0.clone()),
                    WitnessVec::Ring(vr) => try_as_base_scalars::<R>(vr.as_ref()).map(Arc::new),
                }
            } else {
                None
            };

            let tau_cc = mats_const;
            let mtau_cc = mats_const && mtau0_arc.is_some();
            let f_cc = mats_const && f0_arc.is_some();
            let h_cc = false; // h is derived and not const-coeff in SP1 regime

            mles.push(StreamingMleEnum::BaseScalarArc {
                evals: tau0_arc.clone(),
                num_vars: nvars,
                square: false,
            });

            let f_arc_ring: Option<Arc<Vec<R>>> = inst.f.as_ring_arc();

            if mtau_cc {
                mles.push(StreamingMleEnum::BaseScalarArc {
                    evals: mtau0_arc.as_ref().unwrap().clone(),
                    num_vars: nvars,
                    square: false,
                });
            } else {
                match &inst.m_tau {
                    crate::rgchk::MonomialVec::Dense(v) => {
                        let mle = StreamingMleEnum::DenseArc {
                            evals: v.clone(),
                            num_vars: nvars,
                        };
                        if cm_lazy > 0 {
                            mles.push(StreamingMleEnum::LazyFixed {
                                inner: Box::new(mle),
                                num_vars: nvars,
                                fixed: Vec::new(),
                                weights: vec![R::BaseRing::ONE],
                                max_lazy: cm_lazy,
                            });
                        } else {
                            mles.push(mle);
                        }
                    }
                    crate::rgchk::MonomialVec::Digits { digits, exp_table } => {
                        let mle = StreamingMleEnum::MonomialDigitsArc {
                            digits: digits.clone(),
                            exp_table: exp_table.clone(),
                            num_vars: nvars,
                        };
                        if cm_lazy > 0 {
                            mles.push(StreamingMleEnum::LazyFixed {
                                inner: Box::new(mle),
                                num_vars: nvars,
                                fixed: Vec::new(),
                                weights: vec![R::BaseRing::ONE],
                                max_lazy: cm_lazy,
                            });
                        } else {
                            mles.push(mle);
                        }
                    }
                }
            }

            if f_cc {
                mles.push(StreamingMleEnum::BaseScalarArc {
                    evals: f0_arc.as_ref().unwrap().clone(),
                    num_vars: nvars,
                    square: false,
                });
            } else {
                mles.push(StreamingMleEnum::DenseArc {
                    evals: f_arc_ring
                        .as_ref()
                        .expect("Ring witness required when f_cc is false")
                        .clone(),
                    num_vars: nvars,
                });
            }

            // h (on-demand)
            let h_mle = (*h_mles_full[i]).clone();
            if cm_lazy > 0 {
                mles.push(StreamingMleEnum::LazyFixed {
                    inner: Box::new(h_mle),
                    num_vars: nvars,
                    fixed: Vec::new(),
                    weights: vec![R::BaseRing::ONE],
                    max_lazy: cm_lazy,
                });
            } else {
                mles.push(h_mle);
            }

            if profile {
                println!(
                    "[LF+ Cm::sumchecker_streaming] const-coeff mat-vec flags (L_idx={}): mats_const={} tau_cc={} mtau_cc={} f_cc={} h_cc={}",
                    i, mats_const, tau_cc, mtau_cc, f_cc, h_cc
                );
            }

            // Immutable witness MLE for M*h
            let h_witness_arc = h_mles_full[i].clone();

            for m in m_arcs {
                if tau_cc {
                    mles.push(StreamingMleEnum::SparseMatVecConstCoeff {
                        matrix: m.clone(),
                        witness0: tau0_arc.clone(),
                        num_vars: nvars,
                    });
                } else {
                    // unreachable in this path
                    panic!("hmle path expects const-coeff matrices for tau");
                }

                if mtau_cc {
                    mles.push(StreamingMleEnum::SparseMatVecConstCoeff {
                        matrix: m.clone(),
                        witness0: mtau0_arc.as_ref().unwrap().clone(),
                        num_vars: nvars,
                    });
                } else {
                    match &inst.m_tau {
                        crate::rgchk::MonomialVec::Dense(v) => {
                            let mle = StreamingMleEnum::SparseMatVec {
                                matrix: m.clone(),
                                witness: v.clone(),
                                num_vars: nvars,
                            };
                            if cm_lazy > 0 {
                                mles.push(StreamingMleEnum::LazyFixed {
                                    inner: Box::new(mle),
                                    num_vars: nvars,
                                    fixed: Vec::new(),
                                    weights: vec![R::BaseRing::ONE],
                                    max_lazy: cm_lazy,
                                });
                            } else {
                                mles.push(mle);
                            }
                        }
                        crate::rgchk::MonomialVec::Digits { digits, exp_table } => {
                            let mle = StreamingMleEnum::SparseMatVecMonomialDigits {
                                matrix: m.clone(),
                                digits: digits.clone(),
                                exp_table: exp_table.clone(),
                                num_vars: nvars,
                            };
                            if cm_lazy > 0 {
                                mles.push(StreamingMleEnum::LazyFixed {
                                    inner: Box::new(mle),
                                    num_vars: nvars,
                                    fixed: Vec::new(),
                                    weights: vec![R::BaseRing::ONE],
                                    max_lazy: cm_lazy,
                                });
                            } else {
                                mles.push(mle);
                            }
                        }
                    }
                }

                if f_cc {
                    mles.push(StreamingMleEnum::SparseMatVecConstCoeff {
                        matrix: m.clone(),
                        witness0: f0_arc.as_ref().unwrap().clone(),
                        num_vars: nvars,
                    });
                } else {
                    mles.push(StreamingMleEnum::SparseMatVec {
                        matrix: m.clone(),
                        witness: f_arc_ring
                            .as_ref()
                            .expect("Ring witness required when f_cc is false")
                            .clone(),
                        num_vars: nvars,
                    });
                }

                // M * h from on-demand witness
                let mle = StreamingMleEnum::SparseMatVecFromMle {
                    matrix: m.clone(),
                    witness_mle: h_witness_arc.clone(),
                    num_vars: nvars,
                };
                if cm_lazy > 0 {
                    mles.push(StreamingMleEnum::LazyFixed {
                        inner: Box::new(mle),
                        num_vars: nvars,
                        fixed: Vec::new(),
                        weights: vec![R::BaseRing::ONE],
                        max_lazy: cm_lazy,
                    });
                } else {
                    mles.push(mle);
                }
            }
        }

        if profile {
            println!(
                "[LF+ Cm::sumchecker_streaming] build mles: {:?} (mles={})",
                t_build_mles.elapsed(),
                mles.len()
            );
        }

        mles.push(t0_mle.clone());
        mles.push(t1_mle.clone());

        // Pre-compute random-combinator powers (same scheme as sumchecker_streaming)
        let t_rcps = Instant::now();
        let mut rcps = vec![];
        let mut rcp = R::BaseRing::ONE;
        for _ in 0..L {
            for _ in 0..4 {
                rcps.push(rcp);
                rcp *= rc;
            }
            for _ in 0..Mlen {
                for _ in 0..4 {
                    rcps.push(rcp);
                    rcp *= rc;
                }
            }
        }
        rcps.push(rcp); // t(0)
        rcp *= rc;
        rcps.push(rcp); // t(1)
        if profile {
            println!(
                "[LF+ Cm::sumchecker_streaming] build rc powers: {:?} (len={})",
                t_rcps.elapsed(),
                rcps.len()
            );
        }

        // Must match `sumchecker_streaming` combiner
        let comb_fn = |vals: &[R]| -> R {
            (0..L)
                .map(|l| {
                    let l_idx = 1 + l * (4 + 4 * Mlen);
                    vals[0] * ( // eq
                        vals[l_idx] * rcps[l_idx - 1]  // tau
                        + vals[l_idx + 1] * rcps[l_idx] // m_tau
                        + vals[l_idx + 2] * rcps[l_idx + 1] // f
                        + vals[l_idx + 3] * rcps[l_idx + 2] // h
                        + (0..Mlen).map(|i| {
                            let idx = l_idx + 4 + i * 4;
                            vals[idx] * rcps[idx - 1] // M_i * tau
                            + vals[idx + 1] * rcps[idx] // M_i * m_tau
                            + vals[idx + 2] * rcps[idx + 1] // M_i * f
                            + vals[idx + 3] * rcps[idx + 2] // M_i * h
                         }).sum::<R>()
                    )
                    + (vals[l_idx] * vals[vals.len()-2]) * rcps[vals.len() - 3] // t(0)
                    + (vals[l_idx] * vals[vals.len()-1]) * rcps[vals.len() - 2] // t(1)
                })
                .sum::<R>()
        };

        let t_sc = Instant::now();
        let (sumcheck_proof, randomness, final_vals) =
            StreamingSumcheck::prove_as_subprotocol(transcript, mles, nvars, 2, comb_fn);
        if profile {
            println!(
                "[LF+ Cm::sumchecker_streaming] streaming sumcheck: {:?}",
                t_sc.elapsed()
            );
        }

        let ro = randomness.into_iter().map(|x| x.into()).collect::<Vec<R>>();

        let t_evals = Instant::now();
        let evals = (0..L)
            .map(|l| {
                let mut e = Vec::with_capacity(1 + Mlen);
                let l_idx = 1 + l * (4 + 4 * Mlen);
                e.push([
                    final_vals[l_idx],
                    final_vals[l_idx + 1],
                    final_vals[l_idx + 2],
                    final_vals[l_idx + 3],
                ]);
                for i in 0..Mlen {
                    let idx = l_idx + 4 + i * 4;
                    e.push([
                        final_vals[idx],
                        final_vals[idx + 1],
                        final_vals[idx + 2],
                        final_vals[idx + 3],
                    ]);
                }
                InstanceEvals(e)
            })
            .collect::<Vec<_>>();
        if profile {
            println!(
                "[LF+ Cm::sumchecker_streaming] build evals structs: {:?}",
                t_evals.elapsed()
            );
        }

        let t_absorb = Instant::now();
        absorb_evaluations(&evals, transcript);
        if profile {
            println!(
                "[LF+ Cm::sumchecker_streaming] absorb evals: {:?}",
                t_absorb.elapsed()
            );
        }

        if profile {
            println!(
                "[LF+ Cm::sumchecker_streaming] sumcheck+evals: {:?} (mles={}, L={}, Mlen={})",
                t_sumcheck.elapsed(),
                final_vals.len(),
                L,
                Mlen
            );
        }

        (sumcheck_proof, evals, ro)
    }
    fn sumchecker_streaming_base(
        &self,
        dcom: &Dcom<R>,
        h: &[Arc<Vec<R>>],
        t0_mle: &StreamingMleEnum<R>,
        t1_mle: &StreamingMleEnum<R>,
        m_arcs0: &[Arc<SparseMatrix<R::BaseRing>>],
        mats_const: bool,
        transcript: &mut impl Transcript<R>,
        profile: bool,
    ) -> (Proof<R>, Vec<InstanceEvals<R>>, Vec<R>) {
        let t_sumcheck = Instant::now();
        let nvars = self.rg.nvars;

        let rc = transcript.get_challenge();
        let L = self.rg.instances.len();

        let mut mles = Vec::with_capacity(
            1 + L * (4 + 4 * m_arcs0.len()) + 2,
        );

        let r0 = dcom.out.r.clone();
        let one_minus_r0 = r0.iter().copied().map(|x| R::BaseRing::ONE - x).collect();
        mles.push(StreamingMleEnum::EqBase {
            scale: R::BaseRing::ONE,
            r: r0,
            one_minus_r: one_minus_r0,
        });

        let t_build_mles = Instant::now();
        let cm_lazy: usize = std::env::var("LF_PLUS_CM_LAZY_FIX")
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(4);

        for (i, inst) in self.rg.instances.iter().enumerate() {
            let tau0_arc: Arc<Vec<R::BaseRing>> = inst.tau.clone();
            let mtau0_arc: Option<Arc<Vec<R::BaseRing>>> = if mats_const {
                inst.m_tau
                    .as_dense_arc()
                    .and_then(|v| try_as_base_scalars::<R>(v.as_ref()).map(Arc::new))
            } else {
                None
            };
            let f0_arc: Option<Arc<Vec<R::BaseRing>>> = if mats_const {
                match &inst.f {
                    WitnessVec::ConstCoeffBase { values: v0, .. } => Some(v0.clone()),
                    WitnessVec::Ring(vr) => try_as_base_scalars::<R>(vr.as_ref()).map(Arc::new),
                }
            } else {
                None
            };
            let h0_arc: Option<Arc<Vec<R::BaseRing>>> =
                if mats_const { try_as_base_scalars::<R>(h[i].as_ref()).map(Arc::new) } else { None };

            let tau_cc = mats_const;
            let mtau_cc = mats_const && mtau0_arc.is_some();
            let f_cc = mats_const && f0_arc.is_some();
            let h_cc = mats_const && h0_arc.is_some();

            mles.push(StreamingMleEnum::BaseScalarArc {
                evals: tau0_arc.clone(),
                num_vars: nvars,
                square: false,
            });

            let f_arc_ring: Option<Arc<Vec<R>>> = inst.f.as_ring_arc();
            let h_arc_ring: Arc<Vec<R>> = h[i].clone();

            if mtau_cc {
                mles.push(StreamingMleEnum::BaseScalarArc {
                    evals: mtau0_arc.as_ref().unwrap().clone(),
                    num_vars: nvars,
                    square: false,
                });
            } else {
                match &inst.m_tau {
                    crate::rgchk::MonomialVec::Dense(v) => {
                        let mle = StreamingMleEnum::DenseArc {
                            evals: v.clone(),
                            num_vars: nvars,
                        };
                        if cm_lazy > 0 {
                            mles.push(StreamingMleEnum::LazyFixed {
                                inner: Box::new(mle),
                                num_vars: nvars,
                                fixed: Vec::new(),
                                weights: vec![R::BaseRing::ONE],
                                max_lazy: cm_lazy,
                            });
                        } else {
                            mles.push(mle);
                        }
                    }
                    crate::rgchk::MonomialVec::Digits { digits, exp_table } => {
                        let mle = StreamingMleEnum::MonomialDigitsArc {
                            digits: digits.clone(),
                            exp_table: exp_table.clone(),
                            num_vars: nvars,
                        };
                        if cm_lazy > 0 {
                            mles.push(StreamingMleEnum::LazyFixed {
                                inner: Box::new(mle),
                                num_vars: nvars,
                                fixed: Vec::new(),
                                weights: vec![R::BaseRing::ONE],
                                max_lazy: cm_lazy,
                            });
                        } else {
                            mles.push(mle);
                        }
                    }
                }
            }

            if f_cc {
                mles.push(StreamingMleEnum::BaseScalarArc {
                    evals: f0_arc.as_ref().unwrap().clone(),
                    num_vars: nvars,
                    square: false,
                });
            } else {
                mles.push(StreamingMleEnum::DenseArc {
                    evals: f_arc_ring
                        .as_ref()
                        .expect("Ring witness required when f_cc is false")
                        .clone(),
                    num_vars: nvars,
                });
            }

            if h_cc {
                mles.push(StreamingMleEnum::BaseScalarArc {
                    evals: h0_arc.as_ref().unwrap().clone(),
                    num_vars: nvars,
                    square: false,
                });
            } else {
                let mle = StreamingMleEnum::DenseArc {
                    evals: h_arc_ring.clone(),
                    num_vars: nvars,
                };
                if cm_lazy > 0 {
                    mles.push(StreamingMleEnum::LazyFixed {
                        inner: Box::new(mle),
                        num_vars: nvars,
                        fixed: Vec::new(),
                        weights: vec![R::BaseRing::ONE],
                        max_lazy: cm_lazy,
                    });
                } else {
                    mles.push(mle);
                }
            }

            if profile {
                println!(
                    "[LF+ Cm::sumchecker_streaming] const-coeff mat-vec flags (L_idx={}): mats_const={} tau_cc={} mtau_cc={} f_cc={} h_cc={}",
                    i, mats_const, tau_cc, mtau_cc, f_cc, h_cc
                );
            }

            for m0 in m_arcs0 {
                // M * tau : base mat + base witness => base-scalar output
                mles.push(StreamingMleEnum::SparseMatVecConstCoeffBase {
                    matrix: m0.clone(),
                    witness0: tau0_arc.clone(),
                    num_vars: nvars,
                });

                // M * m_tau
                if mtau_cc {
                    mles.push(StreamingMleEnum::SparseMatVecConstCoeffBase {
                        matrix: m0.clone(),
                        witness0: mtau0_arc.as_ref().unwrap().clone(),
                        num_vars: nvars,
                    });
                } else {
                    match &inst.m_tau {
                        crate::rgchk::MonomialVec::Dense(v) => {
                            let mle = StreamingMleEnum::SparseMatVecBaseCoeffRing {
                                matrix0: m0.clone(),
                                witness: v.clone(),
                                num_vars: nvars,
                            };
                            if cm_lazy > 0 {
                                mles.push(StreamingMleEnum::LazyFixed {
                                    inner: Box::new(mle),
                                    num_vars: nvars,
                                    fixed: Vec::new(),
                                    weights: vec![R::BaseRing::ONE],
                                    max_lazy: cm_lazy,
                                });
                            } else {
                                mles.push(mle);
                            }
                        }
                        crate::rgchk::MonomialVec::Digits { digits, exp_table } => {
                            let mle = StreamingMleEnum::SparseMatVecBaseCoeffMonomialDigits {
                                matrix0: m0.clone(),
                                digits: digits.clone(),
                                exp_table: exp_table.clone(),
                                num_vars: nvars,
                            };
                            if cm_lazy > 0 {
                                mles.push(StreamingMleEnum::LazyFixed {
                                    inner: Box::new(mle),
                                    num_vars: nvars,
                                    fixed: Vec::new(),
                                    weights: vec![R::BaseRing::ONE],
                                    max_lazy: cm_lazy,
                                });
                            } else {
                                mles.push(mle);
                            }
                        }
                    }
                }

                // M * f
                if f_cc {
                    mles.push(StreamingMleEnum::SparseMatVecConstCoeffBase {
                        matrix: m0.clone(),
                        witness0: f0_arc.as_ref().unwrap().clone(),
                        num_vars: nvars,
                    });
                } else {
                    let mle = StreamingMleEnum::SparseMatVecBaseCoeffRing {
                        matrix0: m0.clone(),
                        witness: f_arc_ring
                            .as_ref()
                            .expect("Ring witness required when f_cc is false")
                            .clone(),
                        num_vars: nvars,
                    };
                    mles.push(mle);
                }

                // M * h
                if h_cc {
                    mles.push(StreamingMleEnum::SparseMatVecConstCoeffBase {
                        matrix: m0.clone(),
                        witness0: h0_arc.as_ref().unwrap().clone(),
                        num_vars: nvars,
                    });
                } else {
                    let mle = StreamingMleEnum::SparseMatVecBaseCoeffRing {
                        matrix0: m0.clone(),
                        witness: h_arc_ring.clone(),
                        num_vars: nvars,
                    };
                    if cm_lazy > 0 {
                        mles.push(StreamingMleEnum::LazyFixed {
                            inner: Box::new(mle),
                            num_vars: nvars,
                            fixed: Vec::new(),
                            weights: vec![R::BaseRing::ONE],
                            max_lazy: cm_lazy,
                        });
                    } else {
                        mles.push(mle);
                    }
                }
            }
        }

        if profile {
            println!(
                "[LF+ Cm::sumchecker_streaming] build mles: {:?} (mles={})",
                t_build_mles.elapsed(),
                mles.len()
            );
        }

        mles.push(t0_mle.clone());
        mles.push(t1_mle.clone());

        let Mlen = m_arcs0.len();

        let t_rcps = Instant::now();
        let mut rcps = vec![];
        let mut rcp = R::BaseRing::ONE;
        for _ in 0..L {
            for _ in 0..4 {
                rcps.push(rcp);
                rcp *= rc;
            }
            for _ in 0..Mlen {
                for _ in 0..4 {
                    rcps.push(rcp);
                    rcp *= rc;
                }
            }
        }
        rcps.push(rcp);
        rcp *= rc;
        rcps.push(rcp);
        if profile {
            println!(
                "[LF+ Cm::sumchecker_streaming] build rc powers: {:?} (len={})",
                t_rcps.elapsed(),
                rcps.len()
            );
        }

        // MUST match `sumchecker_streaming` combiner exactly (including the `tau * t(z)` terms),
        // otherwise the verifier's sumcheck will fail.
        let comb_fn = |vals: &[R]| -> R {
            (0..L)
                .map(|l| {
                    let l_idx = 1 + l * (4 + 4 * Mlen);
                    vals[0] * ( // eq
                        vals[l_idx] * rcps[l_idx - 1]  // tau
                        + vals[l_idx + 1] * rcps[l_idx] // m_tau
                        + vals[l_idx + 2] * rcps[l_idx + 1] // f
                        + vals[l_idx + 3] * rcps[l_idx + 2] // h
                        + (0..Mlen).map(|i| {
                            let idx = l_idx + 4 + i * 4;
                            vals[idx] * rcps[idx - 1] // M_i * tau
                            + vals[idx + 1] * rcps[idx] // M_i * m_tau
                            + vals[idx + 2] * rcps[idx + 1] // M_i * f
                            + vals[idx + 3] * rcps[idx + 2] // M_i * h
                         }).sum::<R>()
                    )
                    + (vals[l_idx] * vals[vals.len()-2]) * rcps[vals.len() - 3] // t(0)
                    + (vals[l_idx] * vals[vals.len()-1]) * rcps[vals.len() - 2] // t(1)
                })
                .sum::<R>()
        };

        let (sumcheck_proof, randomness, final_vals) =
            StreamingSumcheck::prove_as_subprotocol(transcript, mles, nvars, 2, comb_fn);

        let ro = randomness.into_iter().map(|x| x.into()).collect::<Vec<R>>();

        let t_evals = Instant::now();
        let evals = (0..L)
            .map(|l| {
                let mut e = Vec::with_capacity(1 + Mlen);
                let l_idx = 1 + l * (4 + 4 * Mlen);
                e.push([
                    final_vals[l_idx],
                    final_vals[l_idx + 1],
                    final_vals[l_idx + 2],
                    final_vals[l_idx + 3],
                ]);
                for i in 0..Mlen {
                    let idx = l_idx + 4 + i * 4;
                    e.push([
                        final_vals[idx],
                        final_vals[idx + 1],
                        final_vals[idx + 2],
                        final_vals[idx + 3],
                    ]);
                }
                InstanceEvals(e)
            })
            .collect::<Vec<_>>();
        if profile {
            println!(
                "[LF+ Cm::sumchecker_streaming] build evals structs: {:?}",
                t_evals.elapsed()
            );
        }

        let t_absorb = Instant::now();
        absorb_evaluations(&evals, transcript);
        if profile {
            println!(
                "[LF+ Cm::sumchecker_streaming] absorb evals: {:?}",
                t_absorb.elapsed()
            );
        }

        if profile {
            println!(
                "[LF+ Cm::sumchecker_streaming] sumcheck+evals: {:?} (mles={}, L={}, Mlen={})",
                t_sumcheck.elapsed(),
                final_vals.len(),
                L,
                Mlen
            );
        }

        (sumcheck_proof, evals, ro)
    }
}

impl<R: CoeffRing> CmProof<R>
where
    R::BaseRing: Zq,
{
    #[inline]
    pub fn verify_with_mlen(
        &self,
        mlen: usize,
        transcript: &mut impl Transcript<R>,
    ) -> Result<ComX<R>, SumCheckError<R>> {
        let k = self.dcom.dparams.k;
        let d = R::dimension();
        let nvars = self.dcom.out.nvars;
        let L = self.evals.0.len();

        self.dcom.verify(transcript).unwrap();

        let s = (0..3)
            .map(|_| short_challenge(128, transcript))
            .collect::<Vec<R>>();

        let s_prime = (0..k)
            .map(|_| {
                (0..d)
                    .map(|_| short_challenge(128, transcript))
                    .collect::<Vec<R>>()
            })
            .collect::<Vec<_>>();
        let s_prime_flat = s_prime.clone().into_iter().flatten().collect::<Vec<R>>();

        absorb_comh(&self.comh, transcript);

        let kappa = self.comh[0].len();
        let log_kappa = log2(kappa) as usize;

        let c = (0..2)
            .map(|_| {
                transcript
                    .get_challenges(log_kappa)
                    .into_iter()
                    .map(|x| x.into())
                    .collect::<Vec<R>>()
            })
            .collect::<Vec<_>>();

        let u: Vec<Vec<R>> = (0..L)
            .map(|l| {
                self.dcom
                    .out
                    .e
                    .iter()
                    .map(|e_i| {
                        e_i.iter()
                            .skip(l * k)
                            .take(k)
                            .flatten()
                            .zip(s_prime_flat.iter())
                            .map(|(u_ij, s_ij)| *u_ij * *s_ij)
                            .sum()
                    })
                    .collect::<Vec<R>>()
            })
            .collect();

        let tensor_c0 = tensor(&c[0]);
        let tensor_c1 = tensor(&c[1]);
        let tcch0 = self
            .comh
            .iter()
            .map(|com| {
                tensor_c0
                    .iter()
                    .zip(com)
                    .map(|(&t_i, ch_i)| t_i * ch_i)
                    .sum::<R>()
            })
            .collect::<Vec<R>>();
        let tcch1 = self
            .comh
            .iter()
            .map(|com| {
                tensor_c1
                    .iter()
                    .zip(com)
                    .map(|(&t_i, ch_i)| t_i * ch_i)
                    .sum::<R>()
            })
            .collect::<Vec<R>>();

        let dp = R::dimension() / 2;
        let l = self.dcom.dparams.l;
        let dpp = (0..l)
            .map(|i| R::from(R::BaseRing::from(dp as u128).pow([i as u64])))
            .collect::<Vec<_>>();
        let xp = (0..d).map(|i| unit_monomial::<R>(i)).collect::<Vec<_>>();

        let mut verify_sumcheck =
            |sumcheck_proof: &Proof<R>, evals: &[InstanceEvals<R>]| -> Result<Vec<R>, ()> {
                let rc: R = transcript.get_challenge().into();

                let z_idx = L * (4 + 4 * mlen);

                let claimed_sum = self
                    .dcom
                    .evals
                    .iter()
                    .enumerate()
                    .map(|(l, eval)| {
                        let l_idx = l * (4 + 4 * mlen);

                        R::from(eval.a[0]) * rc.pow([l_idx as u64])
                            + eval.b[0] * rc.pow([l_idx as u64 + 1])
                            + eval.c[0] * rc.pow([l_idx as u64 + 2])
                            + u[l][0] * rc.pow([l_idx as u64 + 3])
                            + (0..mlen)
                                .map(|i| {
                                    let idx = l_idx + 4 + i * 4;
                                    R::from(eval.a[1 + i]) * rc.pow([idx as u64])
                                        + eval.b[1 + i] * rc.pow([idx as u64 + 1])
                                        + eval.c[1 + i] * rc.pow([idx as u64 + 2])
                                        + u[l][1 + i] * rc.pow([idx as u64 + 3])
                                })
                                .sum::<R>()
                            + tcch0[l] * rc.pow([z_idx as u64])
                            + tcch1[l] * rc.pow([z_idx as u64 + 1])
                    })
                    .sum::<R>();

                let subclaim = MLSumcheck::verify_as_subprotocol(
                    transcript,
                    nvars,
                    2,
                    claimed_sum,
                    sumcheck_proof,
                )
                .unwrap();

                let r: Vec<R> = self.dcom.out.r.iter().map(|x| R::from(*x)).collect();
                let ro: Vec<R> = subclaim.point.into_iter().map(|x| x.into()).collect();

                // OPTIMIZED: Use tensor structure for O(small) evaluation instead of O(n)
                // The tensor product t(z) = tensor(c_z) ⊗ s' ⊗ d_powers ⊗ x_powers
                // can be evaluated factor-by-factor in O(κ + k*d + ℓ + d) time.
                use crate::tensor_eval::eval_t_z_optimized;
                let t0_ro = eval_t_z_optimized(&c[0], &s_prime_flat, &dpp, &xp, &ro);
                let t1_ro = eval_t_z_optimized(&c[1], &s_prime_flat, &dpp, &xp, &ro);

                let expected_eval = subclaim.expected_evaluation;

                absorb_evaluations(evals, transcript);

                let eq = eq_eval(&r, &ro).unwrap();

                let eval = evals
                    .iter()
                    .enumerate()
                    .map(|(l, el)| {
                        let el = &el.0;
                        let l_idx = l * (4 + 4 * mlen);
                        eq * (el[0][0] * rc.pow([l_idx as u64])
                            + el[0][1] * rc.pow([l_idx as u64 + 1])
                            + el[0][2] * rc.pow([l_idx as u64 + 2])
                            + el[0][3] * rc.pow([l_idx as u64 + 3])
                            + (0..mlen)
                                .map(|i| {
                                    // M_i
                                    let M_evals = el[i + 1];
                                    let idx = l_idx + 4 + i * 4;
                                    M_evals[0] * rc.pow([idx as u64])
                                        + M_evals[1] * rc.pow([idx as u64 + 1])
                                        + M_evals[2] * rc.pow([idx as u64 + 2])
                                        + M_evals[3] * rc.pow([idx as u64 + 3])
                                })
                                .sum::<R>())
                            + (t0_ro * el[0][0]) * rc.pow([z_idx as u64])
                            + (t1_ro * el[0][0]) * rc.pow([z_idx as u64 + 1])
                    })
                    .sum::<R>();

                assert_eq!(expected_eval, eval);

                Ok(ro)
            };

        let ro0 = verify_sumcheck(&self.sumcheck_proofs.0, &self.evals.0).unwrap();
        let ro1 = verify_sumcheck(&self.sumcheck_proofs.1, &self.evals.1).unwrap();

        let ro = ro0.into_iter().zip(ro1).collect::<Vec<_>>();

        // Step 6
        Ok(self.x(&s, ro))
    }

    pub fn verify(
        &self,
        M: &[Arc<SparseMatrix<R>>],
        transcript: &mut impl Transcript<R>,
    ) -> Result<ComX<R>, SumCheckError<R>> {
        self.verify_with_mlen(M.len(), transcript)
    }

    pub fn x(&self, s: &[R], ro: Vec<(R, R)>) -> ComX<R> {
        let L = self.dcom.fcoms.len();

        // TODO needs more folding challenges `s` for the L instances
        let cm_g = self
            .dcom
            .fcoms
            .iter()
            .enumerate()
            .map(|(l, cmc)| {
                cmc.C_Mf
                    .iter()
                    .zip(&cmc.cm_mtau)
                    .zip(&cmc.cm_f)
                    .zip(&self.comh[l])
                    .map(|(((r_Mf, r_mtau), r_f), r_comh)| {
                        s[0] * r_Mf + s[1] * r_mtau + s[2] * r_f + r_comh
                    })
                    .collect::<Vec<R>>()
            })
            .collect::<Vec<_>>();

        let vo = (0..L)
            .map(|l| {
                let e0l = &self.evals.0[l].0;
                let e1l = &self.evals.1[l].0;
                e0l.iter()
                    .zip(e1l.iter())
                    .map(|(e0li, e1li)| {
                        (
                            (s[0] * e0li[0]) + (s[1] * e0li[1]) + (s[2] * e0li[2]) + e0li[3],
                            (s[0] * e1li[0]) + (s[1] * e1li[1]) + (s[2] * e1li[2]) + e1li[3],
                        )
                    })
                    .collect::<Vec<(R, R)>>()
            })
            .collect::<Vec<Vec<_>>>();

        ComX { cm_g, ro, vo }
    }
}

fn absorb_comh<R: OverField>(comh: &[Vec<R>], transcript: &mut impl Transcript<R>) {
    comh.iter().for_each(|ci| transcript.absorb_slice(ci));
}

fn absorb_evaluations<R: OverField>(
    evals: &[InstanceEvals<R>],
    transcript: &mut impl Transcript<R>,
) {
    evals.iter().for_each(|ieval| {
        ieval.0.iter().for_each(|vals| {
            transcript.absorb_slice(vals);
        });
    });
}

/// t(z) = tensor(c(z)) ⊗ s' ⊗ (1, d', ..., d'^(ℓ-1)) ⊗ (1, X, ..., X^(d-1))
#[allow(dead_code)]
// Dense reference implementation (debugging / cross-checking).
// Hot paths use streaming `Tensor4Padded` (prover) and `eval_t_z_optimized` (verifier).
fn calculate_t_z<T>(c_z: &[T], s_prime: &[T], d_prime_powers: &[T], x_powers: &[T]) -> Vec<T>
where
    T: Clone + One + Sub<Output = T> + Mul<Output = T>,
{
    let tensor_c_z = tensor(c_z);
    let part1 = tensor_product(&tensor_c_z, s_prime);
    let part2 = tensor_product(&part1, d_prime_powers);
    tensor_product(&part2, x_powers)
}

#[cfg(test)]
mod tests {
    use ark_ff::PrimeField;
    use ark_std::Zero;
    use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
    use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as R;
    use stark_rings_linalg::{Matrix, SparseMatrix};
    use std::sync::Arc;

    use super::*;
    use crate::{
        rgchk::{DecompParameters, RgInstance},
        transcript::PoseidonTranscript,
    };

    #[test]
    fn test_com() {
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
        let M: Vec<Arc<SparseMatrix<R>>> = vec![Arc::new(m)];

        let kappa = 2;
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
            dparams: DecompParameters { b, k, l },
        };

        let cm = Cm { rg };

        let mut ts = PoseidonTranscript::empty::<PC>();
        let (_com, proof) = cm.prove(&M, &mut ts);

        let mut ts = PoseidonTranscript::empty::<PC>();
        proof.verify(&M, &mut ts).unwrap();
    }
}
