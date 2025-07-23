use ark_std::{
    iter::once,
    log2,
    ops::{Mul, Sub},
    One, Zero,
};
use latticefold::{
    transcript::Transcript,
    utils::sumcheck::{
        utils::{build_eq_x_r, eq_eval},
        MLSumcheck, Proof, SumCheckError,
    },
};
use stark_rings::{
    balanced_decomposition::{
        convertible_ring::ConvertibleRing, Decompose, DecomposeToVec, GadgetDecompose,
    },
    exp, psi, psi_range_check, unit_monomial, CoeffRing, OverField, PolyRing, Ring, Zq,
};
use stark_rings_linalg::{ops::Transpose, Matrix, SparseMatrix};
use stark_rings_poly::mle::{DenseMultilinearExtension, SparseMultilinearExtension};
use thiserror::Error;

use crate::{
    rgchk::{Dcom, DecompParameters, Rg},
    setchk::{In, MonomialSet, Out},
    utils::{tensor, tensor_product},
};

#[derive(Clone, Debug)]
pub struct Cm<R: PolyRing> {
    pub rg: Rg<R>,
    pub coms: Vec<CmComs<R>>,
}

#[derive(Clone, Debug)]
pub struct CmComs<R> {
    pub cm_f: Vec<R>,
    pub C_Mf: Vec<R>,
    pub cm_mtau: Vec<R>,
}

// eval over r_o of [tau (a), m_tau (b), f (c), h (u)] over 1 + n_lin
#[derive(Clone, Debug)]
pub struct InstanceEvals<R>(Vec<[R; 4]>);

#[derive(Clone, Debug)]
pub struct CmProof<R: PolyRing> {
    pub dcom: Dcom<R>,
    pub comh: Vec<Vec<R>>,
    pub sumcheck_proofs: (Proof<R>, Proof<R>),
    pub evals: (Vec<InstanceEvals<R>>, Vec<InstanceEvals<R>>),

    pub cm_coms: Vec<CmComs<R>>,
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
    pub fn prove(&self, transcript: &mut impl Transcript<R>) -> (Com<R>, CmProof<R>) {
        let k = self.rg.dparams.k;
        let d = R::dimension();
        let dp = R::dimension() / 2;
        let l = self.rg.dparams.l;
        let n = self.rg.instances[0].tau.len();

        let dcom = self.rg.range_check(transcript);

        let s: Vec<R> = transcript
            .get_challenges(3)
            .into_iter()
            .map(|x| x.into())
            .collect();

        let s_prime = (0..k)
            .map(|_| {
                transcript
                    .get_challenges(d)
                    .into_iter()
                    .map(|x| x.into())
                    .collect::<Vec<R>>()
            })
            .collect::<Vec<_>>();
        let s_prime_flat = s_prime.clone().into_iter().flatten().collect::<Vec<R>>();

        let h: Vec<Vec<R>> = self
            .rg
            .instances
            .iter()
            .map(|inst| {
                let n = 1 << self.rg.nvars;
                let h_vectors: Vec<Vec<R>> = inst
                    .M_f
                    .iter()
                    .zip(s_prime.iter())
                    .map(|(M, s_i)| M.try_mul_vec(s_i).unwrap())
                    .collect();

                let mut h = vec![R::zero(); n];
                for v in h_vectors {
                    for (i, val) in v.iter().enumerate() {
                        h[i] += *val;
                    }
                }
                h
            })
            .collect();

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

        let kappa = comh.len();
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

        let dpp = (0..l)
            .map(|i| R::from(R::BaseRing::from(dp as u128).pow([i as u64])))
            .collect::<Vec<_>>();
        let xp = (0..d).map(|i| unit_monomial::<R>(i)).collect::<Vec<_>>();

        let mut t0 = calculate_t_z(&c[0], &s_prime_flat, &dpp, &xp);
        if t0.len() <= n {
            t0.resize(n, R::zero()); // pad
        } else {
            panic!("t0 too large!");
        };

        let mut t1 = calculate_t_z(&c[1], &s_prime_flat, &dpp, &xp);
        if t1.len() <= n {
            t1.resize(n, R::zero()); // pad
        } else {
            panic!("t1 too large!");
        };

        let (proof_a, evals_a, ro_a) =
            self.sumchecker(&dcom, &h, (t0.clone(), t1.clone()), transcript);
        let (proof_b, evals_b, ro_b) = self.sumchecker(&dcom, &h, (t0, t1), transcript);

        // Step 7
        // TODO needs more folding challenges `s` for the L instances
        let g = self
            .rg
            .instances
            .iter()
            .enumerate()
            .map(|(i, inst)| {
                inst.tau
                    .iter()
                    .zip(&inst.m_tau)
                    .zip(&inst.f)
                    .zip(&h[i])
                    .map(|(((r_tau, r_mtau), r_f), r_h)| {
                        (s[0] * R::from(*r_tau)) + (s[1] * r_mtau) + (s[2] * r_f) + r_h
                    })
                    .collect::<Vec<R>>()
            })
            .collect::<Vec<_>>();

        let proof = CmProof {
            dcom,
            comh,
            sumcheck_proofs: (proof_a, proof_b),
            evals: (evals_a, evals_b),
            cm_coms: self.coms.clone(),
        };

        let ro = ro_a
            .into_iter()
            .zip(ro_b.into_iter())
            .map(|rr| rr)
            .collect::<Vec<_>>();

        let x = proof.x(&s, ro);

        let com = Com { g, x };

        (com, proof)
    }

    fn sumchecker(
        &self,
        dcom: &Dcom<R>,
        h: &[Vec<R>],
        t: (Vec<R>, Vec<R>),
        transcript: &mut impl Transcript<R>,
    ) -> (Proof<R>, Vec<InstanceEvals<R>>, Vec<R>) {
        let nvars = self.rg.nvars;
        let r: Vec<R> = dcom.out.r.iter().map(|x| R::from(*x)).collect();

        let rc: R = transcript.get_challenge().into();

        let L = self.rg.instances.len();

        let mut mles = Vec::with_capacity(
            1 // eq
            + L * (
                4  // [tau, m_tau, f, h]
                + 4 * self.rg.M.len() // M * [tau, ...]
            )
            + 2, // t(z)
        );
        let eq = build_eq_x_r(&r).unwrap();
        mles.push(eq);

        for (i, inst) in self.rg.instances.iter().enumerate() {
            let rtau = inst.tau.iter().map(|z| R::from(*z)).collect::<Vec<_>>();
            let tau_mle = DenseMultilinearExtension::from_evaluations_vec(nvars, rtau.clone());
            let m_tau_mle =
                DenseMultilinearExtension::from_evaluations_vec(nvars, inst.m_tau.clone());
            let f_mle = DenseMultilinearExtension::from_evaluations_vec(nvars, inst.f.clone());
            let h_mle = DenseMultilinearExtension::from_evaluations_vec(nvars, h[i].clone());

            mles.push(tau_mle);
            mles.push(m_tau_mle);
            mles.push(f_mle);
            mles.push(h_mle);

            for m in &self.rg.M {
                let Mtau = m.try_mul_vec(&rtau).unwrap();
                mles.push(DenseMultilinearExtension::from_evaluations_vec(nvars, Mtau));

                let Mm_tau = m.try_mul_vec(&inst.m_tau).unwrap();
                mles.push(DenseMultilinearExtension::from_evaluations_vec(
                    nvars, Mm_tau,
                ));

                let Mf = m.try_mul_vec(&inst.f).unwrap();
                mles.push(DenseMultilinearExtension::from_evaluations_vec(nvars, Mf));

                let Mh = m.try_mul_vec(&h[i]).unwrap();
                mles.push(DenseMultilinearExtension::from_evaluations_vec(nvars, Mh));
            }
        }

        let t0_mle = DenseMultilinearExtension::from_evaluations_vec(nvars, t.0.clone());
        let t1_mle = DenseMultilinearExtension::from_evaluations_vec(nvars, t.1.clone());
        mles.push(t0_mle);
        mles.push(t1_mle);

        let Mlen = self.rg.M.len();
        let comb_fn = |vals: &[R]| -> R {
            (0..L)
                .map(|l| {
                    let l_idx = 1 + l * (4 + 4 * Mlen);
                    vals[0] * ( // eq
                    vals[l_idx] * rc.pow([l_idx as u64 - 1])  // tau
                    + vals[l_idx + 1] * rc.pow([l_idx as u64]) // m_tau
                    + vals[l_idx + 2] * rc.pow([l_idx as u64 + 1]) // f
                    + vals[l_idx + 3] * rc.pow([l_idx as u64 + 2]) // h
                    + (0..Mlen).map(|i| {
                        let idx = l_idx + 4 + i * 4;
                        vals[idx] * rc.pow([idx as u64 - 1]) // M_i * tau
                        + vals[idx + 1] * rc.pow([idx as u64]) // M_i * m_tau
                        + vals[idx + 2] * rc.pow([idx as u64 + 1]) // M_i * f
                        + vals[idx + 3] * rc.pow([idx as u64 + 2]) // M_i * h
                     }).sum::<R>()
                )
            + (vals[l_idx] * vals[vals.len()-2]) * rc.pow([vals.len() as u64 - 3]) // t(0)
            + (vals[l_idx] * vals[vals.len()-1]) * rc.pow([vals.len() as u64 - 2])
                    // t(1)
                })
                .sum::<R>()
        };

        let (sumcheck_proof, prover_state) =
            MLSumcheck::prove_as_subprotocol(transcript, mles.clone(), nvars, 2, comb_fn);
        let ro = prover_state
            .randomness
            .into_iter()
            .map(|x| x.into())
            .collect::<Vec<R>>();

        let evals = (0..L)
            .map(|l| {
                let mut e = Vec::with_capacity(1 + Mlen);
                let l_idx = 1 + l * (4 + 4 * Mlen);
                e.push([
                    mles[l_idx].evaluate(&ro).unwrap(),
                    mles[l_idx + 1].evaluate(&ro).unwrap(),
                    mles[l_idx + 2].evaluate(&ro).unwrap(),
                    mles[l_idx + 3].evaluate(&ro).unwrap(),
                ]);
                for i in 0..Mlen {
                    let idx = l_idx + 4 + i * 4;
                    e.push([
                        mles[idx].evaluate(&ro).unwrap(),
                        mles[idx + 1].evaluate(&ro).unwrap(),
                        mles[idx + 2].evaluate(&ro).unwrap(),
                        mles[idx + 3].evaluate(&ro).unwrap(),
                    ]);
                }
                InstanceEvals(e)
            })
            .collect::<Vec<_>>();

        (sumcheck_proof, evals, ro)
    }
}

impl<R: CoeffRing> CmProof<R>
where
    R::BaseRing: Zq,
{
    pub fn verify(&self, transcript: &mut impl Transcript<R>) -> Result<ComX<R>, SumCheckError<R>> {
        let k = self.dcom.dparams.k;
        let d = R::dimension();
        let nvars = self.dcom.out.nvars;
        let M = &self.dcom.out.M;
        let L = self.evals.0.len();

        self.dcom.verify(transcript).unwrap();

        let s: Vec<R> = transcript
            .get_challenges(3)
            .into_iter()
            .map(|x| x.into())
            .collect();

        let s_prime: Vec<Vec<R>> = (0..k)
            .map(|_| {
                transcript
                    .get_challenges(d)
                    .into_iter()
                    .map(|x| x.into())
                    .collect::<Vec<R>>()
            })
            .collect::<Vec<_>>();
        let s_prime_flat = s_prime.clone().into_iter().flatten().collect::<Vec<R>>();

        let kappa = self.comh.len();
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

                let z_idx = L * (4 + 4 * M.len());

                let claimed_sum = self
                    .dcom
                    .evals
                    .iter()
                    .enumerate()
                    .map(|(l, eval)| {
                        let l_idx = l * (4 + 4 * M.len());

                        R::from(eval.a[0]) * rc.pow([l_idx as u64])
                            + eval.b[0] * rc.pow([l_idx as u64 + 1])
                            + eval.c[0] * rc.pow([l_idx as u64 + 2])
                            + u[l][0] * rc.pow([l_idx as u64 + 3])
                            + (0..M.len())
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
                    &sumcheck_proof,
                )
                .unwrap();

                let r: Vec<R> = self.dcom.out.r.iter().map(|x| R::from(*x)).collect();
                let ro: Vec<R> = subclaim.point.into_iter().map(|x| x.into()).collect();
                let t0 = DenseMultilinearExtension::from_evaluations_vec(
                    nvars,
                    calculate_t_z(&c[0], &s_prime_flat, &dpp, &xp),
                );
                let t0_ro = t0.evaluate(&ro).unwrap();
                let t1 = DenseMultilinearExtension::from_evaluations_vec(
                    nvars,
                    calculate_t_z(&c[1], &s_prime_flat, &dpp, &xp),
                );
                let t1_ro = t1.evaluate(&ro).unwrap();

                let expected_eval = subclaim.expected_evaluation;
                let eq = eq_eval(&r, &ro).unwrap();

                let eval = evals
                    .iter()
                    .enumerate()
                    .map(|(l, el)| {
                        let el = &el.0;
                        let l_idx = l * (4 + 4 * M.len());
                        eq * (el[0][0] * rc.pow([l_idx as u64])
                            + el[0][1] * rc.pow([l_idx as u64 + 1])
                            + el[0][2] * rc.pow([l_idx as u64 + 2])
                            + el[0][3] * rc.pow([l_idx as u64 + 3])
                            + (0..M.len())
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

        let ro = ro0
            .into_iter()
            .zip(ro1.into_iter())
            .map(|rr| rr)
            .collect::<Vec<_>>();

        // Step 6
        Ok(self.x(&s, ro))
    }

    pub fn x(&self, s: &[R], ro: Vec<(R, R)>) -> ComX<R> {
        let L = self.cm_coms.len();

        // TODO needs more folding challenges `s` for the L instances
        let cm_g = self
            .cm_coms
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

/// t(z) = tensor(c(z)) ⊗ s' ⊗ (1, d', ..., d'^(ℓ-1)) ⊗ (1, X, ..., X^(d-1))
fn calculate_t_z<T>(c_z: &[T], s_prime: &[T], d_prime_powers: &[T], x_powers: &[T]) -> Vec<T>
where
    T: Clone + One + Sub<Output = T> + Mul<Output = T>,
{
    let tensor_c_z = tensor(c_z);

    let part1 = tensor_product(&tensor_c_z, s_prime);
    let part2 = tensor_product(&part1, d_prime_powers);
    let t_z = tensor_product(&part2, x_powers);

    t_z
}

#[cfg(test)]
mod tests {
    use ark_ff::PrimeField;
    use ark_std::{One, Zero};
    use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
    use latticefold::transcript::poseidon::PoseidonTS;
    use stark_rings::{
        balanced_decomposition::DecomposeToVec, cyclotomic_ring::models::frog_ring::RqPoly as R,
    };
    use stark_rings_linalg::SparseMatrix;

    use super::*;
    use crate::{rgchk::RgInstance, utils::split};

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
        let M = vec![m];

        let kappa = 2;
        let b = (R::dimension() / 2) as u128;
        let k = 2;
        // log_d' (q)
        let l = ((<<R as PolyRing>::BaseRing>::MODULUS.0[0] as f64).ln()
            / ((R::dimension() / 2) as f64).ln())
        .ceil() as usize;

        let A = Matrix::<R>::rand(&mut rand::thread_rng(), kappa, n);

        let dparams = DecompParameters { b, k, l };
        let instance = RgInstance::from_f(f.clone(), &A, &dparams);

        let cm_f = A.try_mul_vec(&instance.f).unwrap();
        let C_Mf = A
            .try_mul_vec(&instance.tau.iter().map(|z| R::from(*z)).collect::<Vec<R>>())
            .unwrap();
        let cm_mtau = A.try_mul_vec(&instance.m_tau).unwrap();

        let rg = Rg {
            nvars: log2(n) as usize,
            instances: vec![instance],
            M,
            dparams: DecompParameters { b, k, l },
        };

        let cm = Cm {
            rg,
            coms: vec![CmComs {
                cm_f,
                C_Mf,
                cm_mtau,
            }],
        };

        let mut ts = PoseidonTS::default::<PC>();
        let (_com, proof) = cm.prove(&mut ts);

        let mut ts = PoseidonTS::default::<PC>();
        proof.verify(&mut ts).unwrap();
    }
}
