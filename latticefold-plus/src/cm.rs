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
    rgchk::{Dcom, Rg},
    setchk::{In, MonomialSet, Out},
    utils::{tensor, tensor_product},
};

pub struct Cm<R: PolyRing> {
    pub rg: Rg<R>,
    pub cm_f: Vec<R>,
    pub C_Mf: Vec<R>,
    pub cm_mtau: Vec<R>,
}

pub struct CmProof<R: PolyRing> {
    pub dcom: Dcom<R>,
    pub comh: Vec<R>,
    pub sumcheck_proofs: (Proof<R>, Proof<R>),
    pub evals: (Vec<[R; 4]>, Vec<[R; 4]>), // eval over r0 of [tau (a), m_tau (b), f (c), h (u)],
    // over 1 + n_lin, for two sumchecks
    pub cm_f: Vec<R>,
    pub C_Mf: Vec<R>,
    pub cm_mtau: Vec<R>,
}

pub struct Com<R: PolyRing> {
    pub ro: Vec<R>,
    pub g: Vec<R>,
}

impl<R: CoeffRing> Cm<R>
where
    R::BaseRing: Zq,
{
    pub fn prove(&self, transcript: &mut impl Transcript<R>) -> CmProof<R> {
        let k = self.rg.comM_f.len();
        let d = R::dimension();
        let dp = R::dimension() / 2;
        let l = self.rg.l;
        let n = self.rg.tau.len();

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

        let h: Vec<R> = {
            let n = 1 << self.rg.nvars;
            let h_vectors: Vec<Vec<R>> = self
                .rg
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
        };

        let comh: Vec<R> = {
            let comh_vectors = self
                .rg
                .comM_f
                .iter()
                .zip(s_prime.iter())
                .map(|(comM_f_i, s_i)| comM_f_i.try_mul_vec(s_i).unwrap())
                .collect::<Vec<_>>();

            let mut comh = vec![R::zero(); self.rg.comM_f[0].nrows];
            for v in comh_vectors {
                for (i, val) in v.iter().enumerate() {
                    comh[i] += *val;
                }
            }
            comh
        };

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

        let (proof_a, evals_a) =
            self.sumchecker(&dcom, &h, &s_prime, (t0.clone(), t1.clone()), transcript);
        let (proof_b, evals_b) = self.sumchecker(&dcom, &h, &s_prime, (t0, t1), transcript);

        // Step 7
        let cm_g = self
            .rg
            .tau
            .iter()
            .zip(&self.rg.m_tau)
            .zip(&self.rg.f)
            .zip(&h)
            .map(|(((r_tau, r_mtau), r_f), r_h)| {
                (s[0] * R::from(*r_tau)) + (s[1] * r_mtau) + (s[2] * r_f) + r_h
            })
            .collect::<Vec<R>>();

        CmProof {
            dcom,
            comh,
            sumcheck_proofs: (proof_a, proof_b),
            evals: (evals_a, evals_b),

            cm_f: self.cm_f.clone(),
            C_Mf: self.C_Mf.clone(),
            cm_mtau: self.cm_mtau.clone(),
        }
    }

    fn sumchecker(
        &self,
        dcom: &Dcom<R>,
        h: &Vec<R>,
        s_prime: &Vec<Vec<R>>,
        t: (Vec<R>, Vec<R>),
        transcript: &mut impl Transcript<R>,
    ) -> (Proof<R>, Vec<[R; 4]>) {
        let nvars = self.rg.nvars;
        let r: Vec<R> = dcom.out.r.iter().map(|x| R::from(*x)).collect();

        let eq = build_eq_x_r(&r).unwrap();
        let rc: R = transcript.get_challenge().into();

        let rtau = self.rg.tau.iter().map(|z| R::from(*z)).collect::<Vec<_>>();
        let tau_mle = DenseMultilinearExtension::from_evaluations_vec(nvars, rtau.clone());
        let m_tau_mle =
            DenseMultilinearExtension::from_evaluations_vec(nvars, self.rg.m_tau.clone());
        let f_mle = DenseMultilinearExtension::from_evaluations_vec(nvars, self.rg.f.clone());
        let h_mle = DenseMultilinearExtension::from_evaluations_vec(nvars, h.clone());
        let t0_mle = DenseMultilinearExtension::from_evaluations_vec(nvars, t.0.clone());
        let t1_mle = DenseMultilinearExtension::from_evaluations_vec(nvars, t.1.clone());

        let mut mles = Vec::with_capacity(7 + 4 * self.rg.M.len());
        mles.push(eq);
        mles.push(tau_mle);
        mles.push(m_tau_mle);
        mles.push(f_mle);
        mles.push(h_mle);
        mles.push(t0_mle);
        mles.push(t1_mle);

        for m in &self.rg.M {
            let Mtau = m.try_mul_vec(&rtau).unwrap();
            mles.push(DenseMultilinearExtension::from_evaluations_vec(nvars, Mtau));

            let Mm_tau = m.try_mul_vec(&self.rg.m_tau).unwrap();
            mles.push(DenseMultilinearExtension::from_evaluations_vec(
                nvars, Mm_tau,
            ));

            let Mf = m.try_mul_vec(&self.rg.f).unwrap();
            mles.push(DenseMultilinearExtension::from_evaluations_vec(nvars, Mf));

            let Mh = m.try_mul_vec(&h).unwrap();
            mles.push(DenseMultilinearExtension::from_evaluations_vec(nvars, Mh));
        }

        let comb_fn = |vals: &[R]| -> R {
            vals[0] // eq
                * (
                    vals[1]  // tau
                    + vals[2] * rc // m_tau
                    + vals[3] * rc.pow([2]) // f
                    + vals[4] * rc.pow([3]) // h
                    + (0..self.rg.M.len()).map(|i| {
                        let idx = 7 + i * 4;
                        vals[idx] * rc.pow([6 + (i * 4) as u64]) // M_i * tau
                        + vals[idx + 1] * rc.pow([6 + (i * 4 + 1) as u64]) // M_i * m_tau
                        + vals[idx + 2] * rc.pow([6 + (i * 4 + 2) as u64]) // M_i * f
                        + vals[idx + 3] * rc.pow([6 + (i * 4 + 3) as u64]) // M_i * h
                     }).sum::<R>()
                )
             + (vals[1] * vals[5]) * rc.pow([4]) // t(0)
             + (vals[1] * vals[6]) * rc.pow([5]) // t(1)
        };

        let (sumcheck_proof, prover_state) =
            MLSumcheck::prove_as_subprotocol(transcript, mles.clone(), nvars, 2, comb_fn);
        let ro = prover_state
            .randomness
            .into_iter()
            .map(|x| x.into())
            .collect::<Vec<R>>();

        let mut evals = Vec::with_capacity(1 + self.rg.M.len());
        evals.push([
            mles[1].evaluate(&ro).unwrap(),
            mles[2].evaluate(&ro).unwrap(),
            mles[3].evaluate(&ro).unwrap(),
            mles[4].evaluate(&ro).unwrap(),
        ]);
        for i in 0..self.rg.M.len() {
            evals.push([
                mles[7 + i * 4].evaluate(&ro).unwrap(),
                mles[7 + 1 + i * 4].evaluate(&ro).unwrap(),
                mles[7 + 2 + i * 4].evaluate(&ro).unwrap(),
                mles[7 + 3 + i * 4].evaluate(&ro).unwrap(),
            ]);
        }

        (sumcheck_proof, evals)
    }

    pub fn verify(
        &self,
        proof: &CmProof<R>,
        transcript: &mut impl Transcript<R>,
    ) -> Result<(), SumCheckError<R>> {
        let k = self.rg.comM_f.len();
        let d = R::dimension();

        proof.dcom.verify(transcript).unwrap();

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

        let kappa = proof.comh.len();
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

        let u: Vec<R> = proof
            .dcom
            .out
            .e
            .iter()
            .map(|e_i| {
                e_i.iter()
                    .flatten()
                    .zip(s_prime.iter().flatten())
                    .map(|(u_ij, s_ij)| *u_ij * *s_ij)
                    .sum()
            })
            .collect();

        let tcch0 = tensor(&c[0])
            .iter()
            .zip(&proof.comh)
            .map(|(&t_i, ch_i)| t_i * ch_i)
            .sum::<R>();
        let tcch1 = tensor(&c[1])
            .iter()
            .zip(&proof.comh)
            .map(|(&t_i, ch_i)| t_i * ch_i)
            .sum::<R>();

        let dp = R::dimension() / 2;
        let l = self.rg.l;
        let dpp = (0..l)
            .map(|i| R::from(R::BaseRing::from(dp as u128).pow([i as u64])))
            .collect::<Vec<_>>();
        let xp = (0..d).map(|i| unit_monomial::<R>(i)).collect::<Vec<_>>();

        let mut verify_sumcheck = |sumcheck_proof: &Proof<R>, evals: &[[R; 4]]| {
            let rc: R = transcript.get_challenge().into();

            let claimed_sum = R::from(proof.dcom.a[0])
                + proof.dcom.b[0] * rc
                + proof.dcom.c[0] * rc.pow([2])
                + u[0] * rc.pow([3])
                + tcch0 * rc.pow([4])
                + tcch1 * rc.pow([5])
                + (0..self.rg.M.len())
                    .map(|i| {
                        R::from(proof.dcom.a[i + 1]) * rc.pow([6 + (i * 4) as u64])
                            + proof.dcom.b[i + 1] * rc.pow([6 + (i * 4 + 1) as u64])
                            + proof.dcom.c[i + 1] * rc.pow([6 + (i * 4 + 2) as u64])
                            + u[i + 1] * rc.pow([6 + (i * 4 + 3) as u64])
                    })
                    .sum::<R>();

            let subclaim = MLSumcheck::verify_as_subprotocol(
                transcript,
                self.rg.nvars,
                2,
                claimed_sum,
                &sumcheck_proof,
            )
            .unwrap();

            let r: Vec<R> = proof.dcom.out.r.iter().map(|x| R::from(*x)).collect();
            let ro: Vec<R> = subclaim.point.into_iter().map(|x| x.into()).collect();
            let t0 = DenseMultilinearExtension::from_evaluations_vec(
                self.rg.nvars,
                calculate_t_z(&c[0], &s_prime_flat, &dpp, &xp),
            );
            let t0_ro = t0.evaluate(&ro).unwrap();
            let t1 = DenseMultilinearExtension::from_evaluations_vec(
                self.rg.nvars,
                calculate_t_z(&c[1], &s_prime_flat, &dpp, &xp),
            );
            let t1_ro = t1.evaluate(&ro).unwrap();

            let expected_eval = subclaim.expected_evaluation;
            let e = eq_eval(&r, &ro).unwrap();

            assert_eq!(
                expected_eval,
                e * (evals[0][0] + rc * evals[0][1] + rc.pow([2]) * evals[0][2] + rc.pow([3]) * evals[0][3] // base
                    + (0..self.rg.M.len()).map(|i| { // M_i
                        let M_evals = evals[i + 1];
                        M_evals[0] * rc.pow([6 + (i * 4) as u64])
                        + M_evals[1] * rc.pow([6 + (i * 4 + 1) as u64])
                        + M_evals[2] * rc.pow([6 + (i * 4 + 2) as u64])
                        + M_evals[3] * rc.pow([6 + (i * 4 + 3) as u64])
                    }).sum::<R>())
                    + rc.pow([4]) * (t0_ro * evals[0][0])
                    + rc.pow([5]) * (t1_ro * evals[0][0])
            );
        };

        verify_sumcheck(&proof.sumcheck_proofs.0, &proof.evals.0);
        verify_sumcheck(&proof.sumcheck_proofs.1, &proof.evals.1);

        // Step 6
        let cm_g = proof
            .C_Mf
            .iter()
            .zip(&proof.cm_mtau)
            .zip(&proof.cm_f)
            .zip(&proof.comh)
            .map(|(((r_Mf, r_mtau), r_f), r_comh)| {
                s[0] * r_Mf + s[1] * r_mtau + s[2] * r_f + r_comh
            })
            .collect::<Vec<R>>();

        let v0 = once(&proof.evals.0)
            .chain(once(&proof.evals.1))
            .map(|evals| {
                let evals = evals[0];
                (s[0] * evals[0]) + (s[1] * evals[1]) + (s[2] * evals[2]) + evals[3]
            })
            .collect::<Vec<R>>();

        Ok(())
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
    use crate::utils::split;

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

        let cfs: Matrix<_> = f
            .iter()
            .map(|r| r.coeffs().to_vec())
            .collect::<Vec<Vec<_>>>()
            .into();
        let dec = cfs
            .vals
            .iter()
            .map(|row| row.decompose_to_vec(b, k as usize))
            .collect::<Vec<_>>();

        let mut D_f = vec![Matrix::zero(n, R::dimension()); k as usize];

        // map dec: (Z n x d x k) to D_f: (Z n x d, k matrices)
        dec.iter().enumerate().for_each(|(n_i, drow)| {
            drow.iter().enumerate().for_each(|(d_i, coeffs)| {
                coeffs.iter().enumerate().for_each(|(k_i, coeff)| {
                    D_f[k_i].vals[n_i][d_i] = *coeff;
                });
            });
        });

        let M_f: Vec<Matrix<R>> = D_f
            .iter()
            .map(|m| {
                m.vals
                    .iter()
                    .map(|row| {
                        row.iter()
                            .map(|c| exp::<R>(*c).unwrap())
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>()
                    .into()
            })
            .collect::<Vec<_>>();

        let A = Matrix::<R>::rand(&mut rand::thread_rng(), kappa, n);

        let comM_f = M_f
            .iter()
            .map(|M| A.try_mul_mat(M).unwrap())
            .collect::<Vec<_>>();
        let com = Matrix::hconcat(&comM_f).unwrap();

        let tau = split(&com, n, (R::dimension() / 2) as u128, l);

        let m_tau = tau
            .iter()
            .map(|c| exp::<R>(*c).unwrap())
            .collect::<Vec<_>>();

        let cm_f = A.try_mul_vec(&f).unwrap();
        let C_Mf = A
            .try_mul_vec(&tau.iter().map(|z| R::from(*z)).collect::<Vec<R>>())
            .unwrap();
        let cm_mtau = A.try_mul_vec(&m_tau).unwrap();

        let rg = Rg {
            nvars: log2(n) as usize,
            M_f,
            f,
            tau,
            m_tau,
            comM_f,
            M,
            b,
            k,
            l,
        };

        let cm = Cm {
            rg,
            cm_f,
            C_Mf,
            cm_mtau,
        };

        let mut ts = PoseidonTS::default::<PC>();
        let proof = cm.prove(&mut ts);

        let mut ts = PoseidonTS::default::<PC>();
        cm.verify(&proof, &mut ts).unwrap();
    }
}
