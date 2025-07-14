use ark_std::{log2, One, Zero};
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
    utils::{from_bits, split, to_bits},
};

pub struct Cm<R: PolyRing> {
    pub rg: Rg<R>,
}

pub struct CmProof<R: PolyRing> {
    pub dcom: Dcom<R>,
    pub comh: Vec<R>,
    pub sumcheck_proof: Proof<R>,
    pub evals: [R; 4], // eval over r0 of tau (a), m_tau (b), f (c), h (u)
}

impl<R: CoeffRing> Cm<R>
where
    R::BaseRing: Zq,
{
    pub fn prove(&self, transcript: &mut impl Transcript<R>) -> CmProof<R> {
        let k = self.rg.comM_f.len();
        let d = R::dimension();

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

        let (proof, evals) = self.sumchecker(&dcom, &h, &s_prime, transcript);

        CmProof {
            dcom,
            comh,
            sumcheck_proof: proof,
            evals,
        }
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

        let rc1: R = transcript.get_challenge().into();

        let u: R = proof
            .dcom
            .out
            .e
            .iter()
            .flatten()
            .zip(s_prime.iter().flatten())
            .map(|(u_ij, s_ij)| *u_ij * *s_ij)
            .sum();
        let claimed_sum = R::from(proof.dcom.a)
            + proof.dcom.b * rc1
            + proof.dcom.c * rc1.pow([2])
            + u * rc1.pow([3]);

        let subclaim = MLSumcheck::verify_as_subprotocol(
            transcript,
            self.rg.nvars,
            2,
            claimed_sum,
            &proof.sumcheck_proof,
        )
        .unwrap();

        let r: Vec<R> = proof.dcom.out.r.iter().map(|x| R::from(*x)).collect();
        let ro: Vec<R> = subclaim.point.into_iter().map(|x| x.into()).collect();
        let s = subclaim.expected_evaluation;
        let e = eq_eval(&r, &ro).unwrap();

        assert_eq!(
            s,
            e * (proof.evals[0]
                + rc1 * proof.evals[1]
                + rc1.pow([2]) * proof.evals[2]
                + rc1.pow([3]) * proof.evals[3])
        );

        Ok(())
    }

    fn sumchecker(
        &self,
        dcom: &Dcom<R>,
        h: &Vec<R>,
        s_prime: &Vec<Vec<R>>,
        transcript: &mut impl Transcript<R>,
    ) -> (Proof<R>, [R; 4]) {
        let nvars = self.rg.nvars;
        let r: Vec<R> = dcom.out.r.iter().map(|x| R::from(*x)).collect();

        let eq = build_eq_x_r(&r).unwrap();
        let rc1: R = transcript.get_challenge().into();

        let tau_mle = DenseMultilinearExtension::from_evaluations_vec(
            nvars,
            self.rg.tau.iter().map(|z| R::from(*z)).collect::<Vec<_>>(),
        );
        let m_tau_mle =
            DenseMultilinearExtension::from_evaluations_vec(nvars, self.rg.m_tau.clone());
        let f_mle = DenseMultilinearExtension::from_evaluations_vec(nvars, self.rg.f.clone());
        let h_mle = DenseMultilinearExtension::from_evaluations_vec(nvars, h.clone());

        let u: R = dcom
            .out
            .e
            .iter()
            .flatten()
            .zip(s_prime.iter().flatten())
            .map(|(u_ij, s_ij)| *u_ij * *s_ij)
            .sum();

        let mles = vec![
            eq,
            tau_mle.clone(),
            m_tau_mle.clone(),
            f_mle.clone(),
            h_mle.clone(),
        ];

        let comb_fn = |vals: &[R]| -> R {
            vals[0] // eq
                * (
                    vals[1]  // tau
                    + vals[2] * rc1 // m_tau
                    + vals[3] * (rc1 * rc1) // f
                    + vals[4] * (rc1 * rc1 * rc1) // h
                )
        };

        let (sumcheck_proof, prover_state) =
            MLSumcheck::prove_as_subprotocol(transcript, mles, nvars, 2, comb_fn);
        let ro = prover_state
            .randomness
            .into_iter()
            .map(|x| x.into())
            .collect::<Vec<R>>();

        let va = tau_mle.evaluate(&ro).unwrap();
        let vb = m_tau_mle.evaluate(&ro).unwrap();
        let vc = f_mle.evaluate(&ro).unwrap();
        let vh = h_mle.evaluate(&ro).unwrap();

        (sumcheck_proof, [va, vb, vc, vh])
    }
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
        let n = 1 << 17;
        let mut f = vec![R::zero(); n];
        f[0].coeffs_mut()[0] = 2u128.into();
        f[0].coeffs_mut()[1] = 5u128.into();
        f[1].coeffs_mut()[0] = 4u128.into();
        f[1].coeffs_mut()[2] = 1u128.into();

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

        let rg = Rg {
            nvars: log2(n) as usize,
            M_f,
            f,
            tau,
            m_tau,
            comM_f,
            b,
            k,
            l,
        };

        let cm = Cm { rg };

        let mut ts = PoseidonTS::default::<PC>();
        let proof = cm.prove(&mut ts);

        let mut ts = PoseidonTS::default::<PC>();
        cm.verify(&proof, &mut ts).unwrap();
    }
}
