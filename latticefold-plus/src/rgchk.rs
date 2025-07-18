use ark_std::{iter::once, log2, One, Zero};
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
    exp, psi, psi_range_check, CoeffRing, OverField, PolyRing, Ring, Zq,
};
use stark_rings_linalg::{ops::Transpose, Matrix, SparseMatrix};
use stark_rings_poly::mle::{DenseMultilinearExtension, SparseMultilinearExtension};
use thiserror::Error;

use crate::{
    setchk::{In, MonomialSet, Out},
    utils::split,
};

// D_f: decomposed cf(f), Z n x dk
// M_f: EXP(D_f)

#[derive(Clone, Debug)]
pub struct DecompParameters {
    pub b: u128,
    pub k: usize,
    pub l: usize,
}

#[derive(Debug)]
pub struct Rg<R: PolyRing> {
    pub nvars: usize,
    pub M_f: Vec<Matrix<R>>,   // n x d, k matrices, monomials
    pub tau: Vec<R::BaseRing>, // n
    pub m_tau: Vec<R>,         // n, monomials
    pub f: Vec<R>,             // n
    pub comM_f: Vec<Matrix<R>>,
    pub M: Vec<SparseMatrix<R>>, // n_lin matrices, n x n
    pub dparams: DecompParameters,
}

#[derive(Debug)]
pub struct Dcom<R: PolyRing> {
    pub v: Vec<R::BaseRing>, // eval over M_f
    pub a: Vec<R::BaseRing>, // eval over tau
    pub b: Vec<R>,           // eval over m_tau
    pub c: Vec<R>,           // eval over f
    pub out: Out<R>,         // set checks
    pub dparams: DecompParameters,
}

impl<R: CoeffRing> Rg<R>
where
    R::BaseRing: Zq,
{
    /// Range checks
    pub fn range_check(&self, transcript: &mut impl Transcript<R>) -> Dcom<R> {
        let sets = self
            .M_f
            .iter()
            .map(|m| MonomialSet::Matrix(SparseMatrix::<R>::from_dense(m)))
            .chain(once(MonomialSet::Vector(self.m_tau.clone())))
            .collect::<Vec<_>>();

        let in_rel = In {
            sets,
            nvars: self.nvars,
            M: self.M.clone(),
        };
        let out_rel = in_rel.set_check(transcript);

        let cfs = self
            .f
            .iter()
            .map(|r| r.coeffs().to_vec())
            .collect::<Vec<_>>()
            .transpose();
        let v = cfs
            .into_iter()
            .map(|evals| {
                let mle = DenseMultilinearExtension::from_evaluations_vec(self.nvars, evals);
                mle.evaluate(&out_rel.r).unwrap()
            })
            .collect::<Vec<_>>();

        let r = out_rel.r.iter().map(|z| R::from(*z)).collect::<Vec<_>>();

        let mut a = Vec::with_capacity(1 + self.M.len());
        let mut b = Vec::with_capacity(1 + self.M.len());
        // Let `c` be the evaluation of `f` over r
        let mut c = Vec::with_capacity(1 + self.M.len());

        a.push(
            DenseMultilinearExtension::from_evaluations_vec(self.nvars, self.tau.clone())
                .evaluate(&out_rel.r)
                .unwrap(),
        );

        b.push(out_rel.b[0]);

        c.push(
            DenseMultilinearExtension::from_evaluations_vec(self.nvars, self.f.clone())
                .evaluate(&out_rel.r.iter().map(|z| R::from(*z)).collect::<Vec<_>>())
                .unwrap(),
        );

        self.M.iter().for_each(|m| {
            let Mtau = m
                .try_mul_vec(&self.tau.iter().map(|z| R::from(*z)).collect::<Vec<R>>())
                .unwrap();
            a.push(
                DenseMultilinearExtension::from_evaluations_vec(self.nvars, Mtau)
                    .evaluate(&r)
                    .unwrap()
                    .ct(),
            );

            let Mm_tau = m
                .try_mul_vec(&self.m_tau.iter().map(|z| R::from(*z)).collect::<Vec<R>>())
                .unwrap();
            b.push(
                DenseMultilinearExtension::from_evaluations_vec(self.nvars, Mm_tau)
                    .evaluate(&r)
                    .unwrap(),
            );

            let Mf = m.try_mul_vec(&self.f).unwrap();
            c.push(
                DenseMultilinearExtension::from_evaluations_vec(self.nvars, Mf)
                    .evaluate(&r)
                    .unwrap(),
            );
        });

        Dcom {
            v,
            a,
            b,
            c,
            out: out_rel,
            dparams: self.dparams.clone(),
        }
    }
}

impl<R: CoeffRing> Dcom<R>
where
    R::BaseRing: Zq,
{
    pub fn verify(&self, transcript: &mut impl Transcript<R>) -> Result<(), ()> {
        self.out.verify(transcript).unwrap(); //.map_err(|_| ())?;

        // ct(psi b) =? a
        for (&a_i, b_i) in self.a.iter().zip(self.b.iter()) {
            ((psi::<R>() * b_i).ct() == a_i)
                .then(|| ())
                .ok_or(())
                .unwrap();
        }

        let d = R::dimension();
        let d_prime = d / 2;
        let u_comb = self.out.e[0].iter().take(self.dparams.k).enumerate().fold(
            vec![R::zero(); d],
            |mut acc, (i, u_i)| {
                let d_ppow = R::BaseRing::from(d_prime as u128).pow([i as u64]);
                u_i.iter()
                    .zip(acc.iter_mut())
                    .for_each(|(u_ij, a_j)| *a_j += *u_ij * d_ppow);
                acc
            },
        );

        // ct(psi (sum d^i u_i)) =? v
        let v_rec = u_comb
            .iter()
            .map(|uc| (psi::<R>() * uc).ct())
            .collect::<Vec<_>>();

        (self.v == v_rec).then(|| ()).ok_or(())?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::PrimeField;
    use ark_std::Zero;
    use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
    use latticefold::transcript::poseidon::PoseidonTS;
    use stark_rings::{
        balanced_decomposition::DecomposeToVec, cyclotomic_ring::models::frog_ring::RqPoly as R,
    };

    use super::*;

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
            M: vec![],
            dparams: DecompParameters { b, k, l },
        };

        let mut ts = PoseidonTS::default::<PC>();
        let dcom = rg.range_check(&mut ts);

        let mut ts = PoseidonTS::default::<PC>();
        dcom.verify(&mut ts).unwrap();
    }
}
