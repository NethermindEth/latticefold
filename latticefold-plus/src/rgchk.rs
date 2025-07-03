use ark_std::{log2, Zero};
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

use crate::setchk::{In, MonomialSet, Out};

// D_f: decomposed cf(f), Z n x dk
// M_f: EXP(D_f)

#[derive(Debug)]
pub struct Rg<R: CoeffRing>
where
    R::BaseRing: Zq,
{
    pub nvars: usize,
    pub M_f: Vec<Matrix<R>>,   // n x d, k matrices, monomials
    pub f: Vec<R>,             // n
    pub tau: Vec<R::BaseRing>, // n
    pub m_tau: Vec<R>,         // n, monomials

    // decomposition
    pub b: u128,
    pub k: usize,
}

#[derive(Debug)]
pub struct Dcom<R: CoeffRing>
where
    R::BaseRing: Zq,
{
    pub v: Vec<R::BaseRing>,
    pub a: R::BaseRing,
    pub b: R,
    pub u: Vec<Vec<R>>,
}

pub fn split<R: PolyRing>(com: Matrix<R>, n: usize, b: u128, k: usize) -> Vec<R::BaseRing>
where
    R: Decompose,
{
    let M_prime = com.gadget_decompose(b, k);
    let M_dprime = M_prime.vals.into_iter().fold(vec![], |mut acc, row| {
        // TODO pre-alloc
        acc.extend(row);
        acc
    });
    let mut tau = M_dprime
        .iter()
        .map(|r| r.coeffs().to_vec())
        .into_iter()
        .fold(vec![], |mut acc, row| {
            // TODO pre-alloc
            acc.extend(row);
            acc
        });
    if tau.len() < n {
        // TODO handle when opposite
        tau.resize(n, R::BaseRing::zero());
    }
    tau
}

impl<R: CoeffRing> Rg<R>
where
    R::BaseRing: Zq,
{
    pub fn range_check(&self, transcript: &mut impl Transcript<R>) -> Dcom<R> {
        let in_rel = In {
            sets: self
                .M_f
                .iter()
                .map(|m| MonomialSet::Matrix(SparseMatrix::<R>::from_dense(m)))
                .chain(
                    ark_std::iter::once(&self.m_tau)
                        .map(|m_tau| MonomialSet::Vector(m_tau.clone())),
                )
                .collect::<Vec<_>>(),

            nvars: self.nvars,
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

        let a = DenseMultilinearExtension::from_evaluations_vec(self.nvars, self.tau.clone())
            .evaluate(&out_rel.r)
            .unwrap();

        let b = out_rel.b[0];

        Dcom {
            v,
            a,
            b,
            u: out_rel.e,
        }
    }
}

impl<R: CoeffRing> Dcom<R>
where
    R::BaseRing: Zq,
{
    pub fn verify(&self, transcript: &mut impl Transcript<R>) -> Result<(), ()> {
        // TODO Run set checks

        ((psi::<R>() * self.b).ct() == self.a)
            .then(|| ())
            .ok_or(())?;

        let d = R::dimension();
        let d_prime = d / 2;
        let u_comb = self
            .u
            .iter()
            .enumerate()
            .fold(vec![R::zero(); d], |mut acc, (i, u_i)| {
                let d_ppow = R::BaseRing::from(d_prime as u128).pow([i as u64]);
                u_i.iter()
                    .zip(acc.iter_mut())
                    .for_each(|(u_ij, a_j)| *a_j += *u_ij * d_ppow);
                acc
            });

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
    use ark_std::{One, Zero};
    use cyclotomic_rings::rings::GoldilocksPoseidonConfig as PC;
    use latticefold::transcript::poseidon::PoseidonTS;
    use stark_rings::{
        balanced_decomposition::DecomposeToVec, cyclotomic_ring::models::goldilocks::RqPoly as R,
    };
    use stark_rings_linalg::SparseMatrix;

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

        let coms = M_f
            .iter()
            .map(|M| A.try_mul_mat(M).unwrap())
            .collect::<Vec<_>>();
        let com = Matrix::hconcat(&coms).unwrap();

        let tau = split(com, n, (R::dimension() / 2) as u128, l);

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
            b,
            k,
        };

        let mut ts = PoseidonTS::default::<PC>();
        let dcom = rg.range_check(&mut ts);

        let mut ts = PoseidonTS::default::<PC>();
        dcom.verify(&mut ts).unwrap();
    }
}
