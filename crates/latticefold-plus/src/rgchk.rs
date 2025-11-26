use ark_std::iter::once;
use latticefold::transcript::Transcript;
use stark_rings::{
    balanced_decomposition::{Decompose, DecomposeToVec},
    exp, psi, CoeffRing, OverField, PolyRing, Ring, Zq,
};
use stark_rings_linalg::{ops::Transpose, Matrix, SparseMatrix};
use stark_rings_poly::mle::DenseMultilinearExtension;
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

#[derive(Clone, Debug)]
pub struct FComs<R> {
    pub cm_f: Vec<R>,
    pub C_Mf: Vec<R>,
    pub cm_mtau: Vec<R>,
}

#[derive(Clone, Debug)]
pub struct Rg<R: PolyRing> {
    pub nvars: usize,
    pub instances: Vec<RgInstance<R>>, // L instances
    pub dparams: DecompParameters,
}

#[derive(Clone, Debug)]
pub struct RgInstance<R: PolyRing> {
    pub M_f: Vec<Matrix<R>>,   // n x d, k matrices, monomials
    pub tau: Vec<R::BaseRing>, // n
    pub m_tau: Vec<R>,         // n, monomials
    pub f: Vec<R>,             // n
    pub comM_f: Vec<Matrix<R>>,
    pub fcoms: FComs<R>,
}

#[derive(Clone, Debug)]
pub struct Dcom<R: PolyRing> {
    pub evals: Vec<DcomEvals<R>>, // L evals
    pub fcoms: Vec<FComs<R>>,     // L commitments
    pub out: Out<R>,              // set checks
    pub dparams: DecompParameters,
}

#[derive(Clone, Debug)]
pub struct DcomEvals<R: PolyRing> {
    pub v: Vec<R::BaseRing>, // eval over M_f
    pub a: Vec<R::BaseRing>, // eval over tau
    pub b: Vec<R>,           // eval over m_tau
    pub c: Vec<R>,           // eval over f
}

#[derive(Debug, Error)]
pub enum RangeCheckError<R: PolyRing> {
    #[error("Psi check failed: a = {0}, b = {1}")]
    PsiCheckAB(R::BaseRing, R),
    #[error("Psi check failed: v = {0}, u-comb = {1}")]
    PsiCheckVU(Vec<R::BaseRing>, Vec<R>),
}

impl<R: CoeffRing> Rg<R>
where
    R::BaseRing: Zq,
{
    /// Range checks
    ///
    /// Support for `L` [`RgInstance`]s mapped to the corresponding [`DcomEvals`].
    pub fn range_check(
        &self,
        M: &[SparseMatrix<R>],
        transcript: &mut impl Transcript<R>,
    ) -> Dcom<R> {
        let mut sets = Vec::with_capacity(self.instances.len() * (self.instances[0].M_f.len() + 1));
        for inst in &self.instances {
            inst.M_f.iter().for_each(|m| {
                sets.push(MonomialSet::Matrix(SparseMatrix::<R>::from_dense(m)));
            });
        }
        for inst in &self.instances {
            sets.push(MonomialSet::Vector(inst.m_tau.clone()));
        }

        let in_rel = In {
            sets,
            nvars: self.nvars,
        };
        let out_rel = in_rel.set_check(M, transcript);

        let evals = self
            .instances
            .iter()
            .enumerate()
            .map(|(l, inst)| {
                let cfs = inst
                    .f
                    .iter()
                    .map(|r| r.coeffs().to_vec())
                    .collect::<Vec<_>>()
                    .transpose();
                let v = cfs
                    .into_iter()
                    .map(|evals| {
                        let mle =
                            DenseMultilinearExtension::from_evaluations_vec(self.nvars, evals);
                        mle.evaluate(&out_rel.r).unwrap()
                    })
                    .collect::<Vec<_>>();

                // TODO v is equal to c[0]

                let r = out_rel.r.iter().map(|z| R::from(*z)).collect::<Vec<_>>();

                let mut a = Vec::with_capacity(1 + M.len());
                let mut b = Vec::with_capacity(1 + M.len());
                // Let `c` be the evaluation of `f` over r
                let mut c = Vec::with_capacity(1 + M.len());

                a.push(
                    DenseMultilinearExtension::from_evaluations_vec(self.nvars, inst.tau.clone())
                        .evaluate(&out_rel.r)
                        .unwrap(),
                );

                b.push(out_rel.b[l]);

                c.push(
                    DenseMultilinearExtension::from_evaluations_vec(self.nvars, inst.f.clone())
                        .evaluate(&out_rel.r.iter().map(|z| R::from(*z)).collect::<Vec<_>>())
                        .unwrap(),
                );

                M.iter().for_each(|m| {
                    let Mtau = m
                        .try_mul_vec(&inst.tau.iter().map(|z| R::from(*z)).collect::<Vec<R>>())
                        .unwrap();
                    a.push(
                        DenseMultilinearExtension::from_evaluations_vec(self.nvars, Mtau)
                            .evaluate(&r)
                            .unwrap()
                            .ct(),
                    );

                    let Mm_tau = m.try_mul_vec(&inst.m_tau).unwrap();
                    b.push(
                        DenseMultilinearExtension::from_evaluations_vec(self.nvars, Mm_tau)
                            .evaluate(&r)
                            .unwrap(),
                    );

                    let Mf = m.try_mul_vec(&inst.f).unwrap();
                    c.push(
                        DenseMultilinearExtension::from_evaluations_vec(self.nvars, Mf)
                            .evaluate(&r)
                            .unwrap(),
                    );
                });
                DcomEvals { v, a, b, c }
            })
            .collect::<Vec<_>>();

        absorb_evaluations(&evals, transcript);

        Dcom {
            evals,
            fcoms: self
                .instances
                .iter()
                .map(|inst| inst.fcoms.clone())
                .collect(),
            out: out_rel,
            dparams: self.dparams.clone(),
        }
    }
}

impl<R: CoeffRing> Dcom<R>
where
    R::BaseRing: Zq,
{
    pub fn verify(&self, transcript: &mut impl Transcript<R>) -> Result<(), RangeCheckError<R>> {
        self.out.verify(transcript).unwrap(); //.map_err(|_| ())?;

        absorb_evaluations(&self.evals, transcript);

        for (l, eval) in self.evals.iter().enumerate() {
            // ct(psi b) =? a
            for (&a_i, b_i) in eval.a.iter().zip(eval.b.iter()) {
                ((psi::<R>() * b_i).ct() == a_i)
                    .then_some(())
                    .ok_or(RangeCheckError::PsiCheckAB(a_i, *b_i))?;
            }

            let d = R::dimension();
            let d_prime = d / 2;
            for (ni, _) in self.out.e.iter().enumerate() {
                let u_comb = self.out.e[ni]
                    .iter()
                    .skip(self.dparams.k * l)
                    .take(self.dparams.k)
                    .enumerate()
                    .fold(vec![R::zero(); d], |mut acc, (i, u_i)| {
                        let d_ppow = R::BaseRing::from(d_prime as u128).pow([i as u64]);
                        u_i.iter()
                            .zip(acc.iter_mut())
                            .for_each(|(u_ij, a_j)| *a_j += *u_ij * d_ppow);
                        acc
                    });

                // ct(psi (sum d^i u_i)) =? v
                let v_rec = u_comb
                    .iter()
                    .map(|&uc| (psi::<R>() * uc).ct())
                    .collect::<Vec<_>>();

                if ni == 0 {
                    (eval.v == v_rec)
                        .then_some(())
                        .ok_or(RangeCheckError::PsiCheckVU(v_rec, u_comb))?;
                } else {
                    (eval.c[ni].coeffs() == v_rec)
                        .then_some(())
                        .ok_or(RangeCheckError::PsiCheckVU(v_rec, u_comb))?;
                }
            }
        }

        Ok(())
    }
}

impl<R: PolyRing> RgInstance<R> {
    /// Construct monomial sets from `M_f` and `m_tau`
    pub fn sets(&self) -> Vec<MonomialSet<R>> {
        self.M_f
            .iter()
            .map(|m| MonomialSet::Matrix(SparseMatrix::<R>::from_dense(m)))
            .chain(once(MonomialSet::Vector(self.m_tau.clone())))
            .collect()
    }
}

impl<R: CoeffRing> RgInstance<R>
where
    R::BaseRing: Decompose + Zq,
    R: Decompose,
{
    pub fn from_f(f: Vec<R>, A: &Matrix<R>, decomp: &DecompParameters) -> Self {
        let n = f.len();

        let cfs: Matrix<_> = f
            .iter()
            .map(|r| r.coeffs().to_vec())
            .collect::<Vec<Vec<_>>>()
            .into();
        let dec = cfs
            .vals
            .iter()
            .map(|row| row.decompose_to_vec(decomp.b, decomp.k))
            .collect::<Vec<_>>();

        let mut D_f = vec![Matrix::zero(n, R::dimension()); decomp.k];

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

        let comM_f = M_f
            .iter()
            .map(|M| A.try_mul_mat(M).unwrap())
            .collect::<Vec<_>>();
        let com = Matrix::hconcat(&comM_f).unwrap();

        let tau = split(&com, n, (R::dimension() / 2) as u128, decomp.l);

        let m_tau = tau
            .iter()
            .map(|c| exp::<R>(*c).unwrap())
            .collect::<Vec<_>>();

        let cm_f = A.try_mul_vec(&f).unwrap();
        let C_Mf = A
            .try_mul_vec(&tau.iter().map(|z| R::from(*z)).collect::<Vec<R>>())
            .unwrap();
        let cm_mtau = A.try_mul_vec(&m_tau).unwrap();
        let fcoms = FComs {
            cm_f,
            C_Mf,
            cm_mtau,
        };

        Self {
            M_f,
            tau,
            m_tau,
            f,
            comM_f,
            fcoms,
        }
    }
}

fn absorb_evaluations<R: OverField>(evals: &[DcomEvals<R>], transcript: &mut impl Transcript<R>) {
    evals.iter().for_each(|eval| {
        transcript.absorb_slice(&eval.a.iter().map(|z| R::from(*z)).collect::<Vec<R>>());
        transcript.absorb_slice(&eval.c);
    });
}

#[cfg(test)]
mod tests {
    use ark_ff::PrimeField;
    use ark_std::{log2, Zero};
    use latticefold_rings::rings::FrogPoseidonConfig as PC;
    use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as R;

    use super::*;
    use crate::transcript::PoseidonTranscript;

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

        let A = Matrix::<R>::rand(&mut ark_std::test_rng(), kappa, n);

        let dparams = DecompParameters { b, k, l };
        let instance = RgInstance::from_f(f.clone(), &A, &dparams);

        let rg = Rg {
            nvars: log2(n) as usize,
            instances: vec![instance],
            dparams,
        };

        let mut ts = PoseidonTranscript::empty::<PC>();
        let dcom = rg.range_check(&[], &mut ts);

        let mut ts = PoseidonTranscript::empty::<PC>();
        dcom.verify(&mut ts).unwrap();
    }

    #[test]
    fn test_range_check_mm() {
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

        let kappa = 1;
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
            dparams,
        };

        let mut ts = PoseidonTranscript::empty::<PC>();
        let dcom = rg.range_check(&M, &mut ts);

        let mut ts = PoseidonTranscript::empty::<PC>();
        dcom.verify(&mut ts).unwrap();
    }
}
