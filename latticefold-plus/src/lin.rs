use ark_std::log2;
use latticefold::{
    arith::r1cs::R1CS,
    transcript::Transcript,
    utils::sumcheck::{
        utils::{build_eq_x_r, eq_eval},
        MLSumcheck,
    },
};
use stark_rings::{
    balanced_decomposition::{
        convertible_ring::ConvertibleRing, Decompose, DecomposeToVec, GadgetDecompose,
    },
    exp, psi, psi_range_check, unit_monomial, CoeffRing, OverField, PolyRing, Ring, Zq,
};
use stark_rings_linalg::{ops::Transpose, Matrix, SparseMatrix};
use stark_rings_poly::mle::DenseMultilinearExtension;

use crate::{
    cm::{Cm, CmProof},
    rgchk::{DecompParameters, Rg},
    utils::split,
};

pub trait Linearize<R: OverField> {
    type Proof: Verify<R>;
    fn linearize(&self, transcript: &mut impl Transcript<R>) -> Self::Proof;
}

pub trait Verify<R: OverField> {
    fn verify(&self, transcript: &mut impl Transcript<R>) -> bool;
}

#[derive(Debug)]
pub struct Parameters {
    pub kappa: usize,
    pub decomp: DecompParameters,
}

#[derive(Debug)]
pub struct Lin<R> {
    pub M: Vec<SparseMatrix<R>>,
    pub v: Vec<R>,
    pub f: Vec<R>,

    pub params: Parameters,
}

impl<R: CoeffRing> Lin<R>
where
    R::BaseRing: ConvertibleRing + Decompose + Zq,
    R: Decompose,
{
    /// Π'cm protocol
    pub fn cm(&self, transcript: &mut impl Transcript<R>) -> CmProof<R> {
        let n = self.f.len();

        let cfs: Matrix<_> = self
            .f
            .iter()
            .map(|r| r.coeffs().to_vec())
            .collect::<Vec<Vec<_>>>()
            .into();
        let dec = cfs
            .vals
            .iter()
            .map(|row| row.decompose_to_vec(self.params.decomp.b, self.params.decomp.k))
            .collect::<Vec<_>>();

        let mut D_f = vec![Matrix::zero(n, R::dimension()); self.params.decomp.k];

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

        let A = Matrix::<R>::rand(&mut rand::thread_rng(), self.params.kappa, n);

        let comM_f = M_f
            .iter()
            .map(|M| A.try_mul_mat(M).unwrap())
            .collect::<Vec<_>>();
        let com = Matrix::hconcat(&comM_f).unwrap();

        let tau = split(&com, n, (R::dimension() / 2) as u128, self.params.decomp.l);

        let m_tau = tau
            .iter()
            .map(|c| exp::<R>(*c).unwrap())
            .collect::<Vec<_>>();

        let cm_f = A.try_mul_vec(&self.f).unwrap();
        let C_Mf = A
            .try_mul_vec(&tau.iter().map(|z| R::from(*z)).collect::<Vec<R>>())
            .unwrap();
        let cm_mtau = A.try_mul_vec(&m_tau).unwrap();

        let rg = Rg {
            nvars: log2(n) as usize,
            M_f,
            f: self.f.clone(),
            tau,
            m_tau,
            comM_f,
            M: self.M.clone(),
            dparams: self.params.decomp.clone(),
        };

        let cm = Cm {
            rg,
            cm_f,
            C_Mf,
            cm_mtau,
        };

        cm.prove(transcript)
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
    use crate::{r1cs::CommittedR1CS, utils::split};

    #[test]
    fn test_lin() {
        let n = 1 << 15;
        let z = vec![R::one(); n];

        let mut r1cs = R1CS::<R> {
            l: 1,
            A: SparseMatrix::identity(n),
            B: SparseMatrix::identity(n),
            C: SparseMatrix::identity(n),
        };

        r1cs.A.coeffs[0][0].0 = 2u128.into();
        r1cs.C.coeffs[0][0].0 = 2u128.into();

        r1cs.A = r1cs.A.gadget_decompose(2, 4);
        r1cs.B = r1cs.B.gadget_decompose(2, 4);
        r1cs.C = r1cs.C.gadget_decompose(2, 4);

        let f = z.gadget_decompose(2, 4);
        r1cs.check_relation(&f).unwrap();

        let cr1cs = CommittedR1CS {
            r1cs,
            f,
            x: vec![1u128.into()],
            cm: vec![],
        };

        let kappa = 2;
        let b = (R::dimension() / 2) as u128;
        let k = 2;
        // log_d' (q)
        let l = ((<<R as PolyRing>::BaseRing>::MODULUS.0[0] as f64).ln()
            / ((R::dimension() / 2) as f64).ln())
        .ceil() as usize;

        let mut ts = PoseidonTS::default::<PC>();
        let lproof = cr1cs.linearize(&mut ts);

        let lin = Lin {
            M: vec![cr1cs.r1cs.A, cr1cs.r1cs.B, cr1cs.r1cs.C],
            v: vec![lproof.va, lproof.vb, lproof.vc],
            f: cr1cs.f,
            params: Parameters {
                kappa,
                decomp: DecompParameters { b, k, l },
            },
        };

        let cmproof = lin.cm(&mut ts);

        let mut ts = PoseidonTS::default::<PC>();
        lproof.verify(&mut ts);
        cmproof.verify(&mut ts).unwrap();
    }
}
