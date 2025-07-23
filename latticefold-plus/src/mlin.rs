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
    cm::{Cm, CmComs, CmProof},
    lin::{Lin, Parameters},
    rgchk::{Dcom, Rg, RgInstance},
    utils::{tensor, tensor_product},
};

#[derive(Clone, Debug)]
pub struct Mlin<R> {
    pub lins: Vec<Lin<R>>,
    pub A: Matrix<R>,
    pub params: Parameters,
}

#[derive(Clone, Debug)]
pub struct LinB2X<R> {
    pub cm_g: Vec<R>,
    pub ro: Vec<(R, R)>,
    pub vo: Vec<(R, R)>,
}

#[derive(Clone, Debug)]
pub struct LinB2<R> {
    pub g: Vec<R>,
    pub x: LinB2X<R>,
}

impl<R: CoeffRing> Mlin<R>
where
    R::BaseRing: ConvertibleRing + Decompose + Zq,
    R: Decompose,
{
    /// Πmlin protocol
    pub fn mlin(&self, transcript: &mut impl Transcript<R>) -> (LinB2<R>, CmProof<R>) {
        let n = self.lins[0].f.len();

        let instances = self
            .lins
            .iter()
            .map(|lin| RgInstance::from_f(lin.f.clone(), &self.A, &self.params.decomp))
            .collect::<Vec<_>>();

        let coms = instances
            .iter()
            .map(|inst| {
                let cm_f = self.A.try_mul_vec(&inst.f).unwrap();
                let C_Mf = self
                    .A
                    .try_mul_vec(&inst.tau.iter().map(|z| R::from(*z)).collect::<Vec<R>>())
                    .unwrap();
                let cm_mtau = self.A.try_mul_vec(&inst.m_tau).unwrap();

                CmComs {
                    cm_f,
                    C_Mf,
                    cm_mtau,
                }
            })
            .collect::<Vec<_>>();

        let rg = Rg {
            nvars: log2(n) as usize,
            instances,
            M: self.lins[0].M.clone(),
            dparams: self.params.decomp.clone(),
        };

        let cm = Cm { rg, coms };

        let (com, proof) = cm.prove(transcript);

        let cm_g = com
            .x
            .cm_g
            .iter()
            .fold(vec![R::zero(); self.params.kappa], |mut acc, cm| {
                acc.iter_mut().zip(cm.iter()).for_each(|(acc_r, cm_r)| {
                    *acc_r += cm_r;
                });
                acc
            });

        let nlin = com.x.vo[0].len();
        let vo = com
            .x
            .vo
            .iter()
            .fold(vec![(R::zero(), R::zero()); nlin], |mut acc, v| {
                v.iter().enumerate().for_each(|(i, v)| {
                    acc[i].0 += v.0;
                    acc[i].1 += v.1;
                });
                acc
            });

        let x = LinB2X {
            cm_g,
            ro: com.x.ro,
            vo,
        };

        let g = com.g.iter().fold(vec![R::zero(); n], |mut acc, gi| {
            acc.iter_mut().zip(gi.iter()).for_each(|(acc_r, gi_r)| {
                *acc_r += gi_r;
            });
            acc
        });
        let linb2 = LinB2 { g, x };

        (linb2, proof)
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::PrimeField;
    use ark_std::{One, Zero};
    use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
    use latticefold::{arith::r1cs::R1CS, transcript::poseidon::PoseidonTS};
    use stark_rings::{
        balanced_decomposition::DecomposeToVec, cyclotomic_ring::models::frog_ring::RqPoly as R,
    };
    use stark_rings_linalg::SparseMatrix;

    use super::*;
    use crate::{
        lin::{Linearize, Verify},
        r1cs::CommittedR1CS,
        rgchk::DecompParameters,
        utils::split,
    };

    #[test]
    fn test_mlin() {
        let n = 1 << 15;
        let k = 2;
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

        let cr1cs0 = CommittedR1CS {
            r1cs: r1cs.clone(),
            f: f0,
            x: vec![1u128.into()],
            cm: vec![],
        };
        let cr1cs1 = CommittedR1CS {
            r1cs,
            f: f1,
            x: vec![1u128.into()],
            cm: vec![],
        };

        let kappa = 2;
        let b = (R::dimension() / 2) as u128;
        // log_d' (q)
        let l = ((<<R as PolyRing>::BaseRing>::MODULUS.0[0] as f64).ln()
            / ((R::dimension() / 2) as f64).ln())
        .ceil() as usize;

        let mut ts = PoseidonTS::default::<PC>();
        let lproof0 = cr1cs0.linearize(&mut ts);
        let lproof1 = cr1cs1.linearize(&mut ts);

        let params = Parameters {
            kappa,
            decomp: DecompParameters { b, k, l },
        };

        let lin0 = Lin {
            M: vec![
                cr1cs0.r1cs.A.clone(),
                cr1cs0.r1cs.B.clone(),
                cr1cs0.r1cs.C.clone(),
            ],
            v: vec![lproof0.va, lproof0.vb, lproof0.vc],
            f: cr1cs0.f.clone(),
            params: params.clone(),
        };

        let lin1 = Lin {
            M: vec![cr1cs1.r1cs.A, cr1cs1.r1cs.B, cr1cs1.r1cs.C],
            v: vec![lproof1.va, lproof1.vb, lproof1.vc],
            f: cr1cs1.f,
            params: params.clone(),
        };

        let A = Matrix::<R>::rand(&mut rand::thread_rng(), params.kappa, n);

        let mlin = Mlin {
            lins: vec![lin0, lin1],
            params,
            A,
        };

        let (_linb2, cmproof) = mlin.mlin(&mut ts);

        let mut ts = PoseidonTS::default::<PC>();
        lproof0.verify(&mut ts);
        lproof1.verify(&mut ts);
        cmproof.verify(&mut ts).unwrap();
    }
}
