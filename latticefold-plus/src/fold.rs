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

#[derive(Debug)]
pub struct Mlin<R> {
    pub lins: Vec<Lin<R>>,
    pub params: Parameters,
}

impl<R: CoeffRing> Mlin<R>
where
    R::BaseRing: ConvertibleRing + Decompose + Zq,
    R: Decompose,
{
    /// Πmlin protocol
    pub fn mlin(&self, transcript: &mut impl Transcript<R>) -> CmProof<R> {
        let n = self.lins[0].f.len();

        let A = Matrix::<R>::rand(&mut rand::thread_rng(), self.params.kappa, n);

        let instances = self
            .lins
            .iter()
            .map(|lin| RgInstance::from_f(lin.f.clone(), &A, &self.params.decomp))
            .collect::<Vec<_>>();

        let coms = instances
            .iter()
            .map(|inst| {
                let cm_f = A.try_mul_vec(&inst.f).unwrap();
                let C_Mf = A
                    .try_mul_vec(&inst.tau.iter().map(|z| R::from(*z)).collect::<Vec<R>>())
                    .unwrap();
                let cm_mtau = A.try_mul_vec(&inst.m_tau).unwrap();

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

        cm.prove(transcript)
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

        let f0 = z.gadget_decompose(2, 4);
        let mut f1 = z.gadget_decompose(2, 4);
        f1[0] = R::from(0u128);
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
        let k = 2;
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

        let mlin = Mlin {
            lins: vec![lin0, lin1],
            params,
        };

        let cmproof = mlin.mlin(&mut ts);

        let mut ts = PoseidonTS::default::<PC>();
        lproof0.verify(&mut ts);
        lproof1.verify(&mut ts);
        cmproof.verify(&mut ts).unwrap();
    }
}
