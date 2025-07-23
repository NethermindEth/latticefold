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
        convertible_ring::ConvertibleRing, recompose, Decompose, DecomposeToVec, GadgetDecompose,
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
pub struct Decomp<R> {
    pub f: Vec<R>,
    pub B: u128,
    pub r: Vec<(R, R)>,
    pub A: Matrix<R>,
    pub M: Vec<SparseMatrix<R>>,
}

#[derive(Clone, Debug)]
pub struct DecompProof<R> {
    /// C = com(F)
    pub C: (Vec<R>, Vec<R>), // kappa x 2
    pub B: u128,
    pub v: (Vec<(R, R)>, Vec<(R, R)>), // (v(0), v(1))
}

//#[derive(Clone, Debug)]
//pub struct LinBX<R> {
//    pub C: Vec<R>,
//    pub r: Vec<R>,
//    pub v:
//}

impl<R: PolyRing> Decomp<R>
where
    R: Decompose,
{
    pub fn decompose(&self) -> DecompProof<R> {
        let nvars = log2(self.A.ncols) as usize;
        let F = self.f.decompose_to_vec(self.B, 2).transpose();

        let vi_calc = |Fi: &[R]| -> Vec<(R, R)> {
            let r_a = self.r.iter().map(|rr| rr.0).collect::<Vec<_>>();
            let r_b = self.r.iter().map(|rr| rr.1).collect::<Vec<_>>();
            let fv = (
                DenseMultilinearExtension::from_evaluations_vec(nvars, Fi.to_vec())
                    .evaluate(&r_a)
                    .unwrap(),
                DenseMultilinearExtension::from_evaluations_vec(nvars, Fi.to_vec())
                    .evaluate(&r_b)
                    .unwrap(),
            );
            let mut vi = vec![fv];
            self.M.iter().for_each(|M_i| {
                let vj = (
                    DenseMultilinearExtension::from_evaluations_vec(
                        nvars,
                        M_i.try_mul_vec(&Fi).unwrap(),
                    )
                    .evaluate(&r_a)
                    .unwrap(),
                    DenseMultilinearExtension::from_evaluations_vec(
                        nvars,
                        M_i.try_mul_vec(&Fi).unwrap(),
                    )
                    .evaluate(&r_b)
                    .unwrap(),
                );
                vi.push(vj);
            });
            vi
        };

        let v0 = vi_calc(&F[0]);
        let v1 = vi_calc(&F[1]);

        let C = (
            self.A.try_mul_vec(&F[0]).unwrap(),
            self.A.try_mul_vec(&F[1]).unwrap(),
        );

        DecompProof {
            C,
            B: self.B,
            v: (v0, v1),
        }
    }
}

impl<R: PolyRing> DecompProof<R> {
    pub fn verify(&self, cm_f: &[R], v: &[(R, R)]) {
        let rec_cm = self
            .C
            .0
            .iter()
            .zip(self.C.1.iter())
            .map(|(&r0, &r1)| recompose(&[r0, r1], R::from(self.B)))
            .collect::<Vec<R>>();

        let rec_v = self
            .v
            .0
            .iter()
            .zip(self.v.1.iter())
            .map(|(v0, v1)| {
                (
                    recompose(&[v0.0, v1.0], R::from(self.B)),
                    recompose(&[v0.1, v1.1], R::from(self.B)),
                )
            })
            .collect::<Vec<(R, R)>>();

        assert_eq!(rec_cm, cm_f);
        assert_eq!(rec_v, v);
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::PrimeField;
    use ark_std::One;
    use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
    use latticefold::{arith::r1cs::R1CS, transcript::poseidon::PoseidonTS};
    use stark_rings::{
        balanced_decomposition::GadgetDecompose, cyclotomic_ring::models::frog_ring::RqPoly as R,
    };
    use stark_rings_linalg::SparseMatrix;

    use super::*;
    use crate::{
        lin::{Lin, Linearize, Parameters, Verify},
        mlin::Mlin,
        r1cs::CommittedR1CS,
    };

    fn identity_cs(n: usize) -> (R1CS<R>, Vec<R>) {
        let r1cs = R1CS::<R> {
            l: 1,
            A: SparseMatrix::identity(n),
            B: SparseMatrix::identity(n),
            C: SparseMatrix::identity(n),
        };
        let z = vec![R::one(); n];
        (r1cs, z)
    }

    #[test]
    fn test_decomp_r1cs() {
        let B = 50u128;
        let kappa = 2;
        let n = 1 << 15;
        let k = 4;

        let (mut r1cs, z) = identity_cs(n / 4);

        r1cs.A.coeffs[0][0].0 = 2u128.into();
        r1cs.C.coeffs[0][0].0 = 2u128.into();

        r1cs.A = r1cs.A.gadget_decompose(2, k);
        r1cs.B = r1cs.B.gadget_decompose(2, k);
        r1cs.C = r1cs.C.gadget_decompose(2, k);
        r1cs.A.pad_rows(n);
        r1cs.B.pad_rows(n);
        r1cs.C.pad_rows(n);

        let f = z.gadget_decompose(2, k);
        let A = Matrix::<R>::rand(&mut rand::thread_rng(), kappa, n);
        let cm_f = A.try_mul_vec(&f).unwrap();

        r1cs.check_relation(&f).unwrap();

        let cr1cs = CommittedR1CS {
            r1cs,
            f,
            x: vec![1u128.into()],
            cm: vec![],
        };

        let mut ts = PoseidonTS::default::<PC>();
        let lin = cr1cs.linearize(&mut ts);

        let mut ts = PoseidonTS::default::<PC>();
        lin.verify(&mut ts);

        let r = lin
            .r
            .iter()
            .map(|r| (r.clone(), r.clone()))
            .collect::<Vec<_>>();

        let decomp = Decomp {
            f: cr1cs.f,
            B,
            r,
            A,
            M: vec![cr1cs.r1cs.A, cr1cs.r1cs.B, cr1cs.r1cs.C],
        };

        let proof = decomp.decompose();

        let v = vec![
            (lin.v.clone(), lin.v),
            (lin.va.clone(), lin.va),
            (lin.vb.clone(), lin.vb),
            (lin.vc.clone(), lin.vc),
        ];
        proof.verify(&cm_f, &v);
    }

    #[test]
    fn test_decomp_g() {
        let B = (<<R as PolyRing>::BaseRing>::MODULUS.0[0] as f64)
            .sqrt()
            .ceil() as u128
            + 1;
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

        let M = vec![
            cr1cs0.r1cs.A.clone(),
            cr1cs0.r1cs.B.clone(),
            cr1cs0.r1cs.C.clone(),
        ];

        let lin0 = Lin {
            M: M.clone(),
            v: vec![lproof0.va, lproof0.vb, lproof0.vc],
            f: cr1cs0.f.clone(),
            params: params.clone(),
        };

        let lin1 = Lin {
            M: M.clone(),
            v: vec![lproof1.va, lproof1.vb, lproof1.vc],
            f: cr1cs1.f,
            params: params.clone(),
        };

        let A = Matrix::<R>::rand(&mut rand::thread_rng(), params.kappa, n);

        let mlin = Mlin {
            lins: vec![lin0, lin1],
            params,
            A: A.clone(),
        };

        let (linb2, cmproof) = mlin.mlin(&mut ts);

        let mut ts = PoseidonTS::default::<PC>();
        lproof0.verify(&mut ts);
        lproof1.verify(&mut ts);
        cmproof.verify(&mut ts).unwrap();

        let decomp = Decomp {
            f: linb2.g,
            B,
            r: linb2.x.ro,
            A,
            M,
        };

        let proof = decomp.decompose();

        proof.verify(&linb2.x.cm_g, &linb2.x.vo);
    }
}
