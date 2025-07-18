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

#[derive(Debug)]
pub struct Decomp<R> {
    pub f: Vec<R>,
    pub B: u128,
    pub r: Vec<R>,
    pub A: Matrix<R>,
    pub M: Vec<SparseMatrix<R>>,
}

#[derive(Debug)]
pub struct DecompProof<R> {
    /// C = com(F)
    pub C: (Vec<R>, Vec<R>), // kappa x 2
    pub B: u128,
    pub M: Vec<SparseMatrix<R>>,
    pub vs: Vec<(R, R)>,
}

impl<R: PolyRing> Decomp<R>
where
    R: Decompose,
{
    pub fn decompose(&self) -> DecompProof<R> {
        let nvars = log2(self.A.ncols) as usize;
        let F = self.f.decompose_to_vec(self.B, 2).transpose();

        let vs = self
            .M
            .iter()
            .map(|M_i| {
                (
                    DenseMultilinearExtension::from_evaluations_vec(
                        nvars,
                        M_i.try_mul_vec(&F[0]).unwrap(),
                    )
                    .evaluate(&self.r)
                    .unwrap(),
                    DenseMultilinearExtension::from_evaluations_vec(
                        nvars,
                        M_i.try_mul_vec(&F[1]).unwrap(),
                    )
                    .evaluate(&self.r)
                    .unwrap(),
                )
            })
            .collect::<Vec<_>>();
        let C = (
            self.A.try_mul_vec(&F[0]).unwrap(),
            self.A.try_mul_vec(&F[1]).unwrap(),
        );

        DecompProof {
            C,
            B: self.B,
            M: self.M.clone(),
            vs,
        }
    }
}

impl<R: PolyRing> DecompProof<R> {
    pub fn verify(&self, cm_f: &[R], v: &[R]) {
        let rec_cm = self
            .C
            .0
            .iter()
            .zip(self.C.1.iter())
            .map(|(&r0, &r1)| recompose(&[r0, r1], R::from(self.B)))
            .collect::<Vec<R>>();

        let rec_v = self
            .vs
            .iter()
            .map(|&(v0, v1)| recompose(&[v0, v1], R::from(self.B)))
            .collect::<Vec<R>>();

        assert_eq!(rec_cm, cm_f);
        assert_eq!(rec_v, v);
    }
}

#[cfg(test)]
mod tests {
    use ark_std::One;
    use cyclotomic_rings::rings::GoldilocksPoseidonConfig as PC;
    use latticefold::{arith::r1cs::R1CS, transcript::poseidon::PoseidonTS};
    use stark_rings::{
        balanced_decomposition::GadgetDecompose, cyclotomic_ring::models::goldilocks::RqPoly as R,
    };
    use stark_rings_linalg::SparseMatrix;

    use super::*;
    use crate::{
        lin::{Linearize, Verify},
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
    fn test_decomp() {
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

        let decomp = Decomp {
            f: cr1cs.f,
            B,
            r: lin.r,
            A,
            M: vec![cr1cs.r1cs.A, cr1cs.r1cs.B, cr1cs.r1cs.C],
        };

        let proof = decomp.decompose();

        let v = vec![lin.va, lin.vb, lin.vc];
        proof.verify(&cm_f, &v);
    }
}
