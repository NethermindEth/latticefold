use ark_std::log2;
use stark_rings::{
    balanced_decomposition::{recompose, Decompose, DecomposeToVec},
    PolyRing, Zq,
};
use stark_rings_linalg::{ops::Transpose, Matrix, SparseMatrix};
use stark_rings_poly::mle::DenseMultilinearExtension;

use crate::lin::{LinB, LinBX};

#[derive(Clone, Debug)]
pub struct Decomp<R> {
    pub f: Vec<R>,
    pub r: Vec<(R, R)>,
    pub M: Vec<SparseMatrix<R>>,
}

#[derive(Clone, Debug)]
pub struct DecompProof<R> {
    /// C = com(F)
    pub C: (Vec<R>, Vec<R>), // kappa x 2
    pub v: (Vec<(R, R)>, Vec<(R, R)>), // (v(0), v(1))
}

impl<R: PolyRing> Decomp<R>
where
    R: Decompose,
    R::BaseRing: Zq,
{
    pub fn decompose(&self, A: &Matrix<R>, B: u128) -> ((LinB<R>, LinB<R>), DecompProof<R>) {
        let nvars = log2(A.ncols) as usize;
        let mut F = self.f.decompose_to_vec(B, 2).transpose().into_iter();
        let F0 = F.next().unwrap();
        let F1 = F.next().unwrap();

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
                        M_i.try_mul_vec(Fi).unwrap(),
                    )
                    .evaluate(&r_a)
                    .unwrap(),
                    DenseMultilinearExtension::from_evaluations_vec(
                        nvars,
                        M_i.try_mul_vec(Fi).unwrap(),
                    )
                    .evaluate(&r_b)
                    .unwrap(),
                );
                vi.push(vj);
            });
            vi
        };

        let v0 = vi_calc(&F0);
        let v1 = vi_calc(&F1);

        let C0 = A.try_mul_vec(&F0).unwrap();
        let C1 = A.try_mul_vec(&F1).unwrap();

        let linb0 = LinB {
            x: LinBX {
                cm_f: C0.clone(),
                r: self.r.clone(),
                v: v0.clone(),
            },
            f: F0,
        };
        let linb1 = LinB {
            x: LinBX {
                cm_f: C1.clone(),
                r: self.r.clone(),
                v: v1.clone(),
            },
            f: F1,
        };
        let proof = DecompProof {
            C: (C0, C1),
            v: (v0, v1),
        };

        ((linb0, linb1), proof)
    }
}

impl<R: PolyRing> DecompProof<R> {
    pub fn verify(&self, cm_f: &[R], v: &[(R, R)], B: u128) {
        let Br = R::from(B);
        let rec_cm = self
            .C
            .0
            .iter()
            .zip(self.C.1.iter())
            .map(|(&r0, &r1)| recompose(&[r0, r1], Br))
            .collect::<Vec<R>>();

        let rec_v = self
            .v
            .0
            .iter()
            .zip(self.v.1.iter())
            .map(|(v0, v1)| (recompose(&[v0.0, v1.0], Br), recompose(&[v0.1, v1.1], Br)))
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
    use latticefold::arith::r1cs::R1CS;
    use stark_rings::{
        balanced_decomposition::GadgetDecompose, cyclotomic_ring::models::frog_ring::RqPoly as R,
    };
    use stark_rings_linalg::SparseMatrix;

    use super::*;
    use crate::{
        lin::{LinParameters, Linearize, Verify},
        mlin::Mlin,
        r1cs::{r1cs_decomposed_square, ComR1CS},
        rgchk::DecompParameters,
        transcript::PoseidonTranscript,
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

        let (mut r1cs, z) = identity_cs(n / k);
        r1cs.A.coeffs[0][0].0 = 2u128.into();
        r1cs.C.coeffs[0][0].0 = 2u128.into();
        let r1cs = r1cs_decomposed_square(r1cs, n, 2, k);

        let A = Matrix::<R>::rand(&mut rand::thread_rng(), kappa, n);

        let cr1cs = ComR1CS::new(r1cs, z, 1, 2, k, &A);

        let mut ts = PoseidonTranscript::empty::<PC>();
        let (linb, lproof) = cr1cs.linearize(&mut ts);

        let mut ts = PoseidonTranscript::empty::<PC>();
        lproof.verify(&mut ts);

        let r = lproof
            .r
            .iter()
            .map(|r| (r.clone(), r.clone()))
            .collect::<Vec<_>>();

        let decomp = Decomp {
            f: cr1cs.f,
            r,
            M: cr1cs.x.matrices(),
        };

        let ((_linb0, _linb1), proof) = decomp.decompose(&A, B);

        proof.verify(&cr1cs.x.cm_f, &linb.x.v, B);
    }

    #[test]
    fn test_decomp_g() {
        let B = (<<R as PolyRing>::BaseRing>::MODULUS.0[0] as f64)
            .sqrt()
            .ceil() as u128
            + 1;
        let n = 1 << 15;
        let k = 2;
        let kappa = 2;
        let b = (R::dimension() / 2) as u128;
        // log_d' (q)
        let l = ((<<R as PolyRing>::BaseRing>::MODULUS.0[0] as f64).ln()
            / ((R::dimension() / 2) as f64).ln())
        .ceil() as usize;

        let params = LinParameters {
            kappa,
            decomp: DecompParameters { b, k, l },
        };

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

        let A = Matrix::<R>::rand(&mut rand::thread_rng(), params.kappa, n);

        let cr1cs0 = ComR1CS::new(r1cs.clone(), z0, 1, B, k, &A);
        let cr1cs1 = ComR1CS::new(r1cs, z1, 1, B, k, &A);

        let mut ts = PoseidonTranscript::empty::<PC>();
        let (linb0, lproof0) = cr1cs0.linearize(&mut ts);
        let (linb1, lproof1) = cr1cs1.linearize(&mut ts);

        let M = cr1cs0.x.matrices();

        let mlin = Mlin {
            lins: vec![linb0, linb1],
            params,
        };

        let (linb2, cmproof) = mlin.mlin(&A, &M, &mut ts);

        let mut ts = PoseidonTranscript::empty::<PC>();
        lproof0.verify(&mut ts);
        lproof1.verify(&mut ts);
        cmproof.verify(&M, &mut ts).unwrap();

        let decomp = Decomp {
            f: linb2.g,
            r: linb2.x.ro,
            M,
        };

        let (_linb, proof) = decomp.decompose(&A, B);

        proof.verify(&linb2.x.cm_g, &linb2.x.vo, B);
    }
}
