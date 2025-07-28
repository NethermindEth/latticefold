use ark_std::log2;
use latticefold::transcript::Transcript;
use stark_rings::{
    balanced_decomposition::{convertible_ring::ConvertibleRing, Decompose},
    CoeffRing, Zq,
};
use stark_rings_linalg::{Matrix, SparseMatrix};

use crate::{
    cm::{Cm, CmComs, CmProof},
    lin::{LinB, LinParameters},
    rgchk::{Rg, RgInstance},
};

#[derive(Clone, Debug)]
pub struct Mlin<R> {
    pub lins: Vec<LinB<R>>,
    pub params: LinParameters,
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
    ///
    /// Folds L `LinB` instances.
    pub fn mlin(
        &self,
        A: &Matrix<R>,
        M: &[SparseMatrix<R>],
        transcript: &mut impl Transcript<R>,
    ) -> (LinB2<R>, CmProof<R>) {
        let n = self.lins[0].f.len();

        let instances = self
            .lins
            .iter()
            .map(|lin| RgInstance::from_f(lin.f.clone(), A, &self.params.decomp))
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
            dparams: self.params.decomp.clone(),
        };

        let cm = Cm { rg, coms };

        let (com, proof) = cm.prove(M, transcript);

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
    use ark_std::One;
    use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
    use latticefold::arith::r1cs::R1CS;
    use stark_rings::{cyclotomic_ring::models::frog_ring::RqPoly as R, PolyRing};
    use stark_rings_linalg::SparseMatrix;

    use super::*;
    use crate::{
        lin::{Linearize, LinearizedVerify},
        r1cs::{r1cs_decomposed_square, ComR1CS},
        rgchk::DecompParameters,
        transcript::PoseidonTranscript,
    };

    #[test]
    fn test_mlin() {
        let n = 1 << 15;
        let k = 2;
        let B = 2;
        let b = (R::dimension() / 2) as u128;
        let kappa = 2;
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

        let mut r1cs = r1cs_decomposed_square(
            R1CS::<R> {
                l: 1,
                A: SparseMatrix::identity(n / k),
                B: SparseMatrix::identity(n / k),
                C: SparseMatrix::identity(n / k),
            },
            n,
            B,
            k,
        );

        r1cs.A.coeffs[0][0].0 = 2u128.into();
        r1cs.C.coeffs[0][0].0 = 2u128.into();

        let A = Matrix::<R>::rand(&mut ark_std::test_rng(), params.kappa, n);

        let cr1cs0 = ComR1CS::new(r1cs.clone(), z0, 1, B, k, &A);
        let cr1cs1 = ComR1CS::new(r1cs, z1, 1, B, k, &A);

        let mut ts = PoseidonTranscript::empty::<PC>();
        let (linb0, lproof0) = cr1cs0.linearize(&mut ts);
        let (linb1, lproof1) = cr1cs1.linearize(&mut ts);

        let M = cr1cs0.x.matrices();

        let A = Matrix::<R>::rand(&mut ark_std::test_rng(), params.kappa, n);

        let mlin = Mlin {
            lins: vec![linb0, linb1],
            params,
        };

        let (_linb2, cmproof) = mlin.mlin(&A, &M, &mut ts);

        let mut ts = PoseidonTranscript::empty::<PC>();
        lproof0.verify(&mut ts);
        lproof1.verify(&mut ts);
        cmproof.verify(&M, &mut ts).unwrap();
    }
}
