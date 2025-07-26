use latticefold::transcript::Transcript;
use stark_rings::{
    balanced_decomposition::{convertible_ring::ConvertibleRing, Decompose},
    CoeffRing, PolyRing, Zq,
};
use stark_rings_linalg::{Matrix, SparseMatrix};

use crate::{
    cm::CmProof,
    decomp::{Decomp, DecompProof},
    lin::{LinB, LinParameters},
    mlin::{LinB2X, Mlin},
};

#[derive(Clone, Debug)]
pub struct PlusProver<R> {
    pub instances: Mlin<R>,
    /// Ajtai matrix
    pub A: Matrix<R>,
    pub params: PlusParameters,
}

#[derive(Clone, Debug)]
pub struct PlusVerifier<R> {
    /// Ajtai matrix
    pub A: Matrix<R>,
    pub params: PlusParameters,
}

#[derive(Clone, Debug)]
pub struct PlusProof<R: PolyRing> {
    pub linb2x: LinB2X<R>,
    pub cmproof: CmProof<R>,
    pub dproof: DecompProof<R>,
}

#[derive(Clone, Debug)]
pub struct PlusParameters {
    pub lin: LinParameters,
    pub B: u128,
}

impl<R: CoeffRing> PlusProver<R>
where
    R::BaseRing: ConvertibleRing + Decompose + Zq,
    R: Decompose,
{
    /// Prove
    pub fn prove(
        &self,
        M: &[SparseMatrix<R>],
        transcript: &mut impl Transcript<R>,
    ) -> ((LinB<R>, LinB<R>), PlusProof<R>) {
        let (linb2, cmproof) = self.instances.mlin(&self.A, M, transcript);
        let decomp = Decomp {
            f: linb2.g,
            r: linb2.x.ro.clone(),
            M: M.to_vec(),
        };
        let (linb, dproof) = decomp.decompose(&self.A, self.params.B);

        let proof = PlusProof {
            linb2x: linb2.x,
            cmproof,
            dproof,
        };

        (linb, proof)
    }
}

impl<R: CoeffRing> PlusVerifier<R>
where
    R::BaseRing: Zq,
{
    /// Verify
    pub fn verify(
        &self,
        proof: &PlusProof<R>,
        M: &[SparseMatrix<R>],
        transcript: &mut impl Transcript<R>,
    ) {
        proof.cmproof.verify(M, transcript).unwrap();
        proof
            .dproof
            .verify(&proof.linb2x.cm_g, &proof.linb2x.vo, self.params.B);
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::PrimeField;
    use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
    use latticefold::arith::r1cs::R1CS;
    use rand::prelude::*;
    use stark_rings::{cyclotomic_ring::models::frog_ring::RqPoly as R, Ring};
    use stark_rings_linalg::SparseMatrix;

    use super::*;
    use crate::{
        lin::{Linearize, Verify},
        r1cs::{r1cs_decomposed_square, ComR1CS},
        rgchk::DecompParameters,
        transcript::PoseidonTranscript,
        utils::estimate_bound,
    };

    #[test]
    fn test_prove() {
        let n = 1 << 15;
        let sop = R::dimension() * 128; // S inf-norm = 128
        let L = 3;
        let k = 2;
        let d = R::dimension();
        let b = (R::dimension() / 2) as u128;
        let B = estimate_bound(sop, L, d, k) + 1;
        let m = n / k;
        let kappa = 2;
        // log_d' (q)
        let l = ((<<R as PolyRing>::BaseRing>::MODULUS.0[0] as f64).ln()
            / ((R::dimension() / 2) as f64).ln())
        .ceil() as usize;
        let params = LinParameters {
            kappa,
            decomp: DecompParameters { b, k, l },
        };

        let mut rng = ark_std::test_rng();
        let pop = [R::ZERO, R::ONE];
        let z0: Vec<R> = (0..m).map(|_| *pop.choose(&mut rng).unwrap()).collect();
        let z1: Vec<R> = (0..m).map(|_| *pop.choose(&mut rng).unwrap()).collect();

        let r1cs = r1cs_decomposed_square(
            R1CS::<R> {
                l: 1,
                A: SparseMatrix::identity(m),
                B: SparseMatrix::identity(m),
                C: SparseMatrix::identity(m),
            },
            n,
            B,
            k,
        );

        let A = Matrix::<R>::rand(&mut ark_std::test_rng(), params.kappa, n);

        let cr1cs0 = ComR1CS::new(r1cs.clone(), z0, 1, B, k, &A);
        let cr1cs1 = ComR1CS::new(r1cs, z1, 1, B, k, &A);

        let M = cr1cs0.x.matrices();

        let mut ts = PoseidonTranscript::empty::<PC>();
        let (linb0, lproof0) = cr1cs0.linearize(&mut ts);
        let (linb1, lproof1) = cr1cs1.linearize(&mut ts);

        let mlin = Mlin {
            lins: vec![linb0, linb1],
            params: params.clone(),
        };

        let prover = PlusProver {
            instances: mlin,
            A: A.clone(),
            params: PlusParameters {
                lin: params.clone(),
                B,
            },
        };

        let (_acc, proof) = prover.prove(&M, &mut ts);

        let verifier = PlusVerifier {
            A,
            params: PlusParameters { lin: params, B },
        };
        let mut ts = PoseidonTranscript::empty::<PC>();
        lproof0.verify(&mut ts);
        lproof1.verify(&mut ts);
        verifier.verify(&proof, &M, &mut ts);
    }

    #[test]
    fn test_prove_multi() {
        let n = 1 << 17;
        let sop = R::dimension() * 128; // S inf-norm = 128
        let L = 3;
        let k = 4;
        let d = R::dimension();
        let b = (R::dimension() / 2) as u128;
        let B = estimate_bound(sop, L, d, k) / 2; // + 1;
        let m = n / k;
        let kappa = 2;
        // log_d' (q)
        let l = ((<<R as PolyRing>::BaseRing>::MODULUS.0[0] as f64).ln()
            / ((R::dimension() / 2) as f64).ln())
        .ceil() as usize;
        let params = LinParameters {
            kappa,
            decomp: DecompParameters { b, k, l },
        };

        let mut rng = ark_std::test_rng();
        let pop = [R::ZERO, R::ONE];
        let z: Vec<R> = (0..m).map(|_| *pop.choose(&mut rng).unwrap()).collect();

        let r1cs = r1cs_decomposed_square(
            R1CS::<R> {
                l: 1,
                A: SparseMatrix::identity(m),
                B: SparseMatrix::identity(m),
                C: SparseMatrix::identity(m),
            },
            n,
            B,
            k,
        );

        let A = Matrix::<R>::rand(&mut ark_std::test_rng(), params.kappa, n);

        let cr1cs = ComR1CS::new(r1cs, z, 1, B, k, &A);

        let M = cr1cs.x.matrices();

        let mut ts = PoseidonTranscript::empty::<PC>();
        let (mut linb, mut lproof) = cr1cs.linearize(&mut ts);

        let mlin = Mlin {
            lins: vec![linb],
            params: params.clone(),
        };

        let mut prover = PlusProver {
            instances: mlin,
            A: A.clone(),
            params: PlusParameters {
                lin: params.clone(),
                B,
            },
        };

        let verifier = PlusVerifier {
            A: A.clone(),
            params: PlusParameters {
                lin: params.clone(),
                B,
            },
        };

        let mut ts_v = PoseidonTranscript::empty::<PC>();

        for _ in 0..3 {
            let (acc, proof) = prover.prove(&M, &mut ts);
            lproof.verify(&mut ts_v);
            verifier.verify(&proof, &M, &mut ts_v);

            (linb, lproof) = cr1cs.linearize(&mut ts);
            let mlin = Mlin {
                lins: vec![acc.0, acc.1, linb],
                params: params.clone(),
            };
            prover.instances = mlin;
        }
    }
}
