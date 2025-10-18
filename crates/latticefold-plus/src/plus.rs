use latticefold::transcript::Transcript;
use stark_rings::{
    balanced_decomposition::{convertible_ring::ConvertibleRing, Decompose},
    CoeffRing, OverField, Zq,
};
use stark_rings_linalg::{Matrix, SparseMatrix};

use crate::{
    cm::CmProof,
    decomp::{Decomp, DecompProof},
    lin::{LinParameters, Linearize, LinearizedVerify},
    mlin::{LinB2X, Mlin},
};

#[derive(Clone, Debug)]
pub struct PlusProver<R: OverField, TS: Transcript<R>> {
    pub acc: Mlin<R>,
    /// Ajtai matrix
    pub A: Matrix<R>,
    pub M: Vec<SparseMatrix<R>>,
    pub transcript: TS,
    pub params: PlusParameters,
}

#[derive(Clone, Debug)]
pub struct PlusVerifier<R: OverField, TS: Transcript<R>> {
    /// Ajtai matrix
    pub A: Matrix<R>,
    pub M: Vec<SparseMatrix<R>>,
    pub transcript: TS,
    pub params: PlusParameters,
}

#[derive(Clone, Debug)]
pub struct PlusProof<R: OverField, P: LinearizedVerify<R>> {
    pub linb2x: LinB2X<R>,
    pub lproof: Vec<P>,
    pub cmproof: CmProof<R>,
    pub dproof: DecompProof<R>,
}

#[derive(Clone, Debug)]
pub struct PlusParameters {
    pub lin: LinParameters,
    pub B: u128,
}

impl<R, TS> PlusProver<R, TS>
where
    R::BaseRing: ConvertibleRing + Decompose + Zq,
    R: CoeffRing + Decompose,
    TS: Transcript<R>,
{
    /// Initialize
    pub fn init(
        A: Matrix<R>,
        M: Vec<SparseMatrix<R>>,
        ncomp: usize,
        params: PlusParameters,
        transcript: TS,
    ) -> Self {
        let mlin = Mlin {
            lins: Vec::with_capacity(2 + ncomp),
            params: params.lin.clone(),
        };

        PlusProver {
            acc: mlin,
            A,
            M,
            transcript,
            params,
        }
    }

    /// Prove
    pub fn prove<L>(&mut self, comp: &[L]) -> PlusProof<R, L::Proof>
    where
        L: Linearize<R>,
    {
        let mut lproof = Vec::with_capacity(comp.len());
        comp.iter().for_each(|compi| {
            let (linb, lp) = compi.linearize(&mut self.transcript);
            lproof.push(lp);
            self.acc.lins.push(linb);
        });
        let (linb2, cmproof) = self.acc.mlin(&self.A, &self.M, &mut self.transcript);
        let decomp = Decomp {
            f: linb2.g,
            r: linb2.x.ro.clone(),
            M: self.M.clone(),
        };
        let (linb, dproof) = decomp.decompose(&self.A, self.params.B);

        let proof = PlusProof {
            linb2x: linb2.x,
            lproof,
            cmproof,
            dproof,
        };

        // Keep only accumulated instance
        self.acc.lins.clear();
        self.acc.lins.push(linb.0);
        self.acc.lins.push(linb.1);

        proof
    }
}

impl<R, TS> PlusVerifier<R, TS>
where
    R::BaseRing: Zq,
    R: CoeffRing,
    TS: Transcript<R>,
{
    /// Initialize
    pub fn init(
        A: Matrix<R>,
        M: Vec<SparseMatrix<R>>,
        params: PlusParameters,
        transcript: TS,
    ) -> Self {
        PlusVerifier {
            A,
            M,
            transcript,
            params,
        }
    }

    /// Verify
    pub fn verify<P: LinearizedVerify<R>>(&mut self, proof: &PlusProof<R, P>) -> bool {
        for lp in &proof.lproof {
            lp.verify(&mut self.transcript);
        }
        proof.cmproof.verify(&self.M, &mut self.transcript).unwrap();
        proof
            .dproof
            .verify(&proof.linb2x.cm_g, &proof.linb2x.vo, self.params.B);
        true
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::PrimeField;
    use latticefold::arith::r1cs::R1CS;
    use latticefold_rings::rings::FrogPoseidonConfig as PC;
    use rand::prelude::*;
    use stark_rings::{cyclotomic_ring::models::frog_ring::RqPoly as R, PolyRing, Ring};
    use stark_rings_linalg::SparseMatrix;

    use super::*;
    use crate::{
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

        let transcript = PoseidonTranscript::empty::<PC>();

        let pparams = PlusParameters { lin: params, B };
        let mut prover = PlusProver::init(A.clone(), M.clone(), 1, pparams.clone(), transcript);

        let proof = prover.prove(&[cr1cs0, cr1cs1]);

        let transcript = PoseidonTranscript::empty::<PC>();
        let mut verifier = PlusVerifier::init(A, M, pparams, transcript);
        verifier.verify(&proof);
    }

    #[test]
    fn test_prove_multi() {
        let n = 1 << 16;
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

        let ts = PoseidonTranscript::empty::<PC>();

        let pparams = PlusParameters { lin: params, B };
        let mut prover = PlusProver::init(A.clone(), M.clone(), 1, pparams.clone(), ts);

        let ts_v = PoseidonTranscript::empty::<PC>();
        let mut verifier = PlusVerifier::init(A, M, pparams, ts_v);

        for _ in 0..3 {
            let proof = prover.prove(std::slice::from_ref(&cr1cs));
            verifier.verify(&proof);
        }
    }
}
