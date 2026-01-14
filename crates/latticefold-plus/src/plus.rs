use ark_ff::{BigInteger, PrimeField};
use latticefold::transcript::Transcript;
use latticefold::commitment::AjtaiCommitmentScheme;
use stark_rings::{
    balanced_decomposition::{convertible_ring::ConvertibleRing, Decompose},
    CoeffRing, OverField, PolyRing, Zq,
};
use stark_rings_linalg::{Matrix, SparseMatrix};
use std::sync::Arc;

use crate::{
    cm::CmProof,
    decomp::DecompProof,
    lin::{LinParameters, Linearize, LinearizedVerify},
    mlin::{LinB2X, Mlin},
    utils::maybe_print_rss,
};

#[derive(Clone, Debug)]
pub struct PlusProver<R: OverField, TS: Transcript<R>> {
    pub acc: Mlin<R>,
    /// Ajtai matrix
    pub A: Matrix<R>,
    pub M: Vec<Arc<SparseMatrix<R>>>,
    pub transcript: TS,
    pub params: PlusParameters,
}


/// Prover variant for the const-coeff/SP1 regime where the external matrices are stored over the base ring.
#[derive(Clone, Debug)]
pub struct PlusProverSparseBase<R: OverField + PolyRing, TS: Transcript<R>> {
    pub acc: Mlin<R>,
    pub scheme: AjtaiCommitmentScheme<R>,
    pub M0: Vec<Arc<SparseMatrix<R::BaseRing>>>,
    pub transcript: TS,
    pub params: PlusParameters,
}

#[derive(Clone, Debug)]
pub struct PlusVerifier<R: OverField, TS: Transcript<R>> {
    /// Ajtai matrix
    pub A: Matrix<R>,
    pub M: Vec<Arc<SparseMatrix<R>>>,
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

impl<R, TS> PlusProverSparseBase<R, TS>
where
    R::BaseRing: ConvertibleRing + Decompose + Zq,
    R: CoeffRing + Decompose,
    TS: Transcript<R>,
{
    pub fn init_seeded_base(
        scheme: AjtaiCommitmentScheme<R>,
        M0: Vec<Arc<SparseMatrix<R::BaseRing>>>,
        ncomp: usize,
        params: PlusParameters,
        transcript: TS,
    ) -> Self {
        let mlin = Mlin {
            lins: Vec::with_capacity(2 + ncomp),
            params: params.lin.clone(),
        };
        PlusProverSparseBase {
            acc: mlin,
            scheme,
            M0,
            transcript,
            params,
        }
    }

    pub fn prove_sparse_base<L>(&mut self, comp: &[L]) -> PlusProof<R, L::Proof>
    where
        L: Linearize<R>,
        R::BaseRing: PrimeField,
        <R::BaseRing as PrimeField>::BigInt: BigInteger,
    {
        maybe_print_rss("PlusProverSparseBase::prove_sparse_base (start)");
        let mut lproof = Vec::with_capacity(comp.len());
        comp.iter().for_each(|compi| {
            let (linb, lp) = compi.linearize(&mut self.transcript);
            lproof.push(lp);
            self.acc.lins.push(linb);
        });

        maybe_print_rss("PlusProverSparseBase::prove_sparse_base (after linearize)");
        let (linb2, cmproof) =
            self.acc
                .mlin_seeded_base(&self.scheme, &self.M0, &mut self.transcript);
        maybe_print_rss("PlusProverSparseBase::prove_sparse_base (after mlin_seeded)");

        let decomp = crate::decomp::DecompBase {
            f: linb2.g,
            r: linb2.x.ro.clone(),
            M0: &self.M0,
        };
        let dproof = decomp.decompose_seeded_base_one_shot(&self.scheme, self.params.B);
        maybe_print_rss("PlusProverSparseBase::prove_sparse_base (after decompose_seeded)");

        let proof = PlusProof {
            linb2x: linb2.x,
            lproof,
            cmproof,
            dproof,
        };

        // One-shot proving: don't retain large witnesses in the prover state.
        self.acc.lins.clear();
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
        M: Vec<Arc<SparseMatrix<R>>>,
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

    /// Get reference to transcript (for metrics after verification)
    pub fn transcript(&self) -> &TS {
        &self.transcript
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::PrimeField;
    use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
    use latticefold::arith::r1cs::R1CS;
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
        let n = 1 << 16; // Increased to accommodate l rounded up to power of 2
        let sop = R::dimension() * 128; // S inf-norm = 128
        let L = 3;
        let k = 2;
        let d = R::dimension();
        let b = (R::dimension() / 2) as u128;
        let B = estimate_bound(sop, L, d, k) + 1;
        let m = n / k;
        let kappa = 2;
        // log_d' (q) - round up to power of 2 for optimized tensor evaluation
        let l_raw = ((<<R as PolyRing>::BaseRing>::MODULUS.0[0] as f64).ln()
            / ((R::dimension() / 2) as f64).ln())
        .ceil() as usize;
        let l = l_raw.next_power_of_two();
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

        let M = cr1cs0.x.matrices_arc();

        let transcript = PoseidonTranscript::empty::<PC>();

        let pparams = PlusParameters { lin: params, B };
        let mut prover = PlusProver::init(A.clone(), M.clone(), 1, pparams.clone(), transcript);

        let proof = prover.prove(&[cr1cs0, cr1cs1]);

        // log_kappa for tensor status printing (kappa already defined above)
        let log_kappa = ark_std::log2(kappa) as usize;
        
        let transcript = PoseidonTranscript::empty::<PC>();
        let mut verifier = PlusVerifier::init(A, M, pparams, transcript);
        
        // Time verification
        let start = std::time::Instant::now();
        verifier.verify(&proof);
        let verify_time = start.elapsed();
        
        // Print transcript metrics for DPP cost estimation
        println!("\n=== LF+ Verifier Metrics (n={}) ===", n);
        println!("  Ring dimension d = {}", R::dimension());
        println!("  Decomposition k = {}, l = {} (padded to {})", k, l, l.next_power_of_two());
        println!("  Folding instances L = {}", L);
        println!("  Verification time: {:?}", verify_time);
        verifier.transcript().print_metrics();
        
        // Print tensor optimization status
        use crate::tensor_eval::print_tensor_optimization_status;
        print_tensor_optimization_status(
            log_kappa,
            k * R::dimension(),
            l.next_power_of_two(),
            R::dimension(),
        );
    }

    #[test]
    fn test_prove_multi() {
        let n = 1 << 17; // Increased to accommodate l rounded up to power of 2
        let sop = R::dimension() * 128; // S inf-norm = 128
        let L = 3;
        let k = 4;
        let d = R::dimension();
        let b = (R::dimension() / 2) as u128;
        let B = estimate_bound(sop, L, d, k) / 2; // + 1;
        let m = n / k;
        let kappa = 2;
        // log_d' (q) - round up to power of 2 for optimized tensor evaluation
        let l_raw = ((<<R as PolyRing>::BaseRing>::MODULUS.0[0] as f64).ln()
            / ((R::dimension() / 2) as f64).ln())
        .ceil() as usize;
        let l = l_raw.next_power_of_two();
        let params = LinParameters {
            kappa,
            decomp: DecompParameters { b, k, l },
        };

        let mut rng = ark_std::test_rng();
        let pop = [R::ZERO, R::ONE];
        // Used by downstream prove path; keep name underscore to avoid unused warnings in tests.
        let _z: Vec<R> = (0..m).map(|_| *pop.choose(&mut rng).unwrap()).collect();

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

        let cr1cs: Vec<_> = (0..4)
            .map(|_| {
                let z: Vec<R> = (0..m).map(|_| *pop.choose(&mut rng).unwrap()).collect();
                ComR1CS::new(r1cs.clone(), z, 1, B, k, &A)
            })
            .collect();

        let M = cr1cs[0].x.matrices_arc();

        let transcript = PoseidonTranscript::empty::<PC>();

        let pparams = PlusParameters { lin: params, B };
        let mut prover = PlusProver::init(A.clone(), M.clone(), 1, pparams.clone(), transcript);

        let proof = prover.prove(&cr1cs);

        let transcript = PoseidonTranscript::empty::<PC>();
        let mut verifier = PlusVerifier::init(A, M, pparams, transcript);
        verifier.verify(&proof);
    }

    /// Large-scale test to measure tensor optimization impact
    /// Run with: cargo test --release test_large_scale -- --nocapture --ignored
    #[test]
    #[ignore] // Only run manually due to long runtime
    fn test_large_scale() {
        use crate::tensor_eval::print_tensor_optimization_status;
        use std::time::Instant;
        
        let profile = std::env::var("LF_PLUS_PROFILE").ok().as_deref() == Some("1");
        if profile {
            #[cfg(feature = "parallel")]
            println!(
                "[LF+ test_large_scale] rayon_threads={}",
                rayon::current_num_threads()
            );
            #[cfg(not(feature = "parallel"))]
            println!("[LF+ test_large_scale] rayon_threads=DISABLED(feature=parallel)");
        }

        let n = 1 << 20; // 1M constraints - closer to SP1 scale
        let sop = R::dimension() * 128;
        let L = 3;
        let k = 4;
        let d = R::dimension();
        let b = (R::dimension() / 2) as u128;
        let B = estimate_bound(sop, L, d, k) / 2;
        let m = n / k;
        let kappa = 2;
        let l_raw = ((<<R as PolyRing>::BaseRing>::MODULUS.0[0] as f64).ln()
            / ((R::dimension() / 2) as f64).ln())
        .ceil() as usize;
        let l = l_raw.next_power_of_two();
        let log_kappa = ark_std::log2(kappa) as usize;
        
        println!("\n========== LARGE SCALE BENCHMARK (n={}) ==========", n);
        println!("Parameters: d={}, k={}, l={} (raw {}), kappa={}", d, k, l, l_raw, kappa);
        print_tensor_optimization_status(log_kappa, k * d, l, d);
        
        let params = LinParameters {
            kappa,
            decomp: DecompParameters { b, k, l },
        };

        let mut rng = ark_std::test_rng();
        let pop = [R::ZERO, R::ONE];

        let t = Instant::now();
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
        if profile {
            println!("[LF+ test_large_scale] build r1cs_decomposed_square: {:?}", t.elapsed());
        }

        let t = Instant::now();
        let A = Matrix::<R>::rand(&mut ark_std::test_rng(), params.kappa, n);
        if profile {
            println!("[LF+ test_large_scale] sample Ajtai A (kappa×n): {:?}", t.elapsed());
        }

        let t = Instant::now();
        let z: Vec<R> = (0..m).map(|_| *pop.choose(&mut rng).unwrap()).collect();
        if profile {
            println!("[LF+ test_large_scale] sample witness z (len=m): {:?}", t.elapsed());
        }

        let t = Instant::now();
        let cr1cs = ComR1CS::new(r1cs, z, 1, B, k, &A);
        if profile {
            println!("[LF+ test_large_scale] ComR1CS::new (decomp+commit): {:?}", t.elapsed());
        }

        let t = Instant::now();
        let M = cr1cs.x.matrices_arc();
        if profile {
            println!("[LF+ test_large_scale] extract M matrices: {:?}", t.elapsed());
        }

        let pparams = PlusParameters { lin: params, B };
        
        // Generate proof
        let ts = PoseidonTranscript::empty::<PC>();
        let mut prover = PlusProver::init(A.clone(), M.clone(), 1, pparams.clone(), ts);
        let t = Instant::now();
        let proof = prover.prove(std::slice::from_ref(&cr1cs));
        if profile {
            println!("[LF+ test_large_scale] prover.prove total: {:?}", t.elapsed());
        }
        
        // Verify with timing
        let ts = PoseidonTranscript::empty::<PC>();
        let mut verifier = PlusVerifier::init(A, M, pparams, ts);
        let start = std::time::Instant::now();
        verifier.verify(&proof);
        let verify_time = start.elapsed();
        
        println!("\n=== VERIFICATION BENCHMARK (n={}) ===", n);
        println!("  Verification time: {:?}", verify_time);
        println!("==========================================\n");
        
        verifier.transcript().print_metrics();
    }
}
