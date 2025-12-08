//! Benchmarks for Construction 5.1: Single instance folding protocol.
//!
//! The single instance folding protocol is a special case of multilinear folding
//! with L=1. It reduces one LinB instance (norm bound B) to one LinB2 instance
//! (norm bound B²). This serves as the baseline folding operation.
//!
//! ## Protocol Overview
//!
//! The single instance fold operates as follows:
//! 1. **Linearization**: Takes a committed R1CS instance and produces a LinB instance
//! 2. **Folding**: Applies the lin protocol to fold LinB → LinB2
//! 3. **Commitment transformation**: Uses the commitment transformation protocol
//!    internally for consistency
//!
//! This construction demonstrates the base case before batching optimizations.
//! Comparing with multilinear folding (L>1) shows the amortization benefits
//! of batching multiple instances.
//!
//! ## Benchmarked Operations
//!
//! - **Prover**: Single instance linearization and folding proof generation
//! - **Verifier**: Verification of the folded LinB2 commitment and constraints
//! - **Baseline comparison**: Performance baseline for measuring multilinear batching efficiency
//!
//! ## Paper Reference
//! Construction 5.1, Section 5.1 of the LatticeFold+ paper

#![allow(non_snake_case)]

use criterion::Criterion;
use latticefold_plus::{
    lin::{LinB, LinParameters, Linearize},
    mlin::LinB2,
    r1cs::ComR1CS,
};
use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as R;
use stark_rings_linalg::{Matrix, SparseMatrix};

#[path = "utils/mod.rs"]
mod utils;
use utils::{
    helpers::{
        bench_prover_protocol, bench_rng, bench_verifier_protocol, create_ajtai_matrix,
        create_test_m_matrix, create_transcript, get_validated_decomp_params, ProverBenchmark,
        R1CSBuilder, VerifierBenchmark, WitnessPattern,
    },
    single_instance_fold,
};

// ============================================================================
// Setup Functions
// ============================================================================

/// Creates single instance fold input for prover benchmarks.
///
/// Generates a committed R1CS instance, linearizes it to produce a LinB instance,
/// then prepares it for the lin folding protocol. Uses witness with all ones to
/// ensure valid R1CS satisfaction with norm bound B.
fn setup_input(n: usize, k: usize, kappa: usize, B: usize) -> LinB<R> {
    let mut rng = bench_rng();
    let dparams = get_validated_decomp_params(k, kappa, n);

    let r1cs = R1CSBuilder::new(n, k, B as u128).build_decomposed_square();
    let A = create_ajtai_matrix(kappa, n, &mut rng);
    let z = WitnessPattern::AllOnes.generate(n / k, &mut rng);

    let cr1cs = ComR1CS::new(r1cs, z, 1, B as u128, k, &A);

    let mut ts = create_transcript();
    let (linb, _lproof) = cr1cs.linearize(&mut ts);

    linb
}

/// Generates a valid single instance fold proof for verifier benchmarks.
///
/// Creates a LinB instance, executes the lin folding protocol to generate
/// a `CmProof`, and validates it before returning. This ensures the verifier
/// benchmarks measure only verification time, not error handling overhead.
fn setup_proof(
    n: usize,
    k: usize,
    kappa: usize,
    B: usize,
) -> (LinB<R>, latticefold_plus::cm::CmProof<R>) {
    let linb = setup_input(n, k, kappa, B);
    let mut rng = bench_rng();
    let dparams = get_validated_decomp_params(k, kappa, n);

    let params = LinParameters {
        kappa,
        decomp: dparams,
    };

    let M = create_test_m_matrix(n);
    let A = create_ajtai_matrix(kappa, n, &mut rng);

    let mut ts = create_transcript();
    let (linb2, proof) = linb.lin(&A, &M, &params, &mut ts);

    let mut verify_ts = create_transcript();
    proof
        .verify(&M, &mut verify_ts)
        .expect("Generated single instance folding proof should be valid");

    (linb, proof)
}

// ============================================================================
// Trait Implementations
// ============================================================================

/// Prover benchmark for single instance folding protocol.
///
/// Measures LinB → LinB2 folding time across varying witness sizes (32K-128K)
/// while keeping k=2, κ=2, and B=50 fixed. This is the baseline folding
/// operation before batching optimizations are applied.
struct SingleInstanceFoldProver;

impl ProverBenchmark for SingleInstanceFoldProver {
    type Input = (LinB<R>, Matrix<R>, Vec<SparseMatrix<R>>, LinParameters);
    type Output = (LinB2<R>, latticefold_plus::cm::CmProof<R>);
    type Params = (usize, usize, usize, usize);

    fn group_name() -> &'static str {
        "SingleInstanceFolding-Prover"
    }

    fn setup_input((n, k, kappa, B): Self::Params) -> Self::Input {
        let mut rng = bench_rng();
        let linb = setup_input(n, k, kappa, B);
        let M = create_test_m_matrix(n);
        let A = create_ajtai_matrix(kappa, n, &mut rng);
        let dparams = get_validated_decomp_params(k, kappa, n);
        let params = LinParameters {
            kappa,
            decomp: dparams,
        };

        (linb, A, M, params)
    }

    fn param_label((n, k, kappa, B): Self::Params) -> String {
        format!("n={}_k={}_κ={}_B={}", n, k, kappa, B)
    }

    fn throughput((n, _, _, _): Self::Params) -> u64 {
        n as u64
    }

    fn run_prover((linb, A, M, params): Self::Input) -> Self::Output {
        let mut ts = create_transcript();
        linb.lin(&A, &M, &params, &mut ts)
    }
}

/// Verifier benchmark for single instance folding protocol.
///
/// Measures verification time for single instance folding proofs. Uses
/// sumcheck verification to confirm the LinB → LinB2 transformation
/// preserves R1CS constraint satisfaction.
struct SingleInstanceFoldVerifier;

impl VerifierBenchmark for SingleInstanceFoldVerifier {
    type Input = LinB<R>;
    type Proof = latticefold_plus::cm::CmProof<R>;
    type Params = (usize, usize, usize, usize);

    fn group_name() -> &'static str {
        "SingleInstanceFolding-Verifier"
    }

    fn setup_proof((n, k, kappa, B): Self::Params) -> (Self::Input, Self::Proof) {
        setup_proof(n, k, kappa, B)
    }

    fn param_label((n, k, kappa, B): Self::Params) -> String {
        format!("n={}_k={}_κ={}_B={}", n, k, kappa, B)
    }

    fn throughput((n, _, _, _): Self::Params) -> u64 {
        n as u64
    }

    fn run_verifier(_input: &Self::Input, proof: &Self::Proof) {
        let n = 65536;
        let M = create_test_m_matrix(n);
        let mut ts = create_transcript();
        proof.verify(&M, &mut ts).unwrap();
    }
}

// ============================================================================
// Benchmark Entry Points
// ============================================================================

/// Benchmark entry point for single instance folding prover with witness size scaling.
fn bench_lin_prover(c: &mut Criterion) {
    bench_prover_protocol::<SingleInstanceFoldProver>(c, single_instance_fold::WITNESS_SCALING);
}

/// Benchmark entry point for single instance folding verifier with witness size scaling.
fn bench_lin_verifier(c: &mut Criterion) {
    bench_verifier_protocol::<SingleInstanceFoldVerifier>(c, single_instance_fold::WITNESS_SCALING);
}

criterion::criterion_group!(benches, bench_lin_prover, bench_lin_verifier);
criterion::criterion_main!(benches);
