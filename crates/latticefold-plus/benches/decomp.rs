//! Benchmarks for Construction 5.3: Decomposition protocol.
//!
//! The decomposition protocol splits a LinB2 instance (norm bound B²) into two
//! LinB instances (each with norm bound B). This is critical for IVC/PCD
//! applications to prevent norm explosion across multiple folding rounds.
//!
//! ## Protocol Overview
//!
//! The decomposition operates by witness splitting:
//! 1. **Witness decomposition**: Splits witness f = F^(0) + B·F^(1) where both
//!    ||F^(0)|| and ||F^(1)|| are bounded by B
//! 2. **Dual commitments**: Creates commitments and linearization proofs for
//!    both decomposed components
//! 3. **Norm control**: Ensures each output maintains norm bound B, preventing
//!    exponential norm growth in recursive composition
//!
//! ## Benchmarked Operations
//!
//! - **Prover**: Witness decomposition and dual linearization proof generation
//! - **Verifier**: Verification of both decomposed LinB instances
//! - **Roundtrip**: Complete fold→decompose cycle to measure IVC overhead
//!
//! ## Paper Reference
//! Construction 5.3, Section 5.3 of the LatticeFold+ paper

#![allow(non_snake_case)]

use criterion::{criterion_group, criterion_main, BatchSize, Criterion};
use latticefold_plus::{
    decomp::{Decomp, DecompProof},
    lin::{LinB, Linearize, LinearizedVerify},
    r1cs::{r1cs_decomposed_square, ComR1CS},
};
use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as R;
use stark_rings_linalg::Matrix;

#[path = "utils/mod.rs"]
mod utils;
use utils::{
    decomposition,
    helpers::{
        bench_prover_protocol, bench_rng, bench_verifier_protocol, configure_benchmark_group,
        create_ajtai_matrix, create_transcript, ProverBenchmark, R1CSBuilder, VerifierBenchmark,
        WitnessPattern,
    },
};

// ============================================================================
// Setup Functions
// ============================================================================

/// Creates decomposition input for prover benchmarks.
///
/// Generates a LinB2 instance by first creating a committed R1CS, linearizing
/// it, and then preparing the decomposition input structure. Uses witness with
/// all ones to ensure valid R1CS satisfaction and norm bounds.
fn setup_input(n: usize, k: usize, kappa: usize, _B: usize) -> Decomp<R> {
    let mut rng = bench_rng();
    let r1cs = R1CSBuilder::new(n, k, 2).build_basic();
    let r1cs = r1cs_decomposed_square(r1cs, n, 2, k);

    let z = WitnessPattern::AllOnes.generate(n / k, &mut rng);
    let A = create_ajtai_matrix(kappa, n, &mut rng);
    let cr1cs = ComR1CS::new(r1cs, z, 1, 2, k, &A);

    let mut ts = create_transcript();
    let (linb, lproof) = cr1cs.linearize(&mut ts);

    let mut ts = create_transcript();
    lproof.verify(&mut ts);

    Decomp {
        f: cr1cs.f,
        r: lproof.r.iter().map(|&r| (r, r)).collect::<Vec<_>>(),
        M: cr1cs.x.matrices(),
    }
}

/// Generates a valid decomposition proof for verifier benchmarks.
///
/// Creates a LinB2 instance, executes the decomposition protocol to generate
/// a `DecompProof`, and validates it before returning. This ensures the
/// verifier benchmarks measure only verification time, not error handling.
fn setup_proof(n: usize, k: usize, kappa: usize, B: usize) -> (Decomp<R>, DecompProof<R>) {
    let mut rng = bench_rng();
    let r1cs = R1CSBuilder::new(n, k, 2).build_basic();
    let r1cs = r1cs_decomposed_square(r1cs, n, 2, k);

    let z = WitnessPattern::AllOnes.generate(n / k, &mut rng);
    let A = create_ajtai_matrix(kappa, n, &mut rng);
    let cr1cs = ComR1CS::new(r1cs, z, 1, 2, k, &A);

    let mut ts = create_transcript();
    let (linb, lproof) = cr1cs.linearize(&mut ts);

    let mut ts = create_transcript();
    lproof.verify(&mut ts);

    let decomp = Decomp {
        f: cr1cs.f,
        r: lproof.r.iter().map(|&r| (r, r)).collect::<Vec<_>>(),
        M: cr1cs.x.matrices(),
    };

    let A_decomp = create_ajtai_matrix(kappa, n, &mut rng);
    let (outputs, proof) = decomp.decompose(&A_decomp, B as u128);

    proof.verify(&decomp.f, &decomp.r, B as u128);

    (decomp, proof)
}

// ============================================================================
// Trait Implementations
// ============================================================================

/// Prover benchmark for standard decomposition protocol.
///
/// Measures witness decomposition and dual linearization proof generation
/// time across varying witness sizes (32K-128K) while keeping k=2, κ=2,
/// and B=50 fixed.
struct DecompositionProver;

impl ProverBenchmark for DecompositionProver {
    type Input = (Decomp<R>, Matrix<R>);
    type Output = ((LinB<R>, LinB<R>), DecompProof<R>);
    type Params = (usize, usize, usize, usize);

    fn group_name() -> &'static str {
        "Decomposition-Prover"
    }

    fn setup_input((n, k, kappa, B): Self::Params) -> Self::Input {
        let mut rng = bench_rng();
        let decomp = setup_input(n, k, kappa, B);
        let A = create_ajtai_matrix(kappa, n, &mut rng);
        (decomp, A)
    }

    fn param_label((n, k, kappa, B): Self::Params) -> String {
        format!("n={}_k={}_κ={}_B={}", n, k, kappa, B)
    }

    fn throughput((n, _, _, _): Self::Params) -> u64 {
        n as u64
    }

    fn run_prover((decomp, A): Self::Input) -> Self::Output {
        let B = 2;
        decomp.decompose(&A, B as u128)
    }
}

/// Verifier benchmark for standard decomposition protocol.
///
/// Measures verification time for decomposition proofs. Verifies that both
/// output LinB instances satisfy the norm bound B and are consistent with
/// the original LinB2 witness.
struct DecompositionVerifier;

impl VerifierBenchmark for DecompositionVerifier {
    type Input = Decomp<R>;
    type Proof = DecompProof<R>;
    type Params = (usize, usize, usize, usize);

    fn group_name() -> &'static str {
        "Decomposition-Verifier"
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

    fn run_verifier(input: &Self::Input, proof: &Self::Proof) {
        let B = 2;
        proof.verify(&input.f, &input.r, B as u128)
    }
}

// ============================================================================
// Benchmark Entry Points
// ============================================================================

/// Benchmark entry point for decomposition prover with witness size scaling.
fn bench_decomp_prover(c: &mut Criterion) {
    bench_prover_protocol::<DecompositionProver>(c, decomposition::WITNESS_SCALING);
}

/// Benchmark entry point for decomposition verifier with witness size scaling.
fn bench_decomp_verifier(c: &mut Criterion) {
    bench_verifier_protocol::<DecompositionVerifier>(c, decomposition::WITNESS_SCALING);
}

/// Benchmark entry point for complete fold→decompose roundtrip cycle.
///
/// Measures the full cost of folding followed by decomposition, which is the
/// typical pattern in IVC/PCD applications. This demonstrates the overhead of
/// maintaining consistent norm bounds across recursive composition.
fn bench_decomp_roundtrip(c: &mut Criterion) {
    let mut group = c.benchmark_group("Decomposition-Roundtrip");
    configure_benchmark_group(&mut group);

    let (n, k, kappa, B) = decomposition::WITNESS_SCALING[1];

    group.throughput(criterion::Throughput::Elements((2 * n) as u64));

    group.bench_function(
        format!("fold_decompose_n={}_k={}_κ={}_B={}", n, k, kappa, B),
        |bencher| {
            bencher.iter_batched(
                || {
                    let mut rng = bench_rng();
                    let decomp = setup_input(n, k, kappa, B);
                    let A = create_ajtai_matrix(kappa, n, &mut rng);
                    (decomp, A)
                },
                |(decomp, A)| {
                    let ((linb0, linb1), _proof) = decomp.decompose(&A, B as u128);
                    (linb0, linb1)
                },
                BatchSize::SmallInput,
            );
        },
    );

    group.finish();
}

criterion_group!(
    benches,
    bench_decomp_prover,
    bench_decomp_verifier,
    bench_decomp_roundtrip,
);
criterion_main!(benches);
