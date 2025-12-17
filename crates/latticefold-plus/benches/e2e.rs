//! Benchmarks for end-to-end LatticeFold+ protocol.
//!
//! The end-to-end (E2E) benchmarks measure the complete LatticeFold+ proving and
//! verification stack, including all sub-protocols integrated together. This provides
//! realistic performance measurements for the full system in production scenarios.
//!
//! ## Protocol Overview
//!
//! The complete LatticeFold+ protocol executes these phases:
//! 1. **R1CS commitment**: Creates committed R1CS instances from witnesses
//! 2. **Linearization**: Converts R1CS to LinB instances
//! 3. **Range check**: Verifies witness coefficients via double commitments
//! 4. **Commitment transformation**: Prepares commitments for folding
//! 5. **Multilinear folding**: Batches L instances into one LinB2 instance
//! 6. **Verification**: Verifier checks all sumcheck proofs
//!
//! ## Benchmarked Operations
//!
//! - **Prover**: Full prove operation on L instances
//! - **Verifier**: Full verification of batched proof
//! - **Protocol scaling**: Performance across n ∈ [64K, 128K]
//! - **Folding arity scaling**: Impact of L ∈ [2,3,4,5] on end-to-end cost
//!
//! ## Paper Reference
//! Complete protocol described in Section 5 of the LatticeFold+ paper

#![allow(non_snake_case)]

use criterion::Criterion;
use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
use latticefold_plus::{
    lin::LinParameters,
    plus::{PlusParameters, PlusProof, PlusProver, PlusVerifier},
    r1cs::{ComR1CS, ComR1CSProof},
    rgchk::DecompParameters,
    transcript::PoseidonTranscript,
    utils::estimate_bound,
};
use stark_rings::{cyclotomic_ring::models::frog_ring::RqPoly as R, PolyRing};

#[path = "utils/mod.rs"]
mod utils;
use utils::{
    e2e,
    helpers::{
        bench_prover_protocol, bench_rng, bench_verifier_protocol, create_ajtai_matrix,
        create_transcript, ProverBenchmark, R1CSBuilder, VerifierBenchmark, WitnessPattern,
    },
};

// ============================================================================
// Setup Functions
// ============================================================================

/// Creates end-to-end protocol input for prover benchmarks.
///
/// Generates a complete setup including decomposed R1CS constraints, L identical
/// committed R1CS instances, Ajtai commitment matrix, and constraint matrices.
/// Returns a PlusProver initialized with all parameters and the instances to fold.
fn setup_input(
    n: usize,
    L: usize,
    k: usize,
    kappa: usize,
) -> (PlusProver<R, PoseidonTranscript<R>>, Vec<ComR1CS<R>>) {
    let mut rng = bench_rng();
    let m = n / k;

    let d = R::dimension();
    let b = (d / 2) as u128;
    let sop = d * 128;
    let B = estimate_bound(sop, L, d, k) / 2;

    let l = ((<<R as PolyRing>::BaseRing as ark_ff::PrimeField>::MODULUS.0[0] as f64).ln()
        / ((d / 2) as f64).ln())
    .ceil() as usize;

    let dparams = DecompParameters { b, k, l };
    let params = LinParameters {
        kappa,
        decomp: dparams,
    };

    let r1cs = R1CSBuilder::new(n, k, B).build_decomposed_square();

    let z = WitnessPattern::BinaryChoice.generate(m, &mut rng);

    let A = create_ajtai_matrix(kappa, n, &mut rng);

    let cr1cs = ComR1CS::new(r1cs, z, 1, B, k, &A);

    let M = cr1cs.x.matrices();

    let pparams = PlusParameters { lin: params, B };

    let ts = create_transcript();
    let prover = PlusProver::init(A, M, 1, pparams, ts);

    let instances = vec![cr1cs; L];

    (prover, instances)
}

/// Generates a valid end-to-end proof for verifier benchmarks.
///
/// Creates a complete prover setup, generates a proof by folding L instances,
/// and prepares the verifier with all necessary parameters. Returns the verifier
/// and proof ready for verification benchmarking.
fn setup_proof(
    n: usize,
    L: usize,
    k: usize,
    kappa: usize,
) -> (
    PlusVerifier<R, PoseidonTranscript<R>>,
    PlusProof<R, ComR1CSProof<R>>,
) {
    let (mut prover, instances) = setup_input(n, L, k, kappa);

    let proof = prover.prove(&instances);

    let mut rng = bench_rng();
    let m = n / k;

    let d = R::dimension();
    let b = (d / 2) as u128;
    let sop = d * 128;
    let B = estimate_bound(sop, L, d, k) / 2;

    let l = ((<<R as PolyRing>::BaseRing as ark_ff::PrimeField>::MODULUS.0[0] as f64).ln()
        / ((d / 2) as f64).ln())
    .ceil() as usize;

    let dparams = DecompParameters { b, k, l };
    let params = LinParameters {
        kappa,
        decomp: dparams,
    };
    let pparams = PlusParameters { lin: params, B };

    let A = create_ajtai_matrix(kappa, n, &mut rng);

    let r1cs = R1CSBuilder::new(n, k, B).build_decomposed_square();
    let z = WitnessPattern::BinaryChoice.generate(m, &mut rng);
    let cr1cs = ComR1CS::new(r1cs, z, 1, B, k, &A);
    let M = cr1cs.x.matrices();

    let ts = create_transcript();
    let verifier = PlusVerifier::init(A, M, pparams, ts);

    (verifier, proof)
}

// ============================================================================
// Trait Implementations
// ============================================================================

/// Prover benchmark for end-to-end LatticeFold+ protocol.
///
/// Measures complete proving time including all sub-protocols: linearization,
/// range check, commitment transformation, and multilinear folding. Folds L
/// identical instances to provide consistent benchmarking conditions.
struct E2EProver;

impl ProverBenchmark for E2EProver {
    type Input = (PlusProver<R, PoseidonTranscript<R>>, Vec<ComR1CS<R>>);
    type Output = PlusProof<R, ComR1CSProof<R>>;
    type Params = (usize, usize, usize, usize);

    fn group_name() -> &'static str {
        "E2E-Prover"
    }

    fn setup_input((n, L, k, kappa): Self::Params) -> Self::Input {
        setup_input(n, L, k, kappa)
    }

    fn param_label((n, L, k, kappa): Self::Params) -> String {
        format!("n={}_L={}_k={}_κ={}", n, L, k, kappa)
    }

    fn throughput((n, L, _, _): Self::Params) -> u64 {
        (n * L) as u64
    }

    fn run_prover((mut prover, instances): Self::Input) -> Self::Output {
        prover.prove(&instances)
    }
}

/// Verifier benchmark for end-to-end LatticeFold+ protocol.
///
/// Measures complete verification time including all sumcheck verifications
/// across the protocol stack. Demonstrates sublinear verification cost
/// relative to prover work.
struct E2EVerifier;

impl VerifierBenchmark for E2EVerifier {
    type Input = PlusVerifier<R, PoseidonTranscript<R>>;
    type Proof = PlusProof<R, ComR1CSProof<R>>;
    type Params = (usize, usize, usize, usize);

    fn group_name() -> &'static str {
        "E2E-Verifier"
    }

    fn setup_proof((n, L, k, kappa): Self::Params) -> (Self::Input, Self::Proof) {
        setup_proof(n, L, k, kappa)
    }

    fn param_label((n, L, k, kappa): Self::Params) -> String {
        format!("n={}_L={}_k={}_κ={}", n, L, k, kappa)
    }

    fn throughput((n, L, _, _): Self::Params) -> u64 {
        (n * L) as u64
    }

    fn run_verifier(verifier: &Self::Input, proof: &Self::Proof) {
        let mut v = verifier.clone();
        assert!(v.verify(proof), "Verification should succeed");
    }
}

/// Prover benchmark measuring folding arity (L) scaling.
///
/// Tests how end-to-end prover performance scales with increasing L values.
/// Demonstrates batching efficiency across the complete protocol when more
/// instances are folded together.
struct E2EFoldingArity;

impl ProverBenchmark for E2EFoldingArity {
    type Input = (PlusProver<R, PoseidonTranscript<R>>, Vec<ComR1CS<R>>);
    type Output = PlusProof<R, ComR1CSProof<R>>;
    type Params = (usize, usize, usize, usize);

    fn group_name() -> &'static str {
        "E2E-FoldingArity"
    }

    fn setup_input((n, L, k, kappa): Self::Params) -> Self::Input {
        setup_input(n, L, k, kappa)
    }

    fn param_label((_n, L, _k, _kappa): Self::Params) -> String {
        format!("L={}", L)
    }

    fn throughput((n, L, _, _): Self::Params) -> u64 {
        (n * L) as u64
    }

    fn run_prover((mut prover, instances): Self::Input) -> Self::Output {
        prover.prove(&instances)
    }
}

// ============================================================================
// Benchmark Entry Points
// ============================================================================

/// Benchmark entry point for E2E prover with protocol scaling.
fn bench_e2e_prover(c: &mut Criterion) {
    bench_prover_protocol::<E2EProver>(c, e2e::PROTOCOL_SCALING);
}

/// Benchmark entry point for E2E verifier with protocol scaling.
fn bench_e2e_verifier(c: &mut Criterion) {
    bench_verifier_protocol::<E2EVerifier>(c, e2e::PROTOCOL_SCALING);
}

/// Benchmark entry point for E2E prover with folding arity scaling.
fn bench_e2e_folding_arity(c: &mut Criterion) {
    bench_prover_protocol::<E2EFoldingArity>(c, e2e::FOLDING_ARITY);
}

criterion::criterion_group!(
    benches,
    bench_e2e_prover,
    bench_e2e_verifier,
    bench_e2e_folding_arity,
);
criterion::criterion_main!(benches);
