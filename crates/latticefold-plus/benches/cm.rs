//! Benchmarks for Construction 4.5: Commitment transformation protocol.
//!
//! The commitment transformation protocol converts double commitments produced
//! by the range check protocol into folded commitments suitable for the main
//! LatticeFold+ multilinear folding protocol. This bridges the range check and
//! folding components of the scheme.
//!
//! ## Protocol Overview
//!
//! The commitment transformation operates in two phases:
//! 1. **Commitment conversion**: Transforms L double commitments from range check
//!    into a single folded commitment structure
//! 2. **Sumcheck verification**: Uses sumcheck protocols to verify the transformation
//!    correctness and commitment consistency
//!
//! ## Benchmarked Operations
//!
//! - **Prover**: Commitment transformation and sumcheck proof generation
//! - **Verifier**: Sumcheck-based verification of the transformation
//! - **Batching efficiency**: Performance scaling across L ∈ [2,3,4,5,6,7,8] instances
//!
//! ## Paper Reference
//! Section 4.5 of the LatticeFold+ paper

#![allow(non_snake_case)]

use criterion::Criterion;
use latticefold_plus::{
    cm::{Cm, CmProof},
    rgchk::{Rg, RgInstance},
};
use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as R;
use stark_rings_linalg::{Matrix, SparseMatrix};

#[path = "utils/mod.rs"]
mod utils;
use utils::{
    commitment_transform,
    helpers::{
        bench_prover_protocol, bench_rng, bench_verifier_protocol, create_test_m_matrix,
        create_transcript, get_validated_decomp_params, ProverBenchmark, VerifierBenchmark,
        WitnessPattern,
    },
};

// ============================================================================
// Setup Functions
// ============================================================================

/// Creates commitment transformation input for prover benchmarks.
///
/// Generates L range check instances with random witnesses and commitment
/// matrices. Uses small random coefficients to ensure range check validity.
/// Each instance is prepared for transformation into folded commitment form.
fn setup_input(L: usize, witness_size: usize, k: usize, kappa: usize) -> Cm<R> {
    let mut rng = bench_rng();
    let dparams = get_validated_decomp_params(k, kappa, witness_size);

    let instances: Vec<RgInstance<R>> = (0..L)
        .map(|_| {
            let f = WitnessPattern::SmallRandom.generate(witness_size, &mut rng);
            let A = Matrix::<R>::rand(&mut rng, kappa, witness_size);
            RgInstance::from_f(f, &A, &dparams)
        })
        .collect();

    let nvars = (witness_size as f64).log2().ceil() as usize;
    let rg = Rg {
        nvars,
        instances,
        dparams,
    };

    Cm { rg }
}

/// Generates a valid commitment transformation proof for verifier benchmarks.
///
/// Creates input with L instances, executes the prover to generate a `CmProof`,
/// and validates it before returning. This ensures the verifier benchmarks
/// measure only verification time, not error handling overhead.
fn setup_proof(L: usize, witness_size: usize, k: usize, kappa: usize) -> (Cm<R>, CmProof<R>) {
    let cm = setup_input(L, witness_size, k, kappa);
    let M = create_test_m_matrix(witness_size);

    let mut ts = create_transcript();
    let (_com, proof) = cm.prove(&M, &mut ts);

    let mut verify_ts = create_transcript();
    proof
        .verify(&M, &mut verify_ts)
        .expect("Generated commitment transformation proof should be valid");

    (cm, proof)
}

// ============================================================================
// Trait Implementations
// ============================================================================

/// Prover benchmark for commitment transformation protocol.
///
/// Measures transformation and sumcheck proof generation time across varying
/// folding arity L ∈ [2,3,4,5,6,7,8] while keeping witness_size=65536, k=2,
/// and κ=2 fixed. Demonstrates batching efficiency of the transformation.
struct CommitmentTransformProver;

impl ProverBenchmark for CommitmentTransformProver {
    type Input = (Cm<R>, Vec<SparseMatrix<R>>);
    type Output = (latticefold_plus::cm::Com<R>, CmProof<R>);
    type Params = (usize, usize, usize, usize);

    fn group_name() -> &'static str {
        "CommitmentTransform-Prover"
    }

    fn setup_input((L, witness_size, k, kappa): Self::Params) -> Self::Input {
        let input = setup_input(L, witness_size, k, kappa);
        let M = create_test_m_matrix(witness_size);
        (input, M)
    }

    fn param_label((L, witness_size, k, kappa): Self::Params) -> String {
        format!("L={}_w={}_k={}_κ={}", L, witness_size, k, kappa)
    }

    fn throughput((L, _, _, _): Self::Params) -> u64 {
        L as u64
    }

    fn run_prover((input, M): Self::Input) -> Self::Output {
        let mut ts = create_transcript();
        input.prove(&M, &mut ts)
    }
}

/// Verifier benchmark for commitment transformation protocol.
///
/// Measures sumcheck verification time for commitment transformation proofs.
/// Tests verification scaling across varying L while other parameters remain
/// fixed, showing sublinear verification cost relative to prover work.
struct CommitmentTransformVerifier;

impl VerifierBenchmark for CommitmentTransformVerifier {
    type Input = Cm<R>;
    type Proof = CmProof<R>;
    type Params = (usize, usize, usize, usize);

    fn group_name() -> &'static str {
        "CommitmentTransform-Verifier"
    }

    fn setup_proof((L, witness_size, k, kappa): Self::Params) -> (Self::Input, Self::Proof) {
        setup_proof(L, witness_size, k, kappa)
    }

    fn param_label((L, witness_size, k, kappa): Self::Params) -> String {
        format!("L={}_w={}_k={}_κ={}", L, witness_size, k, kappa)
    }

    fn throughput((L, _, _, _): Self::Params) -> u64 {
        L as u64
    }

    fn run_verifier(_input: &Self::Input, proof: &Self::Proof) {
        let witness_size = 65536;
        let M = create_test_m_matrix(witness_size);
        let mut ts = create_transcript();
        proof.verify(&M, &mut ts).unwrap();
    }
}

// ============================================================================
// Benchmark Entry Points
// ============================================================================

/// Benchmark entry point for commitment transformation prover with L scaling.
fn bench_cm_prover(c: &mut Criterion) {
    bench_prover_protocol::<CommitmentTransformProver>(c, commitment_transform::FOLDING_ARITY);
}

/// Benchmark entry point for commitment transformation verifier with L scaling.
fn bench_cm_verifier(c: &mut Criterion) {
    bench_verifier_protocol::<CommitmentTransformVerifier>(c, commitment_transform::FOLDING_ARITY);
}

criterion::criterion_group!(benches, bench_cm_prover, bench_cm_verifier);
criterion::criterion_main!(benches);
