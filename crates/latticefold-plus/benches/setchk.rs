//! Benchmarks for Construction 4.2: Set check protocol.
//!
//! The set check protocol verifies that matrices contain monomials, meaning
//! exactly one non-zero entry per row and per column. This property is fundamental
//! for verifying R1CS constraint matrices in the LatticeFold+ scheme.
//!
//! ## Protocol Overview
//!
//! The set check operates using sumcheck protocols to efficiently verify the
//! monomial property:
//! 1. **Monomial verification**: Confirms each matrix row/column has exactly
//!    one non-zero coefficient
//! 2. **Batching**: Multiple matrices can be checked simultaneously for improved
//!    amortization
//!
//! ## Benchmarked Operations
//!
//! - **Prover**: Monomial set generation and sumcheck proof construction
//! - **Verifier**: Sumcheck verification of monomial properties
//! - **Batching efficiency**: Performance across varying batch sizes (1-16 sets)
//!
//! ## Paper Reference
//! Section 4.2 of the LatticeFold+ paper

#![allow(non_snake_case)]

use criterion::Criterion;
use latticefold_plus::setchk::{In, MonomialSet, Out};
use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as R;
use stark_rings_linalg::SparseMatrix;

#[path = "utils/mod.rs"]
mod utils;
use utils::{
    helpers::{
        bench_prover_protocol, bench_verifier_protocol, create_transcript, ProverBenchmark,
        VerifierBenchmark,
    },
    set_check,
};

// ============================================================================
// Setup Functions
// ============================================================================

/// Creates set check input for prover benchmarks.
///
/// Generates identity matrices to guarantee monomial properties (exactly one
/// non-zero per row/column). Uses sumcheck for efficient batch verification
/// across multiple sets.
fn setup_input(set_size: usize, num_batches: usize) -> In<R> {
    let mut sets = Vec::with_capacity(num_batches);

    for _ in 0..num_batches {
        let m = SparseMatrix::<R>::identity(set_size);
        sets.push(MonomialSet::Matrix(m));
    }

    let nvars = (set_size as f64).log2().ceil() as usize;

    In { sets, nvars }
}

/// Generates a valid set check proof for verifier benchmarks.
///
/// Creates input, executes the prover to generate an `Out` proof, and
/// validates it before returning. This ensures the verifier benchmarks
/// measure only verification time, not error handling overhead.
fn setup_proof(set_size: usize, num_batches: usize) -> (In<R>, Out<R>) {
    let input = setup_input(set_size, num_batches);
    let mut ts = create_transcript();

    let output = input.set_check(&[], &mut ts);

    let mut verify_ts = create_transcript();
    output
        .verify(&mut verify_ts)
        .expect("Generated set check proof should be valid");

    (input, output)
}

// ============================================================================
// Trait Implementations
// ============================================================================

/// Prover benchmark for standard set check protocol.
///
/// Measures monomial verification and sumcheck proof generation time
/// across varying set sizes (256-1024) and batch counts.
struct SetCheckProver;

impl ProverBenchmark for SetCheckProver {
    type Input = In<R>;
    type Output = Out<R>;
    type Params = (usize, usize);

    fn group_name() -> &'static str {
        "SetCheck-Prover"
    }

    fn setup_input((set_size, num_batches): Self::Params) -> Self::Input {
        setup_input(set_size, num_batches)
    }

    fn param_label((set_size, num_batches): Self::Params) -> String {
        format!("size={}_batches={}", set_size, num_batches)
    }

    fn throughput((set_size, num_batches): Self::Params) -> u64 {
        (set_size * num_batches) as u64
    }

    fn run_prover(input: Self::Input) -> Self::Output {
        let mut ts = create_transcript();
        input.set_check(&[], &mut ts)
    }
}

/// Verifier benchmark for standard set check protocol.
///
/// Measures proof verification time using sumcheck protocols to verify
/// the monomial property across batched matrix sets.
struct SetCheckVerifier;

impl VerifierBenchmark for SetCheckVerifier {
    type Input = In<R>;
    type Proof = Out<R>;
    type Params = (usize, usize);

    fn group_name() -> &'static str {
        "SetCheck-Verifier"
    }

    fn setup_proof((set_size, num_batches): Self::Params) -> (Self::Input, Self::Proof) {
        setup_proof(set_size, num_batches)
    }

    fn param_label((set_size, num_batches): Self::Params) -> String {
        format!("size={}_batches={}", set_size, num_batches)
    }

    fn throughput((set_size, num_batches): Self::Params) -> u64 {
        (set_size * num_batches) as u64
    }

    fn run_verifier(_input: &Self::Input, proof: &Self::Proof) {
        let mut ts = create_transcript();
        proof.verify(&mut ts).unwrap()
    }
}

/// Prover benchmark measuring batching efficiency.
///
/// Tests how prover performance scales with increasing batch counts
/// (1, 2, 4, 8, 16) while keeping set_size=512 fixed. Demonstrates
/// amortization benefits of batched verification.
struct SetCheckBatching;

impl ProverBenchmark for SetCheckBatching {
    type Input = In<R>;
    type Output = Out<R>;
    type Params = usize;

    fn group_name() -> &'static str {
        "SetCheck-Batching"
    }

    fn setup_input(num_batches: Self::Params) -> Self::Input {
        const SET_SIZE: usize = 512;
        setup_input(SET_SIZE, num_batches)
    }

    fn param_label(num_batches: Self::Params) -> String {
        format!("{}", num_batches)
    }

    fn throughput(num_batches: Self::Params) -> u64 {
        const SET_SIZE: usize = 512;
        (SET_SIZE * num_batches) as u64
    }

    fn run_prover(input: Self::Input) -> Self::Output {
        let mut ts = create_transcript();
        input.set_check(&[], &mut ts)
    }
}

// ============================================================================
// Benchmark Entry Points
// ============================================================================

/// Benchmark entry point for set check prover with size and batch scaling.
fn bench_setchk_prover(c: &mut Criterion) {
    bench_prover_protocol::<SetCheckProver>(c, set_check::SET_SIZES);
}

/// Benchmark entry point for set check verifier with size and batch scaling.
fn bench_setchk_verifier(c: &mut Criterion) {
    bench_verifier_protocol::<SetCheckVerifier>(c, set_check::SET_SIZES);
}

/// Benchmark entry point for set check prover with batching efficiency test.
fn bench_setchk_batching(c: &mut Criterion) {
    const BATCH_COUNTS: [usize; 5] = [1, 2, 4, 8, 16];
    bench_prover_protocol::<SetCheckBatching>(c, &BATCH_COUNTS);
}

criterion::criterion_group!(
    benches,
    bench_setchk_prover,
    bench_setchk_verifier,
    bench_setchk_batching,
);
criterion::criterion_main!(benches);
