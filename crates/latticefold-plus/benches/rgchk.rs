//! Benchmarks for Construction 4.3-4.4: Range check protocol.
//!
//! The range check protocol verifies that committed ring element coefficients
//! lie within a specified range [-B, B] using double commitments and gadget
//! decomposition for efficiency.
//!
//! ## Protocol Overview
//!
//! The range check operates in two phases:
//! 1. **Gadget decomposition**: Witness coefficients are decomposed using base-b
//!    representation to reduce the range check problem
//! 2. **Double commitment**: Commitments are computed in a special form that
//!    enables efficient range verification
//!
//! ## Benchmarked Operations
//!
//! - **Prover**: Witness decomposition and double commitment generation
//! - **Verifier**: Range proof verification via sumcheck protocols
//! - **Parameter scaling**: Performance across varying witness sizes, decomposition
//!   widths (k), and security parameters (κ)
//!
//! ## Paper Reference
//! Section 4.3-4.4 of the LatticeFold+ paper

#![allow(non_snake_case)]

use criterion::Criterion;
use latticefold_plus::rgchk::{Dcom, Rg, RgInstance};
use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as R;
use stark_rings_linalg::Matrix;

#[path = "utils/mod.rs"]
mod utils;
use utils::{
    helpers::{
        bench_prover_protocol, bench_rng, bench_verifier_protocol, create_transcript,
        get_validated_decomp_params, ProverBenchmark, VerifierBenchmark, WitnessPattern,
    },
    range_check,
};

// ============================================================================
// Setup Functions
// ============================================================================

/// Creates range check input for prover benchmarks.
///
/// Generates a single `RgInstance` with a random witness satisfying the
/// decomposition parameter constraints. Uses small random coefficients to
/// ensure range check validity.
fn setup_input(witness_size: usize, k: usize, kappa: usize) -> Rg<R> {
    let mut rng = bench_rng();
    let dparams = get_validated_decomp_params(k, kappa, witness_size);

    let f = WitnessPattern::SmallRandom.generate(witness_size, &mut rng);
    let A = Matrix::<R>::rand(&mut rng, kappa, witness_size);

    let instance = RgInstance::from_f(f, &A, &dparams);
    let nvars = (witness_size as f64).log2().ceil() as usize;

    Rg {
        nvars,
        instances: vec![instance],
        dparams,
    }
}

/// Generates a valid range check proof for verifier benchmarks.
///
/// Creates input, executes the prover to generate a `Dcom` proof, and
/// validates it before returning. This ensures the verifier benchmarks
/// measure only verification time, not error handling overhead.
fn setup_proof(witness_size: usize, k: usize, kappa: usize) -> (Rg<R>, Dcom<R>) {
    let rg = setup_input(witness_size, k, kappa);
    let mut ts = create_transcript();
    let dcom = rg.range_check(&[], &mut ts);

    let mut verify_ts = create_transcript();
    dcom.verify(&mut verify_ts)
        .expect("Generated range check proof should be valid");

    (rg, dcom)
}

// ============================================================================
// Trait Implementations
// ============================================================================

/// Prover benchmark for standard range check protocol.
///
/// Measures witness decomposition and double commitment generation time
/// across varying witness sizes while keeping k=2 and κ=2 fixed.
struct RangeCheckProver;

impl ProverBenchmark for RangeCheckProver {
    type Input = Rg<R>;
    type Output = Dcom<R>;
    type Params = (usize, usize, usize);

    fn group_name() -> &'static str {
        "RangeCheck-Prover"
    }

    fn setup_input((witness_size, k, kappa): Self::Params) -> Self::Input {
        setup_input(witness_size, k, kappa)
    }

    fn param_label((witness_size, k, kappa): Self::Params) -> String {
        format!("w={}_k={}_κ={}", witness_size, k, kappa)
    }

    fn throughput((witness_size, _, _): Self::Params) -> u64 {
        witness_size as u64
    }

    fn run_prover(input: Self::Input) -> Self::Output {
        let mut ts = create_transcript();
        input.range_check(&[], &mut ts)
    }
}

/// Verifier benchmark for standard range check protocol.
///
/// Measures proof verification time using sumcheck protocols to verify
/// the double commitment structure and range constraints.
struct RangeCheckVerifier;

impl VerifierBenchmark for RangeCheckVerifier {
    type Input = Rg<R>;
    type Proof = Dcom<R>;
    type Params = (usize, usize, usize);

    fn group_name() -> &'static str {
        "RangeCheck-Verifier"
    }

    fn setup_proof((witness_size, k, kappa): Self::Params) -> (Self::Input, Self::Proof) {
        setup_proof(witness_size, k, kappa)
    }

    fn param_label((witness_size, k, kappa): Self::Params) -> String {
        format!("w={}_k={}_κ={}", witness_size, k, kappa)
    }

    fn throughput((witness_size, _, _): Self::Params) -> u64 {
        witness_size as u64
    }

    fn run_verifier(_input: &Self::Input, proof: &Self::Proof) {
        let mut ts = create_transcript();
        proof.verify(&mut ts).unwrap()
    }
}

/// Prover benchmark measuring decomposition width (k) scaling.
///
/// Tests how prover performance scales with increasing k ∈ [2,3,4,5]
/// while keeping witness_size=65536 and κ=2 fixed.
struct RangeCheckKScaling;

impl ProverBenchmark for RangeCheckKScaling {
    type Input = Rg<R>;
    type Output = Dcom<R>;
    type Params = (usize, usize, usize);

    fn group_name() -> &'static str {
        "RangeCheck-KScaling"
    }

    fn setup_input((witness_size, k, kappa): Self::Params) -> Self::Input {
        setup_input(witness_size, k, kappa)
    }

    fn param_label((_witness_size, k, _kappa): Self::Params) -> String {
        format!("{}", k)
    }

    fn throughput((witness_size, _, _): Self::Params) -> u64 {
        witness_size as u64
    }

    fn run_prover(input: Self::Input) -> Self::Output {
        let mut ts = create_transcript();
        input.range_check(&[], &mut ts)
    }
}

/// Prover benchmark measuring security parameter (κ) scaling.
///
/// Tests how prover performance scales with increasing κ ∈ [2,3,4,5]
/// while keeping witness_size=65536 and k=2 fixed. Higher κ values
/// increase commitment matrix dimensions and security.
struct RangeCheckKappaScaling;

impl ProverBenchmark for RangeCheckKappaScaling {
    type Input = Rg<R>;
    type Output = Dcom<R>;
    type Params = (usize, usize, usize);

    fn group_name() -> &'static str {
        "RangeCheck-KappaScaling"
    }

    fn setup_input((witness_size, k, kappa): Self::Params) -> Self::Input {
        setup_input(witness_size, k, kappa)
    }

    fn param_label((_witness_size, _k, kappa): Self::Params) -> String {
        format!("{}", kappa)
    }

    fn throughput((witness_size, _, _): Self::Params) -> u64 {
        witness_size as u64
    }

    fn run_prover(input: Self::Input) -> Self::Output {
        let mut ts = create_transcript();
        input.range_check(&[], &mut ts)
    }
}

// ============================================================================
// Benchmark Entry Points
// ============================================================================

/// Benchmark entry point for range check prover with witness size scaling.
fn bench_rgchk_prover(c: &mut Criterion) {
    bench_prover_protocol::<RangeCheckProver>(c, range_check::WITNESS_SCALING);
}

/// Benchmark entry point for range check verifier with witness size scaling.
fn bench_rgchk_verifier(c: &mut Criterion) {
    bench_verifier_protocol::<RangeCheckVerifier>(c, range_check::WITNESS_SCALING);
}

/// Benchmark entry point for range check prover with k-scaling.
fn bench_rgchk_k_scaling(c: &mut Criterion) {
    bench_prover_protocol::<RangeCheckKScaling>(c, range_check::K_SCALING);
}

/// Benchmark entry point for range check prover with κ-scaling.
fn bench_rgchk_kappa_scaling(c: &mut Criterion) {
    bench_prover_protocol::<RangeCheckKappaScaling>(c, range_check::KAPPA_SCALING);
}

criterion::criterion_group!(
    benches,
    bench_rgchk_prover,
    bench_rgchk_verifier,
    bench_rgchk_k_scaling,
    bench_rgchk_kappa_scaling,
);
criterion::criterion_main!(benches);
