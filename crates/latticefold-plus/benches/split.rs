//! Benchmarks for Construction 4.1: Split function.
//!
//! The split function performs gadget decomposition on double commitment matrices,
//! converting ring element representations into base ring scalar form. This is a
//! critical step in the range check protocol for efficient constraint verification.
//!
//! ## Protocol Overview
//!
//! The split function operates on commitment matrices (not witness matrices):
//! 1. **Gadget decomposition**: Applies base-b decomposition to commitment matrices
//! 2. **Ring-to-scalar conversion**: Transforms ring elements into base field scalars
//! 3. **Deterministic operation**: No verifier benchmark needed (computation is deterministic)
//!
//! ## Benchmarked Operations
//!
//! - **Prover**: Matrix decomposition with varying k_first and κ parameters
//! - **No verifier**: Split is a deterministic function with no interactive proof
//! - **Parameter scaling**: Performance across decomposition width (k_first) and
//!   security parameter (κ) variations
//!
//! ## Paper Reference
//! Section 4.1 of the LatticeFold+ paper

#![allow(non_snake_case)]

use criterion::Criterion;
use latticefold_plus::utils::split;
use stark_rings::{cyclotomic_ring::models::frog_ring::RqPoly as R, PolyRing};
use stark_rings_linalg::Matrix;

#[path = "utils/mod.rs"]
mod utils;
use utils::{
    helpers::{bench_prover_protocol, bench_rng, DecompParams, ProverBenchmark},
    split as split_params,
};

// ============================================================================
// Setup Functions
// ============================================================================

/// Creates split function input for prover benchmarks.
///
/// Generates a random κ × (k_first·d) commitment matrix where d is the ring
/// dimension. This matrix represents the output of a double commitment scheme
/// that will be decomposed into base ring scalars.
fn setup_input(k_first: usize, kappa: usize) -> Matrix<R> {
    let mut rng = bench_rng();
    let d = R::dimension();
    let cols = k_first * d;
    let mat = Matrix::<R>::rand(&mut rng, kappa, cols);

    assert_eq!(mat.nrows, kappa, "Matrix should have kappa rows");
    assert_eq!(mat.ncols, cols, "Matrix should have k_first * d columns");

    mat
}

// ============================================================================
// Trait Implementations
// ============================================================================

/// Prover benchmark with varying witness size, k_first, and κ parameters.
///
/// Tests combined parameter scaling to measure how split performance varies
/// across different protocol configurations. Each parameter set uses different
/// values for witness size, decomposition width, and security parameter.
struct SplitVaryingParams;

impl ProverBenchmark for SplitVaryingParams {
    type Input = Matrix<R>;
    type Output = Vec<<R as PolyRing>::BaseRing>;
    type Params = (usize, usize, usize);

    fn group_name() -> &'static str {
        "Split-VaryingParams"
    }

    fn setup_input((_witness_size, k_first, kappa): Self::Params) -> Self::Input {
        setup_input(k_first, kappa)
    }

    fn param_label((witness_size, k_first, kappa): Self::Params) -> String {
        format!("w={}_k={}_κ={}", witness_size, k_first, kappa)
    }

    fn throughput((witness_size, _, _): Self::Params) -> u64 {
        witness_size as u64
    }

    fn run_prover(com: Self::Input) -> Self::Output {
        let params = DecompParams::compute();
        let witness_size = 65536;
        split(&com, witness_size, params.b, params.l)
    }
}

/// Prover benchmark measuring first decomposition width (k_first) scaling.
///
/// Tests how performance scales with k_first ∈ [2,4,6,8] while keeping
/// witness_size=131072 and κ=2 fixed. Higher k_first values increase
/// the commitment matrix width, affecting decomposition cost.
struct SplitScalingKFirst;

impl ProverBenchmark for SplitScalingKFirst {
    type Input = Matrix<R>;
    type Output = Vec<<R as PolyRing>::BaseRing>;
    type Params = (usize, usize, usize);

    fn group_name() -> &'static str {
        "Split-Scaling-KFirst"
    }

    fn setup_input((_witness_size, k_first, kappa): Self::Params) -> Self::Input {
        setup_input(k_first, kappa)
    }

    fn param_label((_witness_size, k_first, _kappa): Self::Params) -> String {
        format!("{}", k_first)
    }

    fn throughput((witness_size, _, _): Self::Params) -> u64 {
        witness_size as u64
    }

    fn run_prover(com: Self::Input) -> Self::Output {
        let params = DecompParams::compute();
        let witness_size = 131072;
        split(&com, witness_size, params.b, params.l)
    }
}

/// Prover benchmark measuring security parameter (κ) scaling.
///
/// Tests how performance scales with κ ∈ [1,2,3,4] while keeping
/// witness_size=131072 and k_first=4 fixed. Higher κ values increase
/// the number of commitment matrix rows.
struct SplitScalingKappa;

impl ProverBenchmark for SplitScalingKappa {
    type Input = Matrix<R>;
    type Output = Vec<<R as PolyRing>::BaseRing>;
    type Params = (usize, usize, usize);

    fn group_name() -> &'static str {
        "Split-Scaling-Kappa"
    }

    fn setup_input((_witness_size, k_first, kappa): Self::Params) -> Self::Input {
        setup_input(k_first, kappa)
    }

    fn param_label((_witness_size, _k_first, kappa): Self::Params) -> String {
        format!("{}", kappa)
    }

    fn throughput((witness_size, _, _): Self::Params) -> u64 {
        witness_size as u64
    }

    fn run_prover(com: Self::Input) -> Self::Output {
        let params = DecompParams::compute();
        let witness_size = 131072;
        split(&com, witness_size, params.b, params.l)
    }
}

// ============================================================================
// Benchmark Entry Points
// ============================================================================

/// Benchmark entry point for split function with varying parameters.
fn bench_split_varying_params(c: &mut Criterion) {
    bench_prover_protocol::<SplitVaryingParams>(c, split_params::WITNESS_SCALING);
}

/// Benchmark entry point for split function with k_first scaling.
fn bench_split_scaling_k_first(c: &mut Criterion) {
    bench_prover_protocol::<SplitScalingKFirst>(c, split_params::K_FIRST_SCALING);
}

/// Benchmark entry point for split function with κ scaling.
fn bench_split_scaling_kappa(c: &mut Criterion) {
    bench_prover_protocol::<SplitScalingKappa>(c, split_params::KAPPA_SCALING);
}

criterion::criterion_group!(
    benches,
    bench_split_varying_params,
    bench_split_scaling_k_first,
    bench_split_scaling_kappa,
);
criterion::criterion_main!(benches);
