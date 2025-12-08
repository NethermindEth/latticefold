//! Benchmarks for double commitment protocol.
//!
//! The double commitment protocol creates RgInstance structures from witness vectors,
//! which form the foundation for range check verification in LatticeFold+. Double
//! commitments commit to witness coefficients in a special form that enables efficient
//! range verification via gadget decomposition.
//!
//! ## Protocol Overview
//!
//! The `RgInstance::from_f` operation performs the following:
//! 1. Takes a witness vector `f` of ring elements
//! 2. Applies gadget decomposition using base-b representation
//! 3. Constructs a double commitment structure suitable for range check protocol
//! 4. Prepares commitments for efficient sumcheck verification
//!
//! ## Benchmarked Operations
//!
//! - **Prover**: RgInstance creation from witness vectors
//! - **Witness scaling**: Performance across n ∈ [32K, 64K, 128K]
//! - **Decomposition width (k) scaling**: Impact of k ∈ [2,4]
//!
//! ## Paper Reference
//! Section 4.1-4.3 of the LatticeFold+ paper

#![allow(non_snake_case)]

use criterion::Criterion;
use latticefold_plus::{
    r1cs::ComR1CS, rgchk::{DecompParameters, RgInstance}, utils::estimate_bound,
};
use stark_rings::{cyclotomic_ring::models::frog_ring::RqPoly as R, PolyRing};
use stark_rings_linalg::Matrix;

#[path = "utils/mod.rs"]
mod utils;
use utils::{
    double_commitment,
    helpers::{
        bench_prover_protocol, bench_rng, create_ajtai_matrix, ProverBenchmark, R1CSBuilder,
        WitnessPattern,
    },
};

// ============================================================================
// Setup Functions
// ============================================================================

/// Creates double commitment input for prover benchmarks.
///
/// Generates a ComR1CS instance with decomposed square R1CS constraints,
/// then extracts the witness vector `f` for RgInstance creation. Returns
/// all components needed for the double commitment operation.
fn setup_input(
    n: usize,
    L: usize,
    k: usize,
    kappa: usize,
) -> (Matrix<R>, Vec<R>, DecompParameters) {
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

    let r1cs = R1CSBuilder::new(n, k, B).build_decomposed_square();

    let z = WitnessPattern::BinaryChoice.generate(m, &mut rng);

    let A = create_ajtai_matrix(kappa, n, &mut rng);

    let cr1cs = ComR1CS::new(r1cs, z, 1, B, k, &A);

    (A, cr1cs.f, dparams)
}

// ============================================================================
// Trait Implementations
// ============================================================================

/// Prover benchmark for double commitment with witness size scaling.
///
/// Measures the time to create an RgInstance from a witness vector using
/// the double commitment structure. Evaluates performance across varying
/// witness sizes to understand scaling characteristics.
struct DoubleCommitmentProver;

impl ProverBenchmark for DoubleCommitmentProver {
    type Input = (Matrix<R>, Vec<R>, DecompParameters);
    type Output = RgInstance<R>;
    type Params = (usize, usize, usize, usize);

    fn group_name() -> &'static str {
        "DoubleCommitment-Prover"
    }

    fn setup_input((n, L, k, kappa): Self::Params) -> Self::Input {
        setup_input(n, L, k, kappa)
    }

    fn param_label((n, L, k, kappa): Self::Params) -> String {
        format!("n={}_L={}_k={}_κ={}", n, L, k, kappa)
    }

    fn throughput((n, _, _, _): Self::Params) -> u64 {
        n as u64
    }

    fn run_prover((A, f, dparams): Self::Input) -> Self::Output {
        RgInstance::from_f(f, &A, &dparams)
    }
}

/// Prover benchmark measuring decomposition width (k) scaling.
///
/// Tests how double commitment performance scales with increasing k values.
/// Higher k creates more decomposition limbs which affects both computation
/// and proof size.
struct DoubleCommitmentKScaling;

impl ProverBenchmark for DoubleCommitmentKScaling {
    type Input = (Matrix<R>, Vec<R>, DecompParameters);
    type Output = RgInstance<R>;
    type Params = (usize, usize, usize, usize);

    fn group_name() -> &'static str {
        "DoubleCommitment-KScaling"
    }

    fn setup_input((n, L, k, kappa): Self::Params) -> Self::Input {
        setup_input(n, L, k, kappa)
    }

    fn param_label((_n, _L, k, _kappa): Self::Params) -> String {
        format!("k={}", k)
    }

    fn throughput((n, _, _, _): Self::Params) -> u64 {
        n as u64
    }

    fn run_prover((A, f, dparams): Self::Input) -> Self::Output {
        RgInstance::from_f(f, &A, &dparams)
    }
}

// ============================================================================
// Benchmark Entry Points
// ============================================================================

/// Benchmark entry point for double commitment with witness size scaling.
fn bench_double_commitment_prover(c: &mut Criterion) {
    bench_prover_protocol::<DoubleCommitmentProver>(c, double_commitment::WITNESS_SCALING);
}

/// Benchmark entry point for double commitment with k-scaling.
fn bench_double_commitment_k_scaling(c: &mut Criterion) {
    bench_prover_protocol::<DoubleCommitmentKScaling>(c, double_commitment::K_SCALING);
}

criterion::criterion_group!(
    benches,
    bench_double_commitment_prover,
    bench_double_commitment_k_scaling,
);
criterion::criterion_main!(benches);
