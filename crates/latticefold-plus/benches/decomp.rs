//! Benchmarks for Construction 5.3: Decomposition protocol
//!
//! The decomposition protocol (Π_decomp,B) splits a LinB2 instance with norm
//! bound B² into two LinB instances each with norm bound B. This is critical
//! for IVC/PCD applications to prevent norm explosion across multiple folding
//! rounds.
//!
//! The protocol works by decomposing the witness f = F^(0) + B·F^(1) where
//! both ||F^(0)|| and ||F^(1)|| are bounded by B, then creating commitments
//! and proofs for both components.
//!
//! Paper reference: Construction 5.3, Section 5.3

#![allow(non_snake_case)]

use criterion::{
    criterion_group, criterion_main, BatchSize, BenchmarkGroup, BenchmarkId, Criterion, Throughput,
};
use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as R;
use stark_rings_linalg::Matrix;
use utils::{bench_rng, decomposition, setup_decomp_input, setup_decomp_proof};

#[path = "utils/mod.rs"]
mod utils;

/// Configure benchmark group with benchmark settings
fn configure_benchmark_group(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    group.sample_size(10);
    group.measurement_time(std::time::Duration::from_secs(10));
    group.warm_up_time(std::time::Duration::from_secs(3));
}

/// Benchmark decomposition prover with varying parameters
///
/// Tests performance across different parameter combinations:
/// - n: witness size (length of witness vector)
/// - k: decomposition width
/// - κ (kappa): number of commitment rows
/// - B: norm bound parameter (output norm, input has B²)
///
/// The decomposition splits one witness with norm B² into two witnesses
/// with norm B each, enabling continued folding in IVC applications.
fn bench_decomp_prover(c: &mut Criterion) {
    let mut group = c.benchmark_group("Decomposition-Prover");
    configure_benchmark_group(&mut group);

    for &(n, k, kappa, B) in decomposition::WITNESS_SCALING {
        // Throughput: number of ring elements in the witness being decomposed
        group.throughput(Throughput::Elements(n as u64));

        let param_label = format!("n={}_k={}_κ={}_B={}", n, k, kappa, B);

        group.bench_with_input(
            BenchmarkId::from_parameter(&param_label),
            &(n, k, kappa, B),
            |bencher, &(n, k, kappa, B)| {
                bencher.iter_batched(
                    || {
                        let mut rng = bench_rng();

                        let decomp = setup_decomp_input(n, k, kappa, B);

                        // Create A matrix for decomposition
                        let A = Matrix::<R>::rand(&mut rng, kappa, n);

                        (decomp, A)
                    },
                    |(decomp, A)| decomp.decompose(&A, B as u128),
                    BatchSize::SmallInput,
                );
            },
        );
    }

    group.finish();
}

/// Benchmark decomposition verifier
///
/// Measures verification time for decomposition proofs.
/// The verifier checks that the two output commitments correctly
/// decompose the input commitment: com(f) = com(F^(0)) + B·com(F^(1)).
fn bench_decomp_verifier(c: &mut Criterion) {
    let mut group = c.benchmark_group("Decomposition-Verifier");
    configure_benchmark_group(&mut group);

    for &(n, k, kappa, B) in decomposition::WITNESS_SCALING {
        group.throughput(Throughput::Elements(n as u64));

        let param_label = format!("n={}_k={}_κ={}_B={}", n, k, kappa, B);

        group.bench_with_input(
            BenchmarkId::from_parameter(&param_label),
            &(n, k, kappa, B),
            |bencher, &(n, k, kappa, B)| {
                // Generate proof once outside benchmark loop
                let (decomp, _outputs, proof) = setup_decomp_proof(n, k, kappa, B);

                bencher.iter(|| {
                    // Verify decomposition
                    proof.verify(&decomp.f, &decomp.r, B as u128)
                });
            },
        );
    }

    group.finish();
}

/// Benchmark fold-then-decompose roundtrip
///
/// This benchmark measures the complete IVC cycle:
/// 1. Start with 2 LinB instances (norm B)
/// 2. Fold them to 1 LinB2 instance (norm B²)
/// 3. Decompose back to 2 LinB instances (norm B)
fn bench_decomp_roundtrip(c: &mut Criterion) {
    let mut group = c.benchmark_group("Decomposition-Roundtrip");
    configure_benchmark_group(&mut group);

    // Use middle parameter set for roundtrip
    let (n, k, kappa, B) = decomposition::WITNESS_SCALING[1];

    group.throughput(Throughput::Elements((2 * n) as u64));

    group.bench_function(
        format!("fold_decompose_n={}_k={}_κ={}_B={}", n, k, kappa, B),
        |bencher| {
            bencher.iter_batched(
                || {
                    let mut rng = bench_rng();

                    // Create 2 LinB instances (this includes folding setup)
                    let decomp = setup_decomp_input(n, k, kappa, B);

                    // Create A matrix
                    let A = Matrix::<R>::rand(&mut rng, kappa, n);

                    (decomp, A)
                },
                |(decomp, A)| {
                    // Decompose: LinB2 (norm B²) → 2×LinB (norm B)
                    let ((linb0, linb1), _proof) = decomp.decompose(&A, B as u128);

                    // Note: Verification is tested separately in bench_decomp_verifier
                    // This roundtrip focuses on decomposition performance only

                    // Return both outputs to prevent dead code elimination
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
