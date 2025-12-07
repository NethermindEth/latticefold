//! Benchmarks for Construction 5.1: Single instance folding
//!
//! The single instance folding protocol (Π_lin,B) is a special case of the
//! multilinear folding protocol with L=1. It reduces one LinB instance with
//! norm bound B to one LinB2 instance with norm bound B².
//!
//! This construction serves as a baseline for measuring the amortization
//! benefits of the multilinear folding protocol (Construction 5.2) with L>1.
//!
//! Paper reference: Construction 5.1, Section 5.1

#![allow(non_snake_case)]

use criterion::{
    criterion_group, criterion_main, BatchSize, BenchmarkGroup, BenchmarkId, Criterion, Throughput,
};
use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
use latticefold_plus::{
    lin::LinParameters, mlin::Mlin, rgchk::DecompParameters, transcript::PoseidonTranscript,
};
use stark_rings::{cyclotomic_ring::models::frog_ring::RqPoly as R, PolyRing};
use stark_rings_linalg::{Matrix, SparseMatrix};
use utils::{bench_rng, setup_lin_input, setup_lin_proof, single_instance_fold};

#[path = "utils/mod.rs"]
mod utils;

/// Configure benchmark group with benchmark settings
fn configure_benchmark_group(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    group.sample_size(10);
    group.measurement_time(std::time::Duration::from_secs(10));
    group.warm_up_time(std::time::Duration::from_secs(3));
}

/// Benchmark single instance folding prover with varying parameters
///
/// Tests performance across different parameter combinations:
/// - n: witness size (length of witness vector after decomposition)
/// - k: decomposition width
/// - κ (kappa): number of commitment rows
/// - B: norm bound parameter
///
/// This construction is implemented as a wrapper that calls mlin() with L=1.
fn bench_lin_prover(c: &mut Criterion) {
    let mut group = c.benchmark_group("SingleInstanceFolding-Prover");
    configure_benchmark_group(&mut group);

    for &(n, k, kappa, B) in single_instance_fold::WITNESS_SCALING {
        // Throughput: number of ring elements in the witness
        group.throughput(Throughput::Elements(n as u64));

        let param_label = format!("n={}_k={}_κ={}_B={}", n, k, kappa, B);

        group.bench_with_input(
            BenchmarkId::from_parameter(&param_label),
            &(n, k, kappa, B),
            |bencher, &(n, k, kappa, B)| {
                bencher.iter_batched(
                    || {
                        let mut rng = bench_rng();

                        let linb = setup_lin_input(n, k, kappa, B);

                        // Create M matrix
                        let mut m = SparseMatrix::identity(n);
                        m.coeffs[0][0].0 = 2u128.into();
                        let M = vec![m];

                        // Create A matrix for folding
                        let A = Matrix::<R>::rand(&mut rng, kappa, n);

                        (linb, A, M)
                    },
                    |(linb, A, M)| {

                        let b = (R::dimension() / 2) as u128;
                        let l = ((<<R as stark_rings::PolyRing>::BaseRing as ark_ff::PrimeField>::MODULUS.0[0] as f64).ln()
                            / ((R::dimension() / 2) as f64).ln())
                            .ceil() as usize;

                        let params = LinParameters {
                            kappa,
                            decomp: DecompParameters { b, k, l },
                        };

                        let mut ts = PoseidonTranscript::empty::<PC>();
                        linb.lin(&A, &M, &params, &mut ts)
                    },
                    BatchSize::SmallInput,
                );
            },
        );
    }

    group.finish();
}

/// Benchmark single instance folding verifier
///
/// Measures verification time for single instance folding proofs.
/// The verifier checks the commitment transformation proof.
fn bench_lin_verifier(c: &mut Criterion) {
    let mut group = c.benchmark_group("SingleInstanceFolding-Verifier");
    configure_benchmark_group(&mut group);

    for &(n, k, kappa, B) in single_instance_fold::WITNESS_SCALING {
        group.throughput(Throughput::Elements(n as u64));

        let param_label = format!("n={}_k={}_κ={}_B={}", n, k, kappa, B);

        group.bench_with_input(
            BenchmarkId::from_parameter(&param_label),
            &(n, k, kappa, B),
            |bencher, &(n, k, kappa, B)| {
                // Generate proof once outside benchmark loop
                let (_linb, _linb2, proof) = setup_lin_proof(n, k, kappa, B);

                // Create matrix M (same pattern as prover)
                let mut m = SparseMatrix::identity(n);
                m.coeffs[0][0].0 = 2u128.into();
                let M = vec![m];

                bencher.iter(|| {
                    let mut ts = PoseidonTranscript::empty::<PC>();
                    proof.verify(&M, &mut ts).unwrap()
                });
            },
        );
    }

    group.finish();
}

/// Benchmark comparison: lin() vs mlin() with L=1
///
/// This benchmark verifies that Construction 5.1 (single instance folding)
/// is correctly implemented as a wrapper around Construction 5.2 (multilinear
/// folding) with L=1. Both should have nearly identical performance.
fn bench_lin_vs_mlin_L1(c: &mut Criterion) {
    let mut group = c.benchmark_group("SingleInstanceFolding-Comparison");
    configure_benchmark_group(&mut group);

    // Use middle parameter set for comparison
    let (n, k, kappa, B) = single_instance_fold::WITNESS_SCALING[1];

    group.throughput(Throughput::Elements(n as u64));

    // Benchmark lin() wrapper
    group.bench_function("lin_wrapper", |bencher| {
        bencher.iter_batched(
            || {
                let mut rng = bench_rng();

                let linb = setup_lin_input(n, k, kappa, B);

                let mut m = SparseMatrix::identity(n);
                m.coeffs[0][0].0 = 2u128.into();
                let M = vec![m];

                let A = Matrix::<R>::rand(&mut rng, kappa, n);

                (linb, A, M)
            },
            |(linb, A, M)| {
                let b = (R::dimension() / 2) as u128;
                let l = ((<<R as stark_rings::PolyRing>::BaseRing as ark_ff::PrimeField>::MODULUS.0
                    [0] as f64)
                    .ln()
                    / ((R::dimension() / 2) as f64).ln())
                .ceil() as usize;

                let params = LinParameters {
                    kappa,
                    decomp: DecompParameters { b, k, l },
                };

                let mut ts = PoseidonTranscript::empty::<PC>();
                linb.lin(&A, &M, &params, &mut ts)
            },
            BatchSize::SmallInput,
        );
    });

    // Benchmark mlin() with L=1
    group.bench_function("mlin_L1", |bencher| {
        bencher.iter_batched(
            || {
                let mut rng = bench_rng();

                let linb = setup_lin_input(n, k, kappa, B);

                let mut m = SparseMatrix::identity(n);
                m.coeffs[0][0].0 = 2u128.into();
                let M = vec![m];

                let A = Matrix::<R>::rand(&mut rng, kappa, n);

                (linb, A, M)
            },
            |(linb, A, M)| {
                let b = (R::dimension() / 2) as u128;
                let l = ((<<R as stark_rings::PolyRing>::BaseRing as ark_ff::PrimeField>::MODULUS.0
                    [0] as f64)
                    .ln()
                    / ((R::dimension() / 2) as f64).ln())
                .ceil() as usize;

                let params = LinParameters {
                    kappa,
                    decomp: DecompParameters { b, k, l },
                };

                let mlin = Mlin {
                    lins: vec![linb], // L=1
                    params,
                };

                let mut ts = PoseidonTranscript::empty::<PC>();
                mlin.mlin(&A, &M, &mut ts)
            },
            BatchSize::SmallInput,
        );
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_lin_prover,
    bench_lin_verifier,
    bench_lin_vs_mlin_L1,
);
criterion_main!(benches);
