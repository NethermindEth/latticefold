//! Benchmarks for Construction 5.2: Multilinear folding
//!
//! The multilinear folding protocol (Π_mlin,L,B) folds L linearized R1CS
//! instances into a single folded instance with aggregated commitments and
//! witnesses. This provides amortization benefits as L increases.
//!
//! Paper reference: Construction 5.2, Section 5.2

#![allow(non_snake_case)]

use criterion::{
    criterion_group, criterion_main, BatchSize, BenchmarkGroup, BenchmarkId, Criterion, Throughput,
};
use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
use latticefold_plus::transcript::PoseidonTranscript;
use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as R;
use stark_rings_linalg::{Matrix, SparseMatrix};
use utils::{bench_rng, multilinear_fold, setup_mlin_input, setup_mlin_proof};

#[path = "utils/mod.rs"]
mod utils;

/// Configure benchmark group with benchmark settings
fn configure_benchmark_group(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    group.sample_size(10);
    group.measurement_time(std::time::Duration::from_secs(10));
    group.warm_up_time(std::time::Duration::from_secs(3));
}

/// Benchmark multilinear folding prover
///
/// Fixed: n=65536, k=2, kappa=2, B=50. Varying: L ∈ [2,3,4,5,6,7,8].
fn bench_mlin_prover(c: &mut Criterion) {
    let mut group = c.benchmark_group("MultilinearFolding-Prover");
    configure_benchmark_group(&mut group);

    for &(L, n, k, kappa, B) in multilinear_fold::FOLDING_ARITY {
        group.throughput(Throughput::Elements((L * n) as u64));

        let param_label = format!("L={}_n={}_k={}_κ={}_B={}", L, n, k, kappa, B);

        group.bench_with_input(
            BenchmarkId::from_parameter(&param_label),
            &(L, n, k, kappa, B),
            |bencher, &(L, n, k, kappa, B)| {
                bencher.iter_batched(
                    || {
                        let mut rng = bench_rng();
                        let mlin = setup_mlin_input(L, n, k, kappa, B);
                        let mut m = SparseMatrix::identity(n);
                        m.coeffs[0][0].0 = 2u128.into();
                        let M = vec![m];

                        // Create A matrix for folding
                        let A = Matrix::<R>::rand(&mut rng, kappa, n);
                        (mlin, A, M)
                    },
                    |(mlin, A, M)| {
                        let mut ts = PoseidonTranscript::empty::<PC>();
                        mlin.mlin(&A, &M, &mut ts)
                    },
                    BatchSize::SmallInput,
                );
            },
        );
    }

    group.finish();
}

/// Benchmark multilinear folding verifier
///
/// Fixed: n=65536, k=2, kappa=2, B=50. Varying: L ∈ [2,3,4,5,6,7,8].
fn bench_mlin_verifier(c: &mut Criterion) {
    let mut group = c.benchmark_group("MultilinearFolding-Verifier");
    configure_benchmark_group(&mut group);

    for &(L, n, k, kappa, B) in multilinear_fold::FOLDING_ARITY {
        group.throughput(Throughput::Elements((L * n) as u64));

        let param_label = format!("L={}_n={}_k={}_κ={}_B={}", L, n, k, kappa, B);

        group.bench_with_input(
            BenchmarkId::from_parameter(&param_label),
            &(L, n, k, kappa, B),
            |bencher, &(L, n, k, kappa, B)| {
                // Generate proof once outside benchmark loop
                let (_mlin, _linb2, proof) = setup_mlin_proof(L, n, k, kappa, B);

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

/// Benchmark decomposition width (k) scaling
///
/// Fixed: L=4, kappa=2, B=50. Varying: k ∈ [2,3,4] with adjusted n.
fn bench_mlin_k_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("MultilinearFolding-KScaling");
    configure_benchmark_group(&mut group);

    for &(L, n, k, kappa, B) in multilinear_fold::K_SCALING {
        group.throughput(Throughput::Elements((L * n) as u64));

        let param_label = format!("k={}_n={}", k, n);

        group.bench_with_input(
            BenchmarkId::from_parameter(&param_label),
            &(L, n, k, kappa, B),
            |bencher, &(L, n, k, kappa, B)| {
                bencher.iter_batched(
                    || {
                        let mut rng = bench_rng();
                        let mlin = setup_mlin_input(L, n, k, kappa, B);
                        let mut m = SparseMatrix::identity(n);
                        m.coeffs[0][0].0 = 2u128.into();
                        let M = vec![m];
                        let A = Matrix::<R>::rand(&mut rng, kappa, n);
                        (mlin, A, M)
                    },
                    |(mlin, A, M)| {
                        let mut ts = PoseidonTranscript::empty::<PC>();
                        mlin.mlin(&A, &M, &mut ts)
                    },
                    BatchSize::SmallInput,
                );
            },
        );
    }

    group.finish();
}

/// Benchmark large witness scaling
///
/// Fixed: L=4, k=2, kappa=2, B=50. Varying: n ∈ [128K, 256K, 512K].
fn bench_mlin_large_witness(c: &mut Criterion) {
    let mut group = c.benchmark_group("MultilinearFolding-LargeWitness");
    configure_benchmark_group(&mut group);

    for &(L, n, k, kappa, B) in multilinear_fold::LARGE_WITNESS {
        group.throughput(Throughput::Elements((L * n) as u64));

        let param_label = format!("n={}", n);

        group.bench_with_input(
            BenchmarkId::from_parameter(&param_label),
            &(L, n, k, kappa, B),
            |bencher, &(L, n, k, kappa, B)| {
                bencher.iter_batched(
                    || {
                        let mut rng = bench_rng();
                        let mlin = setup_mlin_input(L, n, k, kappa, B);
                        let mut m = SparseMatrix::identity(n);
                        m.coeffs[0][0].0 = 2u128.into();
                        let M = vec![m];
                        let A = Matrix::<R>::rand(&mut rng, kappa, n);
                        (mlin, A, M)
                    },
                    |(mlin, A, M)| {
                        let mut ts = PoseidonTranscript::empty::<PC>();
                        mlin.mlin(&A, &M, &mut ts)
                    },
                    BatchSize::SmallInput,
                );
            },
        );
    }

    group.finish();
}

/// Benchmark security parameter (kappa) scaling
///
/// Fixed: L=4, n=65536, k=2, B=50. Varying: kappa ∈ [2,3,4,5].
fn bench_mlin_kappa_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("MultilinearFolding-KappaScaling");
    configure_benchmark_group(&mut group);

    for &(L, n, k, kappa, B) in multilinear_fold::KAPPA_SCALING {
        group.throughput(Throughput::Elements((L * n) as u64));

        group.bench_with_input(
            BenchmarkId::from_parameter(kappa),
            &(L, n, k, kappa, B),
            |bencher, &(L, n, k, kappa, B)| {
                bencher.iter_batched(
                    || {
                        let mut rng = bench_rng();
                        let mlin = setup_mlin_input(L, n, k, kappa, B);
                        let mut m = SparseMatrix::identity(n);
                        m.coeffs[0][0].0 = 2u128.into();
                        let M = vec![m];
                        let A = Matrix::<R>::rand(&mut rng, kappa, n);
                        (mlin, A, M)
                    },
                    |(mlin, A, M)| {
                        let mut ts = PoseidonTranscript::empty::<PC>();
                        mlin.mlin(&A, &M, &mut ts)
                    },
                    BatchSize::SmallInput,
                );
            },
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_mlin_prover,
    bench_mlin_verifier,
    bench_mlin_k_scaling,
    bench_mlin_large_witness,
    bench_mlin_kappa_scaling,
);
criterion_main!(benches);
