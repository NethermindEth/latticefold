//! Benchmarks for Construction 5.2: Multilinear folding
//!
//! The multilinear folding protocol (Π_mlin,L,B) folds L linearized R1CS
//! instances into a single folded instance with aggregated commitments and
//! witnesses. This provides amortization benefits as L increases.
//!
//! Paper reference: Construction 5.2, Section 5.2

#![allow(non_snake_case)]

use criterion::{
    criterion_group, criterion_main, BatchSize, BenchmarkId, BenchmarkGroup, Criterion, Throughput,
};
use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
use latticefold_plus::transcript::PoseidonTranscript;
use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as R;
use stark_rings_linalg::{Matrix, SparseMatrix};

#[path = "utils/mod.rs"]
mod utils;

use utils::{bench_rng, quick, setup_mlin_input, setup_mlin_proof};

/// Configure benchmark group with benchmark settings
fn configure_benchmark_group(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    group.sample_size(10);
    group.measurement_time(std::time::Duration::from_secs(10));
    group.warm_up_time(std::time::Duration::from_secs(3));
}

/// Benchmark multilinear folding prover with varying parameters
///
/// Tests performance across different parameter combinations:
/// - L: number of instances to fold (higher L = better amortization)
/// - n: witness size (length of witness vector after decomposition)
/// - k: decomposition width
/// - κ (kappa): number of commitment rows
/// - B: norm bound parameter
fn bench_mlin_prover(c: &mut Criterion) {
    let mut group = c.benchmark_group("MultilinearFolding-Prover");
    configure_benchmark_group(&mut group);

    for &(L, n, k, kappa, B) in quick::MLIN {
        // Throughput: number of ring elements transformed across L instances
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

                        // Create M matrix
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
/// Measures verification time for multilinear folding proofs.
/// The verifier checks the commitment transformation proof.
fn bench_mlin_verifier(c: &mut Criterion) {
    let mut group = c.benchmark_group("MultilinearFolding-Verifier");
    configure_benchmark_group(&mut group);

    for &(L, n, k, kappa, B) in quick::MLIN {
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

/// Benchmark scaling with folding arity
///
/// Measures how performance scales with L (number of instances to fold).
/// Paper recommends L=8 for optimal amortization.
/// Fixed: n=65536, k=2, κ=2, B=50. Varying: L ∈ [2,4,8].
fn bench_mlin_scaling_L(c: &mut Criterion) {
    let mut group = c.benchmark_group("MultilinearFolding-Scaling-L");
    configure_benchmark_group(&mut group);

    const N: usize = 65536;
    const K: usize = 2;
    const KAPPA: usize = 2;
    const B: usize = 50;
    const L_VALUES: [usize; 3] = [2, 4, 8];

    for L in L_VALUES {
        group.throughput(Throughput::Elements((L * N) as u64));

        group.bench_with_input(
            BenchmarkId::from_parameter(L),
            &L,
            |bencher, &L| {
                bencher.iter_batched(
                    || {
                        let mut rng = bench_rng();

                        let mlin = setup_mlin_input(L, N, K, KAPPA, B);

                        // Create M matrix (modified identity)
                        let mut m = SparseMatrix::identity(N);
                        m.coeffs[0][0].0 = 2u128.into();
                        let M = vec![m];

                        // Create A matrix for folding
                        let A = Matrix::<R>::rand(&mut rng, KAPPA, N);

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
    bench_mlin_scaling_L,
);
criterion_main!(benches);
