//! Benchmarks for Construction 4.5: Commitment transformation
//!
//! The commitment transformation protocol converts double commitments from
//! range check into folded commitments suitable for the main LatticeFold+
//! protocol. Uses sumcheck protocols for verification.
//!
//! Paper reference: Section 4.5

#![allow(non_snake_case)]

use criterion::{
    criterion_group, criterion_main, BatchSize, BenchmarkId, BenchmarkGroup, Criterion, Throughput,
};
use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
use latticefold_plus::transcript::PoseidonTranscript;
use stark_rings_linalg::SparseMatrix;
use utils::{commitment_transform, setup_cm_input, setup_cm_proof};

#[path = "utils/mod.rs"]
mod utils;

/// Configure benchmark group with benchmark settings
fn configure_benchmark_group(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    group.sample_size(10);
    group.measurement_time(std::time::Duration::from_secs(10));
    group.warm_up_time(std::time::Duration::from_secs(3));
}

/// Benchmark commitment transformation prover
///
/// Tests performance across different folding arities L.
/// Fixed: witness_size=65536, k=2, kappa=2. Varying: L ∈ [2,3,4,5,6,7,8].
fn bench_cm_prover(c: &mut Criterion) {
    let mut group = c.benchmark_group("CommitmentTransform-Prover");
    configure_benchmark_group(&mut group);

    for &(L, witness_size, k, kappa) in commitment_transform::FOLDING_ARITY {
        group.throughput(Throughput::Elements(L as u64));

        let param_label = format!("L={}_w={}_k={}_κ={}", L, witness_size, k, kappa);

        group.bench_with_input(
            BenchmarkId::from_parameter(&param_label),
            &(L, witness_size, k, kappa),
            |bencher, &(L, witness_size, k, kappa)| {
                bencher.iter_batched(
                    || {
                        let input = setup_cm_input(L, witness_size, k, kappa);
                        // Create modified identity matrix M
                        let mut m = SparseMatrix::identity(witness_size);
                        m.coeffs[0][0].0 = 2u128.into();
                        let M = vec![m];
                        (input, M)
                    },
                    |(input, M)| {
                        let mut ts = PoseidonTranscript::empty::<PC>();
                        input.prove(&M, &mut ts)
                    },
                    BatchSize::SmallInput,
                );
            },
        );
    }

    group.finish();
}

/// Benchmark commitment transformation verifier
///
/// Tests verification performance across different folding arities L.
/// Fixed: witness_size=65536, k=2, kappa=2. Varying: L ∈ [2,3,4,5,6,7,8].
fn bench_cm_verifier(c: &mut Criterion) {
    let mut group = c.benchmark_group("CommitmentTransform-Verifier");
    configure_benchmark_group(&mut group);

    for &(L, witness_size, k, kappa) in commitment_transform::FOLDING_ARITY {
        group.throughput(Throughput::Elements(L as u64));

        let param_label = format!("L={}_w={}_k={}_κ={}", L, witness_size, k, kappa);

        group.bench_with_input(
            BenchmarkId::from_parameter(&param_label),
            &(L, witness_size, k, kappa),
            |bencher, &(L, witness_size, k, kappa)| {
                // Generate proof once outside benchmark loop
                let (_input, proof) = setup_cm_proof(L, witness_size, k, kappa);

                // Create matrix M (same pattern as prover)
                let mut m = SparseMatrix::identity(witness_size);
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

criterion_group!(benches, bench_cm_prover, bench_cm_verifier,);
criterion_main!(benches);
