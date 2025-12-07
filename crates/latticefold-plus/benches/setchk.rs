//! Benchmarks for Construction 4.2: Monomial set check
//!
//! The set check protocol verifies that matrices contain monomials
//! (exactly one non-zero entry per row/column). The protocol uses
//! sumcheck and batching for efficient verification.
//!
//! Paper reference: Section 4.2

#![allow(non_snake_case)]

use criterion::{
    criterion_group, criterion_main, BatchSize, BenchmarkGroup, BenchmarkId, Criterion, Throughput,
};
use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
use latticefold_plus::transcript::PoseidonTranscript;
use utils::{set_check, setup_setchk_input, setup_setchk_proof};

#[path = "utils/mod.rs"]
mod utils;

/// Configure benchmark group with benchmark settings
fn configure_benchmark_group(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    group.sample_size(10);
    group.measurement_time(std::time::Duration::from_secs(10));
    group.warm_up_time(std::time::Duration::from_secs(3));
}

/// Benchmark set check prover with varying parameters
///
/// Tests performance across different set sizes and batch counts:
/// - Set sizes: powers of 2 from 256 to 1024
/// - Batching: 1-4 sets checked simultaneously
fn bench_setchk_prover(c: &mut Criterion) {
    let mut group = c.benchmark_group("SetCheck-Prover");
    configure_benchmark_group(&mut group);

    for &(set_size, num_batches) in set_check::SET_SIZES {
        // Throughput: number of set elements checked
        group.throughput(Throughput::Elements((set_size * num_batches) as u64));

        let param_label = format!("size={}_batches={}", set_size, num_batches);

        group.bench_with_input(
            BenchmarkId::from_parameter(&param_label),
            &(set_size, num_batches),
            |bencher, &(set_size, num_batches)| {
                bencher.iter_batched(
                    || setup_setchk_input(set_size, num_batches),
                    |input| {
                        let mut ts = PoseidonTranscript::empty::<PC>();
                        input.set_check(&[], &mut ts)
                    },
                    BatchSize::SmallInput,
                );
            },
        );
    }

    group.finish();
}

/// Benchmark set check verifier
///
/// Measures verification time for set check proofs.
fn bench_setchk_verifier(c: &mut Criterion) {
    let mut group = c.benchmark_group("SetCheck-Verifier");
    configure_benchmark_group(&mut group);

    for &(set_size, num_batches) in set_check::SET_SIZES {
        group.throughput(Throughput::Elements((set_size * num_batches) as u64));

        let param_label = format!("size={}_batches={}", set_size, num_batches);

        group.bench_with_input(
            BenchmarkId::from_parameter(&param_label),
            &(set_size, num_batches),
            |bencher, &(set_size, num_batches)| {
                // Generate proof once outside benchmark loop
                let (_input, output) = setup_setchk_proof(set_size, num_batches);

                bencher.iter(|| {
                    let mut ts = PoseidonTranscript::empty::<PC>();
                    output.verify(&mut ts).unwrap()
                });
            },
        );
    }

    group.finish();
}

/// Benchmark batching efficiency
///
/// Measures how performance scales with number of batched sets.
/// Fixed: set_size=512. Varying: num_batches ∈ [1,2,4,8,16].
fn bench_setchk_batching(c: &mut Criterion) {
    let mut group = c.benchmark_group("SetCheck-Batching");
    configure_benchmark_group(&mut group);

    const SET_SIZE: usize = 512;
    const BATCH_COUNTS: [usize; 5] = [1, 2, 4, 8, 16];

    for num_batches in BATCH_COUNTS {
        group.throughput(Throughput::Elements((SET_SIZE * num_batches) as u64));

        group.bench_with_input(
            BenchmarkId::from_parameter(num_batches),
            &num_batches,
            |bencher, &num_batches| {
                bencher.iter_batched(
                    || setup_setchk_input(SET_SIZE, num_batches),
                    |input| {
                        let mut ts = PoseidonTranscript::empty::<PC>();
                        input.set_check(&[], &mut ts)
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
    bench_setchk_prover,
    bench_setchk_verifier,
    bench_setchk_batching,
);
criterion_main!(benches);
