//! Benchmarks for Construction 4.3-4.4: Range check
//!
//! The range check protocol verifies that committed ring element coefficients
//! lie within a specified range (-B, B) using double commitments and gadget
//! decomposition for efficiency.
//!
//! Paper reference: Sections 4.3-4.4

#![allow(non_snake_case)]

use criterion::{
    criterion_group, criterion_main, BatchSize, BenchmarkId, BenchmarkGroup, Criterion, Throughput,
};
use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
use latticefold_plus::transcript::PoseidonTranscript;

#[path = "utils/mod.rs"]
mod utils;

use utils::{quick, setup_rgchk_input, setup_rgchk_proof};

/// Configure benchmark group with benchmark settings
fn configure_benchmark_group(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    group.sample_size(10);
    group.measurement_time(std::time::Duration::from_secs(10));
    group.warm_up_time(std::time::Duration::from_secs(3));
}

/// Benchmark range check prover with varying parameters
///
/// Tests performance across different parameter combinations:
/// - witness_size: number of ring elements in witness
/// - k: decomposition width (determines range B = (d/2)^k)
/// - kappa: number of commitment rows (security parameter)
fn bench_rgchk_prover(c: &mut Criterion) {
    let mut group = c.benchmark_group("RangeCheck-Prover");
    configure_benchmark_group(&mut group);

    for &(witness_size, k, kappa) in quick::RGCHK {
        // Throughput: number of ring elements range-checked
        group.throughput(Throughput::Elements(witness_size as u64));

        let param_label = format!("w={}_k={}_κ={}", witness_size, k, kappa);

        group.bench_with_input(
            BenchmarkId::from_parameter(&param_label),
            &(witness_size, k, kappa),
            |bencher, &(witness_size, k, kappa)| {
                bencher.iter_batched(
                    || setup_rgchk_input(witness_size, k, kappa),
                    |input| {
                        let mut ts = PoseidonTranscript::empty::<PC>();
                        input.range_check(&[], &mut ts)
                    },
                    BatchSize::SmallInput,
                );
            },
        );
    }

    group.finish();
}

/// Benchmark range check verifier
///
/// Measures verification time for range check proofs.
fn bench_rgchk_verifier(c: &mut Criterion) {
    let mut group = c.benchmark_group("RangeCheck-Verifier");
    configure_benchmark_group(&mut group);

    for &(witness_size, k, kappa) in quick::RGCHK {
        group.throughput(Throughput::Elements(witness_size as u64));

        let param_label = format!("w={}_k={}_κ={}", witness_size, k, kappa);

        group.bench_with_input(
            BenchmarkId::from_parameter(&param_label),
            &(witness_size, k, kappa),
            |bencher, &(witness_size, k, kappa)| {
                // Generate proof once outside benchmark loop
                let (_input, output) = setup_rgchk_proof(witness_size, k, kappa);

                bencher.iter(|| {
                    let mut ts = PoseidonTranscript::empty::<PC>();
                    output.verify(&mut ts).unwrap()
                });
            },
        );
    }

    group.finish();
}

/// Benchmark scaling with decomposition width
///
/// Measures how performance scales with the decomposition parameter k,
/// which determines the range B = (d/2)^k.
/// witness_size=65536, kappa=2. Varying: k ∈ [2,3,4,5].
fn bench_rgchk_scaling_k(c: &mut Criterion) {
    let mut group = c.benchmark_group("RangeCheck-Scaling-K");
    configure_benchmark_group(&mut group);

    const WITNESS_SIZE: usize = 65536;
    const KAPPA: usize = 2;
    const K_VALUES: [usize; 4] = [2, 3, 4, 5];

    for k in K_VALUES {
        group.throughput(Throughput::Elements(WITNESS_SIZE as u64));

        group.bench_with_input(
            BenchmarkId::from_parameter(k),
            &k,
            |bencher, &k| {
                bencher.iter_batched(
                    || setup_rgchk_input(WITNESS_SIZE, k, KAPPA),
                    |input| {
                        let mut ts = PoseidonTranscript::empty::<PC>();
                        input.range_check(&[], &mut ts)
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
    bench_rgchk_prover,
    bench_rgchk_verifier,
    bench_rgchk_scaling_k,
);
criterion_main!(benches);