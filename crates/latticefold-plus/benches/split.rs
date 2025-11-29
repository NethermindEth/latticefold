//! Benchmarks for Construction 4.1: Split function
//!
//! The split function performs gadget decomposition on double commitments,
//! converting ring element matrices into base ring scalar representations.
//!
//! Paper reference: Section 4.1
//!
//! Benchmarks the second decomposition step where commitment matrices are
//! split into base ring elements (not witness matrices).

#![allow(non_snake_case)]

use ark_ff::PrimeField;
use criterion::{
    criterion_group, criterion_main, BatchSize, BenchmarkId, BenchmarkGroup, Criterion, Throughput,
};
use latticefold_plus::utils::split;
use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as R;
use stark_rings::PolyRing;

#[path = "utils/mod.rs"]
mod utils;

use utils::{quick, setup_split_input};

/// Decomposition parameters for the split operation
struct DecompParams {
    /// Decomposition base: d/2 where d is ring dimension
    b: u128,
    /// Decomposition width: ⌈log_{b}(q)⌉ where q is the modulus
    l: usize,
}

impl DecompParams {
    /// Compute decomposition parameters for the current ring
    fn compute() -> Self {
        let d = R::dimension();
        let b = (d / 2) as u128;

        // Compute l = ⌈log_{b}(q)⌉ where q is the base ring modulus
        // This determines the width of the second decomposition
        let modulus = <<R as PolyRing>::BaseRing>::MODULUS.0[0] as f64;
        let base = b as f64;
        let l = (modulus.ln() / base.ln()).ceil() as usize;

        Self { b, l }
    }
}

/// Configure benchmark group with standard cryptographic benchmark settings
fn configure_benchmark_group(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    group.sample_size(50);
    group.measurement_time(std::time::Duration::from_secs(10));
    group.warm_up_time(std::time::Duration::from_secs(3));
}

/// Benchmark split with varying witness sizes and commitment parameters
///
/// Tests performance across different parameter combinations:
/// - Input: commitment matrices of size kappa × (k_first * d)
/// - Output: witness_size base ring elements (padded)
/// - Parameters: k_first ∈ [2,8], kappa ∈ [1,4], witness_size ∈ [16K,128K]
fn bench_split_varying_params(c: &mut Criterion) {
    let mut group = c.benchmark_group("Split-VaryingParams");
    configure_benchmark_group(&mut group);

    let params = DecompParams::compute();

    for &(witness_size, k_first, kappa) in quick::SPLIT {
        group.throughput(Throughput::Elements(witness_size as u64));

        let param_label = format!("w={}_k={}_κ={}", witness_size, k_first, kappa);

        group.bench_with_input(
            BenchmarkId::from_parameter(&param_label),
            &(witness_size, k_first, kappa),
            |bencher, &(witness_size, k_first, kappa)| {
                bencher.iter_batched(
                    || setup_split_input(k_first, kappa),
                    |com| split(&com, witness_size, params.b, params.l),
                    BatchSize::SmallInput,
                );
            },
        );
    }

    group.finish();
}

/// Benchmark scaling with first decomposition width
///
/// Measures how performance scales with input matrix width (k_first).
/// Fixed: witness_size=131K, kappa=2. Varying: k_first ∈ [2,4,6,8].
fn bench_split_scaling_k_first(c: &mut Criterion) {
    let mut group = c.benchmark_group("Split-Scaling-KFirst");
    configure_benchmark_group(&mut group);

    let params = DecompParams::compute();

    // Fixed parameters
    const WITNESS_SIZE: usize = 131072;
    const KAPPA: usize = 2;

    // Vary k_first to measure scaling with input matrix width
    const K_VALUES: [usize; 4] = [2, 4, 6, 8];

    for k_first in K_VALUES {
        group.throughput(Throughput::Elements(WITNESS_SIZE as u64));

        group.bench_with_input(
            BenchmarkId::from_parameter(k_first),
            &k_first,
            |bencher, &k_first| {
                bencher.iter_batched(
                    || setup_split_input(k_first, KAPPA),
                    |com| split(&com, WITNESS_SIZE, params.b, params.l),
                    BatchSize::SmallInput,
                );
            },
        );
    }

    group.finish();
}

/// Benchmark scaling with commitment rows
///
/// Measures how performance scales with the security parameter (kappa).
/// Fixed: witness_size=131K, k_first=4. Varying: kappa ∈ [1,2,3,4].
fn bench_split_scaling_kappa(c: &mut Criterion) {
    let mut group = c.benchmark_group("Split-Scaling-Kappa");
    configure_benchmark_group(&mut group);

    let params = DecompParams::compute();

    // Fixed parameters
    const WITNESS_SIZE: usize = 131072;
    const K_FIRST: usize = 4;

    // Vary kappa to measure scaling with security parameter
    const KAPPA_VALUES: [usize; 4] = [1, 2, 3, 4];

    for kappa in KAPPA_VALUES {
        group.throughput(Throughput::Elements(WITNESS_SIZE as u64));

        group.bench_with_input(
            BenchmarkId::from_parameter(kappa),
            &kappa,
            |bencher, &kappa| {
                bencher.iter_batched(
                    || setup_split_input(K_FIRST, kappa),
                    |com| split(&com, WITNESS_SIZE, params.b, params.l),
                    BatchSize::SmallInput,
                );
            },
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_split_varying_params,
    bench_split_scaling_k_first,
    bench_split_scaling_kappa,
);
criterion_main!(benches);
