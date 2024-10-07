#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

use ark_std::{time::Duration, UniformRand};
use criterion::{criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion, PlotConfiguration};
use latticefold::parameters::{BabyBearParams, FrogParams, GoldilocksParams, StarkPrimeParams};
use latticefold::{
    commitment::AjtaiCommitmentScheme,
    parameters::{DecompositionParams, DilithiumTestParams, DILITHIUM_PRIME},
};
use std::fmt::Debug;

use cyclotomic_rings::{
    BabyBearRingNTT, FrogRingNTT, GoldilocksRingNTT, StarkRingNTT, SuitableRing,
};
use lattirust_ring::cyclotomic_ring::models::pow2_debug::Pow2CyclotomicPolyRingNTT;
use rand::thread_rng;
/*
fn generalized_ajtai_benchmark<
    const C: usize,
    const W: usize,
    R: Clone + UniformRand + Debug + SuitableRing,
    P: DecompositionParams + Clone,
>(
    c: &mut Criterion,
    ring_name: &str,
) where
    R: for<'a> std::ops::AddAssign<&'a R>,
{
    let mut rng = thread_rng();

    let ajtai_data: AjtaiCommitmentScheme<C, W, R> = AjtaiCommitmentScheme::rand(&mut rng);
    let witness: Vec<R> = (0..W).map(|_| R::rand(&mut rng)).collect();

    let ajtai_data_2 = ajtai_data.clone();
    let witness_2 = witness.clone();

    let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    let mut group = c.benchmark_group(format!("Ajtai {} C={} W={}", ring_name, C, W));
    group.plot_config(plot_config);

    group.bench_with_input(
        "CommitNTT",
        &(ajtai_data, witness),
        |b, (ajtai_data, witness)| {
            b.iter(|| {
                let _ = ajtai_data.commit_ntt(witness);
            })
        },
    );

    group.bench_with_input(
        "DecomposeCommitNTT",
        &(ajtai_data_2, witness_2),
        |b, (ajtai_data, witness)| {
            b.iter(|| {
                let _ = ajtai_data.decompose_and_commit_ntt::<P>(witness);
            })
        },
    );

    group.finish();
}

macro_rules! run_ajtai_benchmarks {
    ($c:expr, $cw: expr, $($w:expr),+) => {
        $(
            // StarkPrime
            generalized_ajtai_benchmark::<$cw, $w, StarkRingNTT, StarkPrimeParams>($c, "StarkPrime");
            // Goldilocks
            generalized_ajtai_benchmark::<$cw, $w, GoldilocksRingNTT, GoldilocksParams>($c, "Goldilocks");
            // BabyBear
            generalized_ajtai_benchmark::<$cw, $w, BabyBearRingNTT, BabyBearParams>($c, "BabyBear");
            // Frog
            generalized_ajtai_benchmark::<$cw, $w, FrogRingNTT, FrogParams>($c, "Frog");
            // Dilithium
            generalized_ajtai_benchmark::<$cw, $w, Pow2CyclotomicPolyRingNTT<DILITHIUM_PRIME, 256>, DilithiumTestParams>($c, "Dilithium");
        )+
    };
}

fn ajtai_benchmarks(c: &mut Criterion) {
    run_ajtai_benchmarks!(c, 5, { 1 << 16 });
    run_ajtai_benchmarks!(c, 6, { 1 << 16 });
    run_ajtai_benchmarks!(c, 7, { 1 << 16 });
    run_ajtai_benchmarks!(c, 8, { 1 << 16 });
    run_ajtai_benchmarks!(c, 9, { 1 << 16 });
    run_ajtai_benchmarks!(c, 10, { 1 << 16 });
    run_ajtai_benchmarks!(c, 11, { 1 << 16 });
    run_ajtai_benchmarks!(c, 12, { 1 << 16 });
    run_ajtai_benchmarks!(c, 13, { 1 << 16 });
    run_ajtai_benchmarks!(c, 14, { 1 << 16 });
    run_ajtai_benchmarks!(c, 15, { 1 << 16 });
}

pub fn benchmarks_main(c: &mut Criterion) {
    ajtai_benchmarks(c);
}*/

fn ajtai_benchmark<
    const C: usize,
    const W: usize,
    R: Clone + UniformRand + Debug + SuitableRing,
    P: DecompositionParams + Clone,
>(
    group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>,
) where
    R: for<'a> std::ops::AddAssign<&'a R>,
{
    let mut rng = thread_rng();

    let ajtai_data: AjtaiCommitmentScheme<C, W, R> = AjtaiCommitmentScheme::rand(&mut rng);
    let witness: Vec<R> = (0..W).map(|_| R::rand(&mut rng)).collect();

    let ajtai_data_2 = ajtai_data.clone();
    let witness_2 = witness.clone();

    group.bench_with_input(
        BenchmarkId::new("CommitNTT",
                         format!("C={},W={}", C, W)),
        &(ajtai_data, witness),
        |b, (ajtai_data, witness)| {
            b.iter(|| {
                let _ = ajtai_data.commit_ntt(witness);
            })
        },
    );

    group.bench_with_input(
        BenchmarkId::new("DecomposeCommitNTT",
                         format!("C={},W={}", C, W)),
        &(ajtai_data_2, witness_2),
        |b, (ajtai_data, witness)| {
            b.iter(|| {
                let _ = ajtai_data.decompose_and_commit_ntt::<P>(witness);
            })
        },
    );
}

fn ajtai_benchmarks(c: &mut Criterion) {
    let plot_config = PlotConfiguration::default()
        .summary_scale(AxisScale::Logarithmic);

    let c_values = [5, 6, 7, 8, 9, 10];
    let w_values = [{ 1 << 16 }];

    // StarkPrime benchmarks
    {
        let mut group = c.benchmark_group("Ajtai StarkPrime");
        group.plot_config(plot_config.clone());
        for &c in &c_values {
            for &w in &w_values {
                match (c, w) {
                    (5, 65536) => ajtai_benchmark::<5, 65536, StarkRingNTT, StarkPrimeParams>(&mut group),
                    (6, 65536) => ajtai_benchmark::<6, 65536, StarkRingNTT, StarkPrimeParams>(&mut group),
                    (7, 65536) => ajtai_benchmark::<7, 65536, StarkRingNTT, StarkPrimeParams>(&mut group),
                    (8, 65536) => ajtai_benchmark::<8, 65536, StarkRingNTT, StarkPrimeParams>(&mut group),
                    (9, 65536) => ajtai_benchmark::<9, 65536, StarkRingNTT, StarkPrimeParams>(&mut group),
                    (10, 65536) => ajtai_benchmark::<10, 65536, StarkRingNTT, StarkPrimeParams>(&mut group),
                    _ => {} // Skip other combinations
                }
            }
        }
        group.finish();
    }

    // Goldilocks benchmarks
    {
        let mut group = c.benchmark_group("Ajtai Goldilocks");
        group.plot_config(plot_config.clone());
        for &c in &c_values {
            for &w in &w_values {
                match (c, w) {
                    (5, 65536) => ajtai_benchmark::<5, 65536, GoldilocksRingNTT, GoldilocksParams>(&mut group),
                    (6, 65536) => ajtai_benchmark::<6, 65536, GoldilocksRingNTT, GoldilocksParams>(&mut group),
                    (7, 65536) => ajtai_benchmark::<7, 65536, GoldilocksRingNTT, GoldilocksParams>(&mut group),
                    (8, 65536) => ajtai_benchmark::<8, 65536, GoldilocksRingNTT, GoldilocksParams>(&mut group),
                    (9, 65536) => ajtai_benchmark::<9, 65536, GoldilocksRingNTT, GoldilocksParams>(&mut group),
                    (10, 65536) => ajtai_benchmark::<10, 65536, GoldilocksRingNTT, GoldilocksParams>(&mut group),
                    _ => {} // Skip other combinations
                }
            }
        }
        group.finish();
    }

    // BabyBear benchmarks
    {
        let mut group = c.benchmark_group("Ajtai BabyBear");
        group.plot_config(plot_config.clone());
        for &c in &c_values {
            for &w in &w_values {
                match (c, w) {
                    (5, 65536) => ajtai_benchmark::<5, 65536, BabyBearRingNTT, BabyBearParams>(&mut group),
                    (6, 65536) => ajtai_benchmark::<6, 65536, BabyBearRingNTT, BabyBearParams>(&mut group),
                    (7, 65536) => ajtai_benchmark::<7, 65536, BabyBearRingNTT, BabyBearParams>(&mut group),
                    (8, 65536) => ajtai_benchmark::<8, 65536, BabyBearRingNTT, BabyBearParams>(&mut group),
                    (9, 65536) => ajtai_benchmark::<9, 65536, BabyBearRingNTT, BabyBearParams>(&mut group),
                    (10, 65536) => ajtai_benchmark::<10, 65536, BabyBearRingNTT, BabyBearParams>(&mut group),
                    _ => {} // Skip other combinations
                }
            }
        }
        group.finish();
    }
}

pub fn benchmarks_main(c: &mut Criterion) {
    ajtai_benchmarks(c);
}

criterion_group!(
    name=benches;
    config = Criterion::default()
            .sample_size(10)
            .measurement_time(Duration::from_secs(50))
            .warm_up_time(Duration::from_secs(1));
    targets = benchmarks_main
);
criterion_main!(benches);
