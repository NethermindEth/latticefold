#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

use ark_std::{time::Duration, UniformRand};
use criterion::{criterion_group, criterion_main, Criterion};
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

    let mut group = c.benchmark_group(format!("Ajtai {} C={} W={}", ring_name, C, W));

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

fn ajtai_commit_benchmarks(c: &mut Criterion) {
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
    ajtai_commit_benchmarks(c);
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
