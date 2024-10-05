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

    // group.bench_with_input(
    //     "DecomposeCommitNTT",
    //     &(ajtai_data_2, witness_2),
    //     |b, (ajtai_data, witness)| {
    //         b.iter(|| {
    //             let _ = ajtai_data.decompose_and_commit_ntt::<P>(witness);
    //         })
    //     },
    // );
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

macro_rules! run_stark_prime_ajtai_benchmarks {
    ($c:expr, $cw: expr, $($w:expr),+) => {
        $(
            // StarkPrime
            generalized_ajtai_benchmark::<$cw, $w, StarkRingNTT, StarkPrimeParams>($c, "StarkPrime");
        )+
    };
}
macro_rules! run_goldilocks_ajtai_benchmarks {
    ($c:expr, $cw: expr, $($w:expr),+) => {
        $(
            // Goldilocks
            generalized_ajtai_benchmark::<$cw, $w, GoldilocksRingNTT, GoldilocksParams>($c, "Goldilocks");
        )+
    };
}
macro_rules! run_babybear_ajtai_benchmarks {
    ($c:expr, $cw: expr, $($w:expr),+) => {
        $(
            // BabyBear
            generalized_ajtai_benchmark::<$cw, $w, BabyBearRingNTT, BabyBearParams>($c, "BabyBear");
        )+
    };
}
fn ajtai_commit_benchmarks(c: &mut Criterion) {
    run_babybear_ajtai_benchmarks!(c, 1, { 1 << 15 });
    run_babybear_ajtai_benchmarks!(c, 1, { 1 << 16 });
    run_babybear_ajtai_benchmarks!(c, 1, { 1 << 17 });
    run_babybear_ajtai_benchmarks!(c, 1, { 1 << 18 });
    run_babybear_ajtai_benchmarks!(c, 1, { 1 << 19 });
    run_babybear_ajtai_benchmarks!(c, 1, { 1 << 20 });
    run_babybear_ajtai_benchmarks!(c, 2, { 1 << 15 });
    run_babybear_ajtai_benchmarks!(c, 2, { 1 << 16 });
    run_babybear_ajtai_benchmarks!(c, 2, { 1 << 17 });
    run_babybear_ajtai_benchmarks!(c, 2, { 1 << 18 });
    run_babybear_ajtai_benchmarks!(c, 2, { 1 << 19 });
    run_babybear_ajtai_benchmarks!(c, 2, { 1 << 20 });
    run_babybear_ajtai_benchmarks!(c, 3, { 1 << 15 });
    run_babybear_ajtai_benchmarks!(c, 3, { 1 << 16 });
    run_babybear_ajtai_benchmarks!(c, 3, { 1 << 17 });
    run_babybear_ajtai_benchmarks!(c, 3, { 1 << 18 });
    run_babybear_ajtai_benchmarks!(c, 3, { 1 << 19 });
    run_babybear_ajtai_benchmarks!(c, 3, { 1 << 20 });
    run_babybear_ajtai_benchmarks!(c, 4, { 1 << 15 });
    run_babybear_ajtai_benchmarks!(c, 4, { 1 << 16 });
    run_babybear_ajtai_benchmarks!(c, 4, { 1 << 17 });
    run_babybear_ajtai_benchmarks!(c, 4, { 1 << 18 });
    run_babybear_ajtai_benchmarks!(c, 4, { 1 << 19 });
    run_babybear_ajtai_benchmarks!(c, 4, { 1 << 20 });
    run_babybear_ajtai_benchmarks!(c, 5, { 1 << 15 });
    run_babybear_ajtai_benchmarks!(c, 5, { 1 << 16 });
    run_babybear_ajtai_benchmarks!(c, 5, { 1 << 17 });
    run_babybear_ajtai_benchmarks!(c, 5, { 1 << 18 });
    run_babybear_ajtai_benchmarks!(c, 5, { 1 << 19 });
    run_babybear_ajtai_benchmarks!(c, 5, { 1 << 20 });
    run_babybear_ajtai_benchmarks!(c, 6, { 1 << 15 });
    run_babybear_ajtai_benchmarks!(c, 6, { 1 << 16 });
    run_babybear_ajtai_benchmarks!(c, 6, { 1 << 17 });
    run_babybear_ajtai_benchmarks!(c, 6, { 1 << 18 });
    run_babybear_ajtai_benchmarks!(c, 6, { 1 << 19 });
    run_babybear_ajtai_benchmarks!(c, 6, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 1, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 1, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 1, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 1, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 1, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 1, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 2, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 2, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 2, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 2, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 2, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 2, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 3, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 3, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 3, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 3, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 3, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 3, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 4, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 4, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 4, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 4, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 4, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 4, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 5, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 5, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 5, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 5, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 5, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 5, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 6, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 6, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 6, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 6, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 6, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 6, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 7, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 7, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 7, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 7, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 7, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 7, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 8, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 8, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 8, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 8, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 8, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 8, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 9, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 9, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 9, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 9, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 9, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 9, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 10, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 10, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 10, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 10, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 10, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 10, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 11, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 11, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 11, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 11, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 11, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 11, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 12, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 12, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 12, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 12, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 12, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 12, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 13, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 13, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 13, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 13, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 13, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 13, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 14, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 14, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 14, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 14, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 14, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 14, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 15, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 15, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 15, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 15, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 15, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 15, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 16, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 16, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 16, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 16, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 16, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 16, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 17, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 17, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 17, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 17, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 17, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 17, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 18, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 18, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 18, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 18, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 18, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 18, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 19, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 19, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 19, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 19, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 19, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 19, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 20, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 20, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 20, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 20, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 20, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 20, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 21, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 21, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 21, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 21, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 21, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 21, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 22, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 22, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 22, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 22, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 22, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 22, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 23, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 23, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 23, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 23, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 23, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 23, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 24, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 24, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 24, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 24, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 24, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 25, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 25, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 25, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 25, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 25, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 26, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 26, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 26, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 26, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 26, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 26, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 27, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 27, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 27, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 27, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 27, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 27, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 28, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 28, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 28, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 28, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 28, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 28, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 29, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 29, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 29, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 29, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 29, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 29, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 30, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 30, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 30, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 30, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 30, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 30, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 31, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 31, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 31, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 31, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 31, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 31, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 32, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 32, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 32, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 32, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 32, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 32, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 33, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 33, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 33, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 33, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 33, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 33, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 34, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 34, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 34, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 34, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 34, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 34, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 35, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 35, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 35, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 35, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 35, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 35, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 36, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 36, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 36, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 36, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 36, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 36, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 37, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 37, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 37, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 37, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 37, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 37, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 38, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 38, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 38, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 38, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 38, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 38, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 39, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 39, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 39, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 39, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 39, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 39, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 40, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 40, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 40, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 40, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 40, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 40, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 41, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 41, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 41, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 41, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 41, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 41, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 42, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 42, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 42, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 42, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 42, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 42, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 43, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 43, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 43, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 43, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 43, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 43, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 44, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 44, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 44, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 44, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 44, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 44, { 1 << 20 });
    run_goldilocks_ajtai_benchmarks!(c, 45, { 1 << 15 });
    run_goldilocks_ajtai_benchmarks!(c, 45, { 1 << 16 });
    run_goldilocks_ajtai_benchmarks!(c, 45, { 1 << 17 });
    run_goldilocks_ajtai_benchmarks!(c, 45, { 1 << 18 });
    run_goldilocks_ajtai_benchmarks!(c, 45, { 1 << 19 });
    run_goldilocks_ajtai_benchmarks!(c, 45, { 1 << 20 });
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
