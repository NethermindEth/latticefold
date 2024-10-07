#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

use ark_ff::{BigInt, Field, PrimeField};
use ark_std::{time::Duration, UniformRand};
use criterion::{criterion_group, criterion_main, Criterion};
use latticefold::parameters::{BabyBearParams, FrogParams, GoldilocksParams, StarkPrimeParams};
use latticefold::{
    commitment::AjtaiCommitmentScheme,
    parameters::{DecompositionParams, DilithiumTestParams, DILITHIUM_PRIME},
};
use lattirust_ring::{PolyRing, WithLinfNorm};
use std::fmt::Debug;

use cyclotomic_rings::{
    BabyBearRingNTT, FrogRingNTT, GoldilocksRingNTT, StarkRingNTT, SuitableRing,
};
use lattirust_ring::cyclotomic_ring::models::pow2_debug::Pow2CyclotomicPolyRingNTT;
use rand::{thread_rng, RngCore};

fn draw_bellow_bound<R: SuitableRing, const B: u64, Rng>(rng: &mut Rng) -> R
where
    Rng: rand::Rng + ?Sized,
{
    let bound = <<R as PolyRing>::BaseRing as Field>::BasePrimeField::from(B - 1);
    let coeffs = vec![bound; 72];
    let mut poly = R::CoefficientRepresentation::from(coeffs);
    while !all_elements_bellow_bound::<R>(&poly, B) {
        poly = R::CoefficientRepresentation::rand(rng);
    }
    R::from(poly)
}

fn all_elements_bellow_bound<R: SuitableRing>(
    poly: &R::CoefficientRepresentation,
    bound: u64,
) -> bool {
    let coeffs = poly.coeffs();
    let big_int = coeffs
        .iter()
        .map(|coeff: &<<R as PolyRing>::BaseRing as Field>::BasePrimeField| coeff.into_bigint())
        .max()
        .unwrap();
    let bound = <<R as PolyRing>::BaseRing as Field>::BasePrimeField::from(bound).into_bigint();
    big_int < bound
}

fn generalized_ajtai_benchmark<
    const C: usize,
    const W: usize,
    const B: u64,
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
    let witness: Vec<R> = (0..W)
        .map(|_| draw_bellow_bound::<R, B, dyn RngCore>(&mut rng))
        .collect();

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
    ($c:expr, $cw: expr, $($w:expr),+; $($b:expr),+) => {
        $(
            // StarkPrime
            generalized_ajtai_benchmark::<$cw, $w, $b,  StarkRingNTT, StarkPrimeParams>($c, "StarkPrime");
            // Goldilocks
            generalized_ajtai_benchmark::<$cw, $w, $b,  GoldilocksRingNTT, GoldilocksParams>($c, "Goldilocks");
            // BabyBear
            generalized_ajtai_benchmark::<$cw, $w, $b,  BabyBearRingNTT, BabyBearParams>($c, "BabyBear");
            // Frog
            generalized_ajtai_benchmark::<$cw, $w, $b,  FrogRingNTT, FrogParams>($c, "Frog");
            // Dilithium
            generalized_ajtai_benchmark::<$cw, $w, $b,  Pow2CyclotomicPolyRingNTT<DILITHIUM_PRIME, 256>, DilithiumTestParams>($c, "Dilithium");
        )+
    };
}

macro_rules! run_stark_prime_ajtai_benchmarks {
    ($c:expr, $cw: expr, $(($w:expr, $b:expr)),+ ) => {
        $(
            // StarkPrime
            generalized_ajtai_benchmark::<$cw, $w, $b, StarkRingNTT, StarkPrimeParams>($c, "StarkPrime");
        )+
    };
}
macro_rules! run_goldilocks_ajtai_benchmarks {
    ($c:expr, $cw: expr, $(($w:expr, $b:expr)),+ ) => {
        $(
            // Goldilocks
            generalized_ajtai_benchmark::<$cw, $w, $b, GoldilocksRingNTT, GoldilocksParams>($c, "Goldilocks");
        )+
    };
}
macro_rules! run_babybear_ajtai_benchmarks {
    ($c:expr, $cw: expr, $(($w:expr, $b:expr)),+ ) => {
        $(
            // BabyBear
            generalized_ajtai_benchmark::<$cw, $w, $b, BabyBearRingNTT, BabyBearParams>($c, "BabyBear");
        )+
    };
}
fn ajtai_commit_benchmarks(c: &mut Criterion) {
    run_babybear_ajtai_benchmarks!(
        c,
        1,
        ({ 1 << 15 }, 2), // B = 1.646
        ({ 1 << 16 }, 2)  // B = 1.164
                          // ({ 1 << 17 }, 2), // B = 0.823
                          // ({ 1 << 18 }, 2), // B < 1 bellow
                          // ({ 1 << 19 }, 2),
                          // ({ 1 << 20 }, 2)
    );
    run_babybear_ajtai_benchmarks!(
        c,
        2,
        ({ 1 << 15 }, 43), // B = 42.260
        ({ 1 << 16 }, 30), // B = 29.882
        ({ 1 << 17 }, 22), // B = 21.130
        ({ 1 << 18 }, 15), // B = 14.941
        ({ 1 << 19 }, 11), // B = 10.565
        ({ 1 << 20 }, 8)   // B = 7.471
    );
    run_babybear_ajtai_benchmarks!(
        c,
        3,
        ({ 1 << 15 }, 510), // B = 509.889
        ({ 1 << 16 }, 361), // B = 360.546
        ({ 1 << 17 }, 255), // B = 254.945
        ({ 1 << 18 }, 181), // B = 180.273
        ({ 1 << 19 }, 128), // B = 127.472
        ({ 1 << 20 }, 91)   // B = 90.137
    );
    run_babybear_ajtai_benchmarks!(
        c,
        4,
        ({ 1 << 15 }, 4162), // B = 4161.599
        ({ 1 << 16 }, 2943), // B = 2942.694
        ({ 1 << 17 }, 2081), // B = 2080.799
        ({ 1 << 18 }, 1472), // B = 1471.347
        ({ 1 << 19 }, 1041), // B = 1040.3996
        ({ 1 << 20 }, 736)   // B = 735.674
    );
    run_babybear_ajtai_benchmarks!(
        c,
        5,
        ({ 1 << 15 }, 26459), // B = 26458.082
        ({ 1 << 16 }, 18709), // B = 18708.689
        ({ 1 << 17 }, 13230), // B = 13229.041
        ({ 1 << 18 }, 9355),  // B = 9354.345
        ({ 1 << 19 }, 6615),  // B = 6614.521
        ({ 1 << 20 }, 4678)   // B = 4677.172
    );
    run_babybear_ajtai_benchmarks!(
        c,
        6,
        ({ 1 << 15 }, 140863), // B = 140862.511
        ({ 1 << 16 }, 99605),  // B = 99604.837
        ({ 1 << 17 }, 70432),  // B = 70431.256
        ({ 1 << 18 }, 49803),  // B = 49802.418
        ({ 1 << 19 }, 35216),  // B = 35215.628
        ({ 1 << 20 }, 24902)   // B = 24901.209
    );
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
