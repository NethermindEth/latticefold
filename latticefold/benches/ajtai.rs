#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

use ark_std::{time::Duration, UniformRand};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use latticefold::parameters::{FrogParams, GoldilocksParams, StarkPrimeParams};
use latticefold::{
    commitment::AjtaiCommitmentScheme,
    parameters::{
        DecompositionParamData, DecompositionParams, DilithiumTestParams, DILITHIUM_PRIME,
    },
};

use cyclotomic_rings::{FrogRingNTT, GoldilocksRingNTT, StarkRingNTT};
use lattirust_ring::cyclotomic_ring::models::pow2_debug::{
    Pow2CyclotomicPolyRing, Pow2CyclotomicPolyRingNTT,
};
use rand::thread_rng;

fn ajtai_benchmark<
    const Q: u64,
    const N: usize,
    const C: usize,
    const W: usize,
    P: DecompositionParams,
>(
    c: &mut Criterion,
    p: P,
) where
    Pow2CyclotomicPolyRingNTT<Q, N>: From<Pow2CyclotomicPolyRing<Q, N>>,
    Pow2CyclotomicPolyRing<Q, N>: From<Pow2CyclotomicPolyRingNTT<Q, N>>,
{
    let ajtai_data: AjtaiCommitmentScheme<C, W, Pow2CyclotomicPolyRingNTT<Q, N>> =
        AjtaiCommitmentScheme::rand(&mut thread_rng());

    let witness: Vec<Pow2CyclotomicPolyRingNTT<Q, N>> = (0..W)
        .map(|_| Pow2CyclotomicPolyRingNTT::rand(&mut thread_rng()))
        .collect();

    let ajtai_data_2 = ajtai_data.clone();
    let witness_2 = witness.clone();
    let p_2 = p.clone();

    c.bench_with_input(
        BenchmarkId::new(
            format!("Ajtai - CommitNTT - Dilithium C={} W={}", C, W),
            DecompositionParamData::from(p),
        ),
        &(ajtai_data, witness),
        |b, (ajtai_data, witness)| b.iter(|| ajtai_data.commit_ntt(witness)),
    );

    c.bench_with_input(
        BenchmarkId::new(
            format!("Ajtai - DecomposeCommitNTT - Dilithium C={} W={}", C, W),
            DecompositionParamData::from(p_2),
        ),
        &(ajtai_data_2, witness_2),
        |b, (ajtai_data, witness)| b.iter(|| ajtai_data.decompose_and_commit_ntt::<P>(witness)),
    );
}

fn ajtai_starkprime_benchmark<const C: usize, const W: usize, P: DecompositionParams>(
    c: &mut Criterion,
    p: P,
) {
    let ajtai_data: AjtaiCommitmentScheme<C, W, StarkRingNTT> =
        AjtaiCommitmentScheme::rand(&mut thread_rng());

    let witness: Vec<StarkRingNTT> = (0..W)
        .map(|_| StarkRingNTT::rand(&mut thread_rng()))
        .collect();

    let ajtai_data_2 = ajtai_data.clone();
    let witness_2 = witness.clone();
    let p_2 = p.clone();

    c.bench_with_input(
        BenchmarkId::new(
            format!("Ajtai - CommitNTT - Starkprime C={} W={}", C, W),
            DecompositionParamData::from(p),
        ),
        &(ajtai_data, witness),
        |b, (ajtai_data, witness)| b.iter(|| ajtai_data.commit_ntt(witness)),
    );

    c.bench_with_input(
        BenchmarkId::new(
            format!("Ajtai - DecomposeCommitNTT - Starkprime C={} W={}", C, W),
            DecompositionParamData::from(p_2),
        ),
        &(ajtai_data_2, witness_2),
        |b, (ajtai_data, witness)| b.iter(|| ajtai_data.decompose_and_commit_ntt::<P>(witness)),
    );
}

fn ajtai_goldilocks_benchmark<const C: usize, const W: usize, P: DecompositionParams>(
    c: &mut Criterion,
    p: P,
) {
    let ajtai_data: AjtaiCommitmentScheme<C, W, GoldilocksRingNTT> =
        AjtaiCommitmentScheme::rand(&mut thread_rng());

    let witness: Vec<GoldilocksRingNTT> = (0..24)
        .map(|_| GoldilocksRingNTT::rand(&mut thread_rng()))
        .collect();

    let ajtai_data_2 = ajtai_data.clone();
    let witness_2 = witness.clone();
    let p_2 = p.clone();

    c.bench_with_input(
        BenchmarkId::new(
            format!("Ajtai - CommitNTT - Goldilocks C={} W={}", C, W),
            DecompositionParamData::from(p),
        ),
        &(ajtai_data, witness),
        |b, (ajtai_data, witness)| b.iter(|| ajtai_data.commit_ntt(witness)),
    );

    c.bench_with_input(
        BenchmarkId::new(
            format!("Ajtai - DecomposeCommitNTT - Goldilocks C={} W={}", C, W),
            DecompositionParamData::from(p_2),
        ),
        &(ajtai_data_2, witness_2),
        |b, (ajtai_data, witness)| b.iter(|| ajtai_data.decompose_and_commit_ntt::<P>(witness)),
    );
}

fn ajtai_frog_benchmark<const C: usize, const W: usize, P: DecompositionParams>(
    c: &mut Criterion,
    p: P,
) {
    let ajtai_data: AjtaiCommitmentScheme<C, W, FrogRingNTT> =
        AjtaiCommitmentScheme::rand(&mut thread_rng());

    let witness: Vec<FrogRingNTT> = (0..W)
        .map(|_| FrogRingNTT::rand(&mut thread_rng()))
        .collect();

    let ajtai_data_2 = ajtai_data.clone();
    let witness_2 = witness.clone();
    let p_2 = p.clone();

    c.bench_with_input(
        BenchmarkId::new(
            format!("Ajtai - CommitNTT - Frog C={} W={}", C, W),
            DecompositionParamData::from(p),
        ),
        &(ajtai_data, witness),
        |b, (ajtai_data, witness)| b.iter(|| ajtai_data.commit_ntt(witness)),
    );

    c.bench_with_input(
        BenchmarkId::new(
            format!("Ajtai - DecomposeCommitNTT - Frog C={} W={}", C, W),
            DecompositionParamData::from(p_2),
        ),
        &(ajtai_data_2, witness_2),
        |b, (ajtai_data, witness)| b.iter(|| ajtai_data.decompose_and_commit_ntt::<P>(witness)),
    );
}

macro_rules! run_ajtai_benchmarks {
    ($c:expr, $($w:expr),+) => {
        $(
            ajtai_starkprime_benchmark::<9, $w, _>($c, StarkPrimeParams);
            ajtai_goldilocks_benchmark::<9, $w, _>($c, GoldilocksParams);
            ajtai_frog_benchmark::<9, $w, _>($c, FrogParams);
            ajtai_benchmark::<DILITHIUM_PRIME, 256, 9, $w, _>($c, DilithiumTestParams);
        )+
    };
}

fn ajtai_commit_benchmarks(c: &mut Criterion) {
    run_ajtai_benchmarks!(
        c,
        { 1 << 10 },
        { 1 << 11 },
        { 1 << 12 },
        { 1 << 13 },
        { 1 << 14 },
        { 1 << 15 },
        { 1 << 16 },
        { 1 << 17 },
        { 1 << 18 },
        { 1 << 19 },
        { 1 << 20 }
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
