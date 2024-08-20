use ark_ff::UniformRand;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use lattirust_arithmetic::ring::{Pow2CyclotomicPolyRing, Pow2CyclotomicPolyRingNTT, Zq};
use rand::thread_rng;
use std::time::Duration;

use latticefold::commitment::{
    AjtaiCommitmentScheme, AjtaiParamData, AjtaiParams, DilithiumTestParams, DILITHIUM_PRIME,
};

fn ajtai_benchmark<const Q: u64, const N: usize, P: AjtaiParams>(c: &mut Criterion, p: P) {
    let ajtai_data: AjtaiCommitmentScheme<
        Pow2CyclotomicPolyRing<Zq<Q>, N>,
        Pow2CyclotomicPolyRingNTT<Q, N>,
        P,
    > = AjtaiCommitmentScheme::rand(&mut thread_rng());

    let input: Vec<Pow2CyclotomicPolyRingNTT<Q, N>> = (0..P::WITNESS_SIZE)
        .map(|_| Pow2CyclotomicPolyRingNTT::rand(&mut thread_rng()))
        .collect();

    c.bench_with_input(
        BenchmarkId::new("Ajtai", AjtaiParamData::from(p)),
        &(ajtai_data, input),
        |b, (ajtai_data, input)| b.iter(|| ajtai_data.commit_ntt(input)),
    );
}

fn ajtai_benchmarks(c: &mut Criterion) {
    ajtai_benchmark::<DILITHIUM_PRIME, 256, _>(c, DilithiumTestParams);

    // TODO: more benchmarks with different params.
}

pub fn benchmarks_main(c: &mut Criterion) {
    ajtai_benchmarks(c);
}

criterion_group!(
    name=benches;
    config = Criterion::default().sample_size(10).measurement_time(Duration::from_secs(50)).warm_up_time(Duration::from_secs(1));
    targets = benchmarks_main);
criterion_main!(benches);
