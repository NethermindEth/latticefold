use ark_ff::UniformRand;
use ark_std::time::Duration;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use lattirust_ring::Pow2CyclotomicPolyRingNTT;
use rand::thread_rng;

use latticefold::{
    commitment::AjtaiCommitmentScheme,
    parameters::{
        DecompositionParamData, DecompositionParams, DilithiumTestParams, DILITHIUM_PRIME,
    },
};

fn ajtai_benchmark<
    const Q: u64,
    const N: usize,
    const C: usize,
    const W: usize,
    P: DecompositionParams,
>(
    c: &mut Criterion,
    p: P,
) {
    let ajtai_data: AjtaiCommitmentScheme<C, W, Pow2CyclotomicPolyRingNTT<Q, N>> =
        AjtaiCommitmentScheme::rand(&mut thread_rng());

    let witness: Vec<Pow2CyclotomicPolyRingNTT<Q, N>> = (0..W)
        .map(|_| Pow2CyclotomicPolyRingNTT::rand(&mut thread_rng()))
        .collect();

    c.bench_with_input(
        BenchmarkId::new("Ajtai", DecompositionParamData::from(p)),
        &(ajtai_data, witness),
        |b, (ajtai_data, witness)| b.iter(|| ajtai_data.commit_ntt(witness)),
    );
}

fn ajtai_benchmarks(c: &mut Criterion) {
    ajtai_benchmark::<DILITHIUM_PRIME, 256, 9, { 1 << 15 }, _>(c, DilithiumTestParams);

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
