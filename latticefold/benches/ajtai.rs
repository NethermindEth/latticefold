use ark_std::time::Duration;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use cyclotomic_rings::SuitableRing;
use lattirust_ring::Pow2CyclotomicPolyRingNTT;
use rand::thread_rng;

use latticefold::{
    commitment::AjtaiCommitmentScheme,
    parameters::{
        DecompositionParamData, DecompositionParams, DilithiumTestParams, Pow2_57TestParams,
        Pow2_59TestParams, POW2_57_PRIME, POW2_59_PRIME, SOME_FERMAT_PRIME,
    },
};

fn ajtai_benchmark<const C: usize, const W: usize, P: DecompositionParams, R: SuitableRing>(
    c: &mut Criterion,
    p: P,
    prime_name: &str,
) {
    let ajtai_data: AjtaiCommitmentScheme<C, W, R> = AjtaiCommitmentScheme::rand(&mut thread_rng());

    let witness: Vec<R> = (0..W).map(|_| R::rand(&mut thread_rng())).collect();

    let bench_name = format!(
        "Ajtai {}x{} for {} with degree {} and bound {}",
        C,
        W,
        prime_name,
        R::dimension(),
        P::B,
    );
    c.bench_with_input(
        BenchmarkId::new(bench_name, DecompositionParamData::from(p)),
        &(ajtai_data, witness),
        |b, (ajtai_data, witness)| b.iter(|| ajtai_data.commit_ntt(witness)),
    );
}

fn ajtai_benchmarks(c: &mut Criterion) {
    ajtai_benchmark::<9, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<SOME_FERMAT_PRIME, 256>>(
        c,
        DilithiumTestParams,
        "Some fermat(2^16 + 1) prime",
    );
    ajtai_benchmark::<9, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<POW2_59_PRIME, 256>>(
        c,
        Pow2_59TestParams,
        "p = 27*(1<<59) + 1",
    );
    ajtai_benchmark::<9, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<POW2_57_PRIME, 256>>(
        c,
        Pow2_57TestParams,
        "p = 71*(1<<57) + 1",
    );

    // TODO: more benchmarks with different params.
    // ajtai_benchmark::<
    //     6,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    // >(c, BabyBearTestParams, "BabyBear");
}

pub fn benchmarks_main(c: &mut Criterion) {
    ajtai_benchmarks(c);
}

criterion_group!(
    name=benches;
    config = Criterion::default().sample_size(10).measurement_time(Duration::from_secs(50)).warm_up_time(Duration::from_secs(1));
    targets = benchmarks_main);
criterion_main!(benches);
