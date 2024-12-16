use ark_std::{time::Duration, UniformRand};
use criterion::{
    criterion_group, criterion_main, AxisScale, BatchSize::SmallInput, BenchmarkId, Criterion,
    PlotConfiguration,
};
use cyclotomic_rings::rings::{BabyBearRingNTT, FrogRingNTT, GoldilocksRingNTT, StarkRingNTT};
use latticefold::commitment::AjtaiCommitmentScheme;
use stark_rings::cyclotomic_ring::{CRT, ICRT};

include!(concat!(env!("OUT_DIR"), "/generated_ajtai_benchmarks.rs"));

fn ajtai_benchmarks(c: &mut Criterion) {
    bench_ajtai_goldilocks(c);
    bench_ajtai_starkprime(c);
    bench_ajtai_babybear(c);
    bench_ajtai_frog(c);
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
