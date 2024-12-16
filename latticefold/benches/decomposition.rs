use std::time::Duration;

use criterion::{criterion_group, criterion_main, AxisScale, Criterion, PlotConfiguration};
use latticefold::decomposition_parameters::DecompositionParams;

mod utils;

include!(concat!(env!("OUT_DIR"), "/generated_benchmarks.rs"));

pub fn benchmarks_main(c: &mut Criterion) {
    bench_goldilocks_decomposition(c);
    bench_goldilocks_non_scalar_decomposition(c);
    bench_goldilocks_degree_three_non_scalar_decomposition(c);

    bench_stark_prime_decomposition(c);
    bench_stark_prime_non_scalar_decomposition(c);
    bench_stark_prime_degree_three_non_scalar_decomposition(c);

    bench_frog_decomposition(c);
    bench_frog_non_scalar_decomposition(c);
    bench_frog_degree_three_non_scalar_decomposition(c);

    bench_single_babybear_decomposition(c);
    bench_single_babybear_non_scalar_decomposition(c);
    bench_single_babybear_degree_three_non_scalar_decomposition(c);
}

criterion_group!(
    name=benches;
    config = Criterion::default().sample_size(10).measurement_time(Duration::from_secs(50)).warm_up_time(Duration::from_secs(1));
    targets = benchmarks_main);
criterion_main!(benches);
