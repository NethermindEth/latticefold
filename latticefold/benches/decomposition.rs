use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use std::time::Duration;

fn decomposition_benchmark() {}

fn decomposition_benchmarks() {}

fn benchmarks_main(c: &mut Criterion) {}

criterion_group!(
    name=benches;
    config = Criterion::default().sample_size(10).measurement_time(Duration::from_secs(50)).warm_up_time(Duration::from_secs(1));
    targets = benchmarks_main);
criterion_main!(benches);
