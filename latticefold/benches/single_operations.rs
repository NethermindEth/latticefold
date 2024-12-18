#![allow(incomplete_features)]
use ark_std::{test_rng, time::Duration, UniformRand};
use criterion::{
    criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion, PlotConfiguration,
};
use cyclotomic_rings::rings::{GoldilocksRingNTT, SuitableRing};
use lattirust_poly::mle::DenseMultilinearExtension;
use lattirust_ring::{
    cyclotomic_ring::{CRT, ICRT},
    PolyRing,
};
use rand::thread_rng;
use std::fmt::Debug;

fn single_op_benchmark<R: Clone + UniformRand + Debug + SuitableRing>(
    group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>,
) {
    // Get random parameters
    let mut rng = test_rng();
    let a = R::rand(&mut rng);
    let b = R::rand(&mut rng);
    let a_field = <R::CoefficientRepresentation as PolyRing>::BaseRing::rand(&mut rng);
    let b_field = <R::CoefficientRepresentation as PolyRing>::BaseRing::rand(&mut rng);

    // group.bench_with_input("Addition Field", &(a_field, b_field), |bench, (a, b)| {
    //     bench.iter(|| {
    //         let _ = *a + *b;
    //     })
    // });

    // group.bench_with_input(
    //     "Substraction Field",
    //     &(a_field, b_field),
    //     |bench, (a, b)| {
    //         bench.iter(|| {
    //             let _ = *a * *b;
    //         })
    //     },
    // );

    // group.bench_with_input(
    //     "Multiplication Field",
    //     &(a_field, b_field),
    //     |bench, (a, b)| {
    //         bench.iter(|| {
    //             let _ = *a * *b;
    //         })
    //     },
    // );

    // group.bench_with_input("Addition NTT", &(a, b), |bench, (a, b)| {
    //     bench.iter(|| {
    //         let _ = *a + *b;
    //     })
    // });

    // group.bench_with_input("Substraction NTT", &(a, b), |bench, (a, b)| {
    //     bench.iter(|| {
    //         let _ = *a - *b;
    //     })
    // });

    // group.bench_with_input("Multiplication NTT", &(a, b), |bench, (a, b)| {
    //     bench.iter(|| {
    //         let _ = *a * *b;
    //     })
    // });

    // group.bench_with_input("CRT", &a.icrt(), |bench, a| {
    //     bench.iter(|| {
    //         let _ = a.crt();
    //     })
    // });

    // group.bench_with_input("ICRT", &a, |bench, a| {
    //     bench.iter(|| {
    //         let _ = a.icrt();
    //     })
    // });

    for nv in 1..20 {
        let vec_ntt_form = (0..(1 << nv))
            .map(|_| R::rand(&mut rng))
            .collect::<Vec<R>>();
        // let vec_coeff_form = (0..(1 << nv))
        //     .map(|_| R::CoefficientRepresentation::rand(&mut rng))
        //     .collect::<Vec<_>>();

        // group.bench_with_input(
        //     BenchmarkId::new(
        //         "Elementwise CRT",
        //         format!("{} W = {}", "Goldilocks", 1 << nv),
        //     ),
        //     &vec_coeff_form,
        //     |b, input| {
        //         b.iter_batched(
        //             || input.clone(),
        //             |input| CRT::elementwise_crt(input),
        //             criterion::BatchSize::SmallInput,
        //         )
        //     },
        // );
        // group.bench_with_input(
        //     BenchmarkId::new(
        //         "Elementwise ICRT",
        //         format!("{} W = {}", "Goldilocks", 1 << nv),
        //     ),
        //     &vec_ntt_form,
        //     |b, input| {
        //         b.iter_batched(
        //             || input.clone(),
        //             |input| ICRT::elementwise_icrt(input),
        //             criterion::BatchSize::SmallInput,
        //         )
        //     },
        // );

        let mle = DenseMultilinearExtension::from_slice(nv, &vec_ntt_form);
        let points = (0..nv)
            .map(|i| R::rand(&mut thread_rng()))
            .collect::<Vec<_>>();
        group.bench_with_input(
            BenchmarkId::new(
                "MLE Evaluation",
                format!("{} W = {}", "Goldilocks", 1 << nv),
            ),
            &vec_ntt_form,
            |b, input| {
                b.iter_batched(
                    || mle.clone(),
                    |mle| mle.evaluate(&points),
                    criterion::BatchSize::SmallInput,
                )
            },
        );
    }
}
fn single_operation_benchmarks(c: &mut Criterion) {
    let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    let mut group = c.benchmark_group("Single Operations Goldilocks");
    group.plot_config(plot_config.clone());

    single_op_benchmark::<GoldilocksRingNTT>(&mut group);
}

pub fn benchmarks_main(c: &mut Criterion) {
    single_operation_benchmarks(c);
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
