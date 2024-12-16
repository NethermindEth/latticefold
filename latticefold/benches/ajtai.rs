#![allow(incomplete_features)]

use ark_std::{time::Duration, UniformRand};
use criterion::{
    criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion, PlotConfiguration,
};
use cyclotomic_rings::rings::{
    BabyBearRingNTT, FrogRingNTT, GoldilocksRingNTT, StarkRingNTT, SuitableRing,
};
use latticefold::commitment::AjtaiCommitmentScheme;
use std::fmt::Debug;

mod macros;

fn ajtai_benchmark<
    const C: usize, // rows
    const W: usize, // columns
    R: Clone + UniformRand + Debug + SuitableRing,
>(
    group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>,
) {
    let mut rng = ark_std::test_rng();
    let witness: Vec<R> = (0..W).map(|_| R::rand(&mut rng)).collect();
    let ajtai_data: AjtaiCommitmentScheme<C, W, R> = AjtaiCommitmentScheme::rand(&mut rng);

    group.bench_with_input(
        BenchmarkId::new("AjtaiCommitNTT", format!("Kappa={}, W={}", C, W)),
        &(ajtai_data.clone(), witness),
        |b, (ajtai_data, witness)| {
            b.iter(|| {
                let _ = ajtai_data.commit_ntt(witness);
            })
        },
    );
}

macro_rules! run_single_starkprime_benchmark {
    ($crit:expr, $kappa:expr, $w:expr) => {
        ajtai_benchmark::<$kappa, $w, StarkRingNTT>($crit);
    };
}
macro_rules! run_single_goldilocks_benchmark {
    ($crit:expr, $kappa:expr, $w:expr) => {
        ajtai_benchmark::<$kappa, $w, GoldilocksRingNTT>($crit);
    };
}

macro_rules! run_single_babybear_benchmark {
    ($crit:expr, $kappa:expr, $w:expr) => {
        ajtai_benchmark::<$kappa, $w, BabyBearRingNTT>($crit);
    };
}

macro_rules! run_single_frog_benchmark {
    ($crit:expr, $kappa:expr, $w:expr) => {
        ajtai_benchmark::<$kappa, $w, FrogRingNTT>($crit);
    };
}

fn ajtai_benchmarks(c: &mut Criterion) {
    // // BabyBear
    // {
    //     let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    //     let mut group = c.benchmark_group("Ajtai BabyBear");
    //     group.plot_config(plot_config.clone());

    //     run_babybear_ajai_benchmarks!(group);

    //     group.finish();
    // }

    // Goldilocks
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Ajtai Goldilocks");
        group.plot_config(plot_config.clone());

        run_goldilocks_ajai_benchmarks!(group);

        group.finish();
    }

    // // StarkPrime
    // {
    //     let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    //     let mut group = c.benchmark_group("Ajtai StarkPrime");
    //     group.plot_config(plot_config.clone());

    //     run_starkprime_ajai_benchmarks!(group);

    //     group.finish();
    // }

    // // Frog
    // {
    //     let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    //     let mut group = c.benchmark_group("Ajtai Frog");
    //     group.plot_config(plot_config.clone());

    //     run_frog_ajai_benchmarks!(group);

    //     group.finish();
    // }
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
