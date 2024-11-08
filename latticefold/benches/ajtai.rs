#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

use ark_std::{time::Duration, UniformRand};
use criterion::{
    criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion, PlotConfiguration,
};
use cyclotomic_rings::rings::{
    BabyBearRingNTT, FrogRingNTT, GoldilocksRingNTT, StarkRingNTT, SuitableRing,
};
use latticefold::commitment::AjtaiCommitmentScheme;
use lattirust_ring::cyclotomic_ring::models::pow2_debug::Pow2CyclotomicPolyRingNTT;
use rand::thread_rng;
use std::fmt::Debug;

fn ajtai_benchmark<
    const C: usize, // rows
    const W: usize, // columns
    R: Clone + UniformRand + Debug + SuitableRing + for<'a> std::ops::AddAssign<&'a R>,
>(
    group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>,
) {
    let mut rng = thread_rng();

    let witness: Vec<R> = (0..W).map(|_| R::rand(&mut rng)).collect();
    let witness_2 = witness.clone();
    let witness_3 = witness.clone();
    let ajtai_data: AjtaiCommitmentScheme<C, W, R> = AjtaiCommitmentScheme::rand(&mut rng);

    group.bench_with_input(
        BenchmarkId::new("CommitNTT", format!("C={}, W={}", C, W)),
        &(ajtai_data.clone(), witness),
        |b, (ajtai_data, witness)| {
            b.iter(|| {
                let _ = ajtai_data.commit_ntt(witness);
            })
        },
    );

    // NTT -> INTT (coefficients)
    group.bench_with_input(
        BenchmarkId::new("NTT->INTT", format!("C={}, W={}", C, W)),
        &(witness_2),
        |b, witness| {
            b.iter(|| {
                let _ = witness
                    .iter()
                    .map(|&x| x.into())
                    .collect::<Vec<R::CoefficientRepresentation>>();
            })
        },
    );

    // INTT -> NTT
    let coeff = witness_3
        .iter()
        .map(|&x| x.into())
        .collect::<Vec<R::CoefficientRepresentation>>();
    group.bench_with_input(
        BenchmarkId::new("INTT->NTT", format!("C={}, W={}", C, W)),
        &(coeff),
        |b, coeff| {
            b.iter(|| {
                let _: Vec<R> = coeff.iter().map(|&x| R::from(x)).collect();
            })
        },
    );
}

macro_rules! run_single_starkprime_benchmark {
    ($crit:expr, $cw:expr, $w:expr) => {
        ajtai_benchmark::<$cw, $w, StarkRingNTT>($crit);
    };
}

macro_rules! run_single_goldilocks_benchmark {
    ($crit:expr, $cw:expr, $w:expr) => {
        ajtai_benchmark::<$cw, $w, GoldilocksRingNTT>($crit);
    };
}

macro_rules! run_single_babybear_benchmark {
    ($crit:expr, $cw:expr, $w:expr) => {
        ajtai_benchmark::<$cw, $w, BabyBearRingNTT>($crit);
    };
}

macro_rules! run_single_frog_benchmark {
    ($crit:expr, $cw:expr, $w:expr) => {
        ajtai_benchmark::<$cw, $w, FrogRingNTT>($crit);
    };
}

fn ajtai_benchmarks(c: &mut Criterion) {
    // Parameters are (C, W) where C is the number of rows and W is the number of columns.
    // BabyBear
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Ajtai BabyBear");
        group.plot_config(plot_config.clone());

        // Existing benchmarks
        run_single_babybear_benchmark!(&mut group, 1, 32768);
        run_single_babybear_benchmark!(&mut group, 2, 32768);
        run_single_babybear_benchmark!(&mut group, 3, 32768);
        run_single_babybear_benchmark!(&mut group, 4, 32768);
        run_single_babybear_benchmark!(&mut group, 5, 32768);
        run_single_babybear_benchmark!(&mut group, 6, 32768);
        run_single_babybear_benchmark!(&mut group, 1, 65536);
        run_single_babybear_benchmark!(&mut group, 2, 65536);
        run_single_babybear_benchmark!(&mut group, 3, 65536);
        run_single_babybear_benchmark!(&mut group, 4, 65536);
        run_single_babybear_benchmark!(&mut group, 5, 65536);
        run_single_babybear_benchmark!(&mut group, 6, 65536);
        run_single_babybear_benchmark!(&mut group, 1, 131072);
        run_single_babybear_benchmark!(&mut group, 2, 131072);
        run_single_babybear_benchmark!(&mut group, 3, 131072);
        run_single_babybear_benchmark!(&mut group, 4, 131072);
        run_single_babybear_benchmark!(&mut group, 5, 131072);
        run_single_babybear_benchmark!(&mut group, 6, 131072);
        run_single_babybear_benchmark!(&mut group, 1, 262144);
        run_single_babybear_benchmark!(&mut group, 2, 262144);
        run_single_babybear_benchmark!(&mut group, 3, 262144);
        run_single_babybear_benchmark!(&mut group, 4, 262144);
        run_single_babybear_benchmark!(&mut group, 5, 262144);
        run_single_babybear_benchmark!(&mut group, 6, 262144);
        run_single_babybear_benchmark!(&mut group, 1, 524288);
        run_single_babybear_benchmark!(&mut group, 2, 524288);
        run_single_babybear_benchmark!(&mut group, 3, 524288);
        run_single_babybear_benchmark!(&mut group, 4, 524288);
        run_single_babybear_benchmark!(&mut group, 5, 524288);
        run_single_babybear_benchmark!(&mut group, 6, 524288);
        run_single_babybear_benchmark!(&mut group, 1, 1048576);
        run_single_babybear_benchmark!(&mut group, 2, 1048576);
        run_single_babybear_benchmark!(&mut group, 3, 1048576);
        run_single_babybear_benchmark!(&mut group, 4, 1048576);
        run_single_babybear_benchmark!(&mut group, 5, 1048576);
        run_single_babybear_benchmark!(&mut group, 6, 1048576);

        // New benchmarks
        run_single_babybear_benchmark!(&mut group, 12, 32768);
        run_single_babybear_benchmark!(&mut group, 13, 65536);
        run_single_babybear_benchmark!(&mut group, 13, 131072);
        run_single_babybear_benchmark!(&mut group, 14, 262144);
        run_single_babybear_benchmark!(&mut group, 14, 524288);
        run_single_babybear_benchmark!(&mut group, 15, 1048576);

        group.finish();
    }

    // Goldilocks
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Ajtai Goldilocks");
        group.plot_config(plot_config.clone());

        // Existing benchmarks
        run_single_goldilocks_benchmark!(&mut group, 1, 32768);
        run_single_goldilocks_benchmark!(&mut group, 2, 32768);
        run_single_goldilocks_benchmark!(&mut group, 3, 32768);
        run_single_goldilocks_benchmark!(&mut group, 4, 32768);
        run_single_goldilocks_benchmark!(&mut group, 5, 32768);
        run_single_goldilocks_benchmark!(&mut group, 6, 32768);
        run_single_goldilocks_benchmark!(&mut group, 1, 65536);
        run_single_goldilocks_benchmark!(&mut group, 2, 65536);
        run_single_goldilocks_benchmark!(&mut group, 3, 65536);
        run_single_goldilocks_benchmark!(&mut group, 4, 65536);
        run_single_goldilocks_benchmark!(&mut group, 5, 65536);
        run_single_goldilocks_benchmark!(&mut group, 6, 65536);
        run_single_goldilocks_benchmark!(&mut group, 1, 131072);
        run_single_goldilocks_benchmark!(&mut group, 2, 131072);
        run_single_goldilocks_benchmark!(&mut group, 3, 131072);
        run_single_goldilocks_benchmark!(&mut group, 4, 131072);
        run_single_goldilocks_benchmark!(&mut group, 5, 131072);
        run_single_goldilocks_benchmark!(&mut group, 6, 131072);
        run_single_goldilocks_benchmark!(&mut group, 1, 262144);
        run_single_goldilocks_benchmark!(&mut group, 2, 262144);
        run_single_goldilocks_benchmark!(&mut group, 3, 262144);
        run_single_goldilocks_benchmark!(&mut group, 4, 262144);
        run_single_goldilocks_benchmark!(&mut group, 5, 262144);
        run_single_goldilocks_benchmark!(&mut group, 6, 262144);
        run_single_goldilocks_benchmark!(&mut group, 1, 524288);
        run_single_goldilocks_benchmark!(&mut group, 2, 524288);
        run_single_goldilocks_benchmark!(&mut group, 3, 524288);
        run_single_goldilocks_benchmark!(&mut group, 4, 524288);
        run_single_goldilocks_benchmark!(&mut group, 5, 524288);
        run_single_goldilocks_benchmark!(&mut group, 6, 524288);
        run_single_goldilocks_benchmark!(&mut group, 1, 1048576);
        run_single_goldilocks_benchmark!(&mut group, 2, 1048576);
        run_single_goldilocks_benchmark!(&mut group, 3, 1048576);
        run_single_goldilocks_benchmark!(&mut group, 4, 1048576);
        run_single_goldilocks_benchmark!(&mut group, 5, 1048576);
        run_single_goldilocks_benchmark!(&mut group, 6, 1048576);

        // New benchmarks
        run_single_goldilocks_benchmark!(&mut group, 17, 32768);
        run_single_goldilocks_benchmark!(&mut group, 17, 65536);
        run_single_goldilocks_benchmark!(&mut group, 18, 131072);
        run_single_goldilocks_benchmark!(&mut group, 19, 262144);
        run_single_goldilocks_benchmark!(&mut group, 19, 524288);
        run_single_goldilocks_benchmark!(&mut group, 20, 1048576);

        group.finish();
    }

    // StarkPrime
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Ajtai StarkPrime");
        group.plot_config(plot_config.clone());

        // Existing benchmarks
        run_single_starkprime_benchmark!(&mut group, 1, 32768);
        run_single_starkprime_benchmark!(&mut group, 2, 32768);
        run_single_starkprime_benchmark!(&mut group, 3, 32768);
        run_single_starkprime_benchmark!(&mut group, 4, 32768);
        run_single_starkprime_benchmark!(&mut group, 5, 32768);
        run_single_starkprime_benchmark!(&mut group, 6, 32768);
        run_single_starkprime_benchmark!(&mut group, 1, 65536);
        run_single_starkprime_benchmark!(&mut group, 2, 65536);
        run_single_starkprime_benchmark!(&mut group, 3, 65536);
        run_single_starkprime_benchmark!(&mut group, 4, 65536);
        run_single_starkprime_benchmark!(&mut group, 5, 65536);
        run_single_starkprime_benchmark!(&mut group, 6, 65536);
        run_single_starkprime_benchmark!(&mut group, 1, 131072);
        run_single_starkprime_benchmark!(&mut group, 2, 131072);
        run_single_starkprime_benchmark!(&mut group, 3, 131072);
        run_single_starkprime_benchmark!(&mut group, 4, 131072);
        run_single_starkprime_benchmark!(&mut group, 5, 131072);
        run_single_starkprime_benchmark!(&mut group, 6, 131072);
        run_single_starkprime_benchmark!(&mut group, 1, 262144);
        run_single_starkprime_benchmark!(&mut group, 2, 262144);
        run_single_starkprime_benchmark!(&mut group, 3, 262144);
        run_single_starkprime_benchmark!(&mut group, 4, 262144);
        run_single_starkprime_benchmark!(&mut group, 5, 262144);
        run_single_starkprime_benchmark!(&mut group, 6, 262144);
        run_single_starkprime_benchmark!(&mut group, 1, 524288);
        run_single_starkprime_benchmark!(&mut group, 2, 524288);
        run_single_starkprime_benchmark!(&mut group, 3, 524288);
        run_single_starkprime_benchmark!(&mut group, 4, 524288);
        run_single_starkprime_benchmark!(&mut group, 5, 524288);
        run_single_starkprime_benchmark!(&mut group, 6, 524288);
        run_single_starkprime_benchmark!(&mut group, 1, 1048576);
        run_single_starkprime_benchmark!(&mut group, 2, 1048576);
        run_single_starkprime_benchmark!(&mut group, 3, 1048576);
        run_single_starkprime_benchmark!(&mut group, 4, 1048576);
        run_single_starkprime_benchmark!(&mut group, 5, 1048576);
        run_single_starkprime_benchmark!(&mut group, 6, 1048576);

        // New benchmarks
        run_single_starkprime_benchmark!(&mut group, 7, 131072);
        run_single_starkprime_benchmark!(&mut group, 7, 262144);
        run_single_starkprime_benchmark!(&mut group, 7, 524288);
        run_single_starkprime_benchmark!(&mut group, 7, 1048576);

        group.finish();
    }

    // Frog
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Ajtai Frog");
        group.plot_config(plot_config.clone());

        // Existing benchmarks
        run_single_frog_benchmark!(&mut group, 1, 32768);
        run_single_frog_benchmark!(&mut group, 2, 32768);
        run_single_frog_benchmark!(&mut group, 3, 32768);
        run_single_frog_benchmark!(&mut group, 4, 32768);
        run_single_frog_benchmark!(&mut group, 5, 32768);
        run_single_frog_benchmark!(&mut group, 6, 32768);
        run_single_frog_benchmark!(&mut group, 1, 65536);
        run_single_frog_benchmark!(&mut group, 2, 65536);
        run_single_frog_benchmark!(&mut group, 3, 65536);
        run_single_frog_benchmark!(&mut group, 4, 65536);
        run_single_frog_benchmark!(&mut group, 5, 65536);
        run_single_frog_benchmark!(&mut group, 6, 65536);
        run_single_frog_benchmark!(&mut group, 1, 131072);
        run_single_frog_benchmark!(&mut group, 2, 131072);
        run_single_frog_benchmark!(&mut group, 3, 131072);
        run_single_frog_benchmark!(&mut group, 4, 131072);
        run_single_frog_benchmark!(&mut group, 5, 131072);
        run_single_frog_benchmark!(&mut group, 6, 131072);
        run_single_frog_benchmark!(&mut group, 1, 262144);
        run_single_frog_benchmark!(&mut group, 2, 262144);
        run_single_frog_benchmark!(&mut group, 3, 262144);
        run_single_frog_benchmark!(&mut group, 4, 262144);
        run_single_frog_benchmark!(&mut group, 5, 262144);
        run_single_frog_benchmark!(&mut group, 6, 262144);
        run_single_frog_benchmark!(&mut group, 1, 524288);
        run_single_frog_benchmark!(&mut group, 2, 524288);
        run_single_frog_benchmark!(&mut group, 3, 524288);
        run_single_frog_benchmark!(&mut group, 4, 524288);
        run_single_frog_benchmark!(&mut group, 5, 524288);
        run_single_frog_benchmark!(&mut group, 6, 524288);
        run_single_frog_benchmark!(&mut group, 1, 1048576);
        run_single_frog_benchmark!(&mut group, 2, 1048576);
        run_single_frog_benchmark!(&mut group, 3, 1048576);
        run_single_frog_benchmark!(&mut group, 4, 1048576);
        run_single_frog_benchmark!(&mut group, 5, 1048576);
        run_single_frog_benchmark!(&mut group, 6, 1048576);

        // New benchmarks
        run_single_frog_benchmark!(&mut group, 17, 32768);
        run_single_frog_benchmark!(&mut group, 18, 65536);
        run_single_frog_benchmark!(&mut group, 19, 131072);
        run_single_frog_benchmark!(&mut group, 19, 262144);
        run_single_frog_benchmark!(&mut group, 20, 524288);
        run_single_frog_benchmark!(&mut group, 21, 1048576);

        group.finish();
    }
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
