#![allow(incomplete_features)]

use ark_std::{time::Duration, UniformRand};
use criterion::{
    criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion,
    PlotConfiguration,
};
#[allow(unused_imports)]
use cyclotomic_rings::rings::{
    BabyBearRingNTT, FrogRingNTT, GoldilocksRingNTT, StarkRingNTT, SuitableRing,
};
use latticefold::commitment::AjtaiCommitmentScheme;
use std::fmt::Debug;

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
        BenchmarkId::new("CommitNTT", format!("C={}, W={}", C, W)),
        &(ajtai_data.clone(), witness),
        |b, (ajtai_data, witness)| {
            b.iter(|| {
                let _ = ajtai_data.commit_ntt(witness);
            })
        },
    );
}

#[allow(unused_macros)]
macro_rules! run_single_starkprime_benchmark {
    ($crit:expr, $cw:expr, $w:expr) => {
        ajtai_benchmark::<$cw, $w, StarkRingNTT>($crit);
    };
}

#[allow(unused_macros)]
macro_rules! run_single_goldilocks_benchmark {
    ($crit:expr, $cw:expr, $w:expr) => {
        ajtai_benchmark::<$cw, $w, GoldilocksRingNTT>($crit);
    };
}

#[allow(unused_macros)]
macro_rules! run_single_babybear_benchmark {
    ($crit:expr, $cw:expr, $w:expr) => {
        ajtai_benchmark::<$cw, $w, BabyBearRingNTT>($crit);
    };
}

#[allow(unused_macros)]
macro_rules! run_single_frog_benchmark {
    ($crit:expr, $cw:expr, $w:expr) => {
        ajtai_benchmark::<$cw, $w, FrogRingNTT>($crit);
    };
}

fn ajtai_benchmarks(c: &mut Criterion) {
    // Parameters are (C, W) where C is the number of rows and W is the number of columns.
    // Goldilocks
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Ajtai Goldilocks");
        group.plot_config(plot_config.clone());

        run_single_goldilocks_benchmark!(&mut group, 1, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 1, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 1, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 1, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 1, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 1, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 1, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 1, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 1, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 1, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 1, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 1, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 1, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 1, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 1, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 1, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 2, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 2, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 2, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 2, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 2, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 2, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 2, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 2, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 2, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 2, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 2, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 2, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 2, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 2, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 2, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 2, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 3, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 3, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 3, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 3, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 3, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 3, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 3, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 3, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 3, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 3, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 3, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 3, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 3, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 3, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 3, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 3, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 4, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 4, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 4, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 4, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 4, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 4, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 4, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 4, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 4, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 4, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 4, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 4, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 4, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 4, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 4, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 4, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 5, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 5, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 5, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 5, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 5, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 5, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 5, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 5, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 5, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 5, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 5, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 5, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 5, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 5, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 5, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 5, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 6, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 6, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 6, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 6, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 6, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 6, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 6, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 6, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 6, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 6, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 6, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 6, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 6, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 6, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 6, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 6, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 7, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 7, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 7, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 7, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 7, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 7, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 7, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 7, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 7, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 7, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 7, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 7, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 7, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 7, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 7, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 7, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 8, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 8, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 8, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 8, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 8, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 8, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 8, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 8, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 8, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 8, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 8, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 8, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 8, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 8, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 8, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 8, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 9, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 9, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 9, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 9, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 9, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 9, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 9, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 9, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 9, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 9, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 9, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 9, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 9, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 9, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 9, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 9, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 10, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 10, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 10, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 10, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 10, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 10, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 10, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 10, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 10, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 10, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 10, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 10, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 10, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 10, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 10, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 10, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 11, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 11, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 11, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 11, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 11, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 11, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 11, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 11, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 11, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 11, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 11, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 11, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 11, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 11, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 11, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 11, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 12, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 12, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 12, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 12, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 12, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 12, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 12, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 12, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 12, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 12, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 12, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 12, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 12, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 12, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 12, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 12, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 13, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 13, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 13, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 13, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 13, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 13, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 13, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 13, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 13, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 13, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 13, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 13, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 13, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 13, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 13, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 13, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 14, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 14, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 14, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 14, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 14, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 14, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 14, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 14, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 14, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 14, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 14, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 14, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 14, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 14, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 14, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 14, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 15, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 15, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 15, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 15, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 15, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 15, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 15, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 15, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 15, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 15, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 15, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 15, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 15, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 15, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 15, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 15, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 16, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 16, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 16, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 16, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 16, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 16, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 16, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 16, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 16, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 16, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 16, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 16, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 16, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 16, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 16, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 16, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 17, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 17, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 17, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 17, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 17, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 17, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 17, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 17, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 17, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 17, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 17, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 17, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 17, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 17, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 17, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 17, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 18, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 18, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 18, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 18, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 18, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 18, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 18, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 18, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 18, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 18, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 18, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 18, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 18, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 18, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 18, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 18, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 19, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 19, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 19, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 19, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 19, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 19, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 19, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 19, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 19, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 19, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 19, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 19, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 19, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 19, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 19, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 19, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 20, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 20, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 20, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 20, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 20, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 20, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 20, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 20, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 20, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 20, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 20, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 20, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 20, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 20, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 20, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 20, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 21, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 21, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 21, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 21, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 21, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 21, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 21, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 21, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 21, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 21, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 21, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 21, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 21, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 21, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 21, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 21, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 22, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 22, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 22, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 22, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 22, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 22, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 22, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 22, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 22, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 22, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 22, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 22, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 22, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 22, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 22, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 22, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 23, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 23, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 23, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 23, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 23, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 23, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 23, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 23, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 23, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 23, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 23, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 23, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 23, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 23, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 23, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 23, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 24, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 24, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 24, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 24, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 24, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 24, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 24, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 24, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 24, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 24, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 24, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 24, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 24, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 24, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 24, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 24, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 25, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 25, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 25, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 25, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 25, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 25, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 25, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 25, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 25, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 25, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 25, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 25, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 25, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 25, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 25, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 25, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 26, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 26, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 26, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 26, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 26, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 26, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 26, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 26, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 26, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 26, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 26, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 26, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 26, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 26, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 26, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 26, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 27, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 27, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 27, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 27, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 27, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 27, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 27, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 27, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 27, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 27, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 27, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 27, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 27, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 27, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 27, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 27, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 28, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 28, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 28, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 28, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 28, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 28, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 28, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 28, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 28, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 28, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 28, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 28, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 28, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 28, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 28, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 28, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 29, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 29, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 29, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 29, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 29, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 29, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 29, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 29, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 29, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 29, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 29, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 29, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 29, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 29, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 29, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 29, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 30, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 30, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 30, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 30, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 30, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 30, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 30, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 30, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 30, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 30, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 30, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 30, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 30, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 30, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 30, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 30, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 31, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 31, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 31, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 31, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 31, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 31, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 31, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 31, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 31, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 31, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 31, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 31, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 31, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 31, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 31, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 31, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 32, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 32, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 32, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 32, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 32, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 32, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 32, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 32, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 32, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 32, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 32, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 32, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 32, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 32, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 32, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 32, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 33, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 33, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 33, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 33, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 33, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 33, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 33, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 33, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 33, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 33, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 33, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 33, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 33, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 33, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 33, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 33, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 34, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 34, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 34, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 34, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 34, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 34, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 34, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 34, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 34, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 34, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 34, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 34, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 34, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 34, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 34, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 34, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 35, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 35, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 35, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 35, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 35, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 35, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 35, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 35, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 35, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 35, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 35, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 35, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 35, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 35, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 35, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 35, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 36, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 36, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 36, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 36, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 36, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 36, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 36, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 36, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 36, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 36, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 36, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 36, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 36, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 36, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 36, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 36, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 37, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 37, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 37, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 37, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 37, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 37, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 37, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 37, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 37, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 37, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 37, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 37, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 37, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 37, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 37, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 37, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 38, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 38, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 38, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 38, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 38, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 38, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 38, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 38, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 38, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 38, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 38, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 38, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 38, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 38, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 38, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 38, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 39, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 39, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 39, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 39, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 39, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 39, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 39, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 39, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 39, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 39, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 39, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 39, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 39, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 39, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 39, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 39, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 40, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 40, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 40, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 40, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 40, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 40, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 40, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 40, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 40, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 40, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 40, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 40, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 40, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 40, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 40, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 40, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 41, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 41, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 41, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 41, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 41, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 41, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 41, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 41, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 41, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 41, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 41, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 41, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 41, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 41, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 41, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 41, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 42, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 42, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 42, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 42, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 42, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 42, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 42, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 42, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 42, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 42, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 42, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 42, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 42, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 42, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 42, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 42, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 43, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 43, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 43, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 43, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 43, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 43, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 43, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 43, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 43, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 43, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 43, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 43, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 43, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 43, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 43, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 43, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 44, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 44, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 44, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 44, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 44, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 44, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 44, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 44, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 44, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 44, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 44, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 44, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 44, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 44, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 44, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 44, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 45, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 45, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 45, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 45, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 45, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 45, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 45, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 45, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 45, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 45, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 45, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 45, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 45, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 45, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 45, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 45, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 46, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 46, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 46, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 46, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 46, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 46, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 46, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 46, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 46, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 46, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 46, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 46, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 46, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 46, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 46, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 46, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 47, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 47, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 47, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 47, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 47, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 47, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 47, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 47, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 47, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 47, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 47, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 47, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 47, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 47, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 47, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 47, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 48, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 48, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 48, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 48, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 48, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 48, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 48, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 48, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 48, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 48, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 48, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 48, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 48, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 48, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 48, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 48, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 49, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 49, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 49, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 49, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 49, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 49, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 49, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 49, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 49, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 49, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 49, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 49, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 49, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 49, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 49, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 49, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 50, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 50, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 50, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 50, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 50, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 50, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 50, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 50, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 50, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 50, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 50, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 50, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 50, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 50, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 50, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 50, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 51, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 51, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 51, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 51, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 51, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 51, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 51, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 51, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 51, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 51, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 51, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 51, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 51, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 51, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 51, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 51, {1 << 24});

        run_single_goldilocks_benchmark!(&mut group, 52, {1 << 9});
        run_single_goldilocks_benchmark!(&mut group, 52, {1 << 10});
        run_single_goldilocks_benchmark!(&mut group, 52, {1 << 11});
        run_single_goldilocks_benchmark!(&mut group, 52, {1 << 12});
        run_single_goldilocks_benchmark!(&mut group, 52, {1 << 13});
        run_single_goldilocks_benchmark!(&mut group, 52, {1 << 14});
        run_single_goldilocks_benchmark!(&mut group, 52, {1 << 15});
        run_single_goldilocks_benchmark!(&mut group, 52, {1 << 16});
        run_single_goldilocks_benchmark!(&mut group, 52, {1 << 17});
        run_single_goldilocks_benchmark!(&mut group, 52, {1 << 18});
        run_single_goldilocks_benchmark!(&mut group, 52, {1 << 19});
        run_single_goldilocks_benchmark!(&mut group, 52, {1 << 20});
        run_single_goldilocks_benchmark!(&mut group, 52, {1 << 21});
        run_single_goldilocks_benchmark!(&mut group, 52, {1 << 22});
        run_single_goldilocks_benchmark!(&mut group, 52, {1 << 23});
        run_single_goldilocks_benchmark!(&mut group, 52, {1 << 24});

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
