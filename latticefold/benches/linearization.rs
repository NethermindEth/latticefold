#![allow(incomplete_features)]
use criterion::{
    criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion, PlotConfiguration,
};
use cyclotomic_rings::{
    challenge_set::LatticefoldChallengeSet,
    rings::{
        BabyBearChallengeSet, BabyBearRingNTT, FrogChallengeSet, FrogRingNTT,
        GoldilocksChallengeSet, GoldilocksRingNTT, StarkChallengeSet, StarkRingNTT, SuitableRing,
    },
};
use std::time::Duration;
use utils::{wit_and_ccs_gen_degree_three_non_scalar, wit_and_ccs_gen_non_scalar};

mod macros;
mod utils;
use crate::utils::wit_and_ccs_gen;
use latticefold::{
    arith::{Witness, CCCS, CCS, LCCCS},
    decomposition_parameters::DecompositionParams,
    nifs::linearization::{
        LFLinearizationProver, LFLinearizationVerifier, LinearizationProof, LinearizationProver,
        LinearizationVerifier,
    },
    transcript::poseidon::PoseidonTranscript,
};

fn prover_linearization_benchmark<
    const C: usize,
    const W: usize,
    P: DecompositionParams,
    R: SuitableRing,
    CS: LatticefoldChallengeSet<R> + Clone,
>(
    c: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>,
    cm_i: &CCCS<C, R>,
    wit: &Witness<R>,
    ccs: &CCS<R>,
) -> (LCCCS<C, R>, LinearizationProof<R>) {
    let mut transcript = PoseidonTranscript::<R, CS>::default();
    let res = LFLinearizationProver::<_, PoseidonTranscript<R, CS>>::prove(
        cm_i,
        wit,
        &mut transcript,
        ccs,
    )
    .expect("Failed to generate linearization proof");

    c.bench_with_input(
        BenchmarkId::new(
            "Linearization Prover",
            format!(
                "Param. Kappa={}, Cols={},  B={}, L={}, B_small={}, K={}",
                C,
                { W / P::L },
                P::B,
                P::L,
                P::B_SMALL,
                P::K
            ),
        ),
        &(cm_i, wit, ccs),
        |b, (cm_i, wit, ccs)| {
            b.iter_batched(
                || transcript.clone(),
                |mut bench_transcript| {
                    let _ = LFLinearizationProver::<_, PoseidonTranscript<R, CS>>::prove(
                        cm_i,
                        wit,
                        &mut bench_transcript,
                        ccs,
                    )
                    .expect("Failed to generate linearization proof");
                },
                criterion::BatchSize::SmallInput,
            );
        },
    );
    res
}

fn verifier_linearization_benchmark<
    const C: usize,
    const W: usize,
    P: DecompositionParams,
    R: SuitableRing,
    CS: LatticefoldChallengeSet<R> + Clone,
>(
    c: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>,
    cm_i: &CCCS<C, R>,
    ccs: &CCS<R>,
    proof: (LCCCS<C, R>, LinearizationProof<R>),
) {
    c.bench_with_input(
        BenchmarkId::new(
            "Linearization Verifier",
            format!(
                "Param. Kappa={}, Cols={},  B={}, L={}, B_small={}, K={}",
                C,
                { W / P::L },
                P::B,
                P::L,
                P::B_SMALL,
                P::K
            ),
        ),
        &(cm_i, proof.1, ccs),
        |b, (cm_i, proof, ccs)| {
            b.iter(|| {
                let mut transcript = PoseidonTranscript::<R, CS>::default();
                let _ = LFLinearizationVerifier::<_, PoseidonTranscript<R, CS>>::verify(
                    cm_i,
                    proof,
                    &mut transcript,
                    ccs,
                )
                .expect("Failed to verify linearization proof");
            })
        },
    );
}

fn linearization_benchmarks_scalar<
    const X_LEN: usize,
    const C: usize,
    const WIT_LEN: usize,
    const W: usize,
    CS: LatticefoldChallengeSet<R> + Clone,
    R: SuitableRing,
    P: DecompositionParams,
>(
    group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>,
) {
    let r1cs_rows = X_LEN + WIT_LEN + 1;
    let (cm_i, wit, ccs, _) = wit_and_ccs_gen::<X_LEN, C, WIT_LEN, W, P, R>(r1cs_rows);

    let proof = prover_linearization_benchmark::<C, W, P, R, CS>(group, &cm_i, &wit, &ccs);

    verifier_linearization_benchmark::<C, W, P, R, CS>(group, &cm_i, &ccs, proof);
}

fn linearization_benchmarks_non_scalar<
    const X_LEN: usize,
    const C: usize,
    const WIT_LEN: usize,
    const W: usize,
    CS: LatticefoldChallengeSet<R> + Clone,
    R: SuitableRing,
    P: DecompositionParams,
>(
    group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>,
) {
    let r1cs_rows = X_LEN + WIT_LEN + 1;
    let (cm_i, wit, ccs, _) = wit_and_ccs_gen_non_scalar::<X_LEN, C, WIT_LEN, W, P, R>(r1cs_rows);

    let proof = prover_linearization_benchmark::<C, W, P, R, CS>(group, &cm_i, &wit, &ccs);

    verifier_linearization_benchmark::<C, W, P, R, CS>(group, &cm_i, &ccs, proof);
}

fn linearization_benchmarks_degree_three_non_scalar<
    const X_LEN: usize,
    const C: usize,
    const WIT_LEN: usize,
    const W: usize,
    CS: LatticefoldChallengeSet<R> + Clone,
    R: SuitableRing,
    P: DecompositionParams,
>(
    group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>,
) {
    let r1cs_rows = X_LEN + WIT_LEN + 1;
    let (cm_i, wit, ccs, _) =
        wit_and_ccs_gen_degree_three_non_scalar::<X_LEN, C, WIT_LEN, W, P, R>(r1cs_rows);

    let proof = prover_linearization_benchmark::<C, W, P, R, CS>(group, &cm_i, &wit, &ccs);

    verifier_linearization_benchmark::<C, W, P, R, CS>(group, &cm_i, &ccs, proof);
}

// Goldilocks macros
#[allow(unused_macros)]
macro_rules! run_single_goldilocks_benchmark {
    ($crit:expr, $io:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            linearization_benchmarks_scalar::<$io, $cw, $w, {$w * $l}, GoldilocksChallengeSet, GoldilocksRingNTT, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit);
        }
    };
}

#[allow(unused_macros)]
macro_rules! run_single_goldilocks_non_scalar_benchmark {
    ($crit:expr, $io:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            linearization_benchmarks_non_scalar::<$io, $cw, $w, {$w * $l}, GoldilocksChallengeSet, GoldilocksRingNTT, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit);
        }
    };
}

#[allow(unused_macros)]
macro_rules! run_single_goldilocks_degree_three_non_scalar_benchmark {
    ($crit:expr, $io:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            linearization_benchmarks_degree_three_non_scalar::<$io, $cw, $w, {$w * $l}, GoldilocksChallengeSet, GoldilocksRingNTT, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit);
        }
    };
}

// Babybear macros
#[allow(unused_macros)]
macro_rules! run_single_babybear_benchmark {
    ($crit:expr, $io:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            linearization_benchmarks_scalar::<$io, $cw, $w, {$w * $l}, BabyBearChallengeSet, BabyBearRingNTT, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit);
        }
    };
}

#[allow(unused_macros)]
macro_rules! run_single_babybear_non_scalar_benchmark {
    ($crit:expr, $io:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            linearization_benchmarks_non_scalar::<$io, $cw, $w, {$w * $l}, BabyBearChallengeSet, BabyBearRingNTT, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit);
        }
    };
}

#[allow(unused_macros)]
macro_rules! run_single_babybear_degree_three_non_scalar_benchmark {
    ($crit:expr, $io:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
           linearization_benchmarks_degree_three_non_scalar::<$io, $cw, $w, {$w * $l}, BabyBearChallengeSet, BabyBearRingNTT, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit);
        }
    };
}

// StarkPrime macros
#[allow(unused_macros)]
macro_rules! run_single_starkprime_benchmark {
    ($crit:expr, $io:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            linearization_benchmarks_scalar::<$io, $cw, $w, {$w * $l}, StarkChallengeSet, StarkRingNTT, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit);
        }
    };
}

#[allow(unused_macros)]
macro_rules! run_single_starkprime_non_scalar_benchmark {
    ($crit:expr, $io:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            linearization_benchmarks_non_scalar::<$io, $cw, $w, {$w * $l}, StarkChallengeSet, StarkRingNTT, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit);
        }
    };
}

#[allow(unused_macros)]
macro_rules! run_single_starkprime_degree_three_non_scalar_benchmark {
    ($crit:expr, $io:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
           linearization_benchmarks_degree_three_non_scalar::<$io, $cw, $w, {$w * $l}, StarkChallengeSet, StarkRingNTT, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit);
        }
    };
}

// Frog macros
#[allow(unused_macros)]
macro_rules! run_single_frog_benchmark {
    ($crit:expr, $io:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            linearization_benchmarks_scalar::<$io, $cw, $w, {$w * $l}, FrogChallengeSet, FrogRingNTT, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit);
        }
    };
}

#[allow(unused_macros)]
macro_rules! run_single_frog_non_scalar_benchmark {
    ($crit:expr, $io:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            linearization_benchmarks_non_scalar::<$io, $cw, $w, {$w * $l}, FrogChallengeSet, FrogRingNTT, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit);
        }
    };
}

#[allow(unused_macros)]
macro_rules! run_single_frog_degree_three_non_scalar_benchmark {
    ($crit:expr, $io:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
           linearization_benchmarks_degree_three_non_scalar::<$io, $cw, $w, {$w * $l}, FrogChallengeSet, FrogRingNTT, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit);
        }
    };
}

fn benchmarks_main(c: &mut Criterion) {
    // Goldilocks
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Linearization Goldilocks");
        group.plot_config(plot_config.clone());
        #[allow(clippy::identity_op)]
        {
            run_goldilocks_benchmarks!(group);
        }
    }

    // Godlilocks non scalar
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Linearization Goldilocks non scalar");
        group.plot_config(plot_config.clone());

        run_goldilocks_non_scalar_benchmarks!(group);
    }

    // Godlilocks degree three non scalar
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Linearization Goldilocks degree three non scalar");
        group.plot_config(plot_config.clone());

        run_goldilocks_degree_three_non_scalar_benchmarks!(group);
    }

    // BabyBear
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Linearization BabyBear");
        group.plot_config(plot_config.clone());
        #[allow(clippy::identity_op)]
        {
            run_babybear_benchmarks!(group);
        }
    }

    // BabyBear non scalar
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Linearization BabyBear non scalar");
        group.plot_config(plot_config.clone());

        run_babybear_non_scalar_benchmarks!(group);
    }

    // BabyBear non scalar degree three
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Linearization BabyBear non scalar");
        group.plot_config(plot_config.clone());

        run_babybear_degree_three_non_scalar_benchmarks!(group);
    }

    // StarkPrime
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Linearization StarkPrime");
        group.plot_config(plot_config.clone());

        #[allow(clippy::identity_op)]
        {
            run_starkprime_benchmarks!(group);
        }
    }

    // StarkPrime non scalar
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Linearization StarkPrime non scalar");
        group.plot_config(plot_config.clone());

        run_starkprime_non_scalar_benchmarks!(group);
    }

    // StarkPrime non scalar degree three
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Linearization StarkPrime degree three non scalar");
        group.plot_config(plot_config.clone());

        run_starkprime_degree_three_non_scalar_benchmarks!(group);
    }
    // Frog
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Linearization Frog");
        group.plot_config(plot_config.clone());
        #[allow(clippy::identity_op)]
        {
            run_frog_benchmarks!(group);
        }
    }

    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Linearization Frog non scalar");
        group.plot_config(plot_config.clone());

        run_frog_non_scalar_benchmarks!(group);
    }

    // Frog degree three non scalar
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Linearization Frog  degree three non scalar");
        group.plot_config(plot_config.clone());

        run_frog_degree_three_non_scalar_benchmarks!(group);
    }
}

criterion_group!(
    name=benches;
    config = Criterion::default().sample_size(10).measurement_time(Duration::from_secs(50)).warm_up_time(Duration::from_secs(1));
    targets = benchmarks_main);
criterion_main!(benches);
