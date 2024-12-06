#![allow(incomplete_features)]
use std::time::Duration;

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
use utils::wit_and_ccs_gen_non_scalar;
mod macros;
mod utils;

use crate::utils::wit_and_ccs_gen;
use latticefold::nifs::linearization::{
    LFLinearizationProver, LFLinearizationVerifier, LinearizationProver, LinearizationVerifier,
};
use latticefold::{
    arith::{Witness, CCCS, CCS},
    commitment::AjtaiCommitmentScheme,
    decomposition_parameters::DecompositionParams,
    nifs::{NIFSProver, NIFSVerifier},
    transcript::poseidon::PoseidonTranscript,
};

fn prover_e2e_benchmark<
    const C: usize,
    const W: usize,
    P: DecompositionParams,
    R: SuitableRing,
    CS: LatticefoldChallengeSet<R> + Clone + 'static,
>(
    c: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>,
    cm_i: &CCCS<C, R>,
    wit: &Witness<R>,
    ccs: &CCS<R>,
    scheme: &AjtaiCommitmentScheme<C, W, R>,
) {
    let mut prover_transcript = PoseidonTranscript::<R, CS>::default();
    let mut verifier_transcript = PoseidonTranscript::<R, CS>::default();

    let (acc_lcccs, linearization_proof) =
        LFLinearizationProver::<_, PoseidonTranscript<R, CS>>::prove(
            cm_i,
            wit,
            &mut prover_transcript,
            ccs,
        )
        .expect("Failed to generate linearization proof");

    let _ = LFLinearizationVerifier::<_, PoseidonTranscript<R, CS>>::verify(
        cm_i,
        &linearization_proof,
        &mut verifier_transcript,
        ccs,
    )
    .expect("Failed to verify linearization");

    c.bench_with_input(
        BenchmarkId::new(
            "E2E Prover",
            format!(
                "Param. Kappa={}, Cols={}, B={}, L={}, B_small={}, K={}",
                C,
                { W / P::L },
                P::B,
                P::L,
                P::B_SMALL,
                P::K
            ),
        ),
        &(
            acc_lcccs.clone(),
            wit.clone(),
            cm_i.clone(),
            ccs.clone(),
            scheme.clone(),
        ),
        |b, (lcccs, wit, cm_i, ccs, scheme)| {
            b.iter_batched(
                || prover_transcript.clone(),
                |mut bench_prover_transcript| {
                    let _ = NIFSProver::<C, W, R, P, PoseidonTranscript<R, CS>>::prove(
                        lcccs,
                        wit,
                        cm_i,
                        wit,
                        &mut bench_prover_transcript,
                        ccs,
                        scheme,
                    )
                    .expect("Failed to generate proof");
                },
                criterion::BatchSize::SmallInput,
            )
        },
    );
}

fn verifier_e2e_benchmark<
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
    scheme: &AjtaiCommitmentScheme<C, W, R>,
) {
    let mut prover_transcript = PoseidonTranscript::<R, CS>::default();
    let mut verifier_transcript = PoseidonTranscript::<R, CS>::default();

    let (prover_lcccs_acc, acc_linearization_proof) = LFLinearizationProver::<
        _,
        PoseidonTranscript<R, CS>,
    >::prove(
        cm_i, wit, &mut prover_transcript, ccs
    )
    .expect("Failed to generate acc linearization proof");

    let verifier_lcccs_acc = LFLinearizationVerifier::<_, PoseidonTranscript<R, CS>>::verify(
        cm_i,
        &acc_linearization_proof,
        &mut verifier_transcript,
        ccs,
    )
    .expect("Failed to verify acc linearization");

    let (_, _, proof) = NIFSProver::<C, W, R, P, PoseidonTranscript<R, CS>>::prove(
        &prover_lcccs_acc,
        wit,
        cm_i,
        wit,
        &mut prover_transcript,
        ccs,
        scheme,
    )
    .expect("Failed to generate proof");

    c.bench_with_input(
        BenchmarkId::new(
            "E2E Verifier",
            format!(
                "Param. Kappa={}, Cols={}, B={}, L={}, B_small={}, K={}",
                C,
                { W / P::L },
                P::B,
                P::L,
                P::B_SMALL,
                P::K
            ),
        ),
        &(
            verifier_lcccs_acc.clone(),
            cm_i.clone(),
            proof.clone(),
            ccs.clone(),
        ),
        |b, (acc, cm_i, proof, ccs)| {
            b.iter_batched(
                || verifier_transcript.clone(),
                |mut bench_verifier_transcript| {
                    let result = NIFSVerifier::<C, R, P, PoseidonTranscript<R, CS>>::verify(
                        acc,
                        cm_i,
                        proof,
                        &mut bench_verifier_transcript,
                        ccs,
                    );
                    assert!(result.is_ok());
                },
                criterion::BatchSize::SmallInput,
            );
        },
    );
}

fn e2e_benchmarks_scalar<
    const X_LEN: usize,
    const C: usize,
    const WIT_LEN: usize,
    const W: usize,
    CS: LatticefoldChallengeSet<R> + Clone + 'static,
    R: SuitableRing,
    P: DecompositionParams,
>(
    group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>,
) {
    let r1cs_rows = X_LEN + WIT_LEN + 1;

    let (cm_i, wit, ccs, scheme) = wit_and_ccs_gen::<X_LEN, C, WIT_LEN, W, P, R>(r1cs_rows);

    prover_e2e_benchmark::<C, W, P, R, CS>(group, &cm_i, &wit, &ccs, &scheme);
    verifier_e2e_benchmark::<C, W, P, R, CS>(group, &cm_i, &wit, &ccs, &scheme);
}

fn e2e_benchmarks_non_scalar<
    const X_LEN: usize,
    const C: usize,
    const WIT_LEN: usize,
    const W: usize,
    CS: LatticefoldChallengeSet<R> + Clone + 'static,
    R: SuitableRing,
    P: DecompositionParams,
>(
    group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>,
) {
    let r1cs_rows = X_LEN + WIT_LEN + 1;

    let (cm_i, wit, ccs, scheme) =
        wit_and_ccs_gen_non_scalar::<X_LEN, C, WIT_LEN, W, P, R>(r1cs_rows);

    prover_e2e_benchmark::<C, W, P, R, CS>(group, &cm_i, &wit, &ccs, &scheme);
    verifier_e2e_benchmark::<C, W, P, R, CS>(group, &cm_i, &wit, &ccs, &scheme);
}

#[allow(unused_macros)]
macro_rules! run_single_goldilocks_benchmark {
    ($crit:expr, $io:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            e2e_benchmarks_scalar::<$io, $cw, $w, {$w * $l}, GoldilocksChallengeSet, GoldilocksRingNTT, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit);
        }
    };
}

#[allow(unused_macros)]
macro_rules! run_single_babybear_benchmark {
    ($crit:expr, $io:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            e2e_benchmarks_scalar::<$io, $cw, $w, {$w * $l}, BabyBearChallengeSet, BabyBearRingNTT, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit);
        }
    };
}

#[allow(unused_macros)]
macro_rules! run_single_starkprime_benchmark {
    ($crit:expr, $io:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            e2e_benchmarks_scalar::<$io, $cw, $w, {$w * $l}, StarkChallengeSet, StarkRingNTT, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit);
        }
    };
}

#[allow(unused_macros)]
macro_rules! run_single_frog_benchmark {
    ($crit:expr, $io:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            e2e_benchmarks_scalar::<$io, $cw, $w, {$w * $l}, FrogChallengeSet, FrogRingNTT, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit);
        }
    };
}

#[allow(unused_macros)]
macro_rules! run_single_goldilocks_non_scalar_benchmark {
    ($crit:expr, $io:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            e2e_benchmarks_non_scalar::<$io, $cw, $w, {$w * $l}, GoldilocksChallengeSet, GoldilocksRingNTT, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit);
        }
    };
}

#[allow(unused_macros)]
macro_rules! run_single_babybear_non_scalar_benchmark {
    ($crit:expr, $io:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            e2e_benchmarks_non_scalar::<$io, $cw, $w, {$w * $l}, BabyBearChallengeSet, BabyBearRingNTT, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit);
        }
    };
}

#[allow(unused_macros)]
macro_rules! run_single_starkprime_non_scalar_benchmark {
    ($crit:expr, $io:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            e2e_benchmarks_non_scalar::<$io, $cw, $w, {$w * $l}, StarkChallengeSet, StarkRingNTT, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit);
        }
    };
}

#[allow(unused_macros)]
macro_rules! run_single_frog_non_scalar_benchmark {
    ($crit:expr, $io:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            e2e_benchmarks_non_scalar::<$io, $cw, $w, {$w * $l}, FrogChallengeSet, FrogRingNTT, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit);
        }
    };
}

fn benchmarks_main(c: &mut Criterion) {
    // Goldilocks
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("E2E Goldilocks");
        group.plot_config(plot_config.clone());
        #[allow(clippy::identity_op)]
        {
            run_goldilocks_benchmarks!(group);
        }
    }

    // Godlilocks non scalar
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("E2E Goldilocks non scalar");
        group.plot_config(plot_config.clone());

        run_goldilocks_non_scalar_benchmarks!(group);
    }

    // BabyBear
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("E2E BabyBear");
        group.plot_config(plot_config.clone());
        #[allow(clippy::identity_op)]
        {
            run_babybear_benchmarks!(group);
        }
    }

    // BabyBear non scalar
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("E2E BabyBear non scalar");
        group.plot_config(plot_config.clone());

        run_babybear_non_scalar_benchmarks!(group);
    }

    // StarkPrime
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("E2E StarkPrime");
        group.plot_config(plot_config.clone());

        #[allow(clippy::identity_op)]
        {
            run_starkprime_benchmarks!(group);
        }
    }

    // StarkPrime non scalar
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("E2E StarkPrime non scalar");
        group.plot_config(plot_config.clone());

        run_starkprime_non_scalar_benchmarks!(group);
    }

    // Frog
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("E2E Frog");
        group.plot_config(plot_config.clone());
        #[allow(clippy::identity_op)]
        {
            run_frog_benchmarks!(group);
        }
    }

    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("E2E Frog non scalar");
        group.plot_config(plot_config.clone());

        run_frog_non_scalar_benchmarks!(group);
    }
}

criterion_group!(
    name=benches;
    config = Criterion::default().sample_size(10).measurement_time(Duration::from_secs(50)).warm_up_time(Duration::from_secs(1));
    targets = benchmarks_main);
criterion_main!(benches);
