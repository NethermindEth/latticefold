#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
use criterion::{
    criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion, PlotConfiguration,
};
use cyclotomic_rings::{
    challenge_set::LatticefoldChallengeSet,
    rings::{
        GoldilocksChallengeSet, GoldilocksRingNTT, StarkChallengeSet, StarkRingNTT, SuitableRing,
    },
};
use latticefold::nifs::decomposition::{
    DecompositionProver, DecompositionVerifier, LFDecompositionProver, LFDecompositionVerifier,
};
use latticefold::nifs::folding::{
    FoldingProver, FoldingVerifier, LFFoldingProver, LFFoldingVerifier,
};
use std::{fmt::Debug, time::Duration};
mod utils;
use ark_std::UniformRand;

use crate::utils::wit_and_ccs_gen;
use latticefold::{
    arith::{Witness, CCCS, CCS},
    commitment::AjtaiCommitmentScheme,
    decomposition_parameters::DecompositionParams,
    nifs::linearization::{
        LFLinearizationProver, LFLinearizationVerifier, LinearizationProver, LinearizationVerifier,
    },
    transcript::poseidon::PoseidonTranscript,
};

fn prover_folding_benchmark<
    const C: usize,
    const W: usize,
    P: DecompositionParams,
    R: Clone + UniformRand + Debug + SuitableRing,
    CS: LatticefoldChallengeSet<R>,
>(
    c: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>,
    cm_i: &CCCS<C, R>,
    wit: &Witness<R>,
    ccs: &CCS<R>,
    scheme: &AjtaiCommitmentScheme<C, W, R>,
) {
    let mut prover_transcript = PoseidonTranscript::<R, CS>::default();
    let mut verifier_transcript = PoseidonTranscript::<R, CS>::default();

    let (_, linearization_proof) = LFLinearizationProver::<_, PoseidonTranscript<R, CS>>::prove(
        cm_i,
        wit,
        &mut prover_transcript,
        ccs,
    )
    .unwrap();

    let lcccs = LFLinearizationVerifier::<_, PoseidonTranscript<R, CS>>::verify(
        cm_i,
        &linearization_proof,
        &mut verifier_transcript,
        ccs,
    )
    .unwrap();

    let (_, wit_vec, decomposition_proof) =
        LFDecompositionProver::<_, PoseidonTranscript<R, CS>>::prove::<W, C, P>(
            &lcccs,
            wit,
            &mut prover_transcript,
            ccs,
            scheme,
        )
        .unwrap();

    let lcccs_vec = LFDecompositionVerifier::<_, PoseidonTranscript<R, CS>>::verify::<C, P>(
        &lcccs,
        &decomposition_proof,
        &mut verifier_transcript,
        ccs,
    )
    .unwrap();

    let (lcccs, wit_s) = {
        let mut lcccs = lcccs_vec.clone();
        let mut lcccs_r = lcccs_vec;
        lcccs.append(&mut lcccs_r);

        let mut wit_s = wit_vec.clone();
        let mut wit_s_r = wit_vec;
        wit_s.append(&mut wit_s_r);

        (lcccs, wit_s)
    };
    c.bench_with_input(
        BenchmarkId::new(
            "Folding Prover",
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
        &(lcccs, wit_s, ccs),
        |b, (lcccs_vec, wit_vec, ccs)| {
            b.iter(|| {
                let _ = LFFoldingProver::<R, PoseidonTranscript<R, CS>>::prove::<C, P>(
                    lcccs_vec,
                    wit_vec,
                    &mut prover_transcript,
                    ccs,
                )
                .unwrap();
            })
        },
    );
}

fn verifier_folding_benchmark<
    const C: usize,
    const W: usize,
    P: DecompositionParams,
    R: Clone + UniformRand + Debug + SuitableRing,
    CS: LatticefoldChallengeSet<R>,
>(
    c: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>,
    cm_i: &CCCS<C, R>,
    wit: &Witness<R>,
    ccs: &CCS<R>,
    scheme: &AjtaiCommitmentScheme<C, W, R>,
) {
    let mut prover_transcript = PoseidonTranscript::<R, CS>::default();
    let mut verifier_transcript = PoseidonTranscript::<R, CS>::default();

    let (_, linearization_proof) = LFLinearizationProver::<_, PoseidonTranscript<R, CS>>::prove(
        cm_i,
        wit,
        &mut prover_transcript,
        ccs,
    )
    .unwrap();

    let lcccs = LFLinearizationVerifier::<_, PoseidonTranscript<R, CS>>::verify(
        cm_i,
        &linearization_proof,
        &mut verifier_transcript,
        ccs,
    )
    .unwrap();

    let (_, wit_vec, decomposition_proof) =
        LFDecompositionProver::<_, PoseidonTranscript<R, CS>>::prove::<W, C, P>(
            &lcccs,
            wit,
            &mut prover_transcript,
            ccs,
            scheme,
        )
        .unwrap();

    let lcccs_vec = LFDecompositionVerifier::<_, PoseidonTranscript<R, CS>>::verify::<C, P>(
        &lcccs,
        &decomposition_proof,
        &mut verifier_transcript,
        ccs,
    )
    .unwrap();

    let (lcccs, wit_s) = {
        let mut lcccs = lcccs_vec.clone();
        let mut lcccs_r = lcccs_vec;
        lcccs.append(&mut lcccs_r);

        let mut wit_s = wit_vec.clone();
        let mut wit_s_r = wit_vec;
        wit_s.append(&mut wit_s_r);

        (lcccs, wit_s)
    };

    let (_, _, folding_proof) = LFFoldingProver::<R, PoseidonTranscript<R, CS>>::prove::<C, P>(
        &lcccs,
        &wit_s,
        &mut prover_transcript,
        ccs,
    )
    .unwrap();

    c.bench_with_input(
        BenchmarkId::new(
            "Folding Verifier",
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
        &(lcccs, folding_proof, ccs),
        |b, (lcccs_vec, proof, ccs)| {
            b.iter(|| {
                let _ = LFFoldingVerifier::<_, PoseidonTranscript<R, CS>>::verify::<C, P>(
                    lcccs_vec,
                    proof,
                    &mut verifier_transcript,
                    ccs,
                );
            })
        },
    );
}

fn folding_benchmarks<
    const X_LEN: usize,
    const C: usize,
    const WIT_LEN: usize,
    const W: usize,
    CS: LatticefoldChallengeSet<R>,
    R: SuitableRing,
    P: DecompositionParams,
>(
    group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>,
) {
    let r1cs_rows = X_LEN + WIT_LEN + 1; // This makes a square matrix but is too much memory;

    let (cm_i, wit, ccs, scheme) = wit_and_ccs_gen::<X_LEN, C, WIT_LEN, W, P, R>(r1cs_rows);

    prover_folding_benchmark::<C, W, P, R, CS>(group, &cm_i, &wit, &ccs, &scheme);

    verifier_folding_benchmark::<C, W, P, R, CS>(group, &cm_i, &wit, &ccs, &scheme);
}

macro_rules! define_params {
    ($w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        paste::paste! {

            #[derive(Clone)]
            struct [<DecompParamsWithB $b W $w b $b_small K $k>];

            impl DecompositionParams for [<DecompParamsWithB $b W $w b $b_small K $k>] {
                const B: u128 = $b;
                const L: usize = $l;
                const B_SMALL: usize = $b_small;
                const K: usize = $k;
            }
        }
    };
}

#[allow(unused_macros)]
macro_rules! run_single_goldilocks_benchmark {
    ($crit:expr, $io:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            folding_benchmarks::<$io, $cw, $w, {$w * $l}, GoldilocksChallengeSet, GoldilocksRingNTT, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit);
        }
    };
}

// Baybear parameters
#[allow(unused_macros)]
macro_rules! run_single_babybear_benchmark {
    ($crit:expr, $io:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            folding_benchmarks::<$io, $cw, $w, {$w * $l}, BabyBearChallengeSet, BabyBearRingNTT, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit);
        }
    };
}

// Stark parameters
macro_rules! run_single_starkprime_benchmark {
    ($crit:expr, $io:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            folding_benchmarks::<$io, $cw, $w, {$w * $l}, StarkChallengeSet, StarkRingNTT, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit);
        }
    };
}

// Frog parameters
#[allow(unused_macros)]
macro_rules! run_single_frog_benchmark {
    ($crit:expr, $io:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            folding_benchmarks::<$io, $cw, $w, {$w * $l}, FrogChallengeSet, FrogRingNTT, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit);
        }
    };
}

fn benchmarks_main(c: &mut Criterion) {
    // Godlilocks
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Folding Godlilocks");
        group.plot_config(plot_config.clone());

        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        // run_single_goldilocks_benchmark!(&mut group, 1, 6, 512, 120, 9, 2, 7);
        // run_single_goldilocks_benchmark!(&mut group, 1, 7, 512, 256, 8, 2, 8);
        // run_single_goldilocks_benchmark!(&mut group, 1, 7, 512, 256, 8, 4, 4);
        // run_single_goldilocks_benchmark!(&mut group, 1, 8, 512, 512, 7, 2, 9);
        // run_single_goldilocks_benchmark!(&mut group, 1, 8, 1024, 512, 7, 2, 9);
        // run_single_goldilocks_benchmark!(&mut group, 1, 8, 2048, 256, 8, 2, 8);
        // run_single_goldilocks_benchmark!(&mut group, 1, 9, 1024, 1024, 7, 2, 10);
        // run_single_goldilocks_benchmark!(&mut group, 1, 9, 2048, 512, 7, 2, 9);
        // run_single_goldilocks_benchmark!(&mut group, 1, 10, 512, 2048, 6, 2, 11);
        // run_single_goldilocks_benchmark!(&mut group, 1, 10, 1024, 2048, 6, 2, 11);
        // run_single_goldilocks_benchmark!(&mut group, 1, 11, 1024, 4096, 6, 2, 12);
        // run_single_goldilocks_benchmark!(&mut group, 1, 11, 2048, 2048, 6, 2, 12);
        // run_single_goldilocks_benchmark!(&mut group, 1, 12, 1024, 8192, 6, 2, 13);
        // run_single_goldilocks_benchmark!(&mut group, 1, 13, 1024, 16384, 5, 2, 14);
        // run_single_goldilocks_benchmark!(&mut group, 1, 13, 2048, 8192, 5, 2, 13);
        // run_single_goldilocks_benchmark!(&mut group, 1, 14, 1024, 32768, 5, 2, 15);
        // run_single_goldilocks_benchmark!(&mut group, 1, 14, 2048, 16384, 5, 2, 14);
        // run_single_goldilocks_benchmark!(&mut group, 1, 15, 2048, 32768, 4, 2, 15);
        run_single_goldilocks_benchmark!(&mut group, 1, 16, 16384, 65536, 4, 2, 16);
    }

    // BabyBear
    // TODO: Fix f_hat and account for field extensions.
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Folding Goldilocks");
        group.plot_config(plot_config.clone());

        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        /*
        run_single_babybear_benchmark!(&mut group, 1, 6, 1024, 512, 4, 2, 9);
        run_single_babybear_benchmark!(&mut group, 1, 7, 1024, 2048, 3, 2, 11);
        run_single_babybear_benchmark!(&mut group, 1, 8, 4096, 2048, 3, 2, 11);
        run_single_babybear_benchmark!(&mut group, 1, 9, 2048, 8192, 3, 2, 13);
        run_single_babybear_benchmark!(&mut group, 1, 10, 4096, 16384, 3, 2, 14);
        */
    }

    // // StarkPrime
    // {
    //     let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    //     let mut group = c.benchmark_group("Folding StarkPrime");
    //     group.plot_config(plot_config.clone());

    //     // Parameters Criterion, X_LEN, C, W, B, L, B_small, K 3052596316
    //     #[allow(clippy::identity_op)]
    //     {
    //         run_single_starkprime_benchmark!(&mut group, 1, 15, 1024, 3052596316u128, 1, 2, 30);
    //         run_single_starkprime_benchmark!(&mut group, 1, 16, 1024, 4294967296u128, 1, 2, 32);
    //         run_single_starkprime_benchmark!(&mut group, 1, 17, 2048, 8589934592u128, 1, 2, 33);
    //         run_single_starkprime_benchmark!(&mut group, 1, 18, 2048, 20833367754u128, 1, 2, 34);
    //         run_single_starkprime_benchmark!(&mut group, 1, 19, 2048, 34359738368u128, 1, 2, 35);
    //     }
    // }

    // Frog
    // TODO: Fix f_hat and account for field extensions.
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Folding Frog");
        group.plot_config(plot_config.clone());

        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        /*
        run_single_frog_benchmark!(&mut group, 1, 5, 512, 8, 23, 2, 3);
        run_single_frog_benchmark!(&mut group, 1, 9, 1024, 128, 10, 2, 7);
        run_single_frog_benchmark!(&mut group, 1, 10, 1024, 256, 9, 2, 8);
        run_single_frog_benchmark!(&mut group, 1, 12, 512, 1024, 7, 2, 10);
        run_single_frog_benchmark!(&mut group, 1, 15, 1024, 4096, 6, 2, 12);
         */
    }
}

criterion_group!(
    name=benches;
    config = Criterion::default().sample_size(10).measurement_time(Duration::from_secs(50)).warm_up_time(Duration::from_secs(1));
    targets = benchmarks_main);
criterion_main!(benches);
