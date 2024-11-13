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
    R: Clone + UniformRand + Debug + SuitableRing + for<'a> std::ops::AddAssign<&'a R>,
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
    R: Clone + UniformRand + Debug + SuitableRing + for<'a> std::ops::AddAssign<&'a R>,
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

    // // Goldilocks L = 1
    // {
    //     let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    //     let mut group = c.benchmark_group("Folding Godlilocks L = 1");
    //     group.plot_config(plot_config.clone());

    //     // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
    //     run_single_goldilocks_benchmark!(&mut group, 1, 3, 32768, 2, 1, 2, 1);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 4, 32768, 6, 1, 2, 2);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 4, 32768, 6, 1, 4, 1);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 4, 65536, 4, 1, 2, 2);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 4, 65536, 4, 1, 4, 1);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 4, 131072, 2, 1, 2, 1);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 4, 262144, 2, 1, 2, 1);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 5, 32768, 18, 1, 2, 4);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 5, 32768, 18, 1, 4, 2);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 5, 32768, 18, 1, 16, 1);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 5, 65536, 14, 1, 2, 3);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 5, 65536, 14, 1, 8, 1);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 5, 65536, 14, 1, 2, 2);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 5, 131072, 8, 1, 2, 3);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 5, 131072, 8, 1, 8, 1);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 5, 131072, 8, 1, 2, 2);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 5, 262144, 6, 1, 2, 2);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 5, 262144, 6, 1, 4, 1);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 5, 524288, 4, 1, 2, 2);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 5, 524288, 4, 1, 4, 1);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 5, 1048576, 2, 1, 2, 1);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 6, 32768, 50, 1, 2, 5);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 6, 32768, 50, 1, 32, 1);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 6, 32768, 50, 1, 2, 4);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 6, 65536, 34, 1, 2, 5);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 6, 65536, 34, 1, 32, 1);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 6, 65536, 34, 1, 2, 4);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 6, 131072, 24, 1, 2, 4);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 6, 131072, 24, 1, 4, 2);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 6, 131072, 24, 1, 16, 1);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 6, 262144, 16, 1, 2, 4);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 6, 262144, 16, 1, 4, 2);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 6, 262144, 16, 1, 16, 1);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 6, 524288, 12, 1, 2, 3);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 6, 524288, 12, 1, 8, 1);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 6, 524288, 12, 1, 2, 2);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 6, 1048576, 8, 1, 2, 3);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 6, 1048576, 8, 1, 8, 1);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 6, 1048576, 8, 1, 2, 2);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 7, 32768, 118, 1, 2, 6);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 7, 32768, 118, 1, 4, 3);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 7, 32768, 118, 1, 8, 2);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 7, 65536, 84, 1, 2, 6);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 7, 65536, 84, 1, 4, 3);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 7, 65536, 84, 1, 8, 2);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 7, 131072, 58, 1, 2, 5);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 7, 131072, 58, 1, 32, 1);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 7, 131072, 58, 1, 2, 4);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 7, 262144, 42, 1, 2, 5);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 7, 262144, 42, 1, 32, 1);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 7, 262144, 42, 1, 2, 4);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 7, 524288, 28, 1, 2, 4);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 7, 524288, 28, 1, 4, 2);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 7, 524288, 28, 1, 16, 1);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 7, 1048576, 20, 1, 2, 4);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 7, 1048576, 20, 1, 4, 2);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 7, 1048576, 20, 1, 16, 1);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 8, 32768, 264, 1, 2, 8);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 8, 32768, 264, 1, 4, 4);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 8, 32768, 264, 1, 16, 2);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 8, 65536, 186, 1, 2, 7);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 8, 65536, 186, 1, 2, 6);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 8, 65536, 186, 1, 4, 3);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 8, 131072, 132, 1, 2, 7);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 8, 131072, 132, 1, 2, 6);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 8, 131072, 132, 1, 4, 3);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 8, 262144, 92, 1, 2, 6);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 8, 262144, 92, 1, 4, 3);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 8, 262144, 92, 1, 8, 2);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 8, 524288, 66, 1, 2, 6);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 8, 524288, 66, 1, 4, 3);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 8, 524288, 66, 1, 8, 2);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 8, 1048576, 46, 1, 2, 5);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 8, 1048576, 46, 1, 32, 1);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 8, 1048576, 46, 1, 2, 4);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 9, 32768, 560, 1, 2, 9);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 9, 32768, 560, 1, 8, 3);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 9, 32768, 560, 1, 2, 8);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 9, 65536, 396, 1, 2, 8);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 9, 65536, 396, 1, 4, 4);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 9, 65536, 396, 1, 16, 2);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 9, 131072, 280, 1, 2, 8);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 9, 131072, 280, 1, 4, 4);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 9, 131072, 280, 1, 16, 2);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 9, 262144, 198, 1, 2, 7);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 9, 262144, 198, 1, 2, 6);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 9, 262144, 198, 1, 4, 3);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 9, 524288, 140, 1, 2, 7);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 9, 524288, 140, 1, 2, 6);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 9, 524288, 140, 1, 4, 3);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 9, 1048576, 98, 1, 2, 6);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 9, 1048576, 98, 1, 4, 3);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 9, 1048576, 98, 1, 8, 2);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 10, 32768, 1140, 1, 2, 10);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 10, 32768, 1140, 1, 4, 5);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 10, 32768, 1140, 1, 32, 2);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 10, 65536, 806, 1, 2, 9);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 10, 65536, 806, 1, 8, 3);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 10, 65536, 806, 1, 2, 8);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 10, 131072, 570, 1, 2, 9);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 10, 131072, 570, 1, 8, 3);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 10, 131072, 570, 1, 2, 8);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 10, 262144, 402, 1, 2, 8);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 10, 262144, 402, 1, 4, 4);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 10, 262144, 402, 1, 16, 2);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 10, 524288, 284, 1, 2, 8);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 10, 524288, 284, 1, 4, 4);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 10, 524288, 284, 1, 16, 2);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 10, 1048576, 200, 1, 2, 7);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 10, 1048576, 200, 1, 2, 6);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 10, 1048576, 200, 1, 4, 3);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 11, 32768, 2242, 1, 2, 11);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 11, 32768, 2242, 1, 2, 10);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 11, 32768, 2242, 1, 4, 5);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 11, 65536, 1584, 1, 2, 10);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 11, 65536, 1584, 1, 4, 5);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 11, 65536, 1584, 1, 32, 2);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 11, 131072, 1120, 1, 2, 10);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 11, 131072, 1120, 1, 4, 5);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 11, 131072, 1120, 1, 32, 2);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 11, 262144, 792, 1, 2, 9);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 11, 262144, 792, 1, 8, 3);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 11, 262144, 792, 1, 2, 8);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 11, 524288, 560, 1, 2, 9);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 11, 524288, 560, 1, 8, 3);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 11, 524288, 560, 1, 2, 8);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 11, 1048576, 396, 1, 2, 8);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 11, 1048576, 396, 1, 4, 4);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 11, 1048576, 396, 1, 16, 2);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 12, 32768, 4272, 1, 2, 12);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 12, 32768, 4272, 1, 4, 6);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 12, 32768, 4272, 1, 8, 4);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 12, 65536, 3020, 1, 2, 11);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 12, 65536, 3020, 1, 2, 10);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 12, 65536, 3020, 1, 4, 5);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 12, 131072, 2136, 1, 2, 11);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 12, 131072, 2136, 1, 2, 10);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 12, 131072, 2136, 1, 4, 5);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 12, 262144, 1510, 1, 2, 10);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 12, 262144, 1510, 1, 4, 5);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 12, 262144, 1510, 1, 32, 2);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 12, 524288, 1068, 1, 2, 10);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 12, 524288, 1068, 1, 4, 5);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 12, 524288, 1068, 1, 32, 2);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 12, 1048576, 754, 1, 2, 9);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 12, 1048576, 754, 1, 8, 3);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 12, 1048576, 754, 1, 2, 8);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 13, 32768, 7930, 1, 2, 12);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 13, 32768, 7930, 1, 4, 6);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 13, 32768, 7930, 1, 8, 4);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 13, 65536, 5608, 1, 2, 12);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 13, 65536, 5608, 1, 4, 6);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 13, 65536, 5608, 1, 8, 4);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 13, 131072, 3964, 1, 2, 11);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 13, 131072, 3964, 1, 2, 10);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 13, 131072, 3964, 1, 4, 5);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 13, 262144, 2804, 1, 2, 11);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 13, 262144, 2804, 1, 2, 10);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 13, 262144, 2804, 1, 4, 5);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 13, 524288, 1982, 1, 2, 10);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 13, 524288, 1982, 1, 4, 5);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 13, 524288, 1982, 1, 32, 2);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 13, 1048576, 1402, 1, 2, 10);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 13, 1048576, 1402, 1, 4, 5);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 13, 1048576, 1402, 1, 32, 2);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 14, 32768, 14382, 1, 2, 13);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 14, 32768, 14382, 1, 2, 12);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 14, 32768, 14382, 1, 4, 6);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 14, 65536, 10170, 1, 2, 13);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 14, 65536, 10170, 1, 2, 12);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 14, 65536, 10170, 1, 4, 6);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 14, 131072, 7190, 1, 2, 12);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 14, 131072, 7190, 1, 4, 6);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 14, 131072, 7190, 1, 8, 4);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 14, 262144, 5084, 1, 2, 12);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 14, 262144, 5084, 1, 4, 6);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 14, 262144, 5084, 1, 8, 4);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 14, 524288, 3594, 1, 2, 11);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 14, 524288, 3594, 1, 2, 10);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 14, 524288, 3594, 1, 4, 5);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 14, 1048576, 2542, 1, 2, 11);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 14, 1048576, 2542, 1, 2, 10);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 14, 1048576, 2542, 1, 4, 5);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 15, 32768, 25540, 1, 2, 14);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 15, 32768, 25540, 1, 4, 7);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 15, 32768, 25540, 1, 2, 13);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 15, 65536, 18060, 1, 2, 14);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 15, 65536, 18060, 1, 4, 7);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 15, 65536, 18060, 1, 2, 13);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 15, 131072, 12770, 1, 2, 13);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 15, 131072, 12770, 1, 2, 12);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 15, 131072, 12770, 1, 4, 6);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 15, 262144, 9030, 1, 2, 13);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 15, 262144, 9030, 1, 2, 12);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 15, 262144, 9030, 1, 4, 6);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 15, 524288, 6384, 1, 2, 12);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 15, 524288, 6384, 1, 4, 6);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 15, 524288, 6384, 1, 8, 4);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 15, 1048576, 4514, 1, 2, 12);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 15, 1048576, 4514, 1, 4, 6);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 15, 1048576, 4514, 1, 8, 4);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 16, 32768, 44508, 1, 2, 15);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 16, 32768, 44508, 1, 8, 5);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 16, 32768, 44508, 1, 32, 3);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 16, 65536, 31472, 1, 2, 14);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 16, 65536, 31472, 1, 4, 7);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 16, 65536, 31472, 1, 2, 13);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 16, 131072, 22254, 1, 2, 14);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 16, 131072, 22254, 1, 4, 7);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 16, 131072, 22254, 1, 2, 13);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 16, 262144, 15736, 1, 2, 13);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 16, 262144, 15736, 1, 2, 12);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 16, 262144, 15736, 1, 4, 6);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 16, 524288, 11126, 1, 2, 13);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 16, 524288, 11126, 1, 2, 12);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 16, 524288, 11126, 1, 4, 6);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 16, 1048576, 7868, 1, 2, 12);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 16, 1048576, 7868, 1, 4, 6);
    //     run_single_goldilocks_benchmark!(&mut group, 1, 16, 1048576, 7868, 1, 8, 4);
    // }

    // Godlilocks
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Decomposition Godlilocks");
        group.plot_config(plot_config.clone());

        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        run_single_goldilocks_benchmark!(&mut group, 1, 6, 512, 120, 9, 2, 7);
        run_single_goldilocks_benchmark!(&mut group, 1, 7, 512, 256, 8, 2, 8);
        run_single_goldilocks_benchmark!(&mut group, 1, 7, 512, 256, 8, 4, 4);
        run_single_goldilocks_benchmark!(&mut group, 1, 8, 512, 512, 7, 2, 9);
        run_single_goldilocks_benchmark!(&mut group, 1, 8, 1024, 512, 7, 2, 9);
        run_single_goldilocks_benchmark!(&mut group, 1, 8, 2048, 256, 8, 2, 8);
        run_single_goldilocks_benchmark!(&mut group, 1, 9, 1024, 1024, 7, 2, 10);
        run_single_goldilocks_benchmark!(&mut group, 1, 9, 2048, 512, 7, 2, 9);
        run_single_goldilocks_benchmark!(&mut group, 1, 10, 512, 2048, 6, 2, 11);
        run_single_goldilocks_benchmark!(&mut group, 1, 10, 1024, 2048, 6, 2, 11);
        run_single_goldilocks_benchmark!(&mut group, 1, 11, 1024, 4096, 6, 2, 12);
        run_single_goldilocks_benchmark!(&mut group, 1, 11, 2048, 2048, 6, 2, 12);
        run_single_goldilocks_benchmark!(&mut group, 1, 12, 1024, 8192, 6, 2, 13);
        run_single_goldilocks_benchmark!(&mut group, 1, 13, 1024, 16384, 5, 2, 14);
        run_single_goldilocks_benchmark!(&mut group, 1, 13, 2048, 8192, 5, 2, 13);
        run_single_goldilocks_benchmark!(&mut group, 1, 14, 1024, 32768, 5, 2, 15);
        run_single_goldilocks_benchmark!(&mut group, 1, 14, 2048, 16384, 5, 2, 14);
        run_single_goldilocks_benchmark!(&mut group, 1, 15, 2048, 32768, 4, 2, 15);
        run_single_goldilocks_benchmark!(&mut group, 1, 16, 2048, 65536, 4, 2, 16);
    }

    // // BabyBear
    // // TODO: Fix f_hat and account for field extensions.
    // {
    //     let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    //     let mut group = c.benchmark_group("Folding Goldilocks");
    //     group.plot_config(plot_config.clone());

    //     // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
    //     /*
    //     run_single_babybear_benchmark!(&mut group, 1, 6, 1024, 512, 4, 2, 9);
    //     run_single_babybear_benchmark!(&mut group, 1, 7, 1024, 2048, 3, 2, 11);
    //     run_single_babybear_benchmark!(&mut group, 1, 8, 4096, 2048, 3, 2, 11);
    //     run_single_babybear_benchmark!(&mut group, 1, 9, 2048, 8192, 3, 2, 13);
    //     run_single_babybear_benchmark!(&mut group, 1, 10, 4096, 16384, 3, 2, 14);
    //     */
    // }

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

    // // Frog
    // // TODO: Fix f_hat and account for field extensions.
    // {
    //     let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    //     let mut group = c.benchmark_group("Folding Frog");
    //     group.plot_config(plot_config.clone());

    //     // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
    //     /*
    //     run_single_frog_benchmark!(&mut group, 1, 5, 512, 8, 23, 2, 3);
    //     run_single_frog_benchmark!(&mut group, 1, 9, 1024, 128, 10, 2, 7);
    //     run_single_frog_benchmark!(&mut group, 1, 10, 1024, 256, 9, 2, 8);
    //     run_single_frog_benchmark!(&mut group, 1, 12, 512, 1024, 7, 2, 10);
    //     run_single_frog_benchmark!(&mut group, 1, 15, 1024, 4096, 6, 2, 12);
    //      */
    // }
}

criterion_group!(
    name=benches;
    config = Criterion::default().sample_size(10).measurement_time(Duration::from_secs(50)).warm_up_time(Duration::from_secs(1));
    targets = benchmarks_main);
criterion_main!(benches);
