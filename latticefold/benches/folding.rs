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
    // Goldilocks L = 1
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Linearization Godlilocks");
        group.plot_config(plot_config.clone());

        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        
        run_single_goldilocks_benchmark!(&mut group, 1, 2, 512, 4, 1, 2, 2);
        run_single_goldilocks_benchmark!(&mut group, 1, 2, 512, 4, 1, 4, 1);
        run_single_goldilocks_benchmark!(&mut group, 1, 3, 512, 16, 1, 2, 4);
        run_single_goldilocks_benchmark!(&mut group, 1, 3, 512, 16, 1, 4, 2);
        run_single_goldilocks_benchmark!(&mut group, 1, 3, 512, 16, 1, 16, 1);
        run_single_goldilocks_benchmark!(&mut group, 1, 3, 1024, 12, 1, 2, 3);
        run_single_goldilocks_benchmark!(&mut group, 1, 3, 1024, 12, 1, 8, 1);
        run_single_goldilocks_benchmark!(&mut group, 1, 3, 1024, 12, 1, 2, 2);
        run_single_goldilocks_benchmark!(&mut group, 1, 3, 2048, 8, 1, 2, 3);
        run_single_goldilocks_benchmark!(&mut group, 1, 3, 2048, 8, 1, 8, 1);
        run_single_goldilocks_benchmark!(&mut group, 1, 3, 2048, 8, 1, 2, 2);
        run_single_goldilocks_benchmark!(&mut group, 1, 3, 4096, 6, 1, 2, 2);
        run_single_goldilocks_benchmark!(&mut group, 1, 3, 4096, 6, 1, 4, 1);
        run_single_goldilocks_benchmark!(&mut group, 1, 4, 512, 56, 1, 2, 5);
        run_single_goldilocks_benchmark!(&mut group, 1, 4, 512, 56, 1, 32, 1);
        run_single_goldilocks_benchmark!(&mut group, 1, 4, 512, 56, 1, 2, 4);
        run_single_goldilocks_benchmark!(&mut group, 1, 4, 1024, 40, 1, 2, 5);
        run_single_goldilocks_benchmark!(&mut group, 1, 4, 1024, 40, 1, 32, 1);
        run_single_goldilocks_benchmark!(&mut group, 1, 4, 1024, 40, 1, 2, 4);
        run_single_goldilocks_benchmark!(&mut group, 1, 4, 2048, 28, 1, 2, 4);
        run_single_goldilocks_benchmark!(&mut group, 1, 4, 2048, 28, 1, 4, 2);
        run_single_goldilocks_benchmark!(&mut group, 1, 4, 2048, 28, 1, 16, 1);
        run_single_goldilocks_benchmark!(&mut group, 1, 4, 4096, 20, 1, 2, 4);
        run_single_goldilocks_benchmark!(&mut group, 1, 4, 4096, 20, 1, 4, 2);
        run_single_goldilocks_benchmark!(&mut group, 1, 4, 4096, 20, 1, 16, 1);
        run_single_goldilocks_benchmark!(&mut group, 1, 5, 512, 158, 1, 2, 7);
        run_single_goldilocks_benchmark!(&mut group, 1, 5, 512, 158, 1, 2, 6);
        run_single_goldilocks_benchmark!(&mut group, 1, 5, 512, 158, 1, 4, 3);
        run_single_goldilocks_benchmark!(&mut group, 1, 5, 1024, 112, 1, 2, 6);
        run_single_goldilocks_benchmark!(&mut group, 1, 5, 1024, 112, 1, 4, 3);
        run_single_goldilocks_benchmark!(&mut group, 1, 5, 1024, 112, 1, 8, 2);
        run_single_goldilocks_benchmark!(&mut group, 1, 5, 2048, 78, 1, 2, 6);
        run_single_goldilocks_benchmark!(&mut group, 1, 5, 2048, 78, 1, 4, 3);
        run_single_goldilocks_benchmark!(&mut group, 1, 5, 2048, 78, 1, 8, 2);
        run_single_goldilocks_benchmark!(&mut group, 1, 5, 4096, 56, 1, 2, 5);
        run_single_goldilocks_benchmark!(&mut group, 1, 5, 4096, 56, 1, 32, 1);
        run_single_goldilocks_benchmark!(&mut group, 1, 5, 4096, 56, 1, 2, 4);
        run_single_goldilocks_benchmark!(&mut group, 1, 6, 512, 404, 1, 2, 8);
        run_single_goldilocks_benchmark!(&mut group, 1, 6, 512, 404, 1, 4, 4);
        run_single_goldilocks_benchmark!(&mut group, 1, 6, 512, 404, 1, 16, 2);
        run_single_goldilocks_benchmark!(&mut group, 1, 6, 1024, 286, 1, 2, 8);
        run_single_goldilocks_benchmark!(&mut group, 1, 6, 1024, 286, 1, 4, 4);
        run_single_goldilocks_benchmark!(&mut group, 1, 6, 1024, 286, 1, 16, 2);
        run_single_goldilocks_benchmark!(&mut group, 1, 6, 2048, 202, 1, 2, 7);
        run_single_goldilocks_benchmark!(&mut group, 1, 6, 2048, 202, 1, 2, 6);
        run_single_goldilocks_benchmark!(&mut group, 1, 6, 2048, 202, 1, 4, 3);
        run_single_goldilocks_benchmark!(&mut group, 1, 6, 4096, 142, 1, 2, 7);
        run_single_goldilocks_benchmark!(&mut group, 1, 6, 4096, 142, 1, 2, 6);
        run_single_goldilocks_benchmark!(&mut group, 1, 6, 4096, 142, 1, 4, 3);
        run_single_goldilocks_benchmark!(&mut group, 1, 7, 512, 954, 1, 2, 9);
        run_single_goldilocks_benchmark!(&mut group, 1, 7, 512, 954, 1, 8, 3);
        run_single_goldilocks_benchmark!(&mut group, 1, 7, 512, 954, 1, 2, 8);
        run_single_goldilocks_benchmark!(&mut group, 1, 7, 1024, 674, 1, 2, 9);
        run_single_goldilocks_benchmark!(&mut group, 1, 7, 1024, 674, 1, 8, 3);
        run_single_goldilocks_benchmark!(&mut group, 1, 7, 1024, 674, 1, 2, 8);
        run_single_goldilocks_benchmark!(&mut group, 1, 7, 2048, 476, 1, 2, 8);
        run_single_goldilocks_benchmark!(&mut group, 1, 7, 2048, 476, 1, 4, 4);
        run_single_goldilocks_benchmark!(&mut group, 1, 7, 2048, 476, 1, 16, 2);
        run_single_goldilocks_benchmark!(&mut group, 1, 7, 4096, 336, 1, 2, 8);
        run_single_goldilocks_benchmark!(&mut group, 1, 7, 4096, 336, 1, 4, 4);
        run_single_goldilocks_benchmark!(&mut group, 1, 7, 4096, 336, 1, 16, 2);
        run_single_goldilocks_benchmark!(&mut group, 1, 8, 512, 2120, 1, 2, 11);
        run_single_goldilocks_benchmark!(&mut group, 1, 8, 512, 2120, 1, 2, 10);
        run_single_goldilocks_benchmark!(&mut group, 1, 8, 512, 2120, 1, 4, 5);
        run_single_goldilocks_benchmark!(&mut group, 1, 8, 1024, 1500, 1, 2, 10);
        run_single_goldilocks_benchmark!(&mut group, 1, 8, 1024, 1500, 1, 4, 5);
        run_single_goldilocks_benchmark!(&mut group, 1, 8, 1024, 1500, 1, 32, 2);
        run_single_goldilocks_benchmark!(&mut group, 1, 8, 2048, 1060, 1, 2, 10);
        run_single_goldilocks_benchmark!(&mut group, 1, 8, 2048, 1060, 1, 4, 5);
        run_single_goldilocks_benchmark!(&mut group, 1, 8, 2048, 1060, 1, 32, 2);
        run_single_goldilocks_benchmark!(&mut group, 1, 8, 4096, 750, 1, 2, 9);
        run_single_goldilocks_benchmark!(&mut group, 1, 8, 4096, 750, 1, 8, 3);
        run_single_goldilocks_benchmark!(&mut group, 1, 8, 4096, 750, 1, 2, 8);
        run_single_goldilocks_benchmark!(&mut group, 1, 9, 512, 4492, 1, 2, 12);
        run_single_goldilocks_benchmark!(&mut group, 1, 9, 512, 4492, 1, 4, 6);
        run_single_goldilocks_benchmark!(&mut group, 1, 9, 512, 4492, 1, 8, 4);
        run_single_goldilocks_benchmark!(&mut group, 1, 9, 1024, 3176, 1, 2, 11);
        run_single_goldilocks_benchmark!(&mut group, 1, 9, 1024, 3176, 1, 2, 10);
        run_single_goldilocks_benchmark!(&mut group, 1, 9, 1024, 3176, 1, 4, 5);
        run_single_goldilocks_benchmark!(&mut group, 1, 9, 2048, 2246, 1, 2, 11);
        run_single_goldilocks_benchmark!(&mut group, 1, 9, 2048, 2246, 1, 2, 10);
        run_single_goldilocks_benchmark!(&mut group, 1, 9, 2048, 2246, 1, 4, 5);
        run_single_goldilocks_benchmark!(&mut group, 1, 9, 4096, 1588, 1, 2, 10);
        run_single_goldilocks_benchmark!(&mut group, 1, 9, 4096, 1588, 1, 4, 5);
        run_single_goldilocks_benchmark!(&mut group, 1, 9, 4096, 1588, 1, 32, 2);
        run_single_goldilocks_benchmark!(&mut group, 1, 10, 512, 9132, 1, 2, 13);
        run_single_goldilocks_benchmark!(&mut group, 1, 10, 512, 9132, 1, 2, 12);
        run_single_goldilocks_benchmark!(&mut group, 1, 10, 512, 9132, 1, 4, 6);
        run_single_goldilocks_benchmark!(&mut group, 1, 10, 1024, 6458, 1, 2, 12);
        run_single_goldilocks_benchmark!(&mut group, 1, 10, 1024, 6458, 1, 4, 6);
        run_single_goldilocks_benchmark!(&mut group, 1, 10, 1024, 6458, 1, 8, 4);
        run_single_goldilocks_benchmark!(&mut group, 1, 10, 2048, 4566, 1, 2, 12);
        run_single_goldilocks_benchmark!(&mut group, 1, 10, 2048, 4566, 1, 4, 6);
        run_single_goldilocks_benchmark!(&mut group, 1, 10, 2048, 4566, 1, 8, 4);
        run_single_goldilocks_benchmark!(&mut group, 1, 10, 4096, 3228, 1, 2, 11);
        run_single_goldilocks_benchmark!(&mut group, 1, 10, 4096, 3228, 1, 2, 10);
        run_single_goldilocks_benchmark!(&mut group, 1, 10, 4096, 3228, 1, 4, 5);
        run_single_goldilocks_benchmark!(&mut group, 1, 11, 512, 17936, 1, 2, 14);
        run_single_goldilocks_benchmark!(&mut group, 1, 11, 512, 17936, 1, 4, 7);
        run_single_goldilocks_benchmark!(&mut group, 1, 11, 512, 17936, 1, 2, 13);
        run_single_goldilocks_benchmark!(&mut group, 1, 11, 1024, 12682, 1, 2, 13);
        run_single_goldilocks_benchmark!(&mut group, 1, 11, 1024, 12682, 1, 2, 12);
        run_single_goldilocks_benchmark!(&mut group, 1, 11, 1024, 12682, 1, 4, 6);
        run_single_goldilocks_benchmark!(&mut group, 1, 11, 2048, 8968, 1, 2, 13);
        run_single_goldilocks_benchmark!(&mut group, 1, 11, 2048, 8968, 1, 2, 12);
        run_single_goldilocks_benchmark!(&mut group, 1, 11, 2048, 8968, 1, 4, 6);
        run_single_goldilocks_benchmark!(&mut group, 1, 11, 4096, 6340, 1, 2, 12);
        run_single_goldilocks_benchmark!(&mut group, 1, 11, 4096, 6340, 1, 4, 6);
        run_single_goldilocks_benchmark!(&mut group, 1, 11, 4096, 6340, 1, 8, 4);
        run_single_goldilocks_benchmark!(&mut group, 1, 12, 512, 34184, 1, 2, 15);
        run_single_goldilocks_benchmark!(&mut group, 1, 12, 512, 34184, 1, 8, 5);
        run_single_goldilocks_benchmark!(&mut group, 1, 12, 512, 34184, 1, 32, 3);
        run_single_goldilocks_benchmark!(&mut group, 1, 12, 1024, 24172, 1, 2, 14);
        run_single_goldilocks_benchmark!(&mut group, 1, 12, 1024, 24172, 1, 4, 7);
        run_single_goldilocks_benchmark!(&mut group, 1, 12, 1024, 24172, 1, 2, 13);
        run_single_goldilocks_benchmark!(&mut group, 1, 12, 2048, 17092, 1, 2, 14);
        run_single_goldilocks_benchmark!(&mut group, 1, 12, 2048, 17092, 1, 4, 7);
        run_single_goldilocks_benchmark!(&mut group, 1, 12, 2048, 17092, 1, 2, 13);
        run_single_goldilocks_benchmark!(&mut group, 1, 12, 4096, 12086, 1, 2, 13);
        run_single_goldilocks_benchmark!(&mut group, 1, 12, 4096, 12086, 1, 2, 12);
        run_single_goldilocks_benchmark!(&mut group, 1, 12, 4096, 12086, 1, 4, 6);
        run_single_goldilocks_benchmark!(&mut group, 1, 13, 512, 63452, 1, 2, 15);
        run_single_goldilocks_benchmark!(&mut group, 1, 13, 512, 63452, 1, 8, 5);
        run_single_goldilocks_benchmark!(&mut group, 1, 13, 512, 63452, 1, 32, 3);
        run_single_goldilocks_benchmark!(&mut group, 1, 13, 1024, 44868, 1, 2, 15);
        run_single_goldilocks_benchmark!(&mut group, 1, 13, 1024, 44868, 1, 8, 5);
        run_single_goldilocks_benchmark!(&mut group, 1, 13, 1024, 44868, 1, 32, 3);
        run_single_goldilocks_benchmark!(&mut group, 1, 13, 2048, 31726, 1, 2, 14);
        run_single_goldilocks_benchmark!(&mut group, 1, 13, 2048, 31726, 1, 4, 7);
        run_single_goldilocks_benchmark!(&mut group, 1, 13, 2048, 31726, 1, 2, 13);
        run_single_goldilocks_benchmark!(&mut group, 1, 13, 4096, 22434, 1, 2, 14);
        run_single_goldilocks_benchmark!(&mut group, 1, 13, 4096, 22434, 1, 4, 7);
        run_single_goldilocks_benchmark!(&mut group, 1, 13, 4096, 22434, 1, 2, 13);
        run_single_goldilocks_benchmark!(&mut group, 1, 14, 512, 115062, 1, 2, 16);
        run_single_goldilocks_benchmark!(&mut group, 1, 14, 512, 115062, 1, 4, 8);
        run_single_goldilocks_benchmark!(&mut group, 1, 14, 512, 115062, 1, 16, 4);
        run_single_goldilocks_benchmark!(&mut group, 1, 14, 1024, 81360, 1, 2, 16);
        run_single_goldilocks_benchmark!(&mut group, 1, 14, 1024, 81360, 1, 4, 8);
        run_single_goldilocks_benchmark!(&mut group, 1, 14, 1024, 81360, 1, 16, 4);
        run_single_goldilocks_benchmark!(&mut group, 1, 14, 2048, 57530, 1, 2, 15);
        run_single_goldilocks_benchmark!(&mut group, 1, 14, 2048, 57530, 1, 8, 5);
        run_single_goldilocks_benchmark!(&mut group, 1, 14, 2048, 57530, 1, 32, 3);
        run_single_goldilocks_benchmark!(&mut group, 1, 14, 4096, 40680, 1, 2, 15);
        run_single_goldilocks_benchmark!(&mut group, 1, 14, 4096, 40680, 1, 8, 5);
        run_single_goldilocks_benchmark!(&mut group, 1, 14, 4096, 40680, 1, 32, 3);
        run_single_goldilocks_benchmark!(&mut group, 1, 15, 512, 204328, 1, 2, 17);
        run_single_goldilocks_benchmark!(&mut group, 1, 15, 512, 204328, 1, 2, 16);
        run_single_goldilocks_benchmark!(&mut group, 1, 15, 512, 204328, 1, 4, 8);
        run_single_goldilocks_benchmark!(&mut group, 1, 15, 1024, 144482, 1, 2, 17);
        run_single_goldilocks_benchmark!(&mut group, 1, 15, 1024, 144482, 1, 2, 16);
        run_single_goldilocks_benchmark!(&mut group, 1, 15, 1024, 144482, 1, 4, 8);
        run_single_goldilocks_benchmark!(&mut group, 1, 15, 2048, 102164, 1, 2, 16);
        run_single_goldilocks_benchmark!(&mut group, 1, 15, 2048, 102164, 1, 4, 8);
        run_single_goldilocks_benchmark!(&mut group, 1, 15, 2048, 102164, 1, 16, 4);
        run_single_goldilocks_benchmark!(&mut group, 1, 15, 4096, 72240, 1, 2, 16);
        run_single_goldilocks_benchmark!(&mut group, 1, 15, 4096, 72240, 1, 4, 8);
        run_single_goldilocks_benchmark!(&mut group, 1, 15, 4096, 72240, 1, 16, 4);
        run_single_goldilocks_benchmark!(&mut group, 1, 16, 512, 356076, 1, 2, 18);
        run_single_goldilocks_benchmark!(&mut group, 1, 16, 512, 356076, 1, 4, 9);
        run_single_goldilocks_benchmark!(&mut group, 1, 16, 512, 356076, 1, 8, 6);
        run_single_goldilocks_benchmark!(&mut group, 1, 16, 1024, 251784, 1, 2, 17);
        run_single_goldilocks_benchmark!(&mut group, 1, 16, 1024, 251784, 1, 2, 16);
        run_single_goldilocks_benchmark!(&mut group, 1, 16, 1024, 251784, 1, 4, 8);
        run_single_goldilocks_benchmark!(&mut group, 1, 16, 2048, 178038, 1, 2, 17);
        run_single_goldilocks_benchmark!(&mut group, 1, 16, 2048, 178038, 1, 2, 16);
        run_single_goldilocks_benchmark!(&mut group, 1, 16, 2048, 178038, 1, 4, 8);
        run_single_goldilocks_benchmark!(&mut group, 1, 16, 4096, 125892, 1, 2, 16);
        run_single_goldilocks_benchmark!(&mut group, 1, 16, 4096, 125892, 1, 4, 8);
        run_single_goldilocks_benchmark!(&mut group, 1, 16, 4096, 125892, 1, 16, 4);
    } 

    // Godlilocks
    {
    let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    let mut group = c.benchmark_group("Decomposition Godlilocks");
    group.plot_config(plot_config.clone());

    // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
    run_single_goldilocks_benchmark!(&mut group, 1, 6, 512, 128, 9, 2, 7);
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
