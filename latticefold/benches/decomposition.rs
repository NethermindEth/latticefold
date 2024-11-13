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
use rand::thread_rng;
use std::{fmt::Debug, time::Duration};
mod utils;
use ark_std::UniformRand;
use latticefold::{
    arith::{r1cs::get_test_dummy_z_split, Arith, Witness, CCCS, CCS},
    commitment::AjtaiCommitmentScheme,
    decomposition_parameters::DecompositionParams,
    nifs::{
        decomposition::{
            DecompositionProver, DecompositionVerifier, LFDecompositionProver,
            LFDecompositionVerifier,
        },
        linearization::{
            LFLinearizationProver, LFLinearizationVerifier, LinearizationProver,
            LinearizationVerifier,
        },
    },
    transcript::poseidon::PoseidonTranscript,
};
use utils::get_test_dummy_ccs;

fn wit_and_ccs_gen<
    const X_LEN: usize,
    const C: usize, // rows
    const WIT_LEN: usize,
    const W: usize, // columns
    P: DecompositionParams,
    R: Clone + UniformRand + Debug + SuitableRing + for<'a> std::ops::AddAssign<&'a R>,
>(
    r1cs_rows: usize,
) -> (
    CCCS<C, R>,
    Witness<R>,
    CCS<R>,
    AjtaiCommitmentScheme<C, W, R>,
) {
    //TODO: Ensure we draw elements below bound
    let ccs = get_test_dummy_ccs::<R, X_LEN, WIT_LEN, W>(r1cs_rows);
    let (one, x_ccs, w_ccs) = get_test_dummy_z_split::<R, X_LEN, WIT_LEN>();
    let mut z = vec![one];
    z.extend(&x_ccs);
    z.extend(&w_ccs);
    match ccs.check_relation(&z) {
        Ok(_) => println!("R1CS valid!"),
        Err(_) => println!("R1CS invalid"),
    }

    let scheme: AjtaiCommitmentScheme<C, W, R> = AjtaiCommitmentScheme::rand(&mut thread_rng());

    let wit: Witness<R> = Witness::from_w_ccs::<P>(&w_ccs);
    let cm_i: CCCS<C, R> = CCCS {
        cm: wit.commit::<C, W, P>(&scheme).unwrap(),
        x_ccs,
    };

    (cm_i, wit, ccs, scheme)
}

fn prover_decomposition_benchmark<
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
) where
    for<'a> <R as cyclotomic_rings::rings::SuitableRing>::CoefficientRepresentation:
        std::ops::MulAssign<&'a u128>,
{
    println!("Proving decomposition");
    println!("transcript");
    let mut prover_transcript = PoseidonTranscript::<R, CS>::default();
    let mut verifier_transcript = PoseidonTranscript::<R, CS>::default();

    println!("prove linearization");
    let (_, linearization_proof) = LFLinearizationProver::<_, PoseidonTranscript<R, CS>>::prove(
        cm_i,
        wit,
        &mut prover_transcript,
        ccs,
    )
    .unwrap();

    println!("verify linearization");
    let lcccs = LFLinearizationVerifier::<_, PoseidonTranscript<R, CS>>::verify(
        cm_i,
        &linearization_proof,
        &mut verifier_transcript,
        ccs,
    )
    .unwrap();

    c.bench_with_input(
        BenchmarkId::new(
            "Decomposition Prover",
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
        &(lcccs, wit, ccs),
        |b, (lcccs, wit, ccs)| {
            b.iter(|| {
                let (_, _, _) = LFDecompositionProver::<_, PoseidonTranscript<R, CS>>::prove::<
                    W,
                    C,
                    P,
                >(lcccs, wit, &mut prover_transcript, ccs, scheme)
                .unwrap();
            })
        },
    );
}

fn verifier_decomposition_benchmark<
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
) where
    for<'a> <R as cyclotomic_rings::rings::SuitableRing>::CoefficientRepresentation:
        std::ops::MulAssign<&'a u128>,
{
    println!("verify decomposition");
    println!("transcript");
    let mut prover_transcript = PoseidonTranscript::<R, CS>::default();
    let mut verifier_transcript = PoseidonTranscript::<R, CS>::default();

    println!("prove linearization");
    let (_, linearization_proof) = LFLinearizationProver::<_, PoseidonTranscript<R, CS>>::prove(
        cm_i,
        wit,
        &mut prover_transcript,
        ccs,
    )
    .unwrap();

    println!("verify linearization");
    let lcccs = LFLinearizationVerifier::<_, PoseidonTranscript<R, CS>>::verify(
        cm_i,
        &linearization_proof,
        &mut verifier_transcript,
        ccs,
    )
    .unwrap();

    println!("prove decomposition");
    let (_, _, decomposition_proof) =
        LFDecompositionProver::<_, PoseidonTranscript<R, CS>>::prove::<W, C, P>(
            &lcccs,
            wit,
            &mut prover_transcript,
            ccs,
            scheme,
        )
        .unwrap();

    println!("verify decomposition");
    c.bench_with_input(
        BenchmarkId::new(
            "Decomposition Verifier",
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
        &(lcccs, decomposition_proof, ccs),
        |b, (lcccs, proof, ccs)| {
            b.iter(|| {
                let _ = LFDecompositionVerifier::<_, PoseidonTranscript<R, CS>>::verify::<C, P>(
                    lcccs,
                    proof,
                    &mut verifier_transcript,
                    ccs,
                );
            })
        },
    );
}

fn decomposition_benchmarks<
    const X_LEN: usize,
    const C: usize,
    const WIT_LEN: usize,
    const W: usize,
    CS: LatticefoldChallengeSet<R>,
    R: Clone + UniformRand + Debug + SuitableRing + for<'a> std::ops::AddAssign<&'a R>,
    P: DecompositionParams + Clone,
>(
    group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>,
) where
    for<'a> <R as cyclotomic_rings::rings::SuitableRing>::CoefficientRepresentation:
        std::ops::MulAssign<&'a u128>,
{
    let r1cs_rows = X_LEN + WIT_LEN + 1;
    println!("Witness generation");
    let (cm_i, wit, ccs, scheme) = wit_and_ccs_gen::<X_LEN, C, WIT_LEN, W, P, R>(r1cs_rows);
    // N/Q = prime / degree

    prover_decomposition_benchmark::<C, W, P, R, CS>(group, &cm_i, &wit, &ccs, &scheme);

    verifier_decomposition_benchmark::<C, W, P, R, CS>(group, &cm_i, &wit, &ccs, &scheme);
}

// Macros
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

macro_rules! run_single_starkprime_benchmark {
    ($crit:expr, $io:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            decomposition_benchmarks::<$io, $cw, $w,{$w * $l}, StarkChallengeSet, StarkRingNTT, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit);
        }
    };
}

#[macro_export]
macro_rules! run_single_goldilocks_benchmark {
    ($crit:expr, $io:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
         define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            decomposition_benchmarks::<$io, $cw, $w, {$w * $l}, GoldilocksChallengeSet, GoldilocksRingNTT, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit);

        }
    };
}
#[macro_export]
macro_rules! run_single_babybear_benchmark {
    ($crit:expr, $io:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            decomposition_benchmarks::<$io, $cw, $w, {$w * $l}, BabyBearChallengeSet, BabyBearRingNTT, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit);

        }
    };
}
#[macro_export]
macro_rules! run_single_frog_benchmark {
    ($crit:expr, $io:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            decomposition_benchmarks::<$io, $cw, $w, {$w * $l}, FrogChallengeSet, FrogRingNTT, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit);

        }
    };
}

fn benchmarks_main(c: &mut Criterion) {
    // // Babybear
    // {
    //     let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    //     let mut group = c.benchmark_group("Decomposition BabyBear");
    //     group.plot_config(plot_config.clone());

    //     // TODO: Update configurations
    //     run_single_babybear_benchmark!(1, &mut group, 6, 1024, 10, 2);
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

    // StarkPrime
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Decomposition StarkPrime");
        group.plot_config(plot_config.clone());

        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        run_single_starkprime_benchmark!(&mut group, 1, 15, 512, 8633754724, 8, 92918, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 15, 512, 8615125000, 8, 2050, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 15, 512, 8540717056, 8, 304, 4);
        run_single_starkprime_benchmark!(&mut group, 1, 15, 1024, 6104921956, 8, 78134, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 15, 1024, 6088387976, 8, 1826, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 15, 1024, 5972816656, 8, 278, 4);
        run_single_starkprime_benchmark!(&mut group, 1, 15, 2048, 4317015616, 8, 65704, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 15, 2048, 4314825152, 8, 1628, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 15, 2048, 4294967296, 8, 256, 4);
        run_single_starkprime_benchmark!(&mut group, 1, 16, 512, 21195283396, 8, 145586, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 16, 512, 21161991096, 8, 2766, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 16, 512, 20851360000, 8, 380, 4);
        run_single_starkprime_benchmark!(&mut group, 1, 16, 1024, 14987635776, 8, 122424, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 16, 1024, 14959673344, 8, 2464, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 16, 1024, 14666178816, 8, 348, 4);
        run_single_starkprime_benchmark!(&mut group, 1, 16, 2048, 10597878916, 8, 102946, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 16, 2048, 10590025536, 8, 2196, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 16, 2048, 10485760000, 8, 320, 4);
        run_single_starkprime_benchmark!(&mut group, 1, 17, 512, 50614200576, 8, 224976, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 17, 512, 50570904392, 8, 3698, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 17, 512, 50479304976, 8, 474, 4);
        run_single_starkprime_benchmark!(&mut group, 1, 17, 1024, 35789072400, 8, 189180, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 17, 1024, 35741336184, 8, 3294, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 17, 1024, 35477982736, 8, 434, 4);
        run_single_starkprime_benchmark!(&mut group, 1, 17, 2048, 25307082724, 8, 159082, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 17, 2048, 25256916504, 8, 2934, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 17, 2048, 25091827216, 8, 398, 4);
        run_single_starkprime_benchmark!(&mut group, 1, 18, 512, 117850770436, 7, 343294, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 18, 512, 117793118808, 7, 4902, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 18, 512, 116319195136, 7, 584, 4);
        run_single_starkprime_benchmark!(&mut group, 1, 18, 1024, 83332678276, 7, 288674, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 18, 1024, 83224499896, 7, 4366, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 18, 1024, 82538991616, 7, 536, 4);
        run_single_starkprime_benchmark!(&mut group, 1, 18, 2048, 58925620516, 7, 242746, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 18, 2048, 58863869000, 7, 3890, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 18, 2048, 58594980096, 7, 492, 4);
        run_single_starkprime_benchmark!(&mut group, 1, 19, 512, 268120982416, 7, 517804, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 19, 512, 268086587392, 7, 6448, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 19, 512, 265764994576, 7, 718, 4);
        run_single_starkprime_benchmark!(&mut group, 1, 19, 1024, 189590576400, 7, 435420, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 19, 1024, 189514870784, 7, 5744, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 19, 1024, 187457825296, 7, 658, 4);
        run_single_starkprime_benchmark!(&mut group, 1, 19, 2048, 134059964164, 7, 366142, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 19, 2048, 134060503032, 7, 5118, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 19, 2048, 133090713856, 7, 604, 4);
        run_single_starkprime_benchmark!(&mut group, 1, 20, 512, 597108561984, 7, 772728, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 20, 512, 596947688000, 7, 8420, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 20, 512, 594262141456, 7, 878, 4);
        run_single_starkprime_benchmark!(&mut group, 1, 20, 1024, 422219246656, 7, 649784, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 20, 1024, 422212590008, 7, 7502, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 20, 1024, 422026932496, 7, 806, 4);
        run_single_starkprime_benchmark!(&mut group, 1, 20, 2048, 298552960000, 7, 546400, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 20, 2048, 298345446568, 7, 6682, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 20, 2048, 296637086736, 7, 738, 4);
        run_single_starkprime_benchmark!(&mut group, 1, 21, 512, 1303720941636, 7, 1141806, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 21, 512, 1303602169024, 7, 10924, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 21, 512, 1301023109376, 7, 1068, 4);
        run_single_starkprime_benchmark!(&mut group, 1, 21, 1024, 921868819600, 7, 960140, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 21, 1024, 921735471168, 7, 9732, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 21, 1024, 914861642256, 7, 978, 4);
        run_single_starkprime_benchmark!(&mut group, 1, 21, 2048, 651859234884, 7, 807378, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 21, 2048, 651714363000, 7, 8670, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 21, 2048, 650287411216, 7, 898, 4);
        run_single_starkprime_benchmark!(&mut group, 1, 22, 512, 2794687879824, 7, 1671732, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 22, 512, 2793688944704, 7, 14084, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 22, 512, 2786442301696, 7, 1292, 4);
        run_single_starkprime_benchmark!(&mut group, 1, 22, 1024, 1976138685504, 7, 1405752, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 22, 1024, 1975711510592, 7, 12548, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 22, 1024, 1965200244736, 7, 1184, 4);
        run_single_starkprime_benchmark!(&mut group, 1, 22, 2048, 1397341496464, 7, 1182092, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 22, 2048, 1396665211752, 7, 11178, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 22, 2048, 1390974924816, 7, 1086, 4);
        run_single_starkprime_benchmark!(&mut group, 1, 23, 512, 5888940837796, 6, 2426714, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 23, 512, 5888557851112, 6, 18058, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 23, 512, 5861899530496, 6, 1556, 4);
        run_single_starkprime_benchmark!(&mut group, 1, 23, 1024, 4164105496996, 6, 2040614, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 23, 1024, 4163956393472, 6, 16088, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 23, 1024, 4158271385856, 6, 1428, 4);
        run_single_starkprime_benchmark!(&mut group, 1, 23, 2048, 2944470674916, 7, 1715946, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 23, 2048, 2943882002368, 7, 14332, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 23, 2048, 2927055626496, 7, 1308, 4);
        run_single_starkprime_benchmark!(&mut group, 1, 24, 512, 12211739920900, 6, 3494530, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 24, 512, 12211490117952, 6, 23028, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 24, 512, 12176079851776, 6, 1868, 4);
        run_single_starkprime_benchmark!(&mut group, 1, 24, 1024, 8635005577444, 6, 2938538, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 24, 1024, 8632787556744, 6, 20514, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 24, 1024, 8630645337616, 6, 1714, 4);
        run_single_starkprime_benchmark!(&mut group, 1, 24, 2048, 6105870652036, 6, 2471006, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 24, 2048, 6104406528576, 6, 18276, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 24, 2048, 6075732010000, 6, 1570, 4);
        run_single_starkprime_benchmark!(&mut group, 1, 25, 512, 24945070206016, 6, 4994504, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 25, 512, 24943158948232, 6, 29218, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 25, 512, 24907645451536, 6, 2234, 4);
        run_single_starkprime_benchmark!(&mut group, 1, 25, 1024, 17638824019600, 6, 4199860, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 25, 1024, 17636910227000, 6, 26030, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 25, 1024, 17592186044416, 6, 16, 11);
        run_single_starkprime_benchmark!(&mut group, 1, 25, 2048, 12472537595904, 6, 3531648, 2);
        run_single_starkprime_benchmark!(&mut group, 1, 25, 2048, 12471027759000, 6, 23190, 3);
        run_single_starkprime_benchmark!(&mut group, 1, 25, 2048, 12438910749456, 6, 1878, 4);
    }
    // // Frog
    // {
    //     let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    //     let mut group = c.benchmark_group("Decomposition Frog");
    //     group.plot_config(plot_config.clone());

    //     // TODO: Update configurations
    //     run_single_frog_benchmark!(1, &mut group, 6, 1024, 10, 2);
    // }
}

criterion_group!(
    name=benches;
    config = Criterion::default().sample_size(10).measurement_time(Duration::from_secs(50)).warm_up_time(Duration::from_secs(1));
    targets = benchmarks_main);
criterion_main!(benches);
