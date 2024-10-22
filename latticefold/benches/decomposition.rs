#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
use criterion::{
    criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion, PlotConfiguration,
};
use cyclotomic_rings::{challenge_set::LatticefoldChallengeSet, SuitableRing};
use rand::thread_rng;
use std::fmt::Debug;
mod utils;
use ark_std::UniformRand;
use cyclotomic_rings::challenge_set::BinarySmallSet;
use cyclotomic_rings::{
    BabyBearChallengeSet, BabyBearRingNTT, FrogChallengeSet, FrogRingNTT, GoldilocksChallengeSet,
    GoldilocksRingNTT, StarkChallengeSet, StarkRingNTT,
};
use latticefold::{
    arith::{Arith, Witness, CCCS, CCS},
    commitment::AjtaiCommitmentScheme,
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
    parameters::{DecompositionParams, DILITHIUM_PRIME},
    transcript::poseidon::PoseidonTranscript,
};
use lattirust_ring::cyclotomic_ring::models::pow2_debug::Pow2CyclotomicPolyRingNTT;
use paste;
use std::time::Duration;
use utils::{get_test_dummy_ccs, get_test_dummy_z_split};

fn wit_and_ccs_gen<
    const IO: usize,
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
    let ccs = get_test_dummy_ccs::<R, IO, WIT_LEN, W>(r1cs_rows);
    let (one, x_ccs, w_ccs) = get_test_dummy_z_split::<R, IO, WIT_LEN>();
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
) {
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
            format!("Param. B={}, L={}", P::B, P::L),
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
) {
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
            format!("Param. B={}, L={}", P::B, P::L),
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
    const IO: usize,
    const C: usize,
    const WIT_LEN: usize,
    const W: usize,
    CS: LatticefoldChallengeSet<R>,
    R: Clone + UniformRand + Debug + SuitableRing + for<'a> std::ops::AddAssign<&'a R>,
    P: DecompositionParams + Clone,
>(
    group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>,
) {
    let r1cs_rows = 5;
    println!("Witness generation");
    let (cm_i, wit, ccs, scheme) = wit_and_ccs_gen::<IO, C, WIT_LEN, W, P, R>(r1cs_rows);
    // N/Q = prime / degree

    prover_decomposition_benchmark::<C, W, P, R, CS>(group, &cm_i, &wit, &ccs, &scheme);

    verifier_decomposition_benchmark::<C, W, P, R, CS>(group, &cm_i, &wit, &ccs, &scheme);
}

// Macros
macro_rules! define_starkprime_params {
    ($w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        paste::paste! {
            #[derive(Clone)]
            struct [<StarkPrimeParamsWithB $b W $w>];

            impl DecompositionParams for [<StarkPrimeParamsWithB $b W $w>] {
                const B: u128 = $b;
                const L: usize = $l;
                const B_SMALL: usize = $b_small;
                const K: usize = $k;
            }
        }
    };
}

macro_rules! run_single_starkprime_benchmark {
    ($io:expr, $crit:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_starkprime_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            decomposition_benchmarks::<$io, $cw, $w,{$w * $l}, StarkChallengeSet, StarkRingNTT, [<StarkPrimeParamsWithB $b W $w>]>($crit);
        }
    };
}

macro_rules! define_goldilocks_params {
    ($w:expr, $b:expr, $l:expr) => {
        paste::paste! {
            #[derive(Clone)]
            struct [<GoldilocksParamsWithB $b W $w>];

            impl DecompositionParams for [<GoldilocksParamsWithB $b W $w>] {
                const B: u128 = $b;
                const L: usize = $l;
                const B_SMALL: usize = 2; // This is not use in decompose or linearization
                const K: usize = 28;// This is not use in decompose or linearization
            }
        }
    };
}

macro_rules! run_single_goldilocks_benchmark {
    ($io:expr, $crit:expr, $cw:expr, $w:expr, $b:expr, $l:expr) => {
        define_goldilocks_params!($w, $b, $l);
        paste::paste! {
            decomposition_benchmarks::<$io, $cw, $w, {$w * $l}, GoldilocksChallengeSet, GoldilocksRingNTT, [<GoldilocksParamsWithB $b W $w>]>($crit);

        }
    };
}

macro_rules! define_babybear_params {
    ($w:expr, $b:expr, $l:expr) => {
        paste::paste! {
            #[derive(Clone)]
            struct [<BabyBearParamsWithB $b W $w>];

            impl DecompositionParams for [<BabyBearParamsWithB $b W $w>] {
                const B: u128 = $b;
                const L: usize = $l;
                const B_SMALL: usize = 2; // This is not use in decompose or linearization
                const K: usize = 28;// This is not use in decompose or linearization
            }
        }
    };
}

macro_rules! run_single_babybear_benchmark {
    ($io:expr, $crit:expr, $cw:expr, $w:expr, $b:expr, $l:expr) => {
        define_babybear_params!($w, $b, $l);
        paste::paste! {
            decomposition_benchmarks::<$io, $cw, $w, {$w * $l}, BabyBearChallengeSet, BabyBearRingNTT, [<BabyBearParamsWithB $b W $w>]>($crit);

        }
    };
}

macro_rules! define_frog_params {
    ($w:expr, $b:expr, $l:expr) => {
        paste::paste! {
            #[derive(Clone)]
            struct [<FrogParamsWithB $b W $w>];

            impl DecompositionParams for [<FrogParamsWithB $b W $w>] {
                const B: u128 = $b;
                const L: usize = $l;
                const B_SMALL: usize = 2; // This is not use in decompose or commit
                const K: usize = 28;// This is not use in decompose or commit
            }
        }
    };
}

macro_rules! run_single_frog_benchmark {
    ($io:expr, $crit:expr, $cw:expr, $w:expr, $b:expr, $l:expr) => {
        define_frog_params!($w, $b, $l);
        paste::paste! {
            decomposition_benchmarks::<$io, $cw, $w, {$w * $l}, FrogChallengeSet, FrogRingNTT, [<FrogParamsWithB $b W $w>]>($crit);

        }
    };
}

macro_rules! define_dilithium_params {
    ($w:expr, $b:expr, $l:expr) => {
        paste::paste! {
            #[derive(Clone)]
            struct [<DilithiumParamsWithB $b W $w>];

            impl DecompositionParams for [<DilithiumParamsWithB $b W $w>] {
                const B: u128 = $b;
                const L: usize = $l;
                const B_SMALL: usize = 2; // This is not use in decompose or commit
                const K: usize = 28;// This is not use in decompose or commit
            }
        }
    };
}

macro_rules! run_single_dilithium_benchmark {
    ($io:expr, $crit:expr, $cw:expr, $w:expr, $b:expr, $l:expr) => {
        define_dilithium_params!($w, $b, $l);
        paste::paste! {
            decomposition_benchmarks::<$io, $cw, $w, {$w * $l}, BinarySmallSet<DILITHIUM_PRIME, 256>, Pow2CyclotomicPolyRingNTT<DILITHIUM_PRIME, 256>, [<DilithiumParamsWithB $b W $w>]>($crit);
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

    // // Godlilocks
    // {
    //     let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    //     let mut group = c.benchmark_group("Decomposition Godlilocks");
    //     group.plot_config(plot_config.clone());

    //     // TODO: Update configurations
    //     run_single_goldilocks_benchmark!(1, &mut group, 6, 1024, 10, 2);
    // }

    // StarkPrime
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Decomposition StarkPrime");
        group.plot_config(plot_config.clone());

        // TODO: Update configurations
        // Needs to update odd numbers of B
        // Kappa values for B ≈ 2^16 (within a margin of 65536):
        run_single_starkprime_benchmark!(1, &mut group, 6, 32768, 45796, 17, 214, 2);
        // run_single_starkprime_benchmark!(1, &mut group, 6, 65536, 32400, 17, 180, 2);
        // run_single_starkprime_benchmark!(1, &mut group, 7, 131072, 91809, 16, 303, 2);
        // run_single_starkprime_benchmark!(1, &mut group, 7, 262144, 64516, 16, 254, 2);
        // run_single_starkprime_benchmark!(1, &mut group, 7, 524288, 45796, 17, 214, 2);
        // run_single_starkprime_benchmark!(1, &mut group, 7, 1048576, 32400, 17, 180, 2);

        // // Calculating largest B for max_kappa where L is an integer for all num_cols:
        // run_single_starkprime_benchmark!(1, &mut group, 11, 32768, 7086244, 11, 2662, 2);
        // run_single_starkprime_benchmark!(1, &mut group, 12, 65536, 5013121, 12, 2239, 2);
        // run_single_starkprime_benchmark!(1, &mut group, 12, 131072, 3545689, 12, 1883, 2);
        // run_single_starkprime_benchmark!(1, &mut group, 12, 262144, 2505889, 12, 1583, 2);
        // run_single_starkprime_benchmark!(1, &mut group, 13, 524288, 1771561, 13, 11, 6);
        // run_single_starkprime_benchmark!(1, &mut group, 13, 1048576, 1252161, 13, 1119, 2);
    }

    // // Frog
    // {
    //     let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    //     let mut group = c.benchmark_group("Decomposition Frog");
    //     group.plot_config(plot_config.clone());

    //     // TODO: Update configurations
    //     run_single_frog_benchmark!(1, &mut group, 6, 1024, 10, 2);
    // }

    // // Dilithium
    // {
    //     let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    //     let mut group = c.benchmark_group("Decomposition Dilithium");
    //     group.plot_config(plot_config.clone());

    //     // TODO: Update configurations
    //     run_single_dilithium_benchmark!(1, &mut group, 6, 1024, 10, 2);
    // }
}

criterion_group!(
    name=benches;
    config = Criterion::default().sample_size(10).measurement_time(Duration::from_secs(50)).warm_up_time(Duration::from_secs(1));
    targets = benchmarks_main);
criterion_main!(benches);
