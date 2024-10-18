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
    let ccs = get_test_dummy_ccs::<R, IO, WIT_LEN>(r1cs_rows);
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

    let (_, _, decomposition_proof) =
        LFDecompositionProver::<_, PoseidonTranscript<R, CS>>::prove::<W, C, P>(
            &lcccs,
            wit,
            &mut prover_transcript,
            ccs,
            scheme,
        )
        .unwrap();

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
    let (cm_i, wit, ccs, scheme) = wit_and_ccs_gen::<IO, C, WIT_LEN, W, P, R>(r1cs_rows);
    // N/Q = prime / degree

    prover_decomposition_benchmark::<C, W, P, R, CS>(group, &cm_i, &wit, &ccs, &scheme);

    verifier_decomposition_benchmark::<C, W, P, R, CS>(group, &cm_i, &wit, &ccs, &scheme);
}

// Macros
macro_rules! define_starkprime_params {
    ($w:expr, $b:expr, $l:expr) => {
        paste::paste! {
            #[derive(Clone)]
            struct [<StarkPrimeParamsWithB $b W $w>];

            impl DecompositionParams for [<StarkPrimeParamsWithB $b W $w>] {
                const B: u128 = $b;
                const L: usize = $l;
                const B_SMALL: usize = 2; // This is not use in decompose or linearization
                const K: usize = 28;// This is not use in decompose or linearization
            }
        }
    };
}

macro_rules! run_single_starkprime_benchmark {
    ($io:expr, $crit:expr, $cw:expr, $w:expr, $b:expr, $l:expr) => {
        define_starkprime_params!($w, $b, $l);
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
    // Babybear
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Decomposition BabyBear");
        group.plot_config(plot_config.clone());

        // TODO: Update configurations
        run_single_babybear_benchmark!(1, &mut group, 6, 1024, 10, 2);
    }

    // Godlilocks
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Decomposition Godlilocks");
        group.plot_config(plot_config.clone());

        // TODO: Update configurations
        run_single_goldilocks_benchmark!(1, &mut group, 6, 1024, 10, 2);
    }

    // StarkPrime
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Decomposition StarkPrime");
        group.plot_config(plot_config.clone());

        // TODO: Update configurations
        run_single_starkprime_benchmark!(1, &mut group, 6, 1024, 10, 2);
    }

    // Frog
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Decomposition Frog");
        group.plot_config(plot_config.clone());

        // TODO: Update configurations
        run_single_frog_benchmark!(1, &mut group, 6, 1024, 10, 2);
    }

    // Dilithium
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Decomposition Dilithium");
        group.plot_config(plot_config.clone());

        // TODO: Update configurations
        run_single_dilithium_benchmark!(1, &mut group, 6, 1024, 10, 2);
    }
}

criterion_group!(
    name=benches;
    config = Criterion::default().sample_size(10).measurement_time(Duration::from_secs(50)).warm_up_time(Duration::from_secs(1));
    targets = benchmarks_main);
criterion_main!(benches);
