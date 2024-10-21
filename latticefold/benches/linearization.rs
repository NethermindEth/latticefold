#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
use criterion::{
    criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion, PlotConfiguration,
};
use cyclotomic_rings::challenge_set::BinarySmallSet;
use cyclotomic_rings::{challenge_set::LatticefoldChallengeSet, SuitableRing};
use cyclotomic_rings::{
    BabyBearChallengeSet, BabyBearRingNTT, FrogChallengeSet, FrogRingNTT, GoldilocksChallengeSet,
    GoldilocksRingNTT, StarkChallengeSet, StarkRingNTT,
};
use lattirust_ring::cyclotomic_ring::models::pow2_debug::Pow2CyclotomicPolyRingNTT;
use rand::thread_rng;
use std::fmt::Debug;
mod utils;
use ark_std::UniformRand;
use std::time::Duration;
use utils::{get_test_dummy_ccs, get_test_dummy_z_split};

use latticefold::{
    arith::{Arith, Witness, CCCS, CCS, LCCCS},
    commitment::AjtaiCommitmentScheme,
    nifs::linearization::{
        LFLinearizationProver, LFLinearizationVerifier, LinearizationProof, LinearizationProver,
        LinearizationVerifier,
    },
    parameters::{DecompositionParamData, DecompositionParams},
    transcript::poseidon::PoseidonTranscript,
};

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
    let ccs: CCS<R> = get_test_dummy_ccs::<R, IO, WIT_LEN>(r1cs_rows);
    let (one, x_ccs, w_ccs) = get_test_dummy_z_split::<R, IO, WIT_LEN>();
    let mut z = vec![one];
    z.extend(&x_ccs);
    z.extend(&w_ccs);
    match ccs.check_relation(&z) {
        Ok(_) => println!("R1CS valid!"),
        Err(e) => println!("R1CS invalid: {:?}", e),
    }

    let scheme: AjtaiCommitmentScheme<C, W, R> = AjtaiCommitmentScheme::rand(&mut thread_rng());

    let wit: Witness<R> = Witness::from_w_ccs::<P>(&w_ccs);
    let cm_i: CCCS<C, R> = CCCS {
        cm: wit.commit::<C, W, P>(&scheme).unwrap(),
        x_ccs,
    };

    (cm_i, wit, ccs, scheme)
}

fn prover_linearization_benchmark<
    const C: usize,
    const W: usize,
    P: DecompositionParams,
    R: SuitableRing,
    CS: LatticefoldChallengeSet<R>,
>(
    c: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>,
    cm_i: &CCCS<C, R>,
    wit: &Witness<R>,
    ccs: &CCS<R>,
) -> (LCCCS<C, R>, LinearizationProof<R>) {
    println!("Proving linearization");
    println!("transcript");
    let mut transcript = PoseidonTranscript::<R, CS>::default();
    let res = LFLinearizationProver::<_, PoseidonTranscript<R, CS>>::prove(
        &cm_i,
        &wit,
        &mut transcript,
        &ccs,
    );
    match res {
        Ok(_) => println!("Linearization proof generated with success"),
        Err(ref e) => println!("Linearization error: {:?}", e),
    }
    c.bench_with_input(
        BenchmarkId::new(
            "Linearization Prover",
            format!("Param. B={}, L={}", P::B, P::L),
        ),
        &(cm_i, wit, ccs),
        |b, (cm_i, wit, ccs)| {
            b.iter(|| {
                let _ = LFLinearizationProver::<_, PoseidonTranscript<R, CS>>::prove(
                    cm_i,
                    wit,
                    &mut transcript,
                    ccs,
                );
            })
        },
    );
    res.unwrap()
}

fn verifier_linearization_benchmark<
    const C: usize,
    const W: usize,
    P: DecompositionParams,
    R: SuitableRing,
    CS: LatticefoldChallengeSet<R>,
>(
    c: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>,
    cm_i: &CCCS<C, R>,
    ccs: &CCS<R>,
    proof: (LCCCS<C, R>, LinearizationProof<R>),
) {
    println!("Verifying linearization");
    c.bench_with_input(
        BenchmarkId::new(
            "Linearization Verifier",
            format!("Param. B={}, L={}", P::B, P::L),
        ),
        &(cm_i, proof.1, ccs),
        |b, (cm_i, proof, ccs)| {
            b.iter(|| {
                let mut transcript = PoseidonTranscript::<R, CS>::default();
                let _ = LFLinearizationVerifier::<_, PoseidonTranscript<R, CS>>::verify(
                    &cm_i,
                    &proof,
                    &mut transcript,
                    &ccs,
                );
            })
        },
    );
}

fn linearization_benchmarks<
    const IO: usize,
    const C: usize,
    const WIT_LEN: usize,
    const W: usize,
    CS: LatticefoldChallengeSet<R>,
    R: SuitableRing,
    P: DecompositionParams,
>(
    group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>,
) {
    let r1cs_rows = WIT_LEN + IO + 1; // This makes a square matrix but is too much memory
    println!("Witness generation");
    let (cm_i, wit, ccs, _) = wit_and_ccs_gen::<IO, C, WIT_LEN, W, P, R>(r1cs_rows);

    let proof = prover_linearization_benchmark::<C, W, P, R, CS>(group, &cm_i, &wit, &ccs);

    verifier_linearization_benchmark::<C, W, P, R, CS>(group, &cm_i, &ccs, proof);
}

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
            linearization_benchmarks::<$io, $cw, $w,{$w * $l}, StarkChallengeSet, StarkRingNTT, [<StarkPrimeParamsWithB $b W $w>]>($crit);
        }
    };
}

fn benchmarks_main(c: &mut Criterion) {
    // StarkPrime
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Decomposition StarkPrime");
        group.plot_config(plot_config.clone());

        // TODO: Update configurations
        // Kappa values for B ≈ 2^16 (within a margin of 65536):
        run_single_starkprime_benchmark!(1, &mut group, 6, 32768, 45914, 17, 214, 2);
        run_single_starkprime_benchmark!(1, &mut group, 6, 65536, 32466, 17, 180, 2);
        run_single_starkprime_benchmark!(1, &mut group, 7, 131072, 91958, 16, 303, 2);
        run_single_starkprime_benchmark!(1, &mut group, 7, 262144, 65024, 16, 254, 2);
        run_single_starkprime_benchmark!(1, &mut group, 7, 524288, 45978, 17, 214, 2);
        run_single_starkprime_benchmark!(1, &mut group, 7, 1048576, 32512, 17, 180, 2);

        // Calculating largest B for max_kappa where L is an integer for all num_cols:
        run_single_starkprime_benchmark!(1, &mut group, 11, 32768, 7091446, 11, 2662, 2);
        run_single_starkprime_benchmark!(1, &mut group, 12, 65536, 5014410, 12, 2239, 2);
        run_single_starkprime_benchmark!(1, &mut group, 12, 131072, 3545724, 12, 1883, 2);
        run_single_starkprime_benchmark!(1, &mut group, 12, 262144, 2507206, 12, 1583, 2);
        run_single_starkprime_benchmark!(1, &mut group, 13, 524288, 1772862, 13, 11, 6);
        run_single_starkprime_benchmark!(1, &mut group, 13, 1048576, 1253602, 13, 1119, 2);
    }
}

criterion_group!(
    name=benches;
    config = Criterion::default().sample_size(10).measurement_time(Duration::from_secs(50)).warm_up_time(Duration::from_secs(1));
    targets = benchmarks_main);
criterion_main!(benches);
