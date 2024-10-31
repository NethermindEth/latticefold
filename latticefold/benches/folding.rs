#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
use criterion::{
    criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion, PlotConfiguration,
};
use cyclotomic_rings::{challenge_set::LatticefoldChallengeSet, SuitableRing};
use cyclotomic_rings::{StarkChallengeSet, StarkRingNTT};
use latticefold::nifs::decomposition::{
    DecompositionProver, DecompositionVerifier, LFDecompositionProver, LFDecompositionVerifier,
};
use latticefold::nifs::folding::{
    FoldingProver, FoldingVerifier, LFFoldingProver, LFFoldingVerifier,
};
use rand::thread_rng;
use std::{fmt::Debug, time::Duration};
mod utils;
use ark_std::UniformRand;
use utils::{get_test_dummy_ccs, get_test_dummy_z_split};

use latticefold::{
    arith::{Arith, Witness, CCCS, CCS},
    commitment::AjtaiCommitmentScheme,
    decomposition_parameters::DecompositionParams,
    nifs::linearization::{
        LFLinearizationProver, LFLinearizationVerifier, LinearizationProver, LinearizationVerifier,
    },
    transcript::poseidon::PoseidonTranscript,
};

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
    let ccs: CCS<R> = get_test_dummy_ccs::<R, X_LEN, WIT_LEN, W>(r1cs_rows);
    let (one, x_ccs, w_ccs) = get_test_dummy_z_split::<R, X_LEN, WIT_LEN>();
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
    println!("prove folding");
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
    let (_, wit_vec, decomposition_proof) =
        LFDecompositionProver::<_, PoseidonTranscript<R, CS>>::prove::<W, C, P>(
            &lcccs,
            wit,
            &mut prover_transcript,
            ccs,
            scheme,
        )
        .unwrap();

    println!("verify decomposition");
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
    println!("verify folding");
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
    let (_, wit_vec, decomposition_proof) =
        LFDecompositionProver::<_, PoseidonTranscript<R, CS>>::prove::<W, C, P>(
            &lcccs,
            wit,
            &mut prover_transcript,
            ccs,
            scheme,
        )
        .unwrap();

    println!("verify decomposition");
    let lcccs_vec = LFDecompositionVerifier::<_, PoseidonTranscript<R, CS>>::verify::<C, P>(
        &lcccs,
        &decomposition_proof,
        &mut verifier_transcript,
        ccs,
    )
    .unwrap();
    println!("prove folding");
    let (_, _, folding_proof) = LFFoldingProver::<R, PoseidonTranscript<R, CS>>::prove::<C, P>(
        &lcccs_vec,
        &wit_vec,
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
        &(lcccs_vec, folding_proof, ccs),
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

fn linearization_benchmarks<
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
    let r1cs_rows = 5; // This makes a square matrix but is too much memory
    println!("Witness generation");
    let (cm_i, wit, ccs, scheme) = wit_and_ccs_gen::<X_LEN, C, WIT_LEN, W, P, R>(r1cs_rows);

    prover_folding_benchmark::<C, W, P, R, CS>(group, &cm_i, &wit, &ccs, &scheme);

    verifier_folding_benchmark::<C, W, P, R, CS>(group, &cm_i, &wit, &ccs, &scheme);
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
#[macro_export]
macro_rules! run_single_starkprime_benchmark {
    ($crit:expr, $io:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
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

        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        // Using this parameters doesn't trigger OOM killer
        run_single_starkprime_benchmark!(&mut group, 1, 15, 512, 3010936384, 8, 38, 6);
    }
}

criterion_group!(
    name=benches;
    config = Criterion::default().sample_size(10).measurement_time(Duration::from_secs(50)).warm_up_time(Duration::from_secs(1));
    targets = benchmarks_main);
criterion_main!(benches);
