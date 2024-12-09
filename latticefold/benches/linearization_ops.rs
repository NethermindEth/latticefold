#![allow(incomplete_features)]
use ark_std::{test_rng, time::Duration, UniformRand};
use criterion::{
    criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion, PlotConfiguration,
};
use cyclotomic_rings::challenge_set::LatticefoldChallengeSet;
use cyclotomic_rings::rings::{GoldilocksChallengeSet, GoldilocksRingNTT, SuitableRing};
use latticefold::arith::{Instance, Witness, CCCS, CCS};
use latticefold::commitment::AjtaiCommitmentScheme;
use latticefold::decomposition_parameters::DecompositionParams;
use latticefold::nifs::error::LinearizationError;
use latticefold::nifs::linearization::LFLinearizationProver;
use latticefold::transcript::poseidon::PoseidonTranscript;
use latticefold::utils::mle_helpers::calculate_Mz_mles;
use std::fmt::Debug;
use utils::wit_and_ccs_gen;

mod macros;
mod utils;

fn setup_test_environment<
    RqNTT: SuitableRing,
    DP: DecompositionParams,
    const C: usize,
    const W: usize,
    const WIT_LEN: usize,
>() -> (
    Witness<RqNTT>,
    CCCS<C, RqNTT>,
    CCS<RqNTT>,
    AjtaiCommitmentScheme<C, W, RqNTT>,
) {
    let mut rng = test_rng();
    let scheme = AjtaiCommitmentScheme::rand(&mut rng);
    let r1cs_rows = 1 + WIT_LEN + 1;
    let (cm_i, wit, ccs, _) = wit_and_ccs_gen::<1, C, WIT_LEN, W, DP, RqNTT>(r1cs_rows);

    (wit, cm_i, ccs, scheme)
}

fn linearization_operations<
    const C: usize,
    const W: usize,
    const WIT_LEN: usize,
    R: Clone + UniformRand + Debug + SuitableRing,
    CS: LatticefoldChallengeSet<R> + Clone,
    DP: DecompositionParams,
>(
    group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>,
) {
    let mut rng = test_rng();
    let (wit, cm_i, ccs, _) = setup_test_environment::<R, DP, C, W, WIT_LEN>();
    let z_ccs = cm_i.get_z_vector(&wit.w_ccs);

    // Prepare MLEs
    let mz_mles = calculate_Mz_mles::<R, LinearizationError<R>>(&ccs, &z_ccs).unwrap();
    let beta_s = (0..ccs.s).map(|_| R::rand(&mut rng)).collect::<Vec<R>>();

    // MZ mles
    group.bench_with_input(
        BenchmarkId::new(
            "Evaluate Mz_MLEs",
            format!(
                "Kappa = {}, W_CCS= {}, W= {}, L= {}, B= {}, B_SMALL= {}, K= {}",
                C,
                WIT_LEN,
                W,
                DP::L,
                DP::B,
                DP::B_SMALL,
                DP::K
            ),
        ),
        &(ccs.clone(), z_ccs.clone()),
        |bench, (ccs, z_ccs)| {
            bench.iter(|| {
                let _ = calculate_Mz_mles::<R, LinearizationError<R>>(ccs, z_ccs);
            })
        },
    );

    // Prepare the main linearization polynomial.
    group.bench_with_input(
        BenchmarkId::new(
            "Construct Sumcheck Poly",
            format!(
                "Kappa = {}, W_CCS= {}, W= {}, L= {}, B= {}, B_SMALL= {}, K= {}",
                C,
                WIT_LEN,
                W,
                DP::L,
                DP::B,
                DP::B_SMALL,
                DP::K
            ),
        ),
        &(ccs.clone(), mz_mles, beta_s),
        |bench, (ccs, mz_mles, beta_s)| {
            bench.iter(|| {
                // Construct the sumcheck polynomial g
                let _ =
                    LFLinearizationProver::<R, PoseidonTranscript<R, CS>>::construct_polynomial_g(
                        &z_ccs,
                        &mut PoseidonTranscript::<R, CS>::default(),
                        ccs,
                    )
                    .unwrap();
            })
        },
    );

    let point_r = (0..ccs.s).map(|_| R::rand(&mut rng)).collect::<Vec<R>>();
    let Mz_mles = calculate_Mz_mles::<R, LinearizationError<R>>(&ccs, &z_ccs).unwrap();
    group.bench_with_input(
        BenchmarkId::new("Evaluate U and V", format!("Kappa = {}, W_CCS= {}, W= {}, L= {}, B= {}, B_SMALL= {}, K= {}", C, WIT_LEN, W, DP::L, DP::B, DP::B_SMALL, DP::K)),
        &(wit.clone(), point_r.clone(), Mz_mles.clone()),
        |bench, (wit, point_r, Mz_mles)| {
            bench.iter(|| {
                let _ = LFLinearizationProver::<R, PoseidonTranscript<R, CS>>::compute_evaluation_vectors(wit, point_r, &Mz_mles).expect("Failed to compute evaluation vectors");
            })
        },
    );
}

macro_rules! run_single_linearization_goldilocks_benchmark {
    ($crit_group:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            linearization_operations::<$cw, {$w * $l}, $w, GoldilocksRingNTT, GoldilocksChallengeSet, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit_group);
        }
    };
}

fn single_operation_benchmarks(c: &mut Criterion) {
    let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    let mut group = c.benchmark_group("Single Linearization Operations Goldilocks");
    group.plot_config(plot_config.clone());

    // Linearization
    // Please note that C is not used until decomposition.
    // The only parameter that we are interested on varying for linearization is W (as it already includes WIT_LEN and DP::L)
    // We explore parameters in the range  W = 2^9-2^16
    run_goldilocks_linearization_benchmarks!(group);
}

pub fn benchmarks_main(c: &mut Criterion) {
    single_operation_benchmarks(c);
}

criterion_group!(
    name=benches;
    config = Criterion::default()
            .sample_size(10)
            .measurement_time(Duration::from_secs(50))
            .warm_up_time(Duration::from_secs(1));
    targets = benchmarks_main
);
criterion_main!(benches);
