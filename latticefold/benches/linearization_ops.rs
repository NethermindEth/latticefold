#![allow(incomplete_features)]
use ark_std::{test_rng, time::Duration, UniformRand};
use criterion::{
    criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion, PlotConfiguration,
};
use cyclotomic_rings::challenge_set::LatticefoldChallengeSet;
use cyclotomic_rings::rings::{GoldilocksChallengeSet, GoldilocksRingNTT, SuitableRing};
use latticefold::arith::{Instance, Witness, CCS};
use latticefold::decomposition_parameters::DecompositionParams;
use latticefold::nifs::error::LinearizationError;
use latticefold::nifs::linearization::utils::linearization_sumcheck_polynomial_comb_fn;
use latticefold::nifs::linearization::LFLinearizationProver;
use latticefold::transcript::poseidon::PoseidonTranscript;
use latticefold::utils::mle_helpers::calculate_Mz_mles;
use lattirust_poly::mle::DenseMultilinearExtension;
use std::fmt::Debug;
use std::sync::Arc;
use utils::{wit_and_ccs_gen, wit_and_ccs_gen_degree_three_non_scalar, wit_and_ccs_gen_non_scalar};

mod macros;
mod utils;

struct LinearizationSetup<R: SuitableRing> {
    ccs: CCS<R>,
    z_ccs: Vec<R>,
    wit: Witness<R>,
    g_mles: Vec<Arc<DenseMultilinearExtension<R>>>,
    g_degree: usize,
    mz_mles: Vec<DenseMultilinearExtension<R>>,
    point_r: Vec<R>,
}

fn linearization_setup<
    const C: usize,
    const W: usize,
    const WIT_LEN: usize,
    R: SuitableRing,
    CS: LatticefoldChallengeSet<R> + Clone,
    DP: DecompositionParams,
>(
    ccs: CCS<R>,
    z_ccs: Vec<R>,
    wit: Witness<R>,
) -> LinearizationSetup<R> {
    let (g_mles, g_degree, mz_mles) =
        LFLinearizationProver::<R, PoseidonTranscript<R, CS>>::construct_polynomial_g(
            &z_ccs,
            &mut PoseidonTranscript::<R, CS>::default(),
            &ccs,
        )
        .unwrap();

    let point_r = (0..ccs.s)
        .map(|_| R::rand(&mut test_rng()))
        .collect::<Vec<R>>();
    LinearizationSetup {
        ccs,
        z_ccs,
        wit,
        g_mles,
        g_degree,
        mz_mles,
        point_r,
    }
}

fn setup_test_environment<
    const C: usize,
    const WIT_LEN: usize,
    const W: usize,
    R: SuitableRing,
    CS: LatticefoldChallengeSet<R> + Clone,
    DP: DecompositionParams,
>() -> LinearizationSetup<R> {
    let r1cs_rows = 1 + WIT_LEN + 1;
    let (cm_i, wit, ccs, _) = wit_and_ccs_gen::<1, C, WIT_LEN, W, DP, R>(r1cs_rows);
    let z_ccs = cm_i.get_z_vector(&wit.w_ccs);

    linearization_setup::<C, W, WIT_LEN, R, CS, DP>(ccs, z_ccs, wit)
}

fn setup_non_scalar_test_environment<
    const C: usize,
    const WIT_LEN: usize,
    const W: usize,
    R: SuitableRing,
    CS: LatticefoldChallengeSet<R> + Clone,
    DP: DecompositionParams,
>() -> LinearizationSetup<R> {
    let r1cs_rows = 1 + WIT_LEN + 1;
    let (cm_i, wit, ccs, _) = wit_and_ccs_gen_non_scalar::<1, C, WIT_LEN, W, DP, R>(r1cs_rows);
    let z_ccs = cm_i.get_z_vector(&wit.w_ccs);

    linearization_setup::<C, W, WIT_LEN, R, CS, DP>(ccs, z_ccs, wit)
}

fn setup_degree_three_non_scalar_test_environment<
    const C: usize,
    const WIT_LEN: usize,
    const W: usize,
    R: SuitableRing,
    CS: LatticefoldChallengeSet<R> + Clone,
    DP: DecompositionParams,
>() -> LinearizationSetup<R> {
    let r1cs_rows = 1 + WIT_LEN + 1;
    let (cm_i, wit, ccs, _) =
        wit_and_ccs_gen_degree_three_non_scalar::<1, C, WIT_LEN, W, DP, R>(r1cs_rows);
    let z_ccs = cm_i.get_z_vector(&wit.w_ccs);

    linearization_setup::<C, W, WIT_LEN, R, CS, DP>(ccs, z_ccs, wit)
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
    setup: LinearizationSetup<R>,
) {
    group.bench_with_input(
        BenchmarkId::new(
            "Evaluate Mz_MLEs",
            format!(
                "Kappa={}, W_CCS={}, W={}, L={}, B={}, B_SMALL={}, K={}",
                C,
                WIT_LEN,
                W,
                DP::L,
                DP::B,
                DP::B_SMALL,
                DP::K
            ),
        ),
        &setup,
        |bench, setup| {
            bench.iter(|| {
                let _ = calculate_Mz_mles::<R, LinearizationError<R>>(&setup.ccs, &setup.z_ccs);
            })
        },
    );

    // Prepare the main linearization polynomial.
    group.bench_with_input(
        BenchmarkId::new(
            "Construct Sumcheck Poly",
            format!(
                "Kappa={}, W_CCS={}, W={}, L={}, B={}, B_SMALL={}, K={}",
                C,
                WIT_LEN,
                W,
                DP::L,
                DP::B,
                DP::B_SMALL,
                DP::K
            ),
        ),
        &setup,
        |bench, setup| {
            bench.iter(|| {
                // Construct the sumcheck polynomial g
                let _ =
                    LFLinearizationProver::<R, PoseidonTranscript<R, CS>>::construct_polynomial_g(
                        &setup.z_ccs,
                        &mut PoseidonTranscript::<R, CS>::default(),
                        &setup.ccs,
                    )
                    .unwrap();
            })
        },
    );

    // Prepare the main linearization polynomial.
    group.bench_with_input(
        BenchmarkId::new(
            "Sumcheck",
            format!(
                "Kappa={}, W_CCS={}, W={}, L={}, B={}, B_SMALL={}, K={}",
                C,
                WIT_LEN,
                W,
                DP::L,
                DP::B,
                DP::B_SMALL,
                DP::K
            ),
        ),
        &setup,
        |bench, setup| {
            bench.iter(|| {
                let comb_fn = |vals: &[R]| -> R {
                    linearization_sumcheck_polynomial_comb_fn(vals, &setup.ccs)
                };

                // Run sumcheck protocol.
                let _ =
                    LFLinearizationProver::<R, PoseidonTranscript<R, CS>>::generate_sumcheck_proof(
                        &mut PoseidonTranscript::<R, CS>::default(),
                        &setup.g_mles,
                        setup.ccs.s,
                        setup.g_degree,
                        comb_fn,
                    )
                    .expect("Failed to generate sumcheck proof");
            })
        },
    );

    group.bench_with_input(
        BenchmarkId::new("Evaluate U and V", format!("Kappa={}, W_CCS={}, W={}, L={}, B={}, B_SMALL={}, K={}", C, WIT_LEN, W, DP::L, DP::B, DP::B_SMALL, DP::K)),
        &setup,
        |bench, setup| {
            bench.iter(|| {
                let _ = LFLinearizationProver::<R, PoseidonTranscript<R, CS>>::compute_evaluation_vectors(
                        &setup.wit,
                        &setup.point_r,
                        &setup.mz_mles,
                    ).expect("Failed to compute evaluation vectors");
            })
        },
    );
}

macro_rules! run_single_linearization_goldilocks_benchmark {
    ($crit_group:expr, $kappa:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            let setup = setup_test_environment::< $kappa, $w, {$w * $l}, GoldilocksRingNTT, GoldilocksChallengeSet, [<DecompParamsWithB $b W $w b $b_small K $k>]>();
            linearization_operations::<$kappa, {$w * $l}, $w, GoldilocksRingNTT, GoldilocksChallengeSet, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit_group, setup);
        }
    };
}

macro_rules! run_single_linearization_non_scalar_goldilocks_benchmark {
    ($crit_group:expr, $kappa:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            let setup = setup_non_scalar_test_environment::< $kappa, $w, {$w * $l}, GoldilocksRingNTT, GoldilocksChallengeSet, [<DecompParamsWithB $b W $w b $b_small K $k>]>();
            linearization_operations::<$kappa, {$w * $l}, $w, GoldilocksRingNTT, GoldilocksChallengeSet, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit_group, setup);
        }
    };
}

macro_rules! run_single_linearization_degree_three_non_scalar_goldilocks_benchmark {
    ($crit_group:expr, $kappa:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            let setup = setup_degree_three_non_scalar_test_environment::< $kappa, $w, {$w * $l}, GoldilocksRingNTT, GoldilocksChallengeSet, [<DecompParamsWithB $b W $w b $b_small K $k>]>();
            linearization_operations::<$kappa, {$w * $l}, $w, GoldilocksRingNTT, GoldilocksChallengeSet, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit_group, setup);
        }
    };
}

fn single_operation_benchmarks(c: &mut Criterion) {
    // Goldilocks
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Single Linearization Operations Goldilocks");
        group.plot_config(plot_config.clone());
        run_goldilocks_linearization_benchmarks!(group);
    }

    // Goldilocks non-scalar
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Single Linearization Operations Goldilocks Non-Scalar");
        group.plot_config(plot_config.clone());
        run_goldilocks_linearization_non_scalar_benchmarks!(group);
    }

    // Goldilocks degree three non-scalar
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group =
            c.benchmark_group("Single Linearization Operations Goldilocks Degree Three Non-Scalar");
        group.plot_config(plot_config.clone());
        run_goldilocks_linearization_degree_three_non_scalar_benchmarks!(group);
    }
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
