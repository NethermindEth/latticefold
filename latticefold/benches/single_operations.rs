#![allow(incomplete_features)]
use ark_std::{test_rng, time::Duration, UniformRand};
use criterion::{
    criterion_group, criterion_main, AxisScale, Criterion,
    PlotConfiguration,
};
use cyclotomic_rings::rings::{GoldilocksRingNTT, SuitableRing};

use std::fmt::Debug;
use latticefold::arith::{Instance, Witness, CCCS, CCS};
use latticefold::arith::r1cs::{get_test_r1cs, get_test_z_split};
use latticefold::commitment::AjtaiCommitmentScheme;
use latticefold::decomposition_parameters::DecompositionParams;
use latticefold::nifs::error::LinearizationError;
use latticefold::utils::mle_helpers::{calculate_Mz_mles};
use ark_std::rand::Rng;
use lattirust_ring::Ring;
use latticefold::nifs::linearization::utils::prepare_lin_sumcheck_polynomial;


#[derive(Clone)]
pub struct GoldilocksDP;
impl DecompositionParams for GoldilocksDP {
    const B: u128 = 1 << 15;
    const L: usize = 5;
    const B_SMALL: usize = 2;
    const K: usize = 15;
}
pub fn get_test_ccs<R: Ring>(w: usize, l: usize) -> CCS<R> {
    let r1cs = get_test_r1cs::<R>();
    CCS::<R>::from_r1cs_padded(r1cs, w, l)
}

fn setup_test_environment<
    RqNTT: SuitableRing,
    DP: DecompositionParams,
    const C: usize,
    const W: usize
>(
    input: Option<usize>,
) -> (
    Witness<RqNTT>,
    CCCS<C, RqNTT>,
    CCS<RqNTT>,
    AjtaiCommitmentScheme<C, W, RqNTT>,
) {
    let ccs = get_test_ccs::<RqNTT>(W, DP::L);
    let mut rng = test_rng();
    let (_, x_ccs, w_ccs) = get_test_z_split::<RqNTT>(input.unwrap_or(rng.gen_range(0..64)));
    let scheme = AjtaiCommitmentScheme::rand(&mut rng);

    let wit = Witness::from_w_ccs::<DP>(w_ccs);
    let cm_i = CCCS {
        cm: wit.commit::<C, W, DP>(&scheme).unwrap(),
        x_ccs,
    };

    (wit, cm_i, ccs, scheme)
}

fn single_op_benchmark<
    R: Clone + UniformRand + Debug + SuitableRing,
>(
    group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>,
) {
    // Get random parameters
    let mut rng = test_rng();
    let a = R::rand(&mut rng);
    let b = R::rand(&mut rng);

    group.bench_with_input(
        "Addition NTT",
        &(a, b),
        |bench, (a, b)| {
            bench.iter(|| {
                let _ = *a + *b;
            })
        },
    );

    group.bench_with_input(
        "Substraction NTT",
        &(a, b),
        |bench, (a, b)| {
            bench.iter(|| {
                let _ = *a - *b;
            })
        },
    );

    group.bench_with_input(
        "Multiplication NTT",
        &(a, b),
        |bench, (a, b)| {
            bench.iter(|| {
                let _ = *a * *b;
            })
        },
    );
}

fn linearization_operations<
    R: Clone + UniformRand + Debug + SuitableRing,
    DP: DecompositionParams,
    const C: usize,
    const W: usize
>(
    group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>,
) {
    let mut rng = test_rng();
    let (wit, cm_i, ccs, _) =
        setup_test_environment::<R, DP, C, W>(None);
    let z_ccs = cm_i.get_z_vector(&wit.w_ccs);

    // Prepare MLEs
    let mz_mles = calculate_Mz_mles::<R, LinearizationError<R>>(&ccs, &z_ccs).unwrap();
    let beta_s = (0..ccs.s).map(|_| R::rand(&mut rng)).collect::<Vec<R>>();

    // MZ mles
    group.bench_with_input(
        "Evaluate Mz_MLEs",
        &(ccs.clone(), z_ccs.clone()),
        |bench, (ccs, z_ccs)| {
            bench.iter(|| {
                let _ = calculate_Mz_mles::<R, LinearizationError<R>>(ccs, z_ccs);
            })
        },
    );

    // Prepare the main linearization polynomial.
    group.bench_with_input(
        "Prepare Sumcheck Poly",
        &(ccs, mz_mles, beta_s),
        |bench, (ccs, mz_mles, beta_s)| {
            bench.iter(|| {
                // Construct the sumcheck polynomial g
                let _ = prepare_lin_sumcheck_polynomial(ccs.s, &ccs.c, mz_mles, &ccs.S, beta_s).unwrap();
            })
        },
    );

}

fn single_operation_benchmarks(c: &mut Criterion) {

    let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    let mut group = c.benchmark_group("Single Operations Goldilocks");
    group.plot_config(plot_config.clone());

    // Individual operations
    single_op_benchmark::<GoldilocksRingNTT>(&mut group);

    // Linearization (W, L)
    // W = 2^9-16
    // W = WITLEN * L
    linearization_operations::<GoldilocksRingNTT, GoldilocksDP, 4, {4 * GoldilocksDP::L}>(&mut group);


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
