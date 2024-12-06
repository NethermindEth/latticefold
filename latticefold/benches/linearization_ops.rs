#![allow(incomplete_features)]
use ark_std::{test_rng, time::Duration, UniformRand};
use criterion::{
    criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion, PlotConfiguration,
};
use cyclotomic_rings::rings::{GoldilocksRingNTT, SuitableRing};
use latticefold::arith::r1cs::{get_test_dummy_r1cs, get_test_dummy_z_split};
use latticefold::arith::{Arith, Instance, Witness, CCCS, CCS};
use latticefold::commitment::AjtaiCommitmentScheme;
use latticefold::decomposition_parameters::DecompositionParams;
use latticefold::nifs::error::LinearizationError;
use latticefold::nifs::linearization::utils::prepare_lin_sumcheck_polynomial;
use latticefold::utils::mle_helpers::calculate_Mz_mles;
use lattirust_ring::{PolyRing, Ring};
use std::fmt::Debug;

mod macros;

pub fn get_test_dummy_ccs<R: Ring, const X_LEN: usize, const WIT_LEN: usize, const W: usize>(
    rows_size: usize,
    l: usize,
) -> CCS<R> {
    let r1cs = get_test_dummy_r1cs::<R, X_LEN, WIT_LEN>(rows_size);
    CCS::<R>::from_r1cs_padded(r1cs, W, l)
}

pub fn wit_and_ccs_gen<
    const X_LEN: usize,
    const C: usize, // rows
    const WIT_LEN: usize,
    const W: usize, // columns
    P: DecompositionParams,
    R: Clone + UniformRand + Debug + SuitableRing,
>(
    r1cs_rows: usize,
) -> (
    CCCS<C, R>,
    Witness<R>,
    CCS<R>,
    AjtaiCommitmentScheme<C, W, R>,
) {
    let mut rng = test_rng();

    let new_r1cs_rows = if P::L == 1 && (WIT_LEN > 0 && (WIT_LEN & (WIT_LEN - 1)) == 0) {
        r1cs_rows - 2
    } else {
        r1cs_rows // This makes a square matrix but is too much memory
    };
    let ccs: CCS<R> = get_test_dummy_ccs::<R, X_LEN, WIT_LEN, W>(new_r1cs_rows, P::L);
    let (one, x_ccs, w_ccs) = get_test_dummy_z_split::<R, X_LEN, WIT_LEN>();
    let mut z = vec![one];
    z.extend(&x_ccs);
    z.extend(&w_ccs);
    ccs.check_relation(&z).expect("R1CS invalid!");

    let scheme: AjtaiCommitmentScheme<C, W, R> = AjtaiCommitmentScheme::rand(&mut rng);
    let wit: Witness<R> = Witness::from_w_ccs::<P>(w_ccs);

    let cm_i: CCCS<C, R> = CCCS {
        cm: wit.commit::<C, W, P>(&scheme).unwrap(),
        x_ccs,
    };

    (cm_i, wit, ccs, scheme)
}

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
        BenchmarkId::new("Evaluate Mz_MLEs", format!("C= {}, W= {}", C, W)),
        &(ccs.clone(), z_ccs.clone()),
        |bench, (ccs, z_ccs)| {
            bench.iter(|| {
                let _ = calculate_Mz_mles::<R, LinearizationError<R>>(ccs, z_ccs);
            })
        },
    );

    // Prepare the main linearization polynomial.
    group.bench_with_input(
        BenchmarkId::new("Prepare Sumcheck Poly", format!("C= {}, W= {}", C, W)),
        &(ccs, mz_mles, beta_s),
        |bench, (ccs, mz_mles, beta_s)| {
            bench.iter(|| {
                // Construct the sumcheck polynomial g
                let _ = prepare_lin_sumcheck_polynomial(&ccs.c, mz_mles, &ccs.S, beta_s).unwrap();
            })
        },
    );
}

macro_rules! run_single_linearization_goldilocks_benchmark {
    ($crit_group:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            linearization_operations::<$cw, {$w * $l}, $w, GoldilocksRingNTT, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit_group);
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
    run_single_linearization_goldilocks_benchmark!(&mut group, 9, 512, 512, 7, 2, 9);
    run_single_linearization_goldilocks_benchmark!(&mut group, 12, 512, 8192, 5, 2, 13);
    run_single_linearization_goldilocks_benchmark!(&mut group, 15, 512, 65536, 4, 2, 16);
    run_single_linearization_goldilocks_benchmark!(&mut group, 22, 512, 2097152, 3, 2, 21);
    run_single_linearization_goldilocks_benchmark!(&mut group, 39, 512, 4294967296, 2, 2, 32);
    run_single_linearization_goldilocks_benchmark!(&mut group, 8, 1024, 256, 8, 2, 8);
    run_single_linearization_goldilocks_benchmark!(&mut group, 9, 1024, 512, 7, 2, 9);
    run_single_linearization_goldilocks_benchmark!(&mut group, 15, 1024, 65536, 4, 2, 16);
    run_single_linearization_goldilocks_benchmark!(&mut group, 23, 1024, 2097152, 3, 2, 21);
    run_single_linearization_goldilocks_benchmark!(&mut group, 40, 1024, 4294967296, 2, 2, 32);
    run_single_linearization_goldilocks_benchmark!(&mut group, 8, 2048, 256, 8, 2, 8);
    run_single_linearization_goldilocks_benchmark!(&mut group, 10, 2048, 512, 7, 2, 9);
    run_single_linearization_goldilocks_benchmark!(&mut group, 11, 2048, 2048, 6, 2, 11);
    run_single_linearization_goldilocks_benchmark!(&mut group, 24, 2048, 2097152, 3, 2, 21);
    run_single_linearization_goldilocks_benchmark!(&mut group, 41, 2048, 4294967296, 2, 2, 32);
    run_single_linearization_goldilocks_benchmark!(&mut group, 8, 4096, 256, 8, 2, 8);
    run_single_linearization_goldilocks_benchmark!(&mut group, 10, 4096, 512, 7, 2, 9);
    run_single_linearization_goldilocks_benchmark!(&mut group, 11, 4096, 2048, 6, 2, 11);
    run_single_linearization_goldilocks_benchmark!(&mut group, 13, 4096, 8192, 5, 2, 13);
    run_single_linearization_goldilocks_benchmark!(&mut group, 17, 4096, 65536, 4, 2, 16);
    run_single_linearization_goldilocks_benchmark!(&mut group, 25, 4096, 2097152, 3, 2, 21);
    run_single_linearization_goldilocks_benchmark!(&mut group, 42, 4096, 4294967296, 2, 2, 32);
    run_single_linearization_goldilocks_benchmark!(&mut group, 9, 8192, 256, 8, 2, 8);
    run_single_linearization_goldilocks_benchmark!(&mut group, 11, 8192, 512, 7, 2, 9);
    run_single_linearization_goldilocks_benchmark!(&mut group, 12, 8192, 2048, 6, 2, 11);
    run_single_linearization_goldilocks_benchmark!(&mut group, 14, 8192, 8192, 5, 2, 13);
    run_single_linearization_goldilocks_benchmark!(&mut group, 17, 8192, 65536, 4, 2, 16);
    run_single_linearization_goldilocks_benchmark!(&mut group, 26, 8192, 2097152, 3, 2, 21);
    run_single_linearization_goldilocks_benchmark!(&mut group, 43, 8192, 4294967296, 2, 2, 32);
    run_single_linearization_goldilocks_benchmark!(&mut group, 9, 16384, 256, 8, 2, 8);
    run_single_linearization_goldilocks_benchmark!(&mut group, 11, 16384, 512, 7, 2, 9);
    run_single_linearization_goldilocks_benchmark!(&mut group, 12, 16384, 2048, 6, 2, 11);
    run_single_linearization_goldilocks_benchmark!(&mut group, 14, 16384, 8192, 5, 2, 13);
    run_single_linearization_goldilocks_benchmark!(&mut group, 18, 16384, 65536, 4, 2, 16);
    run_single_linearization_goldilocks_benchmark!(&mut group, 26, 16384, 2097152, 3, 2, 21);
    run_single_linearization_goldilocks_benchmark!(&mut group, 44, 16384, 4294967296, 2, 2, 32);
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
