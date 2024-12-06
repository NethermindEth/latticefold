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

#[derive(Clone)]
pub struct GoldilocksDP;
impl DecompositionParams for GoldilocksDP {
    const B: u128 = 1 << 15;
    const L: usize = 4;
    const B_SMALL: usize = 2;
    const K: usize = 15;
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

fn single_op_benchmark<R: Clone + UniformRand + Debug + SuitableRing>(
    group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>,
) {
    // Get random parameters
    let mut rng = test_rng();
    let a = R::rand(&mut rng);
    let b = R::rand(&mut rng);
    let a_field = <R::CoefficientRepresentation as PolyRing>::BaseRing::rand(&mut rng);
    let b_field = <R::CoefficientRepresentation as PolyRing>::BaseRing::rand(&mut rng);

    group.bench_with_input("Addition Field", &(a_field, b_field), |bench, (a, b)| {
        bench.iter(|| {
            let _ = *a + *b;
        })
    });

    group.bench_with_input(
        "Substraction Field",
        &(a_field, b_field),
        |bench, (a, b)| {
            bench.iter(|| {
                let _ = *a * *b;
            })
        },
    );

    group.bench_with_input(
        "Multiplication Field",
        &(a_field, b_field),
        |bench, (a, b)| {
            bench.iter(|| {
                let _ = *a * *b;
            })
        },
    );

    group.bench_with_input("Addition NTT", &(a, b), |bench, (a, b)| {
        bench.iter(|| {
            let _ = *a + *b;
        })
    });

    group.bench_with_input("Substraction NTT", &(a, b), |bench, (a, b)| {
        bench.iter(|| {
            let _ = *a - *b;
        })
    });

    group.bench_with_input("Multiplication NTT", &(a, b), |bench, (a, b)| {
        bench.iter(|| {
            let _ = *a * *b;
        })
    });
}

fn linearization_operations<
    R: Clone + UniformRand + Debug + SuitableRing,
    DP: DecompositionParams,
    const C: usize,
    const W: usize,
    const WIT_LEN: usize,
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

fn single_operation_benchmarks(c: &mut Criterion) {
    let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    let mut group = c.benchmark_group("Single Operations Goldilocks");
    group.plot_config(plot_config.clone());

    // Individual operations
    single_op_benchmark::<GoldilocksRingNTT>(&mut group);

    // Linearization
    // Please note that C is not used until decomposition.
    // The only parameter that we are interested on varying for linearization is W (as it already includes WIT_LEN and DP::L)
    // We explore parameters in the range  W = 2^9-2^16
    linearization_operations::<GoldilocksRingNTT, GoldilocksDP, 4, 512, 128>(&mut group);
    linearization_operations::<GoldilocksRingNTT, GoldilocksDP, 4, 1024, 256>(&mut group);
    linearization_operations::<GoldilocksRingNTT, GoldilocksDP, 4, 2048, 512>(&mut group);
    linearization_operations::<GoldilocksRingNTT, GoldilocksDP, 4, 4096, 1024>(&mut group);
    linearization_operations::<GoldilocksRingNTT, GoldilocksDP, 4, 8192, 2048>(&mut group);
    linearization_operations::<GoldilocksRingNTT, GoldilocksDP, 4, 16384, 4096>(&mut group);
    linearization_operations::<GoldilocksRingNTT, GoldilocksDP, 4, 32768, 8192>(&mut group);
    linearization_operations::<GoldilocksRingNTT, GoldilocksDP, 4, 65536, 16384>(&mut group);
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
