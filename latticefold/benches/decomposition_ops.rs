#![allow(incomplete_features)]
use ark_ff::Zero;
use ark_std::cfg_iter;
use ark_std::{test_rng, time::Duration, UniformRand};
use criterion::{
    criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion, PlotConfiguration,
};
use cyclotomic_rings::challenge_set::LatticefoldChallengeSet;
use cyclotomic_rings::rings::{GoldilocksChallengeSet, GoldilocksRingNTT, SuitableRing};
use latticefold::arith::r1cs::get_test_dummy_z_split;
use latticefold::arith::utils::mat_vec_mul;
use latticefold::arith::{Witness, CCS, LCCCS};
use latticefold::commitment::{AjtaiCommitmentScheme, Commitment};
use latticefold::decomposition_parameters::DecompositionParams;
use latticefold::nifs::decomposition::utils::{
    decompose_B_vec_into_k_vec, decompose_big_vec_into_k_vec_and_compose_back,
};
use latticefold::nifs::error::DecompositionError;
use latticefold::nifs::linearization::utils::compute_u;
use latticefold::utils::mle_helpers::{evaluate_mles, to_mles_err};
use lattirust_poly::mle::DenseMultilinearExtension;
use std::fmt::Debug;
use utils::wit_and_ccs_gen;

mod macros;
mod utils;

fn generate_decomposition_args<
    RqNTT,
    CS,
    DP,
    const WIT_LEN: usize,
    const C: usize,
    const W: usize,
>() -> (
    LCCCS<C, RqNTT>,
    CCS<RqNTT>,
    Witness<RqNTT>,
    Vec<Witness<RqNTT>>,
    Vec<RqNTT>,
    Vec<Vec<DenseMultilinearExtension<RqNTT>>>,
    AjtaiCommitmentScheme<C, W, RqNTT>,
)
where
    RqNTT: SuitableRing,
    CS: LatticefoldChallengeSet<RqNTT> + Clone,
    DP: DecompositionParams,
{
    let mut rng = test_rng();
    let r1cs_rows = 1 + WIT_LEN + 1;
    let (cm_i, wit, ccs, scheme) = wit_and_ccs_gen::<1, C, WIT_LEN, W, DP, RqNTT>(r1cs_rows);
    let (one, x_ccs, w_ccs) = get_test_dummy_z_split::<RqNTT, 1, WIT_LEN>();
    let mut z = vec![one];
    z.extend(&x_ccs);
    z.extend(&w_ccs);

    let log_m = ccs.s;
    let r: Vec<RqNTT> = (0..log_m).map(|_| RqNTT::rand(&mut rng)).collect();
    let single_mz_mles: Vec<DenseMultilinearExtension<RqNTT>> = ccs
        .M
        .iter()
        .map(|m| {
            DenseMultilinearExtension::from_slice(
                log_m,
                &mat_vec_mul(m, &z).expect("Matrix-vector multiplication failed"),
            )
        })
        .collect();

    let u = compute_u(&single_mz_mles, &r).expect("Failed to compute u");

    let v = evaluate_mles::<RqNTT, &DenseMultilinearExtension<RqNTT>, _, DecompositionError>(
        &wit.f_hat, &r,
    )
    .expect("Failed to evaluate MLEs");

    let lcccs = LCCCS {
        r,
        v,
        cm: cm_i.cm,
        u,
        x_w: x_ccs,
        h: RqNTT::one(),
    };

    let wit_s = decompose_B_vec_into_k_vec::<RqNTT, DP>(&wit.f_coeff)
        .into_iter()
        .map(Witness::from_f_coeff::<DP>)
        .collect::<Vec<_>>();
    let point_r = (0..ccs.s)
        .map(|_| RqNTT::rand(&mut test_rng()))
        .collect::<Vec<RqNTT>>();

    let x_s = {
        let mut x_ccs = lcccs.x_w.clone();
        x_ccs.push(lcccs.h);
        decompose_big_vec_into_k_vec_and_compose_back::<RqNTT, DP>(x_ccs)
    };
    let mz_mles = cfg_iter!(wit_s)
        .enumerate()
        .map(|(i, wit)| {
            let z: Vec<_> = {
                let mut z = Vec::with_capacity(x_s[i].len() + wit.w_ccs.len());

                z.extend_from_slice(&x_s[i]);
                z.extend_from_slice(&wit.w_ccs);

                z
            };

            let mles = to_mles_err::<_, _, DecompositionError, _>(
                ccs.s,
                cfg_iter!(ccs.M).map(|m| mat_vec_mul(m, &z)),
            )
            .expect("Failed to convert to MLEs");

            mles
        })
        .collect::<Vec<Vec<_>>>();

    (lcccs, ccs, wit, wit_s, point_r, mz_mles, scheme)
}

fn decomposition_operations<
    const C: usize,
    const W: usize,
    const WIT_LEN: usize,
    R: Clone + UniformRand + Debug + SuitableRing,
    CS: LatticefoldChallengeSet<R> + Clone,
    DP: DecompositionParams,
>(
    group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>,
) {
    let (lcccs, _, wit, wit_s, point_r, mz_mles, scheme) =
        generate_decomposition_args::<R, CS, DP, WIT_LEN, C, W>();

    group.bench_with_input(
        BenchmarkId::new(
            "Decompose witness",
            format!(
                "Kappa={}, W_CCS={}, W={}, B={}, L={}, B_small={}, K={}",
                C,
                WIT_LEN,
                W,
                DP::B,
                DP::L,
                DP::B_SMALL,
                DP::K
            ),
        ),
        &wit.f_coeff,
        |bench, f_coeff| {
            bench.iter_batched(
                || f_coeff,
                |f_coeff| decompose_B_vec_into_k_vec::<R, DP>(f_coeff),
                criterion::BatchSize::SmallInput,
            );
        },
    );

    group.bench_with_input(
        BenchmarkId::new(
            "Decompose x",
            format!(
                "Kappa={}, W_CCS={}, W={}, B={}, L={}, B_small={}, K={}",
                C,
                WIT_LEN,
                W,
                DP::B,
                DP::L,
                DP::B_SMALL,
                DP::K
            ),
        ),
        &(lcccs.clone().x_w, lcccs.h),
        |bench, (x_ccs, h)| {
            bench.iter_batched(
                || (x_ccs.clone(), h),
                |(mut x_ccs, &h)| {
                    x_ccs.push(h);
                    let _ = decompose_big_vec_into_k_vec_and_compose_back::<R, DP>(x_ccs);
                },
                criterion::BatchSize::SmallInput,
            );
        },
    );

    group.bench_with_input(
        BenchmarkId::new(
            "Commit witnesses",
            format!(
                "Kappa={}, W_CCS={}, W={}, B={}, L={}, B_small={}, K={}",
                C,
                WIT_LEN,
                W,
                DP::B,
                DP::L,
                DP::B_SMALL,
                DP::K
            ),
        ),
        &(lcccs.clone(), wit_s.clone(), scheme),
        |bench, (lcccs, wit_s, scheme)| {
            bench.iter(|| {
                let b = R::from(DP::B_SMALL as u128);

                let commitments_k1: Vec<_> = cfg_iter!(wit_s[1..])
                    .map(|wit| wit.commit::<C, W, DP>(scheme))
                    .collect::<Result<_, _>>()
                    .unwrap();

                let b_sum = commitments_k1
                    .iter()
                    .rev()
                    .fold(Commitment::zero(), |acc, y_i| (acc + y_i) * b);

                let mut result = Vec::with_capacity(wit_s.len());
                result.push(&lcccs.cm - b_sum);
                result.extend(commitments_k1);
            });
        },
    );

    group.bench_with_input(
        BenchmarkId::new(
            "compute v's",
            format!(
                "Kappa={}, W_CCS={}, W={}, B={}, L={}, B_small={}, K={}",
                C,
                WIT_LEN,
                W,
                DP::B,
                DP::L,
                DP::B_SMALL,
                DP::K
            ),
        ),
        &(wit_s.clone(), point_r.clone()),
        |bench, (wit_s, point_r)| {
            bench.iter(|| {
                cfg_iter!(wit_s)
                    .map(|wit| evaluate_mles::<R, _, _, DecompositionError>(&wit.f_hat, point_r))
                    .collect::<Result<Vec<_>, _>>()
            });
        },
    );

    group.bench_with_input(
        BenchmarkId::new(
            "compute u's",
            format!(
                "Kappa={}, W_CCS={}, W={}, B={}, L={}, B_small={}, K={}",
                C,
                WIT_LEN,
                W,
                DP::B,
                DP::L,
                DP::B_SMALL,
                DP::K
            ),
        ),
        &(mz_mles, point_r),
        |bench, (mz_mles, point_r)| {
            bench.iter(|| {
                cfg_iter!(mz_mles)
                    .map(|mles| {
                        let u_s_for_i = evaluate_mles::<
                            R,
                            &DenseMultilinearExtension<_>,
                            _,
                            DecompositionError,
                        >(mles, point_r)?;

                        Ok(u_s_for_i)
                    })
                    .collect::<Result<Vec<Vec<R>>, DecompositionError>>()
            });
        },
    );
}

macro_rules! run_single_decomposition_goldilocks_benchmark {
    ($crit_group:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            decomposition_operations::<$cw, {$w * $l}, $w, GoldilocksRingNTT, GoldilocksChallengeSet, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit_group);
        }
    };
}

fn single_operation_benchmarks(c: &mut Criterion) {
    let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    let mut group = c.benchmark_group("Single Decomposition Operations Goldilocks");
    group.plot_config(plot_config.clone());

    // Linearization
    // Please note that C is not used until decomposition.
    // The only parameter that we are interested on varying for linearization is W (as it already includes WIT_LEN and DP::L)
    // We explore parameters in the range  W = 2^9-2^16
    run_goldilocks_decomposition_benchmarks!(group);
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
