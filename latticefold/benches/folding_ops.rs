#![allow(incomplete_features)]
use ark_std::{test_rng, time::Duration, UniformRand};
use criterion::{
    criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion, PlotConfiguration,
};
use cyclotomic_rings::challenge_set::LatticefoldChallengeSet;
use cyclotomic_rings::rings::{GoldilocksChallengeSet, GoldilocksRingNTT, SuitableRing};
use latticefold::arith::{Witness, CCCS, CCS, LCCCS};
use latticefold::commitment::AjtaiCommitmentScheme;
use latticefold::decomposition_parameters::DecompositionParams;
use latticefold::nifs::decomposition::{
    DecompositionProver, DecompositionVerifier, LFDecompositionProver, LFDecompositionVerifier,
};
use latticefold::nifs::folding::utils::{
    compute_v0_u0_x0_cm_0, create_sumcheck_polynomial, sumcheck_polynomial_comb_fn,
};
use latticefold::nifs::folding::LFFoldingProver;
use latticefold::nifs::linearization::{
    LFLinearizationProver, LFLinearizationVerifier, LinearizationProver, LinearizationVerifier,
};
use latticefold::transcript::poseidon::PoseidonTranscript;
use latticefold::utils::sumcheck::MLSumcheck;
use lattirust_poly::mle::DenseMultilinearExtension;
use lattirust_ring::cyclotomic_ring::CRT;
use std::fmt::Debug;
use utils::{wit_and_ccs_gen, wit_and_ccs_gen_degree_three_non_scalar, wit_and_ccs_gen_non_scalar};

mod macros;
mod utils;

struct FoldingSetup<const C: usize, R: SuitableRing> {
    cm_i_s: Vec<LCCCS<C, R>>,
    wit_s: Vec<Witness<R>>,
    ccs: CCS<R>,
    mz_mles: Vec<Vec<DenseMultilinearExtension<R>>>,
    alpha_s: Vec<R>,
    beta_s: Vec<R>,
    zeta_s: Vec<R>,
    mu_s: Vec<R>,
    f_hat_mles: Vec<Vec<DenseMultilinearExtension<R>>>,
    ris: Vec<Vec<R>>,
    prechallenged_m_s: (DenseMultilinearExtension<R>, DenseMultilinearExtension<R>),
    g_mles: Vec<DenseMultilinearExtension<R>>,
    g_degree: usize,
    r_0: Vec<R>,
    theta_s: Vec<Vec<R>>,
    eta_s: Vec<Vec<R>>,
    rho_s: Vec<R::CoefficientRepresentation>,
    rho_s_ntt: Vec<R>,
    f_0: Vec<R>,
}

fn folding_setup<R, CS, DP, const C: usize, const W: usize, const WIT_LEN: usize>(
    cm_i: CCCS<C, R>,
    wit: Witness<R>,
    ccs: CCS<R>,
    scheme: AjtaiCommitmentScheme<C, W, R>,
) -> FoldingSetup<C, R>
where
    R: SuitableRing,
    CS: LatticefoldChallengeSet<R>,
    DP: DecompositionParams,
{
    let mut prover_transcript = PoseidonTranscript::<R, CS>::default();
    let mut verifier_transcript = PoseidonTranscript::<R, CS>::default();

    let (_, linearization_proof) = LFLinearizationProver::<_, PoseidonTranscript<R, CS>>::prove(
        &cm_i,
        &wit,
        &mut prover_transcript,
        &ccs,
    )
    .expect("Linearization proof generation failed");

    let lcccs = LFLinearizationVerifier::<_, PoseidonTranscript<R, CS>>::verify(
        &cm_i,
        &linearization_proof,
        &mut verifier_transcript,
        &ccs,
    )
    .expect("Linearization verification failed");

    let (mz_mles, _, wit_vec, decomposition_proof) =
        LFDecompositionProver::<_, PoseidonTranscript<R, CS>>::prove::<W, C, DP>(
            &lcccs,
            &wit,
            &mut prover_transcript,
            &ccs,
            &scheme,
        )
        .expect("Decomposition proof generation failed");

    let lcccs_vec = LFDecompositionVerifier::<_, PoseidonTranscript<R, CS>>::verify::<C, DP>(
        &lcccs,
        &decomposition_proof,
        &mut verifier_transcript,
        &ccs,
    )
    .expect("Decomposition verification failed");

    let (lcccs, mut wit_s, mz_mles) = {
        let mut lcccs = lcccs_vec.clone();
        let mut lcccs_r = lcccs_vec;
        lcccs.append(&mut lcccs_r);

        let mut wit_s = wit_vec.clone();
        let mut wit_s_r = wit_vec;
        wit_s.append(&mut wit_s_r);

        let mut mz_mles_vec = mz_mles.clone();
        let mut mz_mles_r = mz_mles;
        mz_mles_vec.append(&mut mz_mles_r);
        (lcccs, wit_s, mz_mles_vec)
    };

    let mut rng = test_rng();
    let alpha_s = (0..2 * DP::K)
        .map(|_| R::rand(&mut rng))
        .collect::<Vec<R>>();
    let beta_s = (0..ccs.s).map(|_| R::rand(&mut rng)).collect::<Vec<R>>();
    let zeta_s = (0..2 * DP::K)
        .map(|_| R::rand(&mut rng))
        .collect::<Vec<R>>();
    let mu_s = (0..2 * DP::K)
        .map(|_| R::rand(&mut rng))
        .collect::<Vec<R>>();

    let f_hat_mles = LFFoldingProver::<R, PoseidonTranscript<R, CS>>::setup_f_hat_mles(&mut wit_s);
    let ris = LFFoldingProver::<R, PoseidonTranscript<R, CS>>::get_ris(&lcccs);
    let prechallenged_m_s_1 =
        LFFoldingProver::<R, PoseidonTranscript<R, CS>>::calculate_challenged_mz_mle(
            &mz_mles[0..DP::K],
            &zeta_s[0..DP::K],
        )
        .expect("Failed to calculate first prechallenged_m_s");
    let prechallenged_m_s_2 =
        LFFoldingProver::<R, PoseidonTranscript<R, CS>>::calculate_challenged_mz_mle(
            &mz_mles[DP::K..2 * DP::K],
            &zeta_s[DP::K..2 * DP::K],
        )
        .expect("Failed to calculate second prechallenged_m_s");
    let (g_mles, g_degree) = create_sumcheck_polynomial::<_, DP>(
        ccs.s,
        f_hat_mles.clone(),
        &alpha_s,
        &prechallenged_m_s_1,
        &prechallenged_m_s_2,
        &ris,
        &beta_s,
        &mu_s,
    )
    .expect("Failed to create sumcheck polynomial");

    let r_0 = (0..ccs.s).map(|_| R::rand(&mut rng)).collect::<Vec<R>>();

    let theta_s = LFFoldingProver::<R, PoseidonTranscript<R, CS>>::get_thetas(&f_hat_mles, &r_0)
        .expect("Failed to get thetas");

    let eta_s = LFFoldingProver::<R, PoseidonTranscript<R, CS>>::get_etas(&mz_mles, &r_0)
        .expect("Failed to get etas");

    let rho_s: Vec<R::CoefficientRepresentation> = (0..2 * DP::K)
        .map(|_| R::CoefficientRepresentation::rand(&mut test_rng()))
        .collect();
    let rho_s_ntt = CRT::elementwise_crt(rho_s.clone());

    let f_0 = LFFoldingProver::<R, PoseidonTranscript<R, CS>>::compute_f_0(&rho_s_ntt, &wit_s);

    FoldingSetup {
        cm_i_s: lcccs,
        wit_s,
        ccs,
        mz_mles,
        alpha_s,
        beta_s,
        zeta_s,
        mu_s,
        f_hat_mles,
        ris,
        prechallenged_m_s: (prechallenged_m_s_1, prechallenged_m_s_2),
        g_mles,
        g_degree,
        r_0,
        theta_s,
        eta_s,
        rho_s,
        rho_s_ntt,
        f_0,
    }
}

fn setup_folding_test_environment<
    const C: usize,
    const WIT_LEN: usize,
    const W: usize,
    R: Clone + UniformRand + Debug + SuitableRing,
    CS: LatticefoldChallengeSet<R> + Clone,
    DP: DecompositionParams,
>() -> FoldingSetup<C, R> {
    let r1cs_rows = 1 + WIT_LEN + 1;
    let (cm_i, wit, ccs, scheme) = wit_and_ccs_gen::<1, C, WIT_LEN, W, DP, R>(r1cs_rows);

    folding_setup::<R, CS, DP, C, W, WIT_LEN>(cm_i, wit, ccs, scheme)
}

fn setup_folding_non_scalar_test_environment<
    const C: usize,
    const WIT_LEN: usize,
    const W: usize,
    R: Clone + UniformRand + Debug + SuitableRing,
    CS: LatticefoldChallengeSet<R> + Clone,
    DP: DecompositionParams,
>() -> FoldingSetup<C, R> {
    let r1cs_rows = 1 + WIT_LEN + 1;
    let (cm_i, wit, ccs, scheme) = wit_and_ccs_gen_non_scalar::<1, C, WIT_LEN, W, DP, R>(r1cs_rows);

    folding_setup::<R, CS, DP, C, W, WIT_LEN>(cm_i, wit, ccs, scheme)
}
fn setup_folding_degree_three_non_scalar_test_environment<
    const C: usize,
    const WIT_LEN: usize,
    const W: usize,
    R: Clone + UniformRand + Debug + SuitableRing,
    CS: LatticefoldChallengeSet<R> + Clone,
    DP: DecompositionParams,
>() -> FoldingSetup<C, R> {
    let r1cs_rows = 1 + WIT_LEN + 1;
    let (cm_i, wit, ccs, scheme) =
        wit_and_ccs_gen_degree_three_non_scalar::<1, C, WIT_LEN, W, DP, R>(r1cs_rows);

    folding_setup::<R, CS, DP, C, W, WIT_LEN>(cm_i, wit, ccs, scheme)
}
fn folding_operations<
    const C: usize,
    const W: usize,
    const WIT_LEN: usize,
    R: Clone + UniformRand + Debug + SuitableRing,
    CS: LatticefoldChallengeSet<R> + Clone,
    DP: DecompositionParams,
>(
    group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>,
    setup: &FoldingSetup<C, R>,
) {
    group.bench_with_input(
        BenchmarkId::new(
            "Evaluate Mz_MLEs",
            format!(
                "Kappa={}, W_CCS={}, W={}, B={}, L={}, B_SMALL={}, K={}",
                C,
                WIT_LEN,
                W,
                DP::B,
                DP::L,
                DP::B_SMALL,
                DP::K
            ),
        ),
        &setup,
        |bench, setup| {
            bench.iter(|| {
                let _ =
                    LFFoldingProver::<R, PoseidonTranscript<R, CS>>::calculate_challenged_mz_mle(
                        &setup.mz_mles[0..DP::K],
                        &setup.zeta_s[0..DP::K],
                    )
                    .unwrap();
                let _ =
                    LFFoldingProver::<R, PoseidonTranscript<R, CS>>::calculate_challenged_mz_mle(
                        &setup.mz_mles[DP::K..2 * DP::K],
                        &setup.zeta_s[DP::K..2 * DP::K],
                    )
                    .unwrap();
            })
        },
    );

    group.bench_with_input(
        BenchmarkId::new(
            "Create sumcheck polynomial",
            format!(
                "Kappa={}, W_CCS={}, W={}, B={}, L={}, B_SMALL={}, K={}",
                C,
                WIT_LEN,
                W,
                DP::B,
                DP::L,
                DP::B_SMALL,
                DP::K
            ),
        ),
        &setup,
        |bench, setup| {
            bench.iter_batched(
                || setup.f_hat_mles.clone(),
                |f_hat_mles| {
                    let _ = create_sumcheck_polynomial::<_, DP>(
                        setup.ccs.s,
                        f_hat_mles,
                        &setup.alpha_s,
                        &setup.prechallenged_m_s.0,
                        &setup.prechallenged_m_s.1,
                        &setup.ris,
                        &setup.beta_s,
                        &setup.mu_s,
                    )
                    .unwrap();
                },
                criterion::BatchSize::LargeInput,
            )
        },
    );

    let comb_fn = |vals: &[R]| -> R { sumcheck_polynomial_comb_fn::<R, DP>(vals, &setup.mu_s) };
    let mut transcript = PoseidonTranscript::<R, CS>::default();
    
    group.bench_with_input(
        BenchmarkId::new(
            "Sumcheck",
            format!(
                "Kappa={}, W_CCS={}, W={}, B={}, L={}, B_SMALL={}, K={}",
                C, WIT_LEN, W, DP::B, DP::L, DP::B_SMALL, DP::K
            ),
        ),
        &(&setup, &comb_fn),
        |bench, (setup, comb_fn)| {
            bench.iter_batched(
                || setup.g_mles.clone(),
                |g_mles| {
                    let _ = MLSumcheck::prove_as_subprotocol(
                        &mut transcript,
                        g_mles,
                        setup.ccs.s,
                        setup.g_degree,
                        comb_fn,
                    );
                },
                criterion::BatchSize::LargeInput,
            )
        },
    );

    group.bench_with_input(
        BenchmarkId::new(
            "Get theta's",
            format!(
                "Kappa={}, W_CCS={}, W={}, B={}, L={}, B_SMALL={}, K={}",
                C,
                WIT_LEN,
                W,
                DP::B,
                DP::L,
                DP::B_SMALL,
                DP::K
            ),
        ),
        &setup,
        |bench, setup| {
            bench.iter(|| {
                let _ = LFFoldingProver::<R, PoseidonTranscript<R, CS>>::get_thetas(
                    &setup.f_hat_mles,
                    &setup.r_0,
                )
                .unwrap();
            })
        },
    );

    group.bench_with_input(
        BenchmarkId::new(
            "Get eta's",
            format!(
                "Kappa={}, W_CCS={}, W={}, B={}, L={}, B_SMALL={}, K={}",
                C,
                WIT_LEN,
                W,
                DP::B,
                DP::L,
                DP::B_SMALL,
                DP::K
            ),
        ),
        &setup,
        |bench, setup| {
            bench.iter(|| {
                let _ = LFFoldingProver::<R, PoseidonTranscript<R, CS>>::get_etas(
                    &setup.mz_mles,
                    &setup.r_0,
                )
                .unwrap();
            })
        },
    );
    group.bench_with_input(
        BenchmarkId::new(
            "Compute v0, u0, x0, cm0",
            format!(
                "Kappa={}, W_CCS={}, W={}, B={}, L={}, B_SMALL={}, K={}",
                C,
                WIT_LEN,
                W,
                DP::B,
                DP::L,
                DP::B_SMALL,
                DP::K
            ),
        ),
        &setup,
        |bench, setup| {
            bench.iter_batched(
                || (setup.rho_s.clone(), setup.rho_s_ntt.clone()),
                |(rho_s, rho_s_ntt)| {
                    let _ = compute_v0_u0_x0_cm_0(
                        rho_s,
                        rho_s_ntt,
                        &setup.theta_s,
                        &setup.cm_i_s,
                        &setup.eta_s,
                        &setup.ccs,
                    );
                },
                criterion::BatchSize::LargeInput,
            )
        },
    );

    group.bench_with_input(
        BenchmarkId::new(
            "Compute f0",
            format!(
                "Kappa={}, W_CCS={}, W={}, B={}, L={}, B_SMALL={}, K={}",
                C,
                WIT_LEN,
                W,
                DP::B,
                DP::L,
                DP::B_SMALL,
                DP::K
            ),
        ),
        &setup,
        |bench, setup| {
            bench.iter(|| {
                let _ = LFFoldingProver::<R, PoseidonTranscript<R, CS>>::compute_f_0(
                    &setup.rho_s_ntt,
                    &setup.wit_s,
                );
            })
        },
    );
    group.bench_with_input(
        BenchmarkId::new(
            "Compute w0",
            format!(
                "Kappa={}, W_CCS={}, W={}, B={}, L={}, B_SMALL={}, K={}",
                C,
                WIT_LEN,
                W,
                DP::B,
                DP::L,
                DP::B_SMALL,
                DP::K
            ),
        ),
        &setup,
        |bench, setup| {
            bench.iter_batched(
                || setup.f_0.clone(),
                |f| {
                    let _ = Witness::from_f::<DP>(f);
                },
                criterion::BatchSize::LargeInput,
            )
        },
    );
}

macro_rules! run_single_operations_goldilocks_benchmark {
    ($crit_group:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_SMALL:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_SMALL, $k);
        paste::paste! {
            let setup = setup_folding_test_environment::< $cw, $w, {$w * $l}, GoldilocksRingNTT, GoldilocksChallengeSet, [<DecompParamsWithB $b W $w b $b_SMALL K $k>]>();
            folding_operations::<$cw, {$w * $l}, $w, GoldilocksRingNTT, GoldilocksChallengeSet, [<DecompParamsWithB $b W $w b $b_SMALL K $k>]>($crit_group, &setup);
        }
    };
}

macro_rules! run_single_operations_non_scalar_goldilocks_benchmark {
    ($crit_group:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_SMALL:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_SMALL, $k);
        paste::paste! {
            let setup = setup_folding_non_scalar_test_environment::< $cw, $w, {$w * $l}, GoldilocksRingNTT, GoldilocksChallengeSet, [<DecompParamsWithB $b W $w b $b_SMALL K $k>]>();
            folding_operations::<$cw, {$w * $l}, $w, GoldilocksRingNTT, GoldilocksChallengeSet, [<DecompParamsWithB $b W $w b $b_SMALL K $k>]>($crit_group, &setup);
        }
    };
}

macro_rules! run_single_operations_degree_three_non_scalar_goldilocks_benchmark {
    ($crit_group:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_SMALL:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_SMALL, $k);
        paste::paste! {
            let setup = setup_folding_degree_three_non_scalar_test_environment::< $cw, $w, {$w * $l}, GoldilocksRingNTT, GoldilocksChallengeSet, [<DecompParamsWithB $b W $w b $b_SMALL K $k>]>();
            folding_operations::<$cw, {$w * $l}, $w, GoldilocksRingNTT, GoldilocksChallengeSet, [<DecompParamsWithB $b W $w b $b_SMALL K $k>]>($crit_group, &setup);
        }
    };
}

fn single_operation_benchmarks(c: &mut Criterion) {
    // Goldilocks
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Single Folding Operations Goldilocks");
        group.plot_config(plot_config.clone());

        run_goldilocks_operations_benchmarks!(group);
    }

    // Goldilocks non-scalar
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Single Folding Operations Goldilocks Non-Scalar");
        group.plot_config(plot_config.clone());

        run_goldilocks_operations_non_scalar_benchmarks!(group);
    }

    // // Goldilocks degree three non-scalar
    // {
    //     let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    //     let mut group =
    //         c.benchmark_group("Single Folding Operations Goldilocks Degree Three Non-Scalar");
    //     group.plot_config(plot_config.clone());

    //     run_goldilocks_operations_degree_three_non_scalar_benchmarks!(group);
    // }
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
