#![allow(incomplete_features)]
use ark_std::{test_rng, time::Duration, UniformRand};
use criterion::{
    criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion, PlotConfiguration,
};
use cyclotomic_rings::challenge_set::LatticefoldChallengeSet;
use cyclotomic_rings::rings::{GoldilocksChallengeSet, GoldilocksRingNTT, SuitableRing};
use latticefold::arith::{Witness, CCS, LCCCS};
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
use std::fmt::Debug;
use std::sync::Arc;
use utils::wit_and_ccs_gen;

mod macros;
mod utils;

fn setup_test_environment<RqNTT, CS, DP, const C: usize, const W: usize, const WIT_LEN: usize>() -> (
    Vec<LCCCS<C, RqNTT>>,
    Vec<Witness<RqNTT>>,
    CCS<RqNTT>,
    Vec<Vec<DenseMultilinearExtension<RqNTT>>>,
    (Vec<RqNTT>, Vec<RqNTT>, Vec<RqNTT>, Vec<RqNTT>),
    Vec<Vec<DenseMultilinearExtension<RqNTT>>>,
    Vec<Vec<RqNTT>>,
    (
        DenseMultilinearExtension<RqNTT>,
        DenseMultilinearExtension<RqNTT>,
    ),
    (Vec<Arc<DenseMultilinearExtension<RqNTT>>>, usize),
    Vec<RqNTT>,
    (
        Vec<Vec<RqNTT>>,
        Vec<Vec<RqNTT>>,
        Vec<RqNTT::CoefficientRepresentation>,
    ),
    Vec<RqNTT>,
)
where
    RqNTT: SuitableRing,
    CS: LatticefoldChallengeSet<RqNTT>,
    DP: DecompositionParams,
{
    let r1cs_rows = 1 + WIT_LEN + 1;
    let (cm_i, wit, ccs, scheme) = wit_and_ccs_gen::<1, C, WIT_LEN, W, DP, RqNTT>(r1cs_rows);

    let mut prover_transcript = PoseidonTranscript::<RqNTT, CS>::default();
    let mut verifier_transcript = PoseidonTranscript::<RqNTT, CS>::default();

    let (_, linearization_proof) =
        LFLinearizationProver::<_, PoseidonTranscript<RqNTT, CS>>::prove(
            &cm_i,
            &wit,
            &mut prover_transcript,
            &ccs,
        )
        .unwrap();

    let lcccs = LFLinearizationVerifier::<_, PoseidonTranscript<RqNTT, CS>>::verify(
        &cm_i,
        &linearization_proof,
        &mut verifier_transcript,
        &ccs,
    )
    .unwrap();

    let (mz_mles, _, wit_vec, decomposition_proof) =
        LFDecompositionProver::<_, PoseidonTranscript<RqNTT, CS>>::prove::<W, C, DP>(
            &lcccs,
            &wit,
            &mut prover_transcript,
            &ccs,
            &scheme,
        )
        .unwrap();

    let lcccs_vec = LFDecompositionVerifier::<_, PoseidonTranscript<RqNTT, CS>>::verify::<C, DP>(
        &lcccs,
        &decomposition_proof,
        &mut verifier_transcript,
        &ccs,
    )
    .unwrap();

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

    let alpha_s = (0..2 * DP::K)
        .map(|_| RqNTT::rand(&mut test_rng()))
        .collect::<Vec<_>>();
    let beta_s = (0..2 * DP::K)
        .map(|_| RqNTT::rand(&mut test_rng()))
        .collect::<Vec<_>>();
    let zeta_s = (0..2 * DP::K)
        .map(|_| RqNTT::rand(&mut test_rng()))
        .collect::<Vec<_>>();
    let mu_s = (0..2 * DP::K)
        .map(|_| RqNTT::rand(&mut test_rng()))
        .collect::<Vec<_>>();

    let f_hat_mles =
        LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::setup_f_hat_mles(&mut wit_s);
    let ris = LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::get_ris(&lcccs);
    let prechallenged_m_s_1 =
        LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::calculate_challenged_mz_mle(
            &mz_mles[0..DP::K],
            &zeta_s[0..DP::K],
        )
        .unwrap();
    let prechallenged_m_s_2 =
        LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::calculate_challenged_mz_mle(
            &mz_mles[DP::K..2 * DP::K],
            &zeta_s[DP::K..2 * DP::K],
        )
        .unwrap();
    let (g_mles, g_degree) = create_sumcheck_polynomial::<_, DP>(
        ccs.s,
        &f_hat_mles,
        &alpha_s,
        &prechallenged_m_s_1,
        &prechallenged_m_s_2,
        &ris,
        &beta_s,
        &mu_s,
    )
    .unwrap();
    let r_0 = (0..ccs.s)
        .map(|_| RqNTT::rand(&mut test_rng()))
        .collect::<Vec<_>>();
    let theta_s =
        LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::get_thetas(&f_hat_mles, &r_0)
            .unwrap();
    let eta_s =
        LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::get_etas(&mz_mles, &r_0).unwrap();
    let rho_s = (0..ccs.s)
        .map(|_| RqNTT::CoefficientRepresentation::rand(&mut test_rng()))
        .collect::<Vec<_>>();
    let f_0 = LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::compute_f_0(&rho_s, &wit_s);
    (
        lcccs,
        wit_s,
        ccs,
        mz_mles,
        (alpha_s, beta_s, zeta_s, mu_s),
        f_hat_mles,
        ris,
        (prechallenged_m_s_1, prechallenged_m_s_2),
        (g_mles, g_degree),
        r_0,
        (theta_s, eta_s, rho_s),
        f_0,
    )
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
) {
    let (
        cm_i_s,
        wit_s,
        ccs,
        mz_mles,
        (alpha_s, beta_s, zeta_s, mu_s),
        f_hat_mles,
        ris,
        (prechallenged_m_s_1, prechallenged_m_s_2),
        (g_mles, g_degree),
        r_0,
        (theta_s, eta_s, rho_s),
        f_0,
    ) = setup_test_environment::<R, CS, DP, C, W, WIT_LEN>();

    // MZ mles
    group.bench_with_input(
        BenchmarkId::new(
            "Evaluate Folding Mz_MLEs",
            format!(
                "Kappa={}, W={}, B={}, L={}, B_small={}, K={}",
                C,
                W,
                DP::B,
                DP::L,
                DP::B_SMALL,
                DP::K
            ),
        ),
        &(mz_mles.clone(), zeta_s.clone()),
        |bench, (mz_mles, zeta_s)| {
            bench.iter(|| {
                let _ =
                    LFFoldingProver::<R, PoseidonTranscript<R, CS>>::calculate_challenged_mz_mle(
                        &mz_mles[0..DP::K],
                        &zeta_s[0..DP::K],
                    )
                    .unwrap();
                let _ =
                    LFFoldingProver::<R, PoseidonTranscript<R, CS>>::calculate_challenged_mz_mle(
                        &mz_mles[DP::K..2 * DP::K],
                        &zeta_s[DP::K..2 * DP::K],
                    )
                    .unwrap();
            })
        },
    );

    group.bench_with_input(
        BenchmarkId::new(
            "Folding create sumcheck polynomial",
            format!(
                "Kappa={}, W={}, B={}, L={}, B_small={}, K={}",
                C,
                W,
                DP::B,
                DP::L,
                DP::B_SMALL,
                DP::K
            ),
        ),
        &(
            ccs.s,
            f_hat_mles.clone(),
            alpha_s.clone(),
            prechallenged_m_s_1.clone(),
            prechallenged_m_s_2.clone(),
            ris.clone(),
            beta_s.clone(),
            mu_s.clone(),
        ),
        |bench,
         (
            log_m,
            f_hat_mles,
            alpha_s,
            prechallenged_m_s_1,
            prechallenged_m_s_2,
            ris,
            beta_s,
            mu_s,
        )| {
            bench.iter(|| {
                let _ = create_sumcheck_polynomial::<_, DP>(
                    *log_m,
                    f_hat_mles,
                    alpha_s,
                    prechallenged_m_s_1,
                    prechallenged_m_s_2,
                    ris,
                    beta_s,
                    mu_s,
                )
                .unwrap();
            })
        },
    );

    let comb_fn = |vals: &[R]| -> R { sumcheck_polynomial_comb_fn::<R, DP>(vals, &mu_s) };
    group.bench_with_input(
        BenchmarkId::new(
            "Folding sumcheck",
            format!(
                "Kappa={}, W={}, B={}, L={}, B_small={}, K={}",
                C,
                W,
                DP::B,
                DP::L,
                DP::B_SMALL,
                DP::K
            ),
        ),
        &(
            PoseidonTranscript::<R, CS>::default(),
            g_mles.clone(),
            ccs.s,
            g_degree,
            comb_fn,
        ),
        |bench, (transcript, g_mles, log_m, g_degree, comb_fn)| {
            bench.iter_batched(
                || transcript.clone(),
                |mut t| {
                    let _ = MLSumcheck::prove_as_subprotocol(
                        &mut t, &g_mles, *log_m, *g_degree, comb_fn,
                    );
                },
                criterion::BatchSize::LargeInput,
            )
        },
    );

    group.bench_with_input(
        BenchmarkId::new(
            "Folding get theta's",
            format!(
                "Kappa={}, W={}, B={}, L={}, B_small={}, K={}",
                C,
                W,
                DP::B,
                DP::L,
                DP::B_SMALL,
                DP::K
            ),
        ),
        &(f_hat_mles.clone(), r_0.clone()),
        |bench, (f_hat_mles, r_0)| {
            bench.iter(|| {
                let _ =
                    LFFoldingProver::<R, PoseidonTranscript<R, CS>>::get_thetas(f_hat_mles, r_0)
                        .unwrap();
            })
        },
    );
    group.bench_with_input(
        BenchmarkId::new(
            "Folding get eta's",
            format!(
                "Kappa={}, W={}, B={}, L={}, B_small={}, K={}",
                C,
                W,
                DP::B,
                DP::L,
                DP::B_SMALL,
                DP::K
            ),
        ),
        &(mz_mles.clone(), r_0.clone()),
        |bench, (mz_mles, r_0)| {
            bench.iter(|| {
                let _ = LFFoldingProver::<R, PoseidonTranscript<R, CS>>::get_etas(&mz_mles, &r_0)
                    .unwrap();
            })
        },
    );

    group.bench_with_input(
        BenchmarkId::new(
            "Folding compute v0, u0, x0, cm0",
            format!(
                "Kappa={}, W={}, B={}, L={}, B_small={}, K={}",
                C,
                W,
                DP::B,
                DP::L,
                DP::B_SMALL,
                DP::K
            ),
        ),
        &(
            rho_s.clone(),
            theta_s.clone(),
            cm_i_s.clone(),
            eta_s.clone(),
            ccs.clone(),
        ),
        |bench, (rho_s, theta_s, cm_i_s, eta_s, ccs)| {
            bench.iter(|| {
                let _ = compute_v0_u0_x0_cm_0(rho_s, theta_s, cm_i_s, eta_s, ccs);
            })
        },
    );
    group.bench_with_input(
        BenchmarkId::new(
            "Folding compute f0",
            format!(
                "Kappa={}, W={}, B={}, L={}, B_small={}, K={}",
                C,
                W,
                DP::B,
                DP::L,
                DP::B_SMALL,
                DP::K
            ),
        ),
        &(rho_s.clone(), wit_s.clone()),
        |bench, (rho_s, wit_s)| {
            bench.iter(|| {
                let _ = LFFoldingProver::<R, PoseidonTranscript<R, CS>>::compute_f_0(rho_s, wit_s);
            })
        },
    );
    group.bench_with_input(
        BenchmarkId::new(
            "Folding compute w0",
            format!(
                "Kappa={}, W={}, B={}, L={}, B_small={}, K={}",
                C,
                W,
                DP::B,
                DP::L,
                DP::B_SMALL,
                DP::K
            ),
        ),
        &f_0.clone(),
        |bench, f_0| {
            bench.iter_batched(
                || f_0.clone(),
                |f| {
                    let _ = Witness::from_f::<DP>(f);
                },
                criterion::BatchSize::LargeInput,
            )
        },
    );
}
macro_rules! run_single_folding_goldilocks_benchmark {
    ($crit_group:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        define_params!($w, $b, $l, $b_small, $k);
        paste::paste! {
            folding_operations::<$cw, {$w * $l}, $w, GoldilocksRingNTT, GoldilocksChallengeSet, [<DecompParamsWithB $b W $w b $b_small K $k>]>($crit_group);
        }
    };
}

fn single_operation_benchmarks(c: &mut Criterion) {
    let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    let mut group = c.benchmark_group("Single Folding Operations Goldilocks");
    group.plot_config(plot_config.clone());

    // Linearization
    // Please note that C is not used until decomposition.
    // The only parameter that we are interested on varying for linearization is W (as it already includes WIT_LEN and DP::L)
    // We explore parameters in the range  W = 2^9-2^16
    run_goldilocks_folding_benchmarks!(group);
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
