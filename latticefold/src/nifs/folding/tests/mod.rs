use crate::arith::ccs::{
    get_test_degree_three_ccs_padded, get_test_degree_three_z_non_scalar_split,
};
use crate::arith::{CCS, LCCCS};
use crate::ark_base::Vec;
use crate::decomposition_parameters::test_params::{
    BabyBearDP, FrogDP, GoldilocksDP, StarkFoldingDP,
};
use crate::nifs::folding::utils::SqueezeAlphaBetaZetaMu;
use crate::nifs::folding::{
    prepare_public_output,
    utils::{
        compute_v0_u0_x0_cm_0, create_sumcheck_polynomial, get_rhos, sumcheck_polynomial_comb_fn,
    },
    FoldingProver, FoldingVerifier, LFFoldingProver, LFFoldingVerifier,
};
use crate::nifs::FoldingProof;
use crate::transcript::{Transcript, TranscriptWithShortChallenges};
use crate::utils::sumcheck::MLSumcheck;
use crate::{
    arith::{Witness, CCCS},
    commitment::AjtaiCommitmentScheme,
    decomposition_parameters::DecompositionParams,
    nifs::{
        decomposition::{
            DecompositionProver, DecompositionVerifier, LFDecompositionProver,
            LFDecompositionVerifier,
        },
        linearization::{
            LFLinearizationProver, LFLinearizationVerifier, LinearizationProver,
            LinearizationVerifier,
        },
    },
    transcript::poseidon::PoseidonTranscript,
};
use ark_ff::{Field, PrimeField};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Compress, Validate};
use ark_std::{io::Cursor, test_rng};
use cyclotomic_rings::challenge_set::LatticefoldChallengeSet;
use cyclotomic_rings::rings::{
    BabyBearChallengeSet, BabyBearRingNTT, FrogChallengeSet, GoldilocksChallengeSet,
    StarkChallengeSet, SuitableRing,
};
use lattirust_poly::mle::DenseMultilinearExtension;

use lattirust_ring::cyclotomic_ring::models::{
    frog_ring::RqNTT as FrogRqNTT, goldilocks::RqNTT as GoldilocksRqNTT,
    stark_prime::RqNTT as StarkRqNTT,
};
use lattirust_ring::cyclotomic_ring::{CRT, ICRT};
use lattirust_ring::Ring;
use num_traits::{One, Zero};

const C: usize = 4;
const WIT_LEN: usize = 3;

fn setup_test_environment<RqNTT, CS, DP, const C: usize, const W: usize>(
    generate_proof: bool,
) -> (
    Vec<LCCCS<C, RqNTT>>,
    Vec<Witness<RqNTT>>,
    PoseidonTranscript<RqNTT, CS>,
    CCS<RqNTT>,
    Option<FoldingProof<RqNTT>>,
    Vec<Vec<DenseMultilinearExtension<RqNTT>>>,
)
where
    RqNTT: SuitableRing,
    CS: LatticefoldChallengeSet<RqNTT>,
    DP: DecompositionParams,
{
    let ccs = get_test_degree_three_ccs_padded::<RqNTT>(W, DP::L);

    let mut rng = test_rng();
    let (_, x_ccs, w_ccs) = get_test_degree_three_z_non_scalar_split();

    let scheme = AjtaiCommitmentScheme::rand(&mut rng);

    let wit = Witness::from_w_ccs::<DP>(w_ccs);
    let cm_i = CCCS {
        cm: wit.commit::<C, W, DP>(&scheme).unwrap(),
        x_ccs,
    };
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

    let (lcccs, wit_s, mz_mles) = {
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

    let folding_proof = if generate_proof {
        Some(generate_folding_proof::<_, _, C, DP>(
            &ccs,
            &mut prover_transcript,
            &lcccs,
            wit_s.clone(),
            mz_mles.clone(),
        ))
    } else {
        None
    };

    (
        lcccs,
        wit_s,
        verifier_transcript,
        ccs,
        folding_proof,
        mz_mles,
    )
}

fn generate_folding_proof<RqNTT, CS, const C: usize, DP>(
    ccs: &CCS<RqNTT>,
    prover_transcript: &mut PoseidonTranscript<RqNTT, CS>,
    lcccs: &[LCCCS<C, RqNTT>],
    wit_s: Vec<Witness<RqNTT>>,
    mz_mles: Vec<Vec<DenseMultilinearExtension<RqNTT>>>,
) -> FoldingProof<RqNTT>
where
    RqNTT: SuitableRing,
    CS: LatticefoldChallengeSet<RqNTT>,
    DP: DecompositionParams,
{
    let (_, _, folding_proof) = LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::prove::<
        C,
        DP,
    >(lcccs, wit_s, prover_transcript, ccs, &mz_mles)
    .unwrap();
    folding_proof
}

#[test]
fn test_setup_f_hat_mles() {
    type RqNTT = StarkRqNTT;
    type CS = StarkChallengeSet;
    type DP = StarkFoldingDP;
    const W: usize = WIT_LEN * DP::L;

    let (_, wit_s, _, _, _, _) = setup_test_environment::<RqNTT, CS, DP, C, W>(false);
    let f_hat_mles = LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::setup_f_hat_mles(
        &mut wit_s.clone(),
    );

    let expected_f_hat = wit_s
        .iter()
        .map(|w| w.f_hat.clone())
        .collect::<Vec<Vec<DenseMultilinearExtension<RqNTT>>>>();

    // Validate
    assert!(
        !f_hat_mles.is_empty(),
        "F_hat_mles vector should not be empty"
    );
    assert_eq!(
        f_hat_mles.len(),
        wit_s.len(),
        "Mismatch in F_hat_mles length"
    );
    assert_eq!(
        f_hat_mles, expected_f_hat,
        "F_hat_mles vector does not match expected evaluations"
    );
}

#[test]
fn test_get_ris() {
    type RqNTT = BabyBearRingNTT;
    type CS = BabyBearChallengeSet;
    type DP = BabyBearDP;
    const W: usize = WIT_LEN * DP::L;

    let (lccs, wit_s, _, _, _, _) = setup_test_environment::<RqNTT, CS, DP, C, W>(false);

    let ris = LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::get_ris(&lccs);

    // Compute expected output
    let expected_ris = lccs.iter().map(|cm_i| cm_i.r.clone()).collect::<Vec<_>>();

    // Validate
    assert!(!ris.is_empty(), "Ris vector should not be empty");
    assert_eq!(ris.len(), wit_s.len(), "Mismatch in Ris length");
    assert_eq!(
        ris, expected_ris,
        "Ris vector does not match expected evaluations"
    );
}

#[test]
fn test_get_sumcheck_randomness() {
    type RqNTT = BabyBearRingNTT;
    type CS = BabyBearChallengeSet;
    type DP = BabyBearDP;

    const W: usize = WIT_LEN * DP::L;

    let (lccs, mut wit_s, mut transcript, ccs, _, mz_mles) =
        setup_test_environment::<RqNTT, CS, DP, C, W>(false);
    let (alpha_s, beta_s, zeta_s, mu_s) = transcript.squeeze_alpha_beta_zeta_mu::<DP>(ccs.s);
    let f_hat_mles =
        LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::setup_f_hat_mles(&mut wit_s);

    let ris = LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::get_ris(&lccs);

    let prechallenged_Ms_1 =
        LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::calculate_challenged_mz_mle(
            &mz_mles[0..DP::K],
            &zeta_s[0..DP::K],
        )
        .unwrap();
    let prechallenged_Ms_2 =
        LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::calculate_challenged_mz_mle(
            &mz_mles[DP::K..2 * DP::K],
            &zeta_s[DP::K..2 * DP::K],
        )
        .unwrap();
    let (g_mles, g_degree) = create_sumcheck_polynomial::<_, DP>(
        ccs.s,
        &f_hat_mles,
        &alpha_s,
        &prechallenged_Ms_1,
        &prechallenged_Ms_2,
        &ris,
        &beta_s,
        &mu_s,
    )
    .unwrap();

    let comb_fn =
        |vals: &[RqNTT]| -> RqNTT { sumcheck_polynomial_comb_fn::<RqNTT, DP>(vals, &mu_s) };

    // Compute sumcheck proof
    let (_, prover_state) =
        MLSumcheck::prove_as_subprotocol(&mut transcript, &g_mles, ccs.s, g_degree, comb_fn);
    // Derive randomness
    let r_0 = LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::get_sumcheck_randomness(
        prover_state,
    );

    // Validate - Check dimensions
    assert_eq!(r_0.len(), ccs.s, "Randomness r_0 has the wrong length");
}

#[test]
fn test_get_thetas() {
    type RqNTT = StarkRqNTT;
    type CS = StarkChallengeSet;
    type DP = StarkFoldingDP;

    const W: usize = WIT_LEN * DP::L;

    let (lccs, mut wit_s, mut transcript, ccs, _, mz_mles) =
        setup_test_environment::<RqNTT, CS, DP, C, W>(false);
    let (alpha_s, beta_s, zeta_s, mu_s) = transcript.squeeze_alpha_beta_zeta_mu::<DP>(ccs.s);
    let f_hat_mles =
        LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::setup_f_hat_mles(&mut wit_s);

    let ris = LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::get_ris(&lccs);

    let prechallenged_Ms_1 =
        LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::calculate_challenged_mz_mle(
            &mz_mles[0..DP::K],
            &zeta_s[0..DP::K],
        )
        .unwrap();
    let prechallenged_Ms_2 =
        LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::calculate_challenged_mz_mle(
            &mz_mles[DP::K..2 * DP::K],
            &zeta_s[DP::K..2 * DP::K],
        )
        .unwrap();
    let (g_mles, g_degree) = create_sumcheck_polynomial::<_, DP>(
        ccs.s,
        &f_hat_mles,
        &alpha_s,
        &prechallenged_Ms_1,
        &prechallenged_Ms_2,
        &ris,
        &beta_s,
        &mu_s,
    )
    .unwrap();

    let comb_fn =
        |vals: &[RqNTT]| -> RqNTT { sumcheck_polynomial_comb_fn::<RqNTT, DP>(vals, &mu_s) };

    let (_, prover_state) =
        MLSumcheck::prove_as_subprotocol(&mut transcript, &g_mles, ccs.s, g_degree, comb_fn);
    let r_0 = LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::get_sumcheck_randomness(
        prover_state,
    );

    let theta_s =
        LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::get_thetas(&f_hat_mles, &r_0)
            .unwrap();

    let expected_thetas: Vec<Vec<RqNTT>> = f_hat_mles
        .iter()
        .map(|f_hat_row| {
            f_hat_row
                .iter()
                .map(|f_hat_mle| f_hat_mle.evaluate(&r_0).unwrap())
                .collect::<Vec<_>>()
        })
        .collect();

    // Validate
    assert!(!theta_s.is_empty(), "Thetas vector should not be empty");
    assert_eq!(theta_s.len(), f_hat_mles.len(), "Mismatch in Thetas length");
    assert_eq!(
        theta_s, expected_thetas,
        "Thetas vector does not match expected evaluations"
    );
}

#[test]
fn test_get_etas() {
    type RqNTT = StarkRqNTT;
    type CS = StarkChallengeSet;
    type DP = StarkFoldingDP;

    const W: usize = WIT_LEN * DP::L;

    let (lccs, mut wit_s, mut transcript, ccs, _, mz_mles) =
        setup_test_environment::<RqNTT, CS, DP, C, W>(false);
    let (alpha_s, beta_s, zeta_s, mu_s) = transcript.squeeze_alpha_beta_zeta_mu::<DP>(ccs.s);
    let f_hat_mles =
        LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::setup_f_hat_mles(&mut wit_s);

    let ris = LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::get_ris(&lccs);

    let prechallenged_Ms_1 =
        LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::calculate_challenged_mz_mle(
            &mz_mles[0..DP::K],
            &zeta_s[0..DP::K],
        )
        .unwrap();
    let prechallenged_Ms_2 =
        LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::calculate_challenged_mz_mle(
            &mz_mles[DP::K..2 * DP::K],
            &zeta_s[DP::K..2 * DP::K],
        )
        .unwrap();
    let (g_mles, g_degree) = create_sumcheck_polynomial::<_, DP>(
        ccs.s,
        &f_hat_mles,
        &alpha_s,
        &prechallenged_Ms_1,
        &prechallenged_Ms_2,
        &ris,
        &beta_s,
        &mu_s,
    )
    .unwrap();

    let comb_fn =
        |vals: &[RqNTT]| -> RqNTT { sumcheck_polynomial_comb_fn::<RqNTT, DP>(vals, &mu_s) };

    let (_, prover_state) =
        MLSumcheck::prove_as_subprotocol(&mut transcript, &g_mles, ccs.s, g_degree, comb_fn);
    let r_0 = LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::get_sumcheck_randomness(
        prover_state,
    );

    let eta_s =
        LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::get_etas(&mz_mles, &r_0).unwrap();

    let expected_eta_s: Vec<Vec<RqNTT>> = mz_mles
        .iter()
        .map(|Mz_mles| {
            Mz_mles
                .iter()
                .map(|mle| mle.evaluate(r_0.as_slice()).unwrap())
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    // Validate
    assert!(!eta_s.is_empty(), "Etas vector should not be empty");
    assert_eq!(eta_s.len(), mz_mles.len(), "Mismatch in Etas length");
    assert_eq!(
        eta_s, expected_eta_s,
        "Etas vector does not match expected evaluations"
    );
}

#[test]
fn test_get_rhos() {
    type RqNTT = BabyBearRingNTT;
    type CS = BabyBearChallengeSet;
    type DP = BabyBearDP;

    const W: usize = WIT_LEN * DP::L;

    let (_, _, mut transcript, _, _, _) = setup_test_environment::<RqNTT, CS, DP, C, W>(false);
    let mut transcript_clone = transcript.clone();

    let rho_s = get_rhos::<_, _, DP>(&mut transcript);

    // Compute expected result
    transcript_clone.absorb_field_element(&<_>::from_base_prime_field(
        <_>::from_be_bytes_mod_order(b"rho_s"),
    ));
    let mut expected_rhos = transcript_clone.get_small_challenges((2 * DP::K) - 1); // Note that we are missing the first element
    expected_rhos.push(RqNTT::ONE.icrt());

    // Validate
    assert!(!rho_s.is_empty(), "Rhos vector should not be empty");
    assert_eq!(rho_s.len(), 2 * DP::K, "Mismatch in Rhos length");
    assert_eq!(
        rho_s, expected_rhos,
        "Rhosvector does not match expected evaluations"
    );
}

#[test]
fn test_prepare_public_output() {
    type RqNTT = StarkRqNTT;
    type CS = StarkChallengeSet;
    type DP = StarkFoldingDP;

    const W: usize = WIT_LEN * DP::L;

    let (lccs, mut wit_s, mut transcript, ccs, _, mz_mles) =
        setup_test_environment::<RqNTT, CS, DP, C, W>(false);
    let (alpha_s, beta_s, zeta_s, mu_s) = transcript.squeeze_alpha_beta_zeta_mu::<DP>(ccs.s);
    let f_hat_mles =
        LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::setup_f_hat_mles(&mut wit_s);

    let ris = LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::get_ris(&lccs);

    let prechallenged_Ms_1 =
        LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::calculate_challenged_mz_mle(
            &mz_mles[0..DP::K],
            &zeta_s[0..DP::K],
        )
        .unwrap();
    let prechallenged_Ms_2 =
        LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::calculate_challenged_mz_mle(
            &mz_mles[DP::K..2 * DP::K],
            &zeta_s[DP::K..2 * DP::K],
        )
        .unwrap();
    let (g_mles, g_degree) = create_sumcheck_polynomial::<_, DP>(
        ccs.s,
        &f_hat_mles,
        &alpha_s,
        &prechallenged_Ms_1,
        &prechallenged_Ms_2,
        &ris,
        &beta_s,
        &mu_s,
    )
    .unwrap();

    let comb_fn =
        |vals: &[RqNTT]| -> RqNTT { sumcheck_polynomial_comb_fn::<RqNTT, DP>(vals, &mu_s) };

    let (_, prover_state) =
        MLSumcheck::prove_as_subprotocol(&mut transcript, &g_mles, ccs.s, g_degree, comb_fn);
    let r_0 = LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::get_sumcheck_randomness(
        prover_state,
    );
    let theta_s =
        LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::get_thetas(&f_hat_mles, &r_0)
            .unwrap();

    let eta_s =
        LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::get_etas(&mz_mles, &r_0).unwrap();

    theta_s
        .iter()
        .for_each(|thetas| transcript.absorb_slice(thetas));
    eta_s.iter().for_each(|etas| transcript.absorb_slice(etas));

    let rho_s = get_rhos::<_, _, DP>(&mut transcript);
    let (v_0, cm_0, u_0, x_0) = compute_v0_u0_x0_cm_0(&rho_s, &theta_s, &lccs, &eta_s, &ccs);
    let expected_x_0 = x_0[0..x_0.len() - 1].to_vec();
    let h = x_0.last().copied().unwrap();

    let lcccs = prepare_public_output(
        r_0.clone(),
        v_0.clone(),
        cm_0.clone(),
        u_0.clone(),
        x_0.clone(),
        h,
    );

    assert_eq!(lcccs.r, r_0, "Wrong r in LCCCS");
    assert_eq!(lcccs.v, v_0, "Wrong v in LCCCS");
    assert_eq!(lcccs.cm, cm_0, "Wrong commitment cm in LCCCS");
    assert_eq!(lcccs.u, u_0, "Wrong u in LCCCS");
    assert_eq!(lcccs.x_w, expected_x_0, "Wrong x_w in LCCCS");
    assert_eq!(lcccs.h, h, "Wrong h in LCCCS");
}

#[test]
fn test_compute_f_0() {
    type RqNTT = StarkRqNTT;
    type CS = StarkChallengeSet;
    type DP = StarkFoldingDP;

    const W: usize = WIT_LEN * DP::L;

    let (lccs, mut wit_s, mut transcript, ccs, _, Mz_mles) =
        setup_test_environment::<RqNTT, CS, DP, C, W>(false);
    let (alpha_s, beta_s, zeta_s, mu_s) = transcript.squeeze_alpha_beta_zeta_mu::<DP>(ccs.s);
    let f_hat_mles =
        LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::setup_f_hat_mles(&mut wit_s);

    let ris = LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::get_ris(&lccs);
    let prechallenged_Ms_1 =
        LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::calculate_challenged_mz_mle(
            &Mz_mles[0..DP::K],
            &zeta_s[0..DP::K],
        )
        .unwrap();
    let prechallenged_Ms_2 =
        LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::calculate_challenged_mz_mle(
            &Mz_mles[DP::K..2 * DP::K],
            &zeta_s[DP::K..2 * DP::K],
        )
        .unwrap();
    let (g_mles, g_degree) = create_sumcheck_polynomial::<_, DP>(
        ccs.s,
        &f_hat_mles,
        &alpha_s,
        &prechallenged_Ms_1,
        &prechallenged_Ms_2,
        &ris,
        &beta_s,
        &mu_s,
    )
    .unwrap();

    let comb_fn =
        |vals: &[RqNTT]| -> RqNTT { sumcheck_polynomial_comb_fn::<RqNTT, DP>(vals, &mu_s) };

    let (_, prover_state) =
        MLSumcheck::prove_as_subprotocol(&mut transcript, &g_mles, ccs.s, g_degree, comb_fn);
    let r_0 = LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::get_sumcheck_randomness(
        prover_state,
    );
    let theta_s =
        LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::get_thetas(&f_hat_mles, &r_0)
            .unwrap();

    let eta_s =
        LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::get_etas(&Mz_mles, &r_0).unwrap();
    theta_s
        .iter()
        .for_each(|thetas| transcript.absorb_slice(thetas));
    eta_s.iter().for_each(|etas| transcript.absorb_slice(etas));

    let rho_s = get_rhos::<_, _, DP>(&mut transcript);

    let f_0: Vec<RqNTT> =
        LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::compute_f_0(&rho_s, &wit_s);

    let expected_f_0 =
        rho_s
            .iter()
            .zip(&wit_s)
            .fold(vec![RqNTT::ZERO; wit_s[0].f.len()], |acc, (&rho_i, w_i)| {
                let rho_i: RqNTT = rho_i.crt();

                acc.into_iter()
                    .zip(w_i.f.iter())
                    .map(|(acc_j, w_ij)| acc_j + rho_i * w_ij)
                    .collect()
            });

    // Validate
    assert!(!f_0.is_empty(), "F_0 vector should not be empty");
    assert_eq!(f_0.len(), wit_s[0].f.len(), "Mismatch in F_0 length");
    assert_eq!(
        f_0, expected_f_0,
        "F_0 vector does not match expected evaluations"
    );
}

#[test]
fn test_full_prove() {
    type RqNTT = StarkRqNTT;
    type CS = StarkChallengeSet;
    type DP = StarkFoldingDP;
    const W: usize = WIT_LEN * DP::L;

    let (lccs, wit_s, mut transcript, ccs, _, mz_mles) =
        setup_test_environment::<RqNTT, CS, DP, C, W>(false);
    let _ = LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::prove::<C, DP>(
        &lccs,
        wit_s,
        &mut transcript,
        &ccs,
        &mz_mles,
    )
    .unwrap();
}

#[test]
fn test_proof_serialization() {
    type RqNTT = GoldilocksRqNTT;
    type CS = GoldilocksChallengeSet;
    type DP = GoldilocksDP;
    const W: usize = WIT_LEN * DP::L;

    let (_, _, _, _, proof, _) = setup_test_environment::<RqNTT, CS, DP, C, W>(true);

    let proof = proof.unwrap();

    let mut serialized_data = Vec::new();

    proof
        .serialize_with_mode(&mut serialized_data, Compress::Yes)
        .expect("Failed to serialize folding proof");

    let mut cursor = Cursor::new(serialized_data);

    let deserialized_proof =
        FoldingProof::deserialize_with_mode(&mut cursor, Compress::Yes, Validate::Yes)
            .expect("Failed to deserialize folding proof");

    assert_eq!(proof, deserialized_proof);
}

#[test]
fn test_verify_evaluation() {
    type RqNTT = GoldilocksRqNTT;
    type CS = GoldilocksChallengeSet;
    type DP = GoldilocksDP;
    const W: usize = WIT_LEN * DP::L;

    let (lccs_vec, _, mut transcript, ccs, proof, _) =
        setup_test_environment::<RqNTT, CS, DP, C, W>(true);
    let proof = proof.unwrap();

    let (alpha_s, beta_s, zeta_s, mu_s) = transcript.squeeze_alpha_beta_zeta_mu::<DP>(ccs.s);

    let nvars = ccs.s;
    let degree = 2 * DP::B_SMALL;

    let (claim_g1, claim_g3) =
        LFFoldingVerifier::<RqNTT, PoseidonTranscript<RqNTT, CS>>::calculate_claims::<C>(
            &alpha_s, &zeta_s, &lccs_vec,
        );

    let (r_0, expected_evaluation) =
        LFFoldingVerifier::<RqNTT, PoseidonTranscript<RqNTT, CS>>::verify_sumcheck_proof(
            &mut transcript,
            nvars,
            degree,
            claim_g1 + claim_g3,
            &proof,
        )
        .unwrap();

    let result =
        LFFoldingVerifier::<RqNTT, PoseidonTranscript<RqNTT, CS>>::verify_evaluation::<C, DP>(
            &alpha_s,
            &beta_s,
            &mu_s,
            &zeta_s,
            &r_0,
            expected_evaluation,
            &proof,
            &lccs_vec,
        );

    assert!(result.is_ok());
}

#[test]
fn test_calculate_claims() {
    type RqNTT = GoldilocksRqNTT;
    type CS = GoldilocksChallengeSet;
    type DP = GoldilocksDP;
    const W: usize = WIT_LEN * DP::L;

    let (lcccs_vec, _, _, _, _, _) = setup_test_environment::<RqNTT, CS, DP, C, W>(false);

    let alpha_s = vec![RqNTT::one(); 2 * DP::K];
    let zeta_s = vec![RqNTT::one(); 2 * DP::K];

    let (claim_g1, claim_g3) =
        LFFoldingVerifier::<RqNTT, PoseidonTranscript<RqNTT, CS>>::calculate_claims::<C>(
            &alpha_s, &zeta_s, &lcccs_vec,
        );

    assert_ne!(claim_g1, RqNTT::zero());
    assert_ne!(claim_g3, RqNTT::zero());

    let zero_alphas = vec![RqNTT::zero(); 2 * DP::K];
    let zero_zetas = vec![RqNTT::zero(); 2 * DP::K];

    let (zero_claim_g1, zero_claim_g3) =
        LFFoldingVerifier::<RqNTT, PoseidonTranscript<RqNTT, CS>>::calculate_claims::<C>(
            &zero_alphas,
            &zero_zetas,
            &lcccs_vec,
        );

    assert_eq!(zero_claim_g1, RqNTT::zero());
    assert_eq!(zero_claim_g3, RqNTT::zero());
}

#[test]
fn test_verify_sumcheck_proof() {
    type RqNTT = FrogRqNTT;
    type CS = FrogChallengeSet;
    type DP = FrogDP;

    const W: usize = WIT_LEN * DP::L;

    let (lcccs_vec, _, mut transcript, ccs, proof, _) =
        setup_test_environment::<RqNTT, CS, DP, C, W>(true);
    let proof = proof.unwrap();

    let (alpha_s, _, zeta_s, _) = transcript.squeeze_alpha_beta_zeta_mu::<DP>(ccs.s);

    let nvars = ccs.s;
    let degree = 2 * DP::B_SMALL;

    let (claim_g1, claim_g3) =
        LFFoldingVerifier::<RqNTT, PoseidonTranscript<RqNTT, CS>>::calculate_claims::<C>(
            &alpha_s, &zeta_s, &lcccs_vec,
        );

    let result = LFFoldingVerifier::<RqNTT, PoseidonTranscript<RqNTT, CS>>::verify_sumcheck_proof(
        &mut transcript,
        nvars,
        degree,
        claim_g1 + claim_g3,
        &proof,
    );

    match result {
        Ok((r_0, _)) => {
            assert_eq!(r_0.len(), ccs.s);
            // We can add more assertions here if needed
        }
        Err(e) => panic!("Sumcheck verification failed: {:?}", e),
    }
}

#[test]
fn test_full_verify() {
    type RqNTT = GoldilocksRqNTT;
    type CS = GoldilocksChallengeSet;
    type DP = GoldilocksDP;
    const W: usize = WIT_LEN * DP::L;

    let (lccs_vec, _, mut transcript, ccs, proof, _) =
        setup_test_environment::<RqNTT, CS, DP, C, W>(true);
    let proof = proof.unwrap();

    let result = LFFoldingVerifier::<RqNTT, PoseidonTranscript<RqNTT, CS>>::verify::<C, DP>(
        &lccs_vec,
        &proof,
        &mut transcript,
        &ccs,
    );

    assert!(result.is_ok());
}
