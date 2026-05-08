use ark_std::test_rng;
use cyclotomic_rings::{
    challenge_set::LatticefoldChallengeSet,
    rings::{BabyBearChallengeSet, FrogChallengeSet, GoldilocksChallengeSet, StarkChallengeSet},
};
use num_traits::One;
use rand::Rng;
use stark_rings::cyclotomic_ring::models::{
    babybear::RqNTT as BabyBearRqNTT, frog_ring::RqNTT as FrogRqNTT,
    goldilocks::RqNTT as GoldilocksRqNTT, stark_prime::RqNTT as StarkRqNTT,
};

use super::*;
use crate::{
    arith::{r1cs::get_test_z_split, tests::get_test_ccs, utils::mat_vec_mul},
    commitment::AjtaiCommitmentScheme,
    decomposition_parameters::{
        test_params::{dp_babybear, dp_frog, dp_goldilocks, dp_stark},
        DecompositionParams,
    },
    nifs::linearization::utils::{sumcheck_polynomial_comb_fn, SqueezeBeta},
    transcript::poseidon::PoseidonTranscript,
};

const KAPPA: usize = 4;
const WIT_LEN: usize = 4;

fn setup_test_environment<RqNTT: SuitableRing>(
    dparams: &DecompositionParams,
    input: Option<usize>,
    n: usize,
) -> (
    Witness<RqNTT>,
    CCCS<RqNTT>,
    CCS<RqNTT>,
    AjtaiCommitmentScheme<RqNTT>,
) {
    let ccs = get_test_ccs::<RqNTT>(n, dparams.l);
    let mut rng = test_rng();
    let (_, x_ccs, w_ccs) = get_test_z_split::<RqNTT>(input.unwrap_or(rng.gen_range(0..64)));
    let scheme = AjtaiCommitmentScheme::rand(KAPPA, n, &mut rng);

    let wit = Witness::from_w_ccs(w_ccs, dparams.B, dparams.l);
    let cm_i = CCCS {
        cm: wit.commit(&scheme).unwrap(),
        x_ccs,
    };

    (wit, cm_i, ccs, scheme)
}

#[test]
fn test_compute_z_ccs() {
    type RqNTT = StarkRqNTT;
    let dparams = dp_stark();
    let n = WIT_LEN * dparams.l;
    let (wit, cm_i, _, scheme) = setup_test_environment::<RqNTT>(&dparams, None, n);

    let z_ccs = cm_i.get_z_vector(&wit.w_ccs);

    // Check z_ccs structure
    assert_eq!(z_ccs.len(), cm_i.x_ccs.len() + 1 + wit.w_ccs.len());
    assert_eq!(z_ccs[cm_i.x_ccs.len()], RqNTT::one());

    // Check commitment
    assert_eq!(cm_i.cm, wit.commit(&scheme).unwrap());
}

#[test]
fn test_construct_polynomial() {
    type RqNTT = GoldilocksRqNTT;
    type CS = GoldilocksChallengeSet;
    let dparams = dp_goldilocks();
    let n = WIT_LEN * dparams.l;
    let (wit, cm_i, ccs, _) = setup_test_environment::<RqNTT>(&dparams, None, n);

    let z_ccs = cm_i.get_z_vector(&wit.w_ccs);

    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();
    let mut prover = LFLinearizationProver::new(&mut transcript);
    let (_, g_degree, mz_mles) = prover.construct_polynomial_g(&z_ccs, &ccs).unwrap();

    // Check dimensions
    assert_eq!(mz_mles.len(), ccs.t);

    // Check degree of g
    assert!(g_degree <= ccs.q + 1)
}

#[test]
fn test_generate_sumcheck() {
    type RqNTT = FrogRqNTT;
    type CS = FrogChallengeSet;
    let dparams = dp_frog();
    let n = WIT_LEN * dparams.l;
    let (wit, cm_i, ccs, _) = setup_test_environment::<RqNTT>(&dparams, None, n);

    let z_ccs = cm_i.get_z_vector(&wit.w_ccs);

    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();
    let mut prover = LFLinearizationProver::new(&mut transcript);
    let (g_mles, g_degree, _) = prover.construct_polynomial_g(&z_ccs, &ccs).unwrap();

    let comb_fn = |vals: &[RqNTT]| -> RqNTT { sumcheck_polynomial_comb_fn::<RqNTT>(vals, &ccs) };

    let (_, point_r) = prover
        .generate_sumcheck_proof(g_mles, ccs.s, g_degree, comb_fn)
        .unwrap();

    // Check dimensions
    assert_eq!(point_r.len(), ccs.s);
}

fn prepare_test_vectors<RqNTT: SuitableRing, CS: LatticefoldChallengeSet<RqNTT>>(
    wit: &Witness<RqNTT>,
    cm_i: &CCCS<RqNTT>,
    ccs: &CCS<RqNTT>,
) -> (Vec<RqNTT>, Vec<DenseMultilinearExtension<RqNTT>>) {
    let z_ccs = cm_i.get_z_vector(&wit.w_ccs);

    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();
    let mut prover = LFLinearizationProver::new(&mut transcript);
    let (g_mles, g_degree, Mz_mles) = prover.construct_polynomial_g(&z_ccs, ccs).unwrap();

    let comb_fn = |vals: &[RqNTT]| -> RqNTT { sumcheck_polynomial_comb_fn::<RqNTT>(vals, ccs) };

    let (_, point_r) = prover
        .generate_sumcheck_proof(g_mles, ccs.s, g_degree, comb_fn)
        .unwrap();

    (point_r, Mz_mles)
}

#[test]
fn test_compute_v() {
    type RqNTT = BabyBearRqNTT;
    type CS = BabyBearChallengeSet;
    let dparams = dp_babybear();
    let n = WIT_LEN * dparams.l;

    // Setup shared test state
    let (wit, cm_i, ccs, _) = setup_test_environment::<RqNTT>(&dparams, None, n);
    let (point_r, Mz_mles) = prepare_test_vectors::<RqNTT, CS>(&wit, &cm_i, &ccs);

    // Compute actual v vector
    let (point_r, v, _) =
        LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::compute_evaluation_vectors(
            &wit, &point_r, &Mz_mles,
        )
        .unwrap();

    // Compute expected v vector (witness evaluations)
    let expected_v =
        evaluate_mles::<RqNTT, _, _, LinearizationError<RqNTT>>(&wit.f_hat, &point_r).unwrap();

    // Validate
    assert_eq!(point_r.len(), ccs.s, "point_r length mismatch");
    assert!(!v.is_empty(), "v vector should not be empty");
    assert_eq!(
        v, expected_v,
        "v vector doesn't match expected witness evaluations"
    );
}

#[test]
fn test_compute_u() {
    type RqNTT = FrogRqNTT;
    type CS = FrogChallengeSet;
    let dparams = dp_frog();
    let n = WIT_LEN * dparams.l;
    // Setup shared test state
    let (wit, cm_i, ccs, _) = setup_test_environment::<RqNTT>(&dparams, None, n);
    let (point_r, Mz_mles) = prepare_test_vectors::<RqNTT, CS>(&wit, &cm_i, &ccs);

    // Compute actual u vector
    let (point_r, _, u) =
        LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::compute_evaluation_vectors(
            &wit, &point_r, &Mz_mles,
        )
        .unwrap();

    // Compute expected u vector
    let z_ccs = cm_i.get_z_vector(&wit.w_ccs);
    let expected_Mz_mles: Vec<DenseMultilinearExtension<RqNTT>> = ccs
        .M
        .iter()
        .map(|M| {
            DenseMultilinearExtension::from_evaluations_vec(ccs.s, mat_vec_mul(M, &z_ccs).unwrap())
        })
        .collect();
    let expected_u = compute_u(&expected_Mz_mles, &point_r).unwrap();

    // Validate
    assert_eq!(point_r.len(), ccs.s, "point_r length mismatch");
    assert!(!u.is_empty(), "u vector should not be empty");
    assert_eq!(u, expected_u, "u vector doesn't match expected evaluations");
}

#[test]
fn test_full_prove() {
    type RqNTT = GoldilocksRqNTT;
    type CS = GoldilocksChallengeSet;
    let dparams = dp_goldilocks();
    let n = WIT_LEN * dparams.l;
    let (wit, cm_i, ccs, _) = setup_test_environment::<RqNTT>(&dparams, None, n);
    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();
    let mut prover = LFLinearizationProver::new(&mut transcript);

    let (lcccs, proof) = prover.prove(&cm_i, &wit, &ccs).unwrap();

    assert_eq!(lcccs.r.len(), ccs.s);
    assert_eq!(lcccs.v.len(), proof.v.len());
    assert_eq!(lcccs.u.len(), proof.u.len());
}

#[test]
fn test_verify_sumcheck_proof() {
    type RqNTT = StarkRqNTT;
    type CS = StarkChallengeSet;
    let dparams = dp_stark();
    let n = WIT_LEN * dparams.l;
    let (wit, cm_i, ccs, _) = setup_test_environment::<RqNTT>(&dparams, None, n);
    let mut prover_transcript = PoseidonTranscript::<RqNTT, CS>::default();
    let mut prover = LFLinearizationProver::new(&mut prover_transcript);

    // Generate proof
    let (lcccs, proof) = prover.prove(&cm_i, &wit, &ccs).unwrap();

    // We need to recreate the exact same transcript state
    let mut verify_transcript = PoseidonTranscript::<RqNTT, CS>::default();
    let mut verifier = LFLinearizationVerifier::new(&mut verify_transcript);

    // Generate beta challenges to match prover's transcript state
    let _ = verifier.transcript.squeeze_beta_challenges(ccs.s);

    let result = verifier.verify_sumcheck_proof(&proof, &ccs);

    // Instead of unwrapping, handle the result
    match result {
        Ok((point_r, _)) => {
            assert_eq!(point_r.len(), ccs.s);
            // We know that point_r from lcccs is valid
            assert_eq!(point_r, lcccs.r);
        }
        Err(e) => panic!("Sumcheck verification failed: {:?}", e),
    }
}

#[test]
fn test_verify_evaluation_claim() {
    type RqNTT = BabyBearRqNTT;
    type CS = BabyBearChallengeSet;
    let dparams = dp_babybear();
    let n = WIT_LEN * dparams.l;
    let (wit, cm_i, ccs, _) = setup_test_environment::<RqNTT>(&dparams, None, n);
    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();
    let mut prover = LFLinearizationProver::new(&mut transcript);

    // Generate proof
    let (_, proof) = prover.prove(&cm_i, &wit, &ccs).unwrap();

    // Reset transcript and generate verification data
    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();
    let mut verifier = LFLinearizationVerifier::new(&mut transcript);
    let beta_s = verifier.transcript.squeeze_beta_challenges(ccs.s);

    let (point_r, s) = verifier.verify_sumcheck_proof(&proof, &ccs).unwrap();

    // Test the evaluation claim verification
    let result = LFLinearizationVerifier::<_, ()>::verify_evaluation_claim(
        &beta_s, &point_r, s, &proof, &ccs,
    );

    assert!(result.is_ok());
}

#[test]
fn test_prepare_verifier_output() {
    type RqNTT = FrogRqNTT;
    type CS = FrogChallengeSet;
    let dparams = dp_frog();
    let n = WIT_LEN * dparams.l;
    let (wit, cm_i, ccs, _) = setup_test_environment::<RqNTT>(&dparams, None, n);
    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();
    let mut prover = LFLinearizationProver::new(&mut transcript);

    let (_, proof) = prover.prove(&cm_i, &wit, &ccs).unwrap();

    let point_r = vec![RqNTT::one(); ccs.s]; // Example point_r

    let lcccs =
        LFLinearizationVerifier::<RqNTT, PoseidonTranscript<RqNTT, CS>>::prepare_verifier_output(
            &cm_i,
            point_r.clone(),
            &proof,
        );

    // Verify final state structure
    assert_eq!(lcccs.r, point_r);
    assert_eq!(lcccs.v, proof.v);
    assert_eq!(lcccs.u, proof.u);
    assert_eq!(lcccs.cm, cm_i.cm);
    assert_eq!(lcccs.x_w, cm_i.x_ccs);
    assert_eq!(lcccs.h, RqNTT::one());
}

#[test]
fn test_verify_invalid_proof() {
    type RqNTT = GoldilocksRqNTT;
    type CS = GoldilocksChallengeSet;
    let dparams = dp_goldilocks();
    let n = WIT_LEN * dparams.l;
    let (wit, cm_i, ccs, _) = setup_test_environment::<RqNTT>(&dparams, None, n);
    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();
    let mut prover = LFLinearizationProver::new(&mut transcript);

    let (_, mut proof) = prover.prove(&cm_i, &wit, &ccs).unwrap();

    // Corrupt the proof
    if !proof.u.is_empty() {
        proof.u[0] += RqNTT::one();
    }

    // Reset transcript for verification
    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();
    let mut verifier = LFLinearizationVerifier::new(&mut transcript);

    let result = verifier.verify(&cm_i, &proof, &ccs);

    assert!(result.is_err());
}
