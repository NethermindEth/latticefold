use super::*;
use crate::decomposition_parameters::test_params::PP_STARK;
use crate::{
    arith::{r1cs::get_test_z_split, tests::get_test_ccs},
    commitment::AjtaiCommitmentScheme,
    decomposition_parameters::{test_params::PP, DecompositionParams},
    transcript::poseidon::PoseidonTranscript,
};
use cyclotomic_rings::rings::StarkChallengeSet;
use lattirust_ring::cyclotomic_ring::models::stark_prime::RqNTT;
use num_traits::One;
use rand::thread_rng;

const C: usize = 4;
const WIT_LEN: usize = 4;
const W: usize = WIT_LEN * PP::L;
fn setup_test_environment<RqNTT: SuitableRing>() -> (
    Witness<RqNTT>,
    CCCS<4, RqNTT>,
    CCS<RqNTT>,
    AjtaiCommitmentScheme<C, W, RqNTT>,
) {
    let ccs = get_test_ccs::<RqNTT>(W);
    let (_, x_ccs, w_ccs) = get_test_z_split::<RqNTT>(3);
    let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());

    let wit = Witness::from_w_ccs::<PP>(w_ccs);
    let cm_i = CCCS {
        cm: wit.commit::<C, W, PP>(&scheme).unwrap(),
        x_ccs,
    };

    (wit, cm_i, ccs, scheme)
}

#[test]
fn test_compute_z_ccs() {
    let (wit, cm_i, _, scheme) = setup_test_environment::<RqNTT>();

    let z_ccs = cm_i.get_z_vector(&wit.w_ccs);

    // Check z_ccs structure
    assert_eq!(z_ccs.len(), cm_i.x_ccs.len() + 1 + wit.w_ccs.len());
    assert_eq!(z_ccs[cm_i.x_ccs.len()], RqNTT::one());

    // Check commitment
    assert_eq!(cm_i.cm, wit.commit::<C, W, PP_STARK>(&scheme).unwrap());
}

#[test]
fn test_construct_polynomial() {
    let (wit, cm_i, ccs, _scheme) = setup_test_environment::<RqNTT>();

    let z_ccs = cm_i.get_z_vector(&wit.w_ccs);

    let mut transcript = PoseidonTranscript::<RqNTT, StarkChallengeSet>::default();
    let (_g, beta_s) = LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, StarkChallengeSet>>::construct_polynomial_g(
        &z_ccs,
        &mut transcript,
        &ccs
    ).unwrap();

    // Check dimensions
    assert_eq!(beta_s.len(), ccs.s);
    // Check degree of g
    //assert!(g.aux_info.max_degree <=  )
}

#[test]
fn test_generate_sumcheck() {
    let (wit, cm_i, ccs, _scheme) = setup_test_environment::<RqNTT>();

    let z_ccs = cm_i.get_z_vector(&wit.w_ccs);

    let mut transcript = PoseidonTranscript::<RqNTT, StarkChallengeSet>::default();
    let (g, beta_s) = LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, StarkChallengeSet>>::construct_polynomial_g(
        &z_ccs,
        &mut transcript,
        &ccs
    ).unwrap();

    let (_proof, point_r) = LFLinearizationProver::<
        RqNTT,
        PoseidonTranscript<RqNTT, StarkChallengeSet>,
    >::generate_sumcheck_proof(&g, &beta_s, &mut transcript)
    .unwrap();

    assert_eq!(point_r.len(), ccs.s);
}

#[test]
fn test_compute_evaluation_vectors() {
    let (wit, cm_i, ccs, _scheme) = setup_test_environment::<RqNTT>();

    let z_ccs = cm_i.get_z_vector(&wit.w_ccs);

    let mut transcript = PoseidonTranscript::<RqNTT, StarkChallengeSet>::default();
    let (g, beta_s) = LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, StarkChallengeSet>>::construct_polynomial_g(
            &z_ccs,
            &mut transcript,
            &ccs
        ).unwrap();

    let (_, point_r) = LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, StarkChallengeSet>>::generate_sumcheck_proof(
            &g,
            &beta_s,
            &mut transcript
        ).unwrap();

    let (point_r, v, u) = LFLinearizationProver::<
        RqNTT,
        PoseidonTranscript<RqNTT, StarkChallengeSet>,
    >::compute_evaluation_vectors(&wit, &point_r, &ccs, &z_ccs)
    .unwrap();

    assert_eq!(point_r.len(), ccs.s);
    assert!(!v.is_empty());
    assert!(!u.is_empty());
}

#[test]
fn test_full_prove() {
    let (wit, cm_i, ccs, _scheme) = setup_test_environment::<RqNTT>();
    let mut transcript = PoseidonTranscript::<RqNTT, StarkChallengeSet>::default();

    let (lcccs, proof) =
        LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, StarkChallengeSet>>::prove(
            &cm_i,
            &wit,
            &mut transcript,
            &ccs,
        )
        .unwrap();

    assert_eq!(lcccs.r.len(), ccs.s);
    assert_eq!(lcccs.v.len(), proof.v.len());
    assert_eq!(lcccs.u.len(), proof.u.len());
}

#[test]
fn test_verify_sumcheck_proof() {
    let (wit, cm_i, ccs, _scheme) = setup_test_environment::<RqNTT>();
    let mut prove_transcript = PoseidonTranscript::<RqNTT, StarkChallengeSet>::default();

    // Generate proof
    let (lcccs, proof) =
        LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, StarkChallengeSet>>::prove(
            &cm_i,
            &wit,
            &mut prove_transcript,
            &ccs,
        )
        .unwrap();

    // We need to recreate the exact same transcript state
    let mut verify_transcript = PoseidonTranscript::<RqNTT, StarkChallengeSet>::default();

    // Generate beta challenges to match prover's transcript state
    let _beta_s = BetaChallengeGenerator::generate_challenges(&mut verify_transcript, ccs.s);

    let result = LFLinearizationVerifier::<
    RqNTT,
    PoseidonTranscript<RqNTT, StarkChallengeSet>,
    >::verify_sumcheck_proof(&proof, &mut verify_transcript, &ccs);

    // Instead of unwrapping, handle the result
    match result {
        Ok((point_r, _s)) => {
            assert_eq!(point_r.len(), ccs.s);
            // We know that point_r from lcccs is valid
            assert_eq!(point_r, lcccs.r);
        }
        Err(e) => panic!("Sumcheck verification failed: {:?}", e),
    }
}

#[test]
fn test_verify_evaluation_claim() {
    let (wit, cm_i, ccs, _scheme) = setup_test_environment::<RqNTT>();
    let mut transcript = PoseidonTranscript::<RqNTT, StarkChallengeSet>::default();

    // Generate proof
    let (_, proof) =
        LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, StarkChallengeSet>>::prove(
            &cm_i,
            &wit,
            &mut transcript,
            &ccs,
        )
        .unwrap();

    // Reset transcript and generate verification data
    let mut transcript = PoseidonTranscript::<RqNTT, StarkChallengeSet>::default();
    let beta_s = BetaChallengeGenerator::generate_challenges(&mut transcript, ccs.s);

    let (point_r, s) = LFLinearizationVerifier::<
    RqNTT,
    PoseidonTranscript<RqNTT, StarkChallengeSet>,
    >::verify_sumcheck_proof(&proof, &mut transcript, &ccs)
        .unwrap();

    // Test the evaluation claim verification
    let result = LFLinearizationVerifier::<
    RqNTT,
    PoseidonTranscript<RqNTT, StarkChallengeSet>,
    >::verify_evaluation_claim(&beta_s, &point_r, s, &proof, &ccs);

    assert!(result.is_ok());
}

#[test]
fn test_prepare_final_state() {
    let (wit, cm_i, ccs, _scheme) = setup_test_environment::<RqNTT>();
    let mut transcript = PoseidonTranscript::<RqNTT, StarkChallengeSet>::default();

    let (_, proof) =
        LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, StarkChallengeSet>>::prove(
            &cm_i,
            &wit,
            &mut transcript,
            &ccs,
        )
        .unwrap();

    let point_r = vec![RqNTT::one(); ccs.s]; // Example point_r

    let lcccs = LFLinearizationVerifier::<
    RqNTT,
    PoseidonTranscript<RqNTT, StarkChallengeSet>,
    >::prepare_final_state::<C>(&cm_i, point_r.clone(), &proof);

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
    let (wit, cm_i, ccs, _scheme) = setup_test_environment::<RqNTT>();
    let mut transcript = PoseidonTranscript::<RqNTT, StarkChallengeSet>::default();

    let (_, mut proof) =
        LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, StarkChallengeSet>>::prove(
            &cm_i,
            &wit,
            &mut transcript,
            &ccs,
        )
        .unwrap();

    // Corrupt the proof
    if !proof.u.is_empty() {
        proof.u[0] += RqNTT::one();
    }

    // Reset transcript for verification
    let mut transcript = PoseidonTranscript::<RqNTT, StarkChallengeSet>::default();

    let result =
        LFLinearizationVerifier::<RqNTT, PoseidonTranscript<RqNTT, StarkChallengeSet>>::verify::<C>(
            &cm_i,
            &proof,
            &mut transcript,
            &ccs,
        );

    assert!(result.is_err());
}
