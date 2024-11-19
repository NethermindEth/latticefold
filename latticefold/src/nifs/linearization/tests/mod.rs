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

    let z_ccs = LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, StarkChallengeSet>>::compute_z_ccs::<C>(
        &wit,
        &cm_i.x_ccs
    ).unwrap();

    // Check z_ccs structure
    assert_eq!(z_ccs.len(), cm_i.x_ccs.len() + 1 + wit.w_ccs.len());
    assert_eq!(z_ccs[cm_i.x_ccs.len()], RqNTT::one());

    // Check commitment
    assert_eq!(cm_i.cm, wit.commit::<C, W, PP_STARK>(&scheme).unwrap());
}

#[test]
fn test_construct_polynomial() {
    let (wit, cm_i, ccs, scheme) = setup_test_environment::<RqNTT>();

    let z_ccs = LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, StarkChallengeSet>>::compute_z_ccs::<C>(
        &wit,
        &cm_i.x_ccs
    ).unwrap();

    let mut transcript = PoseidonTranscript::<RqNTT, StarkChallengeSet>::default();
    let (g, beta_s) = LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, StarkChallengeSet>>::construct_polynomial_g(
        &z_ccs,
        &mut transcript,
        &ccs
    ).unwrap();

    // Check dimensions
    assert_eq!(beta_s.len(), ccs.s);

    // Check degree of g
    assert!(g.aux_info.max_degree <= ccs.q + 1)
}

#[test]
fn test_generate_sumcheck() {
    let (wit, cm_i, ccs, scheme) = setup_test_environment::<RqNTT>();

    let z_state = LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, StarkChallengeSet>>::compute_z_ccs::<C>(
        &wit,
        &cm_i.x_ccs
    ).unwrap();

    let mut transcript = PoseidonTranscript::<RqNTT, StarkChallengeSet>::default();
    let (g, beta_s) = LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, StarkChallengeSet>>::construct_polynomial_g(
        &z_state,
        &mut transcript,
        &ccs
    ).unwrap();

    let (proof, point_r) = LFLinearizationProver::<
        RqNTT,
        PoseidonTranscript<RqNTT, StarkChallengeSet>,
    >::generate_sumcheck_proof(&g, &beta_s, &mut transcript)
    .unwrap();

    // Check dimensions
    assert_eq!(point_r.len(), ccs.s);
}

#[test]
fn test_compute_evaluation_vectors() {
    let (wit, cm_i, ccs, scheme) = setup_test_environment::<RqNTT>();

    let z_ccs = LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, StarkChallengeSet>>::compute_z_ccs::<C>(
            &wit,
            &cm_i.x_ccs
        ).unwrap();

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

    let eval_state = LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, StarkChallengeSet>>::compute_evaluation_vectors(
            &wit,
            &point_r,
            &ccs,
            &z_ccs
        ).unwrap();

    // Check lengths and non-empty values
    assert_eq!(eval_state.point_r.len(), ccs.s);
    assert!(!eval_state.v.is_empty());
    assert!(!eval_state.u.is_empty());

    // Check v evaluations
    let witness_f_hat: Vec<RqNTT> = cfg_iter!(wit.f_hat)
        .map(|f_hat_row| {
            DenseMultilinearExtension::from_slice(ccs.s, f_hat_row)
                .evaluate(&point_r)
                .expect("cannot end up here, because the sumcheck subroutine must yield a point of the length log m")
        })
        .collect();
    assert_eq!(eval_state.v, witness_f_hat);
    
    // Check u evaluations
    let Mz_mles: Vec<DenseMultilinearExtension<RqNTT>> = ccs
        .M
        .iter()
        .map(|M| {
            DenseMultilinearExtension::from_slice(
                ccs.s,
                &mat_vec_mul(M, &z_ccs).unwrap(),
            )
        })
        .collect();

    let u = compute_u(&Mz_mles, &point_r).unwrap();
    assert_eq!(eval_state.u, u);



}

#[test]
fn test_full_prove() {
    let (wit, cm_i, ccs, scheme) = setup_test_environment::<RqNTT>();
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
