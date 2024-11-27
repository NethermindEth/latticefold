use crate::arith::r1cs::get_test_z_split;
use crate::arith::tests::get_test_ccs;
use crate::arith::utils::mat_vec_mul;
use crate::arith::{Witness, CCS, LCCCS};
use crate::commitment::{AjtaiCommitmentScheme, Commitment};
use crate::decomposition_parameters::test_params::{BabyBearDP, GoldilocksDP, StarkDP};
use crate::decomposition_parameters::DecompositionParams;
use crate::nifs::decomposition::{
    DecompositionProver, DecompositionVerifier, LFDecompositionProver, LFDecompositionVerifier,
};
use crate::nifs::error::DecompositionError;
use crate::nifs::linearization::utils::compute_u;
use crate::nifs::mle_helpers::{evaluate_mles, to_mles};
use crate::transcript::poseidon::PoseidonTranscript;
use ark_std::vec::Vec;
use cyclotomic_rings::challenge_set::LatticefoldChallengeSet;
use cyclotomic_rings::rings::{
    BabyBearChallengeSet, BabyBearRingNTT, GoldilocksChallengeSet, GoldilocksRingNTT,
    StarkChallengeSet, StarkRingNTT, SuitableRing,
};
use lattirust_poly::mle::DenseMultilinearExtension;
use num_traits::One;
use rand::Rng;

fn generate_decomposition_args<RqNTT, CS, DP, const WIT_LEN: usize, const W: usize>() -> (
    LCCCS<4, RqNTT>,
    PoseidonTranscript<RqNTT, CS>,
    PoseidonTranscript<RqNTT, CS>,
    CCS<RqNTT>,
    Witness<RqNTT>,
    AjtaiCommitmentScheme<4, W, RqNTT>,
)
where
    RqNTT: SuitableRing,
    CS: LatticefoldChallengeSet<RqNTT>,
    DP: DecompositionParams,
{
    let mut rng = ark_std::test_rng();
    let input: usize = rng.gen_range(1..101);
    let ccs = get_test_ccs(W, DP::L);
    let log_m = ccs.s;

    let scheme = AjtaiCommitmentScheme::rand(&mut rng);
    let (_, x_ccs, _) = get_test_z_split::<RqNTT>(input);

    let wit = Witness::rand::<_, DP>(&mut rng, WIT_LEN);
    let mut z: Vec<RqNTT> = Vec::with_capacity(x_ccs.len() + WIT_LEN + 1);

    z.extend_from_slice(&x_ccs);
    z.push(RqNTT::one());
    z.extend_from_slice(&wit.w_ccs);

    let cm: Commitment<4, RqNTT> = scheme.commit_ntt(&wit.f).unwrap();

    let r: Vec<RqNTT> = (0..log_m).map(|_| RqNTT::rand(&mut rng)).collect();
    let Mz_mles: Vec<DenseMultilinearExtension<RqNTT>> = ccs
        .M
        .iter()
        .map(|M| DenseMultilinearExtension::from_slice(log_m, &mat_vec_mul(M, &z).unwrap()))
        .collect();

    let u = compute_u(&Mz_mles, &r).unwrap();

    let v = evaluate_mles::<RqNTT, &DenseMultilinearExtension<RqNTT>, _, DecompositionError>(
        &to_mles::<_, _, DecompositionError>(log_m, &wit.f_hat).unwrap(),
        &r,
    )
    .unwrap();

    let lcccs = LCCCS {
        r,
        v,
        cm,
        u,
        x_w: x_ccs,
        h: RqNTT::one(),
    };

    (
        lcccs,
        PoseidonTranscript::<RqNTT, CS>::default(),
        PoseidonTranscript::<RqNTT, CS>::default(),
        ccs,
        wit,
        scheme,
    )
}

#[test]
fn test_recompose_commitment() {
    type CS = GoldilocksChallengeSet;
    type RqNTT = GoldilocksRingNTT;
    type DP = GoldilocksDP;
    type T = PoseidonTranscript<RqNTT, CS>;
    type Verifier = LFDecompositionVerifier<RqNTT, T>;
    const WIT_LEN: usize = 4;
    const W: usize = WIT_LEN * DP::L;
    const C: usize = 4;

    let (lcccs, _, mut prover_transcript, ccs, wit, scheme) =
        generate_decomposition_args::<RqNTT, CS, DP, WIT_LEN, W>();

    let (_, _, proof) = LFDecompositionProver::<_, T>::prove::<W, C, DP>(
        &lcccs,
        &wit,
        &mut prover_transcript,
        &ccs,
        &scheme,
    )
    .unwrap();

    let b_s = Verifier::calculate_b_s::<DP>();

    let should_equal_y0 =
        Verifier::recompose_commitment::<C>(&proof.y_s, &b_s).expect("Recomposing proof failed");

    assert_eq!(should_equal_y0, lcccs.cm);
}

#[test]
fn test_recompose_u() {
    type CS = StarkChallengeSet;
    type RqNTT = StarkRingNTT;
    type DP = StarkDP;
    type T = PoseidonTranscript<RqNTT, CS>;
    type Verifier = LFDecompositionVerifier<RqNTT, T>;
    const WIT_LEN: usize = 4;
    const W: usize = WIT_LEN * DP::L;
    const C: usize = 4;

    let (lcccs, _, mut prover_transcript, ccs, wit, scheme) =
        generate_decomposition_args::<RqNTT, CS, DP, WIT_LEN, W>();

    let (_, _, proof) = LFDecompositionProver::<_, T>::prove::<W, C, DP>(
        &lcccs,
        &wit,
        &mut prover_transcript,
        &ccs,
        &scheme,
    )
    .unwrap();

    let b_s = Verifier::calculate_b_s::<DP>();

    let should_equal_u0 =
        Verifier::recompose_u(&proof.u_s, &b_s).expect("Recomposing proof failed");

    assert_eq!(should_equal_u0, lcccs.u);
}

#[test]
fn test_recompose_v() {
    type CS = BabyBearChallengeSet;
    type RqNTT = BabyBearRingNTT;
    type DP = BabyBearDP;
    type T = PoseidonTranscript<RqNTT, CS>;
    type Verifier = LFDecompositionVerifier<RqNTT, T>;
    const WIT_LEN: usize = 4;
    const W: usize = WIT_LEN * DP::L;
    const C: usize = 4;

    let (lcccs, _, mut prover_transcript, ccs, wit, scheme) =
        generate_decomposition_args::<RqNTT, CS, DP, WIT_LEN, W>();

    let (_, _, proof) = LFDecompositionProver::<_, T>::prove::<W, C, DP>(
        &lcccs,
        &wit,
        &mut prover_transcript,
        &ccs,
        &scheme,
    )
    .unwrap();

    let b_s = Verifier::calculate_b_s::<DP>();

    for (row, &cm_i_value) in lcccs.v.iter().enumerate() {
        let should_equal_v0 = Verifier::recompose_v(&proof.v_s, &b_s, row);

        assert_eq!(should_equal_v0, cm_i_value);
    }
}

#[test]
fn test_recompose_xw_and_h() {
    type CS = GoldilocksChallengeSet;
    type RqNTT = GoldilocksRingNTT;
    type DP = GoldilocksDP;
    type T = PoseidonTranscript<RqNTT, CS>;
    type Verifier = LFDecompositionVerifier<RqNTT, T>;
    const WIT_LEN: usize = 4;
    const W: usize = WIT_LEN * DP::L;
    const C: usize = 4;

    let (lcccs, _, mut prover_transcript, ccs, wit, scheme) =
        generate_decomposition_args::<RqNTT, CS, DP, WIT_LEN, W>();

    let (_, _, proof) = LFDecompositionProver::<_, T>::prove::<W, C, DP>(
        &lcccs,
        &wit,
        &mut prover_transcript,
        &ccs,
        &scheme,
    )
    .unwrap();

    let b_s = Verifier::calculate_b_s::<DP>();

    let (should_equal_xw, should_equal_h) =
        Verifier::recompose_xw_and_h(&proof.x_s, &b_s).expect("Recomposing proof failed");

    assert_eq!(should_equal_h, lcccs.h);
    assert_eq!(should_equal_xw, lcccs.x_w);
}

#[test]
fn test_verify_full() {
    type CS = StarkChallengeSet;
    type RqNTT = StarkRingNTT;
    type DP = StarkDP;
    type T = PoseidonTranscript<RqNTT, CS>;
    type Verifier = LFDecompositionVerifier<RqNTT, T>;
    const WIT_LEN: usize = 4;
    const W: usize = WIT_LEN * DP::L;
    const C: usize = 4;

    let (lcccs, mut verifier_transcript, mut prover_transcript, ccs, wit, scheme) =
        generate_decomposition_args::<RqNTT, CS, DP, WIT_LEN, W>();

    let (_, _, proof) = LFDecompositionProver::<_, T>::prove::<W, C, DP>(
        &lcccs,
        &wit,
        &mut prover_transcript,
        &ccs,
        &scheme,
    )
    .unwrap();

    let _ = Verifier::verify::<C, DP>(&lcccs, &proof, &mut verifier_transcript, &ccs)
        .expect("Failed to verify decomposition proof");
}

#[test]
fn test_verify_invalid_proof() {
    type CS = GoldilocksChallengeSet;
    type RqNTT = GoldilocksRingNTT;
    type DP = GoldilocksDP;
    type T = PoseidonTranscript<RqNTT, CS>;
    type Verifier = LFDecompositionVerifier<RqNTT, T>;
    const WIT_LEN: usize = 4;
    const W: usize = WIT_LEN * DP::L;
    const C: usize = 4;

    let (lcccs, mut verifier_transcript, mut prover_transcript, ccs, wit, scheme) =
        generate_decomposition_args::<RqNTT, CS, DP, WIT_LEN, W>();

    let (_, _, mut proof) = LFDecompositionProver::<_, T>::prove::<W, C, DP>(
        &lcccs,
        &wit,
        &mut prover_transcript,
        &ccs,
        &scheme,
    )
    .unwrap();

    // Make proof components have mismatched lengths
    if !proof.v_s.is_empty() {
        proof.v_s[0][0] += RqNTT::one();
    }

    let result = Verifier::verify::<C, DP>(&lcccs, &proof, &mut verifier_transcript, &ccs);
    assert!(
        result.is_err(),
        "Verification should fail with mismatched proof component lengths"
    );
}
