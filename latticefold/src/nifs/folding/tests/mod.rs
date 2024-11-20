use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Compress, Validate};
use ark_std::io::Cursor;
use rand::thread_rng;

use crate::nifs::folding::FoldingProof;
use crate::{
    arith::{r1cs::get_test_z_split, tests::get_test_ccs, Witness, CCCS},
    commitment::AjtaiCommitmentScheme,
    decomposition_parameters::{test_params::DPL1, DecompositionParams},
    nifs::{
        decomposition::{
            DecompositionProver, DecompositionVerifier, LFDecompositionProver,
            LFDecompositionVerifier,
        },
        folding::{FoldingProver, FoldingVerifier, LFFoldingProver, LFFoldingVerifier},
        linearization::{
            LFLinearizationProver, LFLinearizationVerifier, LinearizationProver,
            LinearizationVerifier,
        },
    },
    transcript::poseidon::PoseidonTranscript,
};
use cyclotomic_rings::rings::{StarkChallengeSet, StarkRingNTT};

// Boilerplate code to generate values needed for testing
type R = StarkRingNTT;
type CS = StarkChallengeSet;
type T = PoseidonTranscript<StarkRingNTT, CS>;

#[test]
fn test_folding() {
    const WIT_LEN: usize = 4; // 4 is the length of witness in this (Vitalik's) example
    const W: usize = WIT_LEN * DPL1::L; // the number of columns of the Ajtai matrix

    let ccs = get_test_ccs::<R>(W);
    let (_, x_ccs, w_ccs) = get_test_z_split::<R>(3);
    let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());
    let wit: Witness<R> = Witness::from_w_ccs::<DPL1>(w_ccs);
    let cm_i: CCCS<4, R> = CCCS {
        cm: wit.commit::<4, 4, DPL1>(&scheme).unwrap(),
        x_ccs,
    };

    let mut prover_transcript = PoseidonTranscript::<R, CS>::default();
    let mut verifier_transcript = PoseidonTranscript::<R, CS>::default();

    let (_, linearization_proof) =
        LFLinearizationProver::<_, T>::prove(&cm_i, &wit, &mut prover_transcript, &ccs)
            .unwrap();

    let lcccs = LFLinearizationVerifier::<_, PoseidonTranscript<R, CS>>::verify(
        &cm_i,
        &linearization_proof,
        &mut verifier_transcript,
        &ccs,
    )
        .unwrap();

    let (_, vec_wit, decomposition_proof) = LFDecompositionProver::<_, T>::prove::<4, 4, DPL1>(
        &lcccs,
        &wit,
        &mut prover_transcript,
        &ccs,
        &scheme,
    )
        .unwrap();

    let vec_lcccs = LFDecompositionVerifier::<_, T>::verify::<4, DPL1>(
        &lcccs,
        &decomposition_proof,
        &mut verifier_transcript,
        &ccs,
    )
        .unwrap();
    let (lcccs, wit_s) = {
        let mut lcccs = vec_lcccs.clone();
        let mut lcccs_r = vec_lcccs;
        lcccs.append(&mut lcccs_r);

        let mut wit_s = vec_wit.clone();
        let mut wit_s_r = vec_wit;
        wit_s.append(&mut wit_s_r);

        (lcccs, wit_s)
    };
    let (lcccs_prover, _, folding_proof) =
        LFFoldingProver::<_, T>::prove::<4, DPL1>(&lcccs, &wit_s, &mut prover_transcript, &ccs)
            .unwrap();

    let lcccs_verifier = LFFoldingVerifier::<_, T>::verify::<4, DPL1>(
        &lcccs,
        &folding_proof,
        &mut verifier_transcript,
        &ccs,
    )
        .unwrap();

    assert_eq!(lcccs_prover, lcccs_verifier);
}

#[test]
fn test_failing_folding_prover() {
    const WIT_LEN: usize = 4; // 4 is the length of witness in this (Vitalik's) example
    const W: usize = WIT_LEN * DPL1::L; // the number of columns of the Ajtai matrix

    let ccs = get_test_ccs::<R>(W);
    let (_, x_ccs, w_ccs) = get_test_z_split::<R>(3);
    let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());
    let wit: Witness<R> = Witness::from_w_ccs::<DPL1>(w_ccs.clone());
    let cm_i: CCCS<4, R> = CCCS {
        cm: wit.commit::<4, 4, DPL1>(&scheme).unwrap(),
        x_ccs,
    };

    let mut prover_transcript = PoseidonTranscript::<R, CS>::default();
    let mut verifier_transcript = PoseidonTranscript::<R, CS>::default();

    let (_, linearization_proof) =
        LFLinearizationProver::<_, T>::prove(&cm_i, &wit, &mut prover_transcript, &ccs)
            .unwrap();

    let lcccs = LFLinearizationVerifier::<_, PoseidonTranscript<R, CS>>::verify(
        &cm_i,
        &linearization_proof,
        &mut verifier_transcript,
        &ccs,
    )
        .unwrap();

    let (_, mut vec_wit, decomposition_proof) =
        LFDecompositionProver::<_, T>::prove::<4, 4, DPL1>(
            &lcccs,
            &wit,
            &mut prover_transcript,
            &ccs,
            &scheme,
        )
            .unwrap();

    let vec_lcccs = LFDecompositionVerifier::<_, T>::verify::<4, DPL1>(
        &lcccs,
        &decomposition_proof,
        &mut verifier_transcript,
        &ccs,
    )
        .unwrap();

    vec_wit[0] = Witness::<R>::from_w_ccs::<DPL1>(w_ccs);

    let res = LFFoldingProver::<_, T>::prove::<4, DPL1>(
        &vec_lcccs,
        &vec_wit,
        &mut prover_transcript,
        &ccs,
    );

    assert!(res.is_err())
}

#[test]
fn test_folding_proof_serialization() {
    const WIT_LEN: usize = 4; // 4 is the length of witness in this (Vitalik's) example
    const W: usize = WIT_LEN * DPL1::L; // the number of columns of the Ajtai matrix

    let ccs = get_test_ccs::<R>(W);
    let (_, x_ccs, w_ccs) = get_test_z_split::<R>(3);
    let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());
    let wit: Witness<R> = Witness::from_w_ccs::<DPL1>(w_ccs);
    let cm_i: CCCS<4, R> = CCCS {
        cm: wit.commit::<4, 4, DPL1>(&scheme).unwrap(),
        x_ccs,
    };

    let mut prover_transcript = PoseidonTranscript::<R, CS>::default();
    let mut verifier_transcript = PoseidonTranscript::<R, CS>::default();

    let (_, linearization_proof) =
        LFLinearizationProver::<_, T>::prove(&cm_i, &wit, &mut prover_transcript, &ccs)
            .unwrap();

    let lcccs = LFLinearizationVerifier::<_, PoseidonTranscript<R, CS>>::verify(
        &cm_i,
        &linearization_proof,
        &mut verifier_transcript,
        &ccs,
    )
        .unwrap();

    let (_, vec_wit, decomposition_proof) = LFDecompositionProver::<_, T>::prove::<4, 4, DPL1>(
        &lcccs,
        &wit,
        &mut prover_transcript,
        &ccs,
        &scheme,
    )
        .unwrap();

    let vec_lcccs = LFDecompositionVerifier::<_, T>::verify::<4, DPL1>(
        &lcccs,
        &decomposition_proof,
        &mut verifier_transcript,
        &ccs,
    )
        .unwrap();
    let (lcccs, wit_s) = {
        let mut lcccs = vec_lcccs.clone();
        let mut lcccs_r = vec_lcccs;
        lcccs.append(&mut lcccs_r);

        let mut wit_s = vec_wit.clone();
        let mut wit_s_r = vec_wit;
        wit_s.append(&mut wit_s_r);

        (lcccs, wit_s)
    };
    let (_, _, folding_proof) =
        LFFoldingProver::<_, T>::prove::<4, DPL1>(&lcccs, &wit_s, &mut prover_transcript, &ccs)
            .unwrap();

    let mut serialized = Vec::new();
    folding_proof
        .serialize_with_mode(&mut serialized, Compress::Yes)
        .expect("Failed to serialize proof");

    let mut cursor = Cursor::new(&serialized);
    assert_eq!(
        folding_proof,
        FoldingProof::deserialize_with_mode(&mut cursor, Compress::Yes, Validate::Yes)
            .expect("Failed to deserialize proof")
    );
}
