use super::*;
use crate::{
    arith::{r1cs::get_test_z_split, tests::get_test_ccs},
    commitment::AjtaiCommitmentScheme,
    decomposition_parameters::{test_params::PP, DecompositionParams},
    transcript::poseidon::PoseidonTranscript,
};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Compress, Validate};
use ark_std::io::Cursor;
use cyclotomic_rings::challenge_set::LatticefoldChallengeSet;
use cyclotomic_rings::rings::{
    BabyBearChallengeSet, FrogChallengeSet, GoldilocksChallengeSet, StarkChallengeSet,
};
use lattirust_poly::polynomials::{build_eq_x_r, VirtualPolynomial};
use lattirust_ring::OverField;
use std::sync::Arc;

fn test_compute_ui<RqNTT: OverField>() {
    let mut mles = Vec::with_capacity(10);
    let mut rng = ark_std::test_rng();

    for _i in 0..10 {
        let evals: Vec<RqNTT> = (0..8).map(|_| RqNTT::rand(&mut rng)).collect();
        mles.push(DenseMultilinearExtension::from_evaluations_slice(3, &evals))
    }

    for b in 0..8_u8 {
        let us: Vec<RqNTT> = compute_u(
            &mles,
            &[
                (b & 0x01).into(),
                ((b & 0x2) >> 1).into(),
                ((b & 0x4) >> 2).into(),
            ],
        )
        .unwrap();

        for (i, &u) in us.iter().enumerate() {
            assert_eq!(u, mles[i].evaluations[b.to_le() as usize]);
        }
    }
}

fn test_linearization_polynomial<RqNTT: OverField>() {
    let mut rng = ark_std::test_rng();

    let n_c = 4;
    let n_r = 4;
    let log_m = 2;
    let s_i = 3;

    let mut g = VirtualPolynomial::<RqNTT>::new(log_m);
    let z: Vec<RqNTT> = (0..n_c).map(|_| RqNTT::rand(&mut rng)).collect();
    let c = RqNTT::rand(&mut rng);
    let beta: Vec<RqNTT> = (0..log_m).map(|_| RqNTT::rand(&mut rng)).collect();
    let mut M_z_mles: Vec<DenseMultilinearExtension<RqNTT>> = Vec::with_capacity(s_i);

    for _ in 0..s_i {
        let mut mle = Vec::new();
        for _ in 0..n_r {
            let mut row = Vec::new();
            for _ in 0..n_c {
                let random_value = RqNTT::rand(&mut rng);
                row.push(random_value);
            }
            let row_z = row
                .iter()
                .zip(&z)
                .map(|(&r_i, z_i)| r_i * z_i)
                .sum::<RqNTT>();
            mle.push(row_z);
        }
        M_z_mles.push(DenseMultilinearExtension::from_slice(log_m, &mle));
    }

    let _ = g.add_mle_list(M_z_mles.clone().into_iter().map(Arc::new), c);
    let eq_b_r = build_eq_x_r(&beta).unwrap();
    let _ = g.mul_by_mle(eq_b_r, RqNTT::one());

    let polynomial =
        prepare_lin_sumcheck_polynomial(log_m, &[c], &M_z_mles, &[vec![0, 1, 2]], &beta).unwrap();

    for _ in 0..20 {
        let point: Vec<RqNTT> = (0..log_m).map(|_| RqNTT::rand(&mut rng)).collect();
        assert_eq!(
            g.evaluate(&point).unwrap(),
            polynomial.evaluate(&point).unwrap()
        )
    }
}

use rand::thread_rng;

fn setup_test_environment<RqNTT: SuitableRing>() -> (Witness<RqNTT>, CCCS<4, RqNTT>, CCS<RqNTT>) {
    const WIT_LEN: usize = 4;
    const W: usize = WIT_LEN * PP::L;

    let ccs = get_test_ccs::<RqNTT>(W);
    let (_, x_ccs, w_ccs) = get_test_z_split::<RqNTT>(3);
    let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());

    let wit = Witness::from_w_ccs::<PP>(w_ccs);
    let cm_i = CCCS {
        cm: wit.commit::<4, W, PP>(&scheme).unwrap(),
        x_ccs,
    };

    (wit, cm_i, ccs)
}

fn test_prover_state<RqNTT: SuitableRing, CS: LatticefoldChallengeSet<RqNTT>>() {
    let (wit, cm_i, ccs) = setup_test_environment::<RqNTT>();
    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();

    let state =
        LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::prepare_prover_state(
            &wit,
            &cm_i,
            &mut transcript,
            &ccs,
        )
        .unwrap();

    // Verify state properties
    assert_eq!(state.beta_s.len(), ccs.s);
    assert!(!state.z_ccs.is_empty());
    assert!(!state.Mz_mles.is_empty());

    // Verify Mz_mles dimensions
    for mle in &state.Mz_mles {
        assert_eq!(mle.evaluations.len(), 1 << ccs.s);
    }
}

fn test_sumcheck_generation<RqNTT: SuitableRing, CS: LatticefoldChallengeSet<RqNTT>>() {
    let (wit, cm_i, ccs) = setup_test_environment::<RqNTT>();
    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();

    let state =
        LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::prepare_prover_state(
            &wit,
            &cm_i,
            &mut transcript,
            &ccs,
        )
        .unwrap();

    let (_, point_r) =
        LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::generate_sumcheck_proof(
            &state,
            &mut transcript,
            &ccs,
        )
        .unwrap();

    // Verify dimensions
    assert_eq!(point_r.len(), ccs.s);
}

fn test_verifier_state<RqNTT: SuitableRing, CS: LatticefoldChallengeSet<RqNTT>>() {
    let (wit, cm_i, ccs) = setup_test_environment::<RqNTT>();
    let mut prover_transcript = PoseidonTranscript::<RqNTT, CS>::default();

    let (_, proof) = LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::prove(
        &cm_i,
        &wit,
        &mut prover_transcript,
        &ccs,
    )
    .unwrap();

    let mut verifier_transcript = PoseidonTranscript::<RqNTT, CS>::default();
    let state =
        LFLinearizationVerifier::<RqNTT, PoseidonTranscript<RqNTT, CS>>::prepare_verifier_state(
            &mut verifier_transcript,
            &proof,
            &ccs,
        )
        .unwrap();

    assert_eq!(state.beta_s.len(), ccs.s);
    assert_eq!(state.point_r.len(), ccs.s);
}

fn generate_test_proof<RqNTT: SuitableRing, CS: LatticefoldChallengeSet<RqNTT>>(
) -> LinearizationProof<RqNTT> {
    let (wit, cm_i, ccs) = setup_test_environment::<RqNTT>();
    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();

    let (_, proof) = LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::prove(
        &cm_i,
        &wit,
        &mut transcript,
        &ccs,
    )
    .unwrap();

    proof
}

fn test_proof_serialization<RqNTT: SuitableRing, CS: LatticefoldChallengeSet<RqNTT>>() {
    let proof = generate_test_proof::<RqNTT, CS>();

    let mut serialized = Vec::new();
    proof
        .serialize_with_mode(&mut serialized, Compress::Yes)
        .expect("Failed to serialize proof");

    let mut cursor = Cursor::new(&serialized);
    assert_eq!(
        proof,
        LinearizationProof::deserialize_with_mode(&mut cursor, Compress::Yes, Validate::Yes)
            .expect("Failed to deserialize proof")
    );
}

fn test_evaluation_claim<RqNTT: SuitableRing, CS: LatticefoldChallengeSet<RqNTT>>() {
    let (wit, cm_i, ccs) = setup_test_environment::<RqNTT>();
    let mut prover_transcript = PoseidonTranscript::<RqNTT, CS>::default();

    let (_, proof) = LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::prove(
        &cm_i,
        &wit,
        &mut prover_transcript,
        &ccs,
    )
    .unwrap();

    let mut verifier_transcript = PoseidonTranscript::<RqNTT, CS>::default();
    let state =
        LFLinearizationVerifier::<RqNTT, PoseidonTranscript<RqNTT, CS>>::prepare_verifier_state(
            &mut verifier_transcript,
            &proof,
            &ccs,
        )
        .unwrap();

    assert!(
        LFLinearizationVerifier::<RqNTT, PoseidonTranscript<RqNTT, CS>>::verify_evaluation_claim(
            &state, &proof, &ccs,
        )
        .is_ok()
    );
}

mod tests_stark {
    use super::*;
    use crate::arith::r1cs::get_test_dummy_z_split;
    use crate::arith::tests::get_test_dummy_ccs;
    use crate::decomposition_parameters::test_params::PP_STARK;
    use crate::utils::security_check::check_ring_modulus_128_bits_security;
    use lattirust_ring::cyclotomic_ring::models::stark_prime::RqNTT;
    use num_bigint::BigUint;

    #[test]
    fn test_prover_state() {
        super::test_prover_state::<RqNTT, StarkChallengeSet>();
    }

    #[test]
    fn test_sumcheck() {
        test_sumcheck_generation::<RqNTT, StarkChallengeSet>();
    }

    #[test]
    fn test_verifier_state() {
        super::test_verifier_state::<RqNTT, StarkChallengeSet>();
    }

    #[test]
    fn test_evaluation() {
        test_evaluation_claim::<RqNTT, StarkChallengeSet>();
    }

    #[test]
    fn test_serialization() {
        test_proof_serialization::<RqNTT, StarkChallengeSet>();
    }

    #[test]
    fn test_linearization_polynomial() {
        super::test_linearization_polynomial::<RqNTT>();
    }

    #[test]
    fn test_compute_ui() {
        super::test_compute_ui::<RqNTT>();
    }

    #[test]
    fn test_dummy_linearization() {
        type R = RqNTT;
        type CS = StarkChallengeSet;
        type T = PoseidonTranscript<R, CS>;

        const C: usize = 16;
        const X_LEN: usize = 1;
        const WIT_LEN: usize = 2048;
        const W: usize = WIT_LEN * PP_STARK::L; // the number of columns of the Ajtai matrix
        let r1cs_rows_size = X_LEN + WIT_LEN + 1; // Let's have a square matrix

        let ccs = get_test_dummy_ccs::<R, X_LEN, WIT_LEN, W>(r1cs_rows_size);
        let (_, x_ccs, w_ccs) = get_test_dummy_z_split::<R, X_LEN, WIT_LEN>();
        let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());

        let wit = Witness::from_w_ccs::<PP_STARK>(w_ccs);

        // Make bound and securitty checks
        let witness_within_bound = wit.within_bound(PP_STARK::B);
        let stark_modulus = BigUint::parse_bytes(
            b"3618502788666131000275863779947924135206266826270938552493006944358698582017",
            10,
        )
        .expect("Failed to parse stark_modulus");

        if check_ring_modulus_128_bits_security(
            &stark_modulus,
            C,
            16,
            W,
            PP_STARK::B,
            PP_STARK::L,
            witness_within_bound,
        ) {
            println!(" Bound condition satisfied for 128 bits security");
        } else {
            println!("Bound condition not satisfied for 128 bits security");
        }

        let cm_i = CCCS {
            cm: wit.commit::<C, W, PP_STARK>(&scheme).unwrap(),
            x_ccs,
        };

        let mut transcript = PoseidonTranscript::<R, CS>::default();

        let res = LFLinearizationProver::<_, T>::prove(&cm_i, &wit, &mut transcript, &ccs);

        let mut transcript = PoseidonTranscript::<R, CS>::default();

        let res = LFLinearizationVerifier::<_, PoseidonTranscript<R, CS>>::verify(
            &cm_i,
            &res.expect("Linearization proof generation error").1,
            &mut transcript,
            &ccs,
        );

        res.expect("Linearization Verification error");
    }
}

mod tests_goldilocks {
    use super::*;
    use lattirust_ring::cyclotomic_ring::models::goldilocks::RqNTT;

    #[test]
    fn test_prover_state() {
        super::test_prover_state::<RqNTT, GoldilocksChallengeSet>();
    }

    #[test]
    fn test_sumcheck() {
        test_sumcheck_generation::<RqNTT, GoldilocksChallengeSet>();
    }

    #[test]
    fn test_verifier_state() {
        super::test_verifier_state::<RqNTT, GoldilocksChallengeSet>();
    }

    #[test]
    fn test_evaluation() {
        test_evaluation_claim::<RqNTT, GoldilocksChallengeSet>();
    }

    #[test]
    fn test_serialization() {
        test_proof_serialization::<RqNTT, GoldilocksChallengeSet>();
    }

    #[test]
    fn test_linearization_polynomial() {
        super::test_linearization_polynomial::<
            lattirust_ring::cyclotomic_ring::models::stark_prime::RqNTT,
        >();
    }

    #[test]
    fn test_compute_ui() {
        super::test_compute_ui::<RqNTT>();
    }
}

mod tests_frog {
    use super::*;
    use lattirust_ring::cyclotomic_ring::models::frog_ring::RqNTT;

    #[test]
    fn test_prover_state() {
        super::test_prover_state::<RqNTT, FrogChallengeSet>();
    }

    #[test]
    fn test_sumcheck() {
        test_sumcheck_generation::<RqNTT, FrogChallengeSet>();
    }

    #[test]
    fn test_verifier_state() {
        super::test_verifier_state::<RqNTT, FrogChallengeSet>();
    }

    #[test]
    fn test_evaluation() {
        test_evaluation_claim::<RqNTT, FrogChallengeSet>();
    }

    #[test]
    fn test_serialization() {
        test_proof_serialization::<RqNTT, FrogChallengeSet>();
    }

    #[test]
    fn test_linearization_polynomial() {
        super::test_linearization_polynomial::<
            lattirust_ring::cyclotomic_ring::models::stark_prime::RqNTT,
        >();
    }

    #[test]
    fn test_compute_ui() {
        super::test_compute_ui::<RqNTT>();
    }
}

mod tests_babybear {
    use super::*;
    use lattirust_ring::cyclotomic_ring::models::babybear::RqNTT;

    #[test]
    fn test_prover_state() {
        super::test_prover_state::<RqNTT, BabyBearChallengeSet>();
    }

    #[test]
    fn test_sumcheck() {
        test_sumcheck_generation::<RqNTT, BabyBearChallengeSet>();
    }

    #[test]
    fn test_verifier_state() {
        super::test_verifier_state::<RqNTT, BabyBearChallengeSet>();
    }

    #[test]
    fn test_evaluation() {
        test_evaluation_claim::<RqNTT, BabyBearChallengeSet>();
    }

    #[test]
    fn test_serialization() {
        test_proof_serialization::<RqNTT, BabyBearChallengeSet>();
    }

    #[test]
    fn test_linearization_polynomial() {
        super::test_linearization_polynomial::<
            lattirust_ring::cyclotomic_ring::models::stark_prime::RqNTT,
        >();
    }

    #[test]
    fn test_compute_ui() {
        super::test_compute_ui::<RqNTT>();
    }
}
