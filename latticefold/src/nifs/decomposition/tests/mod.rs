use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Compress, Validate};
use ark_std::io::Cursor;
use cyclotomic_rings::{challenge_set::LatticefoldChallengeSet, rings::SuitableRing};
use rand::thread_rng;

use crate::{
    arith::{r1cs::get_test_z_split, tests::get_test_ccs, Witness, CCCS, CCS},
    commitment::AjtaiCommitmentScheme,
    decomposition_parameters::{test_params::PP, DecompositionParams},
    nifs::{
        decomposition::{
            DecompositionProver, DecompositionVerifier, LFDecompositionProver,
            LFDecompositionVerifier,
        },
        linearization::{
            LFLinearizationProver, LFLinearizationVerifier, LinearizationProof,
            LinearizationProver, LinearizationVerifier,
        },
    },
    transcript::poseidon::PoseidonTranscript,
};

const WIT_LEN: usize = 4;
const W: usize = WIT_LEN * PP::L;

fn generate_decomposition_proof<RqNTT, CS>() -> (
    LinearizationProof<RqNTT>,
    CCCS<4, RqNTT>,
    PoseidonTranscript<RqNTT, CS>,
    PoseidonTranscript<RqNTT, CS>,
    CCS<RqNTT>,
    Witness<RqNTT>,
    AjtaiCommitmentScheme<4, W, RqNTT>,
)
where
    RqNTT: SuitableRing,
    CS: LatticefoldChallengeSet<RqNTT>,
{
    let ccs = get_test_ccs::<RqNTT>(W);
    let (_, x_ccs, w_ccs) = get_test_z_split::<RqNTT>(3);
    let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());
    let wit: Witness<RqNTT> = Witness::from_w_ccs::<PP>(w_ccs);
    let cm_i: CCCS<4, RqNTT> = CCCS {
        cm: wit.commit::<4, W, PP>(&scheme).unwrap(),
        x_ccs,
    };

    let mut prover_transcript = PoseidonTranscript::<RqNTT, CS>::default();
    let verifier_transcript = PoseidonTranscript::<RqNTT, CS>::default();

    let (_, linearization_proof) =
        LFLinearizationProver::<_, PoseidonTranscript<RqNTT, CS>>::prove(
            &cm_i,
            &wit,
            &mut prover_transcript,
            &ccs,
        )
        .unwrap();

    (
        linearization_proof,
        cm_i,
        verifier_transcript,
        prover_transcript,
        ccs,
        wit,
        scheme,
    )
}

fn test_decomposition<RqNTT, CS>()
where
    RqNTT: SuitableRing,
    CS: LatticefoldChallengeSet<RqNTT>,
{
    let (
        linearization_proof,
        cm_i,
        mut verifier_transcript,
        mut prover_transcript,
        ccs,
        wit,
        scheme,
    ) = generate_decomposition_proof::<RqNTT, CS>();

    let lcccs = LFLinearizationVerifier::<_, PoseidonTranscript<RqNTT, CS>>::verify(
        &cm_i,
        &linearization_proof,
        &mut verifier_transcript,
        &ccs,
    )
    .unwrap();

    let (_, _, decomposition_proof) =
        LFDecompositionProver::<_, PoseidonTranscript<RqNTT, CS>>::prove::<W, 4, PP>(
            &lcccs,
            &wit,
            &mut prover_transcript,
            &ccs,
            &scheme,
        )
        .unwrap();

    let res = LFDecompositionVerifier::<_, PoseidonTranscript<RqNTT, CS>>::verify::<4, PP>(
        &lcccs,
        &decomposition_proof,
        &mut verifier_transcript,
        &ccs,
    );

    assert!(res.is_ok());
}

fn test_decomposition_proof_serialization<RqNTT, CS>()
where
    RqNTT: SuitableRing,
    CS: LatticefoldChallengeSet<RqNTT>,
{
    let proof = generate_decomposition_proof::<RqNTT, CS>().0;

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

mod pow2 {
    use cyclotomic_rings::challenge_set::BinarySmallSet;
    use lattirust_ring::cyclotomic_ring::models::pow2_debug::Pow2CyclotomicPolyRingNTT;

    const Q: u64 = 17;
    const N: usize = 8;
    type RqNTT = Pow2CyclotomicPolyRingNTT<Q, N>;
    type CS = BinarySmallSet<Q, N>;

    #[test]
    fn test_decomposition() {
        super::test_decomposition::<RqNTT, CS>();
    }

    #[test]
    fn test_decomposition_proof_serialization() {
        super::test_decomposition_proof_serialization::<RqNTT, CS>();
    }
}

mod stark {
    use crate::arith::r1cs::get_test_dummy_z_split;
    use crate::arith::tests::get_test_dummy_ccs;
    use crate::arith::{Witness, CCCS};
    use crate::commitment::AjtaiCommitmentScheme;
    use crate::decomposition_parameters::{test_params::PP_STARK, DecompositionParams};
    use crate::nifs::linearization::{
        LFLinearizationProver, LFLinearizationVerifier, LinearizationProver, LinearizationVerifier,
    };
    use crate::transcript::poseidon::PoseidonTranscript;
    use crate::utils::security_check::check_ring_modulus_128_bits_security;
    use cyclotomic_rings::rings::StarkChallengeSet;
    use lattirust_ring::cyclotomic_ring::models::stark_prime::RqNTT;
    use num_bigint::BigUint;
    use rand::thread_rng;

    type CS = StarkChallengeSet;

    #[test]
    fn test_decomposition() {
        super::test_decomposition::<RqNTT, CS>();
    }

    #[test]
    fn test_decomposition_proof_serialization() {
        super::test_decomposition_proof_serialization::<RqNTT, CS>();
    }

    #[test]
    fn test_dummy_decomposition() {
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

        // Make bound and security checks
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

mod goldilocks {
    use cyclotomic_rings::rings::GoldilocksChallengeSet;
    use lattirust_ring::cyclotomic_ring::models::goldilocks::RqNTT;
    type CS = GoldilocksChallengeSet;

    #[test]
    fn test_decomposition() {
        super::test_decomposition::<RqNTT, CS>();
    }

    #[test]
    fn test_decomposition_proof_serialization() {
        super::test_decomposition_proof_serialization::<RqNTT, CS>();
    }
}

mod frog {
    use cyclotomic_rings::rings::FrogChallengeSet;
    use lattirust_ring::cyclotomic_ring::models::frog_ring::RqNTT;
    type CS = FrogChallengeSet;

    #[test]
    fn test_decomposition() {
        super::test_decomposition::<RqNTT, CS>();
    }

    #[test]
    fn test_decomposition_proof_serialization() {
        super::test_decomposition_proof_serialization::<RqNTT, CS>();
    }
}

mod babybear {
    use cyclotomic_rings::rings::BabyBearChallengeSet;
    use lattirust_ring::cyclotomic_ring::models::babybear::RqNTT;
    type CS = BabyBearChallengeSet;

    #[test]
    fn test_decomposition() {
        super::test_decomposition::<RqNTT, CS>();
    }

    #[test]
    fn test_decomposition_proof_serialization() {
        super::test_decomposition_proof_serialization::<RqNTT, CS>();
    }
}
