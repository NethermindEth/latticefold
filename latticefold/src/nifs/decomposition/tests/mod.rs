use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Compress, Validate};
use ark_std::io::Cursor;
use cyclotomic_rings::{challenge_set::LatticefoldChallengeSet, rings::SuitableRing};
use lattirust_ring::{
    balanced_decomposition::{decompose_balanced_vec, recompose},
    cyclotomic_ring::CRT,
    PolyRing,
};
use rand::{rngs::ThreadRng, thread_rng, Rng};

use crate::{
    arith::{r1cs::get_test_z_split, tests::get_test_ccs, Witness, CCCS, CCS},
    commitment::AjtaiCommitmentScheme,
    decomposition_parameters::{test_params::DP, DecompositionParams},
    nifs::{
        decomposition::{
            utils::{decompose_B_vec_into_k_vec, decompose_big_vec_into_k_vec_and_compose_back},
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

fn draw_ring_bellow_bound<RqPoly, const B: u128>(rng: &mut ThreadRng) -> RqPoly
where
    RqPoly: PolyRing + CRT,
{
    let degree = <RqPoly as PolyRing>::dimension();
    let mut coeffs = Vec::with_capacity(degree);
    for _ in 0..degree {
        let random_coeff = rng.gen_range(0..B);
        coeffs.push(<RqPoly as PolyRing>::BaseRing::from(random_coeff));
    }
    RqPoly::from(coeffs)
}

const WIT_LEN: usize = 4;
const W: usize = WIT_LEN * DP::L;

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
    let wit: Witness<RqNTT> = Witness::from_w_ccs::<DP>(w_ccs);
    let cm_i: CCCS<4, RqNTT> = CCCS {
        cm: wit.commit::<4, W, DP>(&scheme).unwrap(),
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
        LFDecompositionProver::<_, PoseidonTranscript<RqNTT, CS>>::prove::<W, 4, DP>(
            &lcccs,
            &wit,
            &mut prover_transcript,
            &ccs,
            &scheme,
        )
        .unwrap();

    let res = LFDecompositionVerifier::<_, PoseidonTranscript<RqNTT, CS>>::verify::<4, DP>(
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

fn test_decompose_B_vec_into_k_vec<RqNTT, RqPoly>()
where
    RqNTT: SuitableRing<CoefficientRepresentation = RqPoly>,
    RqPoly: PolyRing + CRT,
{
    // Create a test vector
    const N: usize = 32;
    let mut rng = thread_rng();
    let test_vector: Vec<RqPoly> = (0..N)
        .map(|_| draw_ring_bellow_bound::<RqPoly, { DP::B }>(&mut rng))
        .collect();

    // Call the function
    let decomposed = decompose_B_vec_into_k_vec::<RqNTT, DP>(&test_vector);

    // Check that we get K vectors back from the decomposition
    assert_eq!(
        decomposed.len(),
        DP::K,
        "Decomposition should output K={} vectors",
        DP::K
    );

    // Check the length of each inner vector
    for vec in &decomposed {
        assert_eq!(vec.len(), N);
    }

    // Check that the decomposition is correct
    for i in 0..N {
        let decomp_i = decomposed.iter().map(|d_j| d_j[i]).collect::<Vec<_>>();
        assert_eq!(
            test_vector[i],
            recompose(&decomp_i, RqPoly::from(DP::B_SMALL as u128))
        );
    }
}

fn recompose_from_k_vec_to_big_vec<NTT: SuitableRing>(
    k_vecs: &[Vec<NTT>],
) -> Vec<NTT::CoefficientRepresentation> {
    let decomposed_in_b: Vec<Vec<NTT::CoefficientRepresentation>> = k_vecs
        .iter()
        .map(|vec| {
            vec.iter()
                .flat_map(|&x| decompose_balanced_vec(&[x.icrt()], DP::B, DP::L))
                .flatten()
                .collect()
        })
        .collect();

    // Transpose the decomposed vectors
    let mut transposed = vec![vec![]; decomposed_in_b[0].len()];
    for row in &decomposed_in_b {
        for (j, &val) in row.iter().enumerate() {
            transposed[j].push(val);
        }
    }

    // Recompose first with B_SMALL, then with B
    transposed
        .iter()
        .map(|vec| {
            recompose(
                vec,
                NTT::CoefficientRepresentation::from(DP::B_SMALL as u128),
            )
        })
        .collect::<Vec<_>>()
        .chunks(DP::L)
        .map(|chunk| recompose(chunk, NTT::CoefficientRepresentation::from(DP::B)))
        .collect()
}

fn test_decompose_big_vec_into_k_vec_and_compose_back<RqNTT, RqPoly>()
where
    RqNTT: SuitableRing<CoefficientRepresentation = RqPoly>,
    RqPoly: PolyRing + CRT,
    Vec<RqNTT>: FromIterator<<RqPoly as CRT>::CRTForm>,
{
    // Create a test vector
    const N: usize = 32;
    let mut rng = thread_rng();
    let test_vector: Vec<RqNTT> = (0..N)
        .map(|_| draw_ring_bellow_bound::<RqPoly, { DP::B }>(&mut rng).crt())
        .collect();
    let decomposed_and_composed_back =
        decompose_big_vec_into_k_vec_and_compose_back::<RqNTT, DP>(test_vector.clone());
    let restore_decomposed =
        recompose_from_k_vec_to_big_vec::<RqNTT>(&decomposed_and_composed_back);

    // Check each entry matches
    for i in 0..N {
        assert_eq!(
            restore_decomposed[i],
            test_vector[i].icrt(),
            "Mismatch at index {}: decomposed_and_composed_back={}, test_vector={}",
            i,
            restore_decomposed[i],
            test_vector[i].icrt()
        );
    }
}

mod stark {
    use crate::arith::r1cs::get_test_dummy_z_split;
    use crate::arith::tests::get_test_dummy_ccs;
    use crate::arith::{Witness, CCCS};
    use crate::commitment::AjtaiCommitmentScheme;
    use crate::decomposition_parameters::{test_params::StarkDP, DecompositionParams};
    use crate::nifs::linearization::{
        LFLinearizationProver, LFLinearizationVerifier, LinearizationProver, LinearizationVerifier,
    };
    use crate::transcript::poseidon::PoseidonTranscript;
    use crate::utils::security_check::check_ring_modulus_128_bits_security;
    use cyclotomic_rings::rings::StarkChallengeSet;
    use lattirust_ring::cyclotomic_ring::models::stark_prime::{RqNTT, RqPoly};
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
    fn test_decompose_B_vec_into_k_vec() {
        super::test_decompose_B_vec_into_k_vec::<RqNTT, RqPoly>();
    }

    #[test]
    fn test_decompose_big_vec_into_k_vec_and_compose_back() {
        super::test_decompose_big_vec_into_k_vec_and_compose_back::<RqNTT, RqPoly>();
    }

    #[test]
    fn test_dummy_decomposition() {
        type R = RqNTT;
        type CS = StarkChallengeSet;
        type T = PoseidonTranscript<R, CS>;

        const C: usize = 16;
        const X_LEN: usize = 1;
        const WIT_LEN: usize = 2048;
        const W: usize = WIT_LEN * StarkDP::L; // the number of columns of the Ajtai matrix
        let r1cs_rows_size = X_LEN + WIT_LEN + 1; // Let's have a square matrix
        let ccs = get_test_dummy_ccs::<R, X_LEN, WIT_LEN, W>(r1cs_rows_size);
        let (_, x_ccs, w_ccs) = get_test_dummy_z_split::<R, X_LEN, WIT_LEN>();
        let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());
        let wit = Witness::from_w_ccs::<StarkDP>(w_ccs);

        // Make bound and security checks
        let witness_within_bound = wit.within_bound(StarkDP::B);
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
            StarkDP::B,
            StarkDP::L,
            witness_within_bound,
        ) {
            println!(" Bound condition satisfied for 128 bits security");
        } else {
            println!("Bound condition not satisfied for 128 bits security");
        }
        let cm_i = CCCS {
            cm: wit.commit::<C, W, StarkDP>(&scheme).unwrap(),
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
    use lattirust_ring::cyclotomic_ring::models::goldilocks::{RqNTT, RqPoly};
    type CS = GoldilocksChallengeSet;

    #[test]
    fn test_decomposition() {
        super::test_decomposition::<RqNTT, CS>();
    }

    #[test]
    fn test_decomposition_proof_serialization() {
        super::test_decomposition_proof_serialization::<RqNTT, CS>();
    }

    #[test]
    fn test_decompose_B_vec_into_k_vec() {
        super::test_decompose_B_vec_into_k_vec::<RqNTT, RqPoly>();
    }

    #[test]
    fn test_decompose_big_vec_into_k_vec_and_compose_back() {
        super::test_decompose_big_vec_into_k_vec_and_compose_back::<RqNTT, RqPoly>();
    }
}

mod frog {
    use cyclotomic_rings::rings::FrogChallengeSet;
    use lattirust_ring::cyclotomic_ring::models::frog_ring::{RqNTT, RqPoly};
    type CS = FrogChallengeSet;

    #[test]
    fn test_decomposition() {
        super::test_decomposition::<RqNTT, CS>();
    }

    #[test]
    fn test_decomposition_proof_serialization() {
        super::test_decomposition_proof_serialization::<RqNTT, CS>();
    }

    #[test]
    fn test_decompose_B_vec_into_k_vec() {
        super::test_decompose_B_vec_into_k_vec::<RqNTT, RqPoly>();
    }

    #[test]
    fn test_decompose_big_vec_into_k_vec_and_compose_back() {
        super::test_decompose_big_vec_into_k_vec_and_compose_back::<RqNTT, RqPoly>();
    }
}

mod babybear {
    use cyclotomic_rings::rings::BabyBearChallengeSet;
    use lattirust_ring::cyclotomic_ring::models::babybear::{RqNTT, RqPoly};
    type CS = BabyBearChallengeSet;

    #[test]
    fn test_decomposition() {
        super::test_decomposition::<RqNTT, CS>();
    }

    #[test]
    fn test_decomposition_proof_serialization() {
        super::test_decomposition_proof_serialization::<RqNTT, CS>();
    }

    #[test]
    fn test_decompose_B_vec_into_k_vec() {
        super::test_decompose_B_vec_into_k_vec::<RqNTT, RqPoly>();
    }

    #[test]
    fn test_decompose_big_vec_into_k_vec_and_compose_back() {
        super::test_decompose_big_vec_into_k_vec_and_compose_back::<RqNTT, RqPoly>();
    }
}
