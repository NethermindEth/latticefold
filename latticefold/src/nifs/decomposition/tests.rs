macro_rules! generate_decomposition_tests {
    ( $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        use crate::{
            arith::{r1cs::get_test_z_split, tests::get_test_ccs, Witness, CCCS},
            commitment::AjtaiCommitmentScheme,
            nifs::{
                decomposition::{
                    structs::{DecompositionProof, LFDecompositionProver, LFDecompositionVerifier},
                    utils::{
                        decompose_B_vec_into_k_vec, decompose_big_vec_into_k_vec_and_compose_back,
                    },
                    DecompositionParams, DecompositionProver, DecompositionVerifier,
                },
                linearization::{
                    LFLinearizationProver, LFLinearizationVerifier, LinearizationProver,
                    LinearizationVerifier,
                },
            },
            transcript::poseidon::PoseidonTranscript,
        };
        use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Compress, Validate};
        use ark_std::io::Cursor;
        use cyclotomic_rings::rings::SuitableRing;
        use lattirust_ring::{
            balanced_decomposition::{decompose_balanced_vec, recompose},
            cyclotomic_ring::{CRT, ICRT},
            PolyRing,
        };
        use rand::thread_rng;
        use rand::{rngs::ThreadRng, Rng};

        type T = PoseidonTranscript<RqNTT, CS>;

        #[derive(Clone)]
        struct PP;
        impl DecompositionParams for PP {
            const B: u128 = $b;
            const L: usize = $l;
            const B_SMALL: usize = $b_small;
            const K: usize = $k;
        }

        fn draw_ring_bellow_bound<const B: u128>(rng: &mut ThreadRng) -> RqPoly {
            let degree = <RqPoly as PolyRing>::dimension();
            let mut coeffs = Vec::with_capacity(degree);
            for _ in 0..degree {
                let random_coeff = rng.gen_range(0..B);
                coeffs.push(<RqPoly as PolyRing>::BaseRing::from(random_coeff));
            }
            RqPoly::from(coeffs)
        }

        #[test]
        fn test_decomposition() {
            const WIT_LEN: usize = 4;
            const W: usize = WIT_LEN * PP::L;

            let ccs = get_test_ccs::<RqNTT>(W);
            let (_, x_ccs, w_ccs) = get_test_z_split::<RqNTT>(3);
            let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());
            let wit: Witness<RqNTT> = Witness::from_w_ccs::<PP>(&w_ccs);
            let cm_i: CCCS<4, RqNTT> = CCCS {
                cm: wit.commit::<4, W, PP>(&scheme).unwrap(),
                x_ccs,
            };

            let mut prover_transcript = PoseidonTranscript::<RqNTT, CS>::default();
            let mut verifier_transcript = PoseidonTranscript::<RqNTT, CS>::default();

            let (_, linearization_proof) =
                LFLinearizationProver::<_, T>::prove(&cm_i, &wit, &mut prover_transcript, &ccs)
                    .unwrap();

            let lcccs = LFLinearizationVerifier::<_, PoseidonTranscript<RqNTT, CS>>::verify(
                &cm_i,
                &linearization_proof,
                &mut verifier_transcript,
                &ccs,
            )
            .unwrap();

            let (_, _, decomposition_proof) = LFDecompositionProver::<_, T>::prove::<W, 4, PP>(
                &lcccs,
                &wit,
                &mut prover_transcript,
                &ccs,
                &scheme,
            )
            .unwrap();

            let res = LFDecompositionVerifier::<_, T>::verify::<4, PP>(
                &lcccs,
                &decomposition_proof,
                &mut verifier_transcript,
                &ccs,
            );

            assert!(res.is_ok());
        }

        #[test]
        fn test_decompose_B_vec_into_k_vec() {
            // Create a test vector
            const N: usize = 32;
            let mut rng = thread_rng();
            let test_vector: Vec<RqPoly> = (0..N)
                .map(|_| draw_ring_bellow_bound::<{ PP::B }>(&mut rng))
                .collect();

            // Call the function
            let decomposed = decompose_B_vec_into_k_vec::<RqNTT, PP>(&test_vector);

            // Check that we get K vectors back from the decomposition
            assert_eq!(
                decomposed.len(),
                PP::K,
                "Decomposition should output K={} vectors",
                PP::K
            );

            // Check the length of each inner vector
            for vec in &decomposed {
                assert_eq!(vec.len(), N);
            }

            // Check that the decomposition is correct
            for i in 0..N {
                let decomp_i = decomposed.iter().map(|d_j| d_j[i]).collect::<Vec<_>>();
                assert_eq!(test_vector[i], recompose(&decomp_i, PP::B_SMALL as u128));
            }
        }

        #[test]
        fn test_decompose_big_vec_into_k_vec_and_compose_back() {
            // Create a test vector
            const N: usize = 32;
            let mut rng = thread_rng();
            let test_vector: Vec<RqNTT> = (0..N)
                .map(|_| draw_ring_bellow_bound::<{ PP::B }>(&mut rng).crt())
                .collect();
            let decomposed_and_composed_back =
                decompose_big_vec_into_k_vec_and_compose_back::<RqNTT, PP>(&test_vector);
            let restore_decomposed =
                recompose_from_k_vec_to_big_vec::<RqNTT, PP>(&decomposed_and_composed_back);

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

        fn recompose_from_k_vec_to_big_vec<NTT: SuitableRing, DP: DecompositionParams>(
            k_vecs: &[Vec<NTT>],
        ) -> Vec<NTT::CoefficientRepresentation> {
            let decomposed_in_b: Vec<Vec<NTT::CoefficientRepresentation>> = k_vecs
                .iter()
                .map(|vec| {
                    vec.iter()
                        .map(|&x| decompose_balanced_vec(&[x.icrt()], DP::B, DP::L))
                        .flatten()
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
                .map(|vec| recompose(vec, DP::B_SMALL as u128))
                .collect::<Vec<_>>()
                .chunks(DP::L)
                .map(|chunk| recompose(chunk, DP::B))
                .collect()
        }

        #[test]
        fn test_decomposition_proof_serialization() {
            const WIT_LEN: usize = 4; // 4 is the length of witness in this (Vitalik's) example
            const W: usize = WIT_LEN * PP::L; // the number of columns of the Ajtai matrix

            let ccs = get_test_ccs::<RqNTT>(W);
            let (_, x_ccs, w_ccs) = get_test_z_split::<RqNTT>(3);
            let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());
            let wit: Witness<RqNTT> = Witness::from_w_ccs::<PP>(&w_ccs);
            let cm_i: CCCS<4, RqNTT> = CCCS {
                cm: wit.commit::<4, W, PP>(&scheme).unwrap(),
                x_ccs,
            };

            let mut prover_transcript = PoseidonTranscript::<RqNTT, CS>::default();
            let mut verifier_transcript = PoseidonTranscript::<RqNTT, CS>::default();

            let (_, linearization_proof) =
                LFLinearizationProver::<_, T>::prove(&cm_i, &wit, &mut prover_transcript, &ccs)
                    .unwrap();

            let lcccs = LFLinearizationVerifier::<_, PoseidonTranscript<RqNTT, CS>>::verify(
                &cm_i,
                &linearization_proof,
                &mut verifier_transcript,
                &ccs,
            )
            .unwrap();

            let (_, _, decomposition_proof) = LFDecompositionProver::<_, T>::prove::<W, 4, PP>(
                &lcccs,
                &wit,
                &mut prover_transcript,
                &ccs,
                &scheme,
            )
            .unwrap();

            let mut serialized = Vec::new();
            decomposition_proof
                .serialize_with_mode(&mut serialized, Compress::Yes)
                .expect("Failed to serialize proof");

            let mut cursor = Cursor::new(&serialized);
            assert_eq!(
                decomposition_proof,
                DecompositionProof::deserialize_with_mode(
                    &mut cursor,
                    Compress::Yes,
                    Validate::Yes
                )
                .expect("Failed to deserialize decomposition proof")
            );
        }
    };
}

#[cfg(test)]
mod tests_pow2 {
    use cyclotomic_rings::challenge_set::BinarySmallSet;
    use lattirust_ring::cyclotomic_ring::models::pow2_debug::{
        Pow2CyclotomicPolyRing, Pow2CyclotomicPolyRingNTT,
    };

    const Q: u64 = 17;
    const N: usize = 8;
    type RqNTT = Pow2CyclotomicPolyRingNTT<Q, N>;
    type RqPoly = Pow2CyclotomicPolyRing<Q, N>;
    type CS = BinarySmallSet<Q, N>;
    generate_decomposition_tests!(1024, 2, 2, 10);
}

#[cfg(test)]
mod tests_stark {

    use cyclotomic_rings::rings::StarkChallengeSet;
    use lattirust_ring::cyclotomic_ring::models::stark_prime::{RqNTT, RqPoly};
    use num_bigint::BigUint;

    use crate::{
        arith::{r1cs::get_test_dummy_z_split, tests::get_test_dummy_ccs},
        utils::security_check::{check_ring_modulus_128_bits_security, check_witness_bound},
    };
    type CS = StarkChallengeSet;
    generate_decomposition_tests!(1024, 2, 2, 10);

    #[test]
    fn test_dummy_decomposition() {
        type R = RqNTT;
        type CS = StarkChallengeSet;
        type T = PoseidonTranscript<R, CS>;

        #[derive(Clone)]
        struct PP;
        impl DecompositionParams for PP {
            const B: u128 = 10485760000;
            const L: usize = 8;
            const B_SMALL: usize = 320;
            const K: usize = 4;
        }

        const C: usize = 16;
        const X_LEN: usize = 1;
        const WIT_LEN: usize = 2048;
        const W: usize = WIT_LEN * PP::L; // the number of columns of the Ajtai matrix
        let r1cs_rows_size = X_LEN + WIT_LEN + 1; // Let's have a square matrix

        let ccs = get_test_dummy_ccs::<R, X_LEN, WIT_LEN, W>(r1cs_rows_size);
        let (_, x_ccs, w_ccs) = get_test_dummy_z_split::<R, X_LEN, WIT_LEN>();
        let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());

        let wit = Witness::from_w_ccs::<PP>(&w_ccs);

        // Make bound and securitty checks
        let witness_within_bound = check_witness_bound(&wit, PP::B);
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
            PP::B,
            PP::L,
            witness_within_bound,
        ) {
            println!(" Bound condition satisfied for 128 bits security");
        } else {
            println!("Bound condition not satisfied for 128 bits security");
        }

        let cm_i = CCCS {
            cm: wit.commit::<C, W, PP>(&scheme).unwrap(),
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

#[cfg(test)]
mod tests_goldilocks {

    use cyclotomic_rings::rings::GoldilocksChallengeSet;
    use lattirust_ring::cyclotomic_ring::models::goldilocks::{RqNTT, RqPoly};
    type CS = GoldilocksChallengeSet;
    generate_decomposition_tests!(1024, 2, 2, 10);
}

#[cfg(test)]
mod tests_frog {
    use cyclotomic_rings::rings::FrogChallengeSet;
    use lattirust_ring::cyclotomic_ring::models::frog_ring::{RqNTT, RqPoly};
    type CS = FrogChallengeSet;
    generate_decomposition_tests!(1024, 2, 2, 10);
}

#[cfg(test)]
mod tests_babybear {

    use cyclotomic_rings::rings::BabyBearChallengeSet;
    use lattirust_ring::cyclotomic_ring::models::babybear::{RqNTT, RqPoly};
    type CS = BabyBearChallengeSet;

    generate_decomposition_tests!(1024, 2, 2, 10);
}
