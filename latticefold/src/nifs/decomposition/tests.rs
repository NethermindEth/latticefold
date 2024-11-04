macro_rules! generate_decomposition_tests {
    ( $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        use crate::{
            arith::{r1cs::get_test_z_split, tests::get_test_ccs, Witness, CCCS},
            commitment::AjtaiCommitmentScheme,
            nifs::{
                decomposition::{
                    structs::{LFDecompositionProver, LFDecompositionVerifier},
                    utils::decompose_B_vec_into_k_vec,
                    DecompositionParams, DecompositionProver, DecompositionVerifier,
                },
                linearization::{
                    LFLinearizationProver, LFLinearizationVerifier, LinearizationProver,
                    LinearizationVerifier,
                },
            },
            transcript::poseidon::PoseidonTranscript,
        };
        use ark_ff::UniformRand;
        use lattirust_ring::balanced_decomposition::recompose;
        use rand::thread_rng;

        type T = PoseidonTranscript<RqNTT, CS>;

        #[derive(Clone)]
        struct PP;
        impl DecompositionParams for PP {
            const B: u128 = $b;
            const L: usize = $l;
            const B_SMALL: usize = $b_small;
            const K: usize = $k;
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
    };
}

#[cfg(test)]
mod tests_pow2 {
    use cyclotomic_rings::{challenge_set::BinarySmallSet, SuitableRing};
    use lattirust_ring::{
        balanced_decomposition::{
            convertible_ring::ConvertibleRing, decompose_balanced, decompose_balanced_vec,
            pad_and_transpose, Decompose,
        }, cyclotomic_ring::models::pow2_debug::Pow2CyclotomicPolyRingNTT, PolyRing
    };
    // Remove previous imports after macro is implemented

    use super::*;
    const Q: u64 = 17;
    const N: usize = 8;
    type RqNTT = Pow2CyclotomicPolyRingNTT<Q, N>;
    type CS = BinarySmallSet<Q, N>;
    generate_decomposition_tests!(1024, 2, 2, 10);

    // decompose_balance_vec in rust passes, so this should too
    // double check this test is correct for other rings
    #[test]
    fn test_decompose_B_vec_into_k_vec() {
        // Create a test vector
        const N: usize = 20;
        let test_vector: Vec<RqNTT> = (0..N * PP::L)
            .map(|_| RqNTT::rand(&mut thread_rng()))
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
            assert_eq!(vec.len(), N * PP::L);
        }

        // Check that the decomposition is correct
        for i in 0..test_vector.len() {
            let decomp_i = decomposed.iter().map(|d_j| d_j[i]).collect::<Vec<_>>();
            assert_eq!(
                test_vector[i],
                recompose(&decomp_i, RqNTT::from(PP::B_SMALL as u128))
            );
        }
    }

    #[test]
    fn test_decompose_big_vec_into_k_vec_and_compose_back() {
        // Create a test vector
        const N: usize = 10;
        let test_vector: Vec<RqNTT> = (0..N).map(|_| RqNTT::rand(&mut thread_rng())).collect();
        let coeff_repr: Vec<<RqNTT as SuitableRing>::CoefficientRepresentation> =
            test_vector.iter().map(|&x| x.into()).collect();

        // radix-B
        let decomposed_in_B: Vec<<RqNTT as SuitableRing>::CoefficientRepresentation> =
            pad_and_transpose(decompose_balanced_vec(&coeff_repr, PP::B, Some(PP::L)))
                .into_iter()
                .flatten()
                .collect();
        let decomposed_in_b_small =
            decompose_balanced_vec(&decomposed_in_B, PP::B_SMALL as u128, Some(PP::K));

        let recomposed_in_l: Vec<Vec<RqNTT>> = decomposed_in_b_small
            .into_iter()
            .map(|vec| {
                vec.chunks(PP::L)
                    .map(|chunk| {
                        recompose(
                            chunk,
                            <RqNTT as SuitableRing>::CoefficientRepresentation::from(PP::B),
                        )
                        .into()
                    })
                    .collect()
            })
            .collect();

        // Decompose and recompose
        let recomposed_in_b_small = recompose_from_k_vec_to_big_vec::<RqNTT, PP>(&recomposed_in_l);

        // Partially working
        assert_eq!(recomposed_in_b_small[0], decomposed_in_B[0]);
    }
    fn recompose_from_k_vec_to_big_vec<NTT: SuitableRing, DP: DecompositionParams>(
        k_vecs: &[Vec<NTT>],
    ) -> Vec<NTT::CoefficientRepresentation> {
        let decomposed_in_b: Vec<Vec<NTT::CoefficientRepresentation>> = k_vecs
            .iter()
            .map(|vec| {
                let decomposed_vec = vec
                    .into_iter()
                    .map(|&x| {
                        let coeff_repr: NTT::CoefficientRepresentation = x.into();
                        decompose_balanced_vec(&[coeff_repr], DP::B, Some(DP::L))
                            .into_iter()
                            .flatten()
                            .collect::<Vec<NTT::CoefficientRepresentation>>()[0]
                            .clone()
                    })
                    .collect::<Vec<NTT::CoefficientRepresentation>>();
                decomposed_vec
            })
            .collect();

        let mut decompose_in_l =
            vec![NTT::CoefficientRepresentation::default(); decomposed_in_b[0].len()];
        for j in 0..decomposed_in_b[0].len() {
            let coeffs: Vec<_> = decomposed_in_b
                .iter()
                .map(|decomposed_vec| decomposed_vec[j])
                .collect();
            decompose_in_l[j] = recompose(
                &coeffs,
                NTT::CoefficientRepresentation::from(DP::B_SMALL as u128),
            );
        }

        decompose_in_l
    }
}

#[cfg(test)]
mod tests_stark {

    use cyclotomic_rings::StarkChallengeSet;
    use lattirust_ring::cyclotomic_ring::models::stark_prime::RqNTT;
    use num_bigint::BigUint;

    use crate::{
        arith::{r1cs::get_test_dummy_z_split, tests::get_test_dummy_ccs},
        utils::security_check::{check_ring_modulus_128_bits_security, check_witness_bound},
    };
    type CS = StarkChallengeSet;
    generate_decomposition_tests!(1024, 2, 2, 10);

    #[test]
    fn test_decompose_B_vec_into_k_vec() {
        // Create a test vector
        const N: usize = 20;
        let test_vector: Vec<RqNTT> = (0..N * PP::L)
            .map(|_| RqNTT::rand(&mut thread_rng()))
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
            assert_eq!(vec.len(), N * PP::L);
        }

        // Check that the decomposition is correct
        for i in 0..test_vector.len() {
            let decomp_i = decomposed.iter().map(|d_j| d_j[i]).collect::<Vec<_>>();
            assert_eq!(
                test_vector[i],
                recompose(&decomp_i, RqNTT::from(PP::B_SMALL as u128))
            );
        }
    }

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

    use cyclotomic_rings::GoldilocksChallengeSet;
    use lattirust_ring::cyclotomic_ring::models::goldilocks::RqNTT;
    type CS = GoldilocksChallengeSet;
    generate_decomposition_tests!(1024, 2, 2, 10);
}

#[cfg(test)]
mod tests_frog {
    use cyclotomic_rings::FrogChallengeSet;
    use lattirust_ring::cyclotomic_ring::models::frog_ring::RqNTT;
    type CS = FrogChallengeSet;
    generate_decomposition_tests!(1024, 2, 2, 10);
}

#[cfg(test)]
mod tests_babybear {

    use cyclotomic_rings::BabyBearChallengeSet;
    use lattirust_ring::cyclotomic_ring::models::babybear::RqNTT;
    type CS = BabyBearChallengeSet;

    generate_decomposition_tests!(1024, 2, 2, 10);
}
