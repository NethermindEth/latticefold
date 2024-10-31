#[cfg(test)]
mod tests_pow2 {
    use ark_ff::UniformRand;
    use cyclotomic_rings::challenge_set::BinarySmallSet;
    use lattirust_ring::{
        balanced_decomposition::recompose,
        cyclotomic_ring::models::pow2_debug::Pow2CyclotomicPolyRingNTT,
    };

    // Boilerplate code to generate values needed for testing
    const Q: u64 = 17; // Replace with an appropriate modulus
    const N: usize = 8;
    type RqNTT = Pow2CyclotomicPolyRingNTT<Q, N>;
    type CS = BinarySmallSet<Q, N>;

    use rand::thread_rng;

    use crate::{
        arith::{r1cs::tests::get_test_z_split, tests::get_test_ccs, Witness, CCCS},
        commitment::AjtaiCommitmentScheme,
        nifs::{
            decomposition::{
                utils::decompose_B_vec_into_k_vec, DecompositionParams, DecompositionProver,
                DecompositionVerifier, LFDecompositionProver, LFDecompositionVerifier,
            },
            linearization::{
                LFLinearizationProver, LFLinearizationVerifier, LinearizationProver,
                LinearizationVerifier,
            },
        },
        transcript::poseidon::PoseidonTranscript,
    };

    #[derive(Clone)]
    struct PP;

    type T = PoseidonTranscript<RqNTT, CS>;
    impl DecompositionParams for PP {
        const B: u128 = 1_024;
        const L: usize = 2;
        const B_SMALL: usize = 2;
        const K: usize = 10;
    }
    #[test]
    // decompose_balance_vec in rust passes, so this should too
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

        // Additional checks can be added here based on expected properties of the decomposition
    }

    // #[test]
    // fn test_decompose_big_vec_into_k_vec_and_compose_back() {
    //     // Create a test vector
    //     const N: usize = 10;
    //     let test_vector: Vec<RqNTT> = (0..N).map(|_| RqNTT::rand(&mut thread_rng())).collect();
    //     let coeff_repr: Vec<<RqNTT as SuitableRing>::CoefficientRepresentation> =
    //         test_vector.iter().map(|&x| x.into()).collect();

    //     // radix-B
    //     let decomposed_in_B: Vec<<RqNTT as SuitableRing>::CoefficientRepresentation> =
    //         pad_and_transpose(decompose_balanced_vec(&coeff_repr, PP::B, Some(PP::L)))
    //             .into_iter()
    //             .flatten()
    //             .collect();

    //     // Decompose and recompose
    //     let decomposed_in_k_recompose_in_l =
    //         decompose_big_vec_into_k_vec_and_compose_back::<RqNTT, PP>(&test_vector);
    //     let recomposed_in_b_small =
    //         recompose_from_k_vec_to_big_vec::<RqNTT, PP>(&decomposed_in_k_recompose_in_l);

    //     let recomposed_in_b_small_repr: Vec<<RqNTT as SuitableRing>::CoefficientRepresentation> =
    //         recomposed_in_b_small.iter().map(|&x| x.into()).collect();

    //     for (i, recomposed) in recomposed_in_b_small_repr.iter().enumerate() {
    //         assert_eq!(
    //             recomposed,
    //             &decomposed_in_B[i],
    //             "Mismatch at index {}: recomposed={:?}, original={:?}",
    //             i,
    //             recomposed,
    //             decomposed_in_B[i]
    //         );
    //     }
    // }
    // fn recompose_from_k_vec_to_big_vec<NTT: SuitableRing, DP: DecompositionParams>(
    //     k_vecs: &[Vec<NTT>],
    // ) -> Vec<NTT> {
    //     let decompose_in_l: Vec<Vec<NTT::CoefficientRepresentation>> = k_vecs
    //         .iter()
    //         .map(|vec| {
    //             vec.iter()
    //                 .map(|&inner_ring| {
    //                     let coeff: NTT::CoefficientRepresentation = inner_ring.into();
    //                     decompose_balanced_vec(&[coeff], DP::B, Some(DP::L))
    //                         .into_iter()
    //                         .flatten()
    //                         .collect::<Vec<NTT::CoefficientRepresentation>>()
    //                 })
    //                 .flatten()
    //                 .collect()
    //         })
    //         .collect();

    //     let recomposed_in_b_small: Vec<NTT> = decompose_in_l
    //         .into_iter()
    //         .map(|vec_k| {
    //             let ntt_repr: Vec<NTT> =
    //                 vec_k.iter().map(|&x| x.into()).collect();
    //             recompose(&ntt_repr, NTT::from(DP::B_SMALL as u128))
    //         })
    //         .collect();

    //     recomposed_in_b_small
    // }

    // Actual Tests
    #[test]
    fn test_decomposition() {
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

        let res = LFDecompositionVerifier::<_, T>::verify::<4, PP>(
            &lcccs,
            &decomposition_proof,
            &mut verifier_transcript,
            &ccs,
        );

        assert!(res.is_ok());
    }

    #[test]
    fn test_failing_decomposition() {
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

        let (_, _, w_ccs) = get_test_z_split::<RqNTT>(100);
        let fake_witness = Witness::<RqNTT>::from_w_ccs::<PP>(&w_ccs);

        let (_, _, decomposition_proof) = LFDecompositionProver::<_, T>::prove::<W, 4, PP>(
            &lcccs,
            &fake_witness,
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

        assert!(res.is_err());
    }
}
