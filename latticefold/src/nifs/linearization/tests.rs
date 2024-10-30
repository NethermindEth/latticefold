#[macro_export]
macro_rules! generate_tests {
    () => {
        use ark_ff::UniformRand;
        use lattirust_poly::mle::DenseMultilinearExtension;
        use rand::thread_rng;
        use $crate::{
            arith::{r1cs::tests::get_test_z_split, tests::get_test_ccs, Witness, CCCS},
            commitment::AjtaiCommitmentScheme,
            decomposition_parameters::DecompositionParams,
            nifs::linearization::{
                structs::LinearizationProver, utils::compute_u, LFLinearizationProver,
                LFLinearizationVerifier, LinearizationVerifier,
            },
            transcript::poseidon::PoseidonTranscript,
        };

        #[test]
        fn test_compute_u() {
            let mut mles = Vec::with_capacity(10);
            let mut rng = ark_std::test_rng();
            // generate evals
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

        // Actual Tests
        #[test]
        fn test_linearization() {
            type R = RqNTT;
            type T = PoseidonTranscript<R, CS>;

            impl DecompositionParams for PP {
                const B: u128 = 1_024;
                const L: usize = 2;
                const B_SMALL: usize = 2;
                const K: usize = 10;
            }

            const WIT_LEN: usize = 4; // 4 is the length of witness in this (Vitalik's) example
            const W: usize = WIT_LEN * PP::L; // the number of columns of the Ajtai matrix

            let ccs = get_test_ccs::<R>(W);
            let (_, x_ccs, w_ccs) = get_test_z_split::<R>(3);
            let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());
            #[derive(Clone)]
            struct PP;

            let wit: Witness<R> = Witness::from_w_ccs::<PP>(&w_ccs);
            let cm_i: CCCS<4, R> = CCCS {
                cm: wit.commit::<4, W, PP>(&scheme).unwrap(),
                x_ccs,
            };
            let mut transcript = PoseidonTranscript::<R, CS>::default();

            let res = LFLinearizationProver::<_, T>::prove(&cm_i, &wit, &mut transcript, &ccs);

            let mut transcript = PoseidonTranscript::<R, CS>::default();

            let res = LFLinearizationVerifier::<_, PoseidonTranscript<R, CS>>::verify(
                &cm_i,
                &res.unwrap().1,
                &mut transcript,
                &ccs,
            );

            res.unwrap();
        }
    };
}

#[cfg(test)]
mod tests_pow2 {

    use cyclotomic_rings::challenge_set::BinarySmallSet;
    use lattirust_ring::cyclotomic_ring::models::pow2_debug::Pow2CyclotomicPolyRingNTT;

    // Boilerplate code to generate values needed for testing
    const Q: u64 = 17; // Replace with an appropriate modulus
    const N: usize = 8;
    type RqNTT = Pow2CyclotomicPolyRingNTT<Q, N>;
    type CS = BinarySmallSet<Q, N>;
    generate_tests!();
}

#[cfg(test)]
mod tests_stark {

    use cyclotomic_rings::StarkChallengeSet;
    use lattirust_ring::cyclotomic_ring::models::stark_prime::RqNTT;
    type CS = StarkChallengeSet;
    generate_tests!();
}

#[cfg(test)]
mod tests_goldilocks {

    use cyclotomic_rings::GoldilocksChallengeSet;
    use lattirust_ring::cyclotomic_ring::models::goldilocks::RqNTT;
    type CS = GoldilocksChallengeSet;
    generate_tests!();
}

#[cfg(test)]
mod tests_frog {
    use cyclotomic_rings::FrogChallengeSet;
    use lattirust_ring::cyclotomic_ring::models::frog_ring::RqNTT;
    type CS = FrogChallengeSet;
    generate_tests!();
}

#[cfg(test)]
mod tests_babybear {

    use cyclotomic_rings::BabyBearChallengeSet;
    use lattirust_ring::cyclotomic_ring::models::babybear::RqNTT;
    type CS = BabyBearChallengeSet;

    generate_tests!();
}
