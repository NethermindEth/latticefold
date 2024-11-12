#[macro_export]
macro_rules! generate_linearization_tests {
    ( $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        use ark_std::sync::Arc;

        use ark_ff::UniformRand;
        use ark_std::One;
        use lattirust_poly::mle::DenseMultilinearExtension;
        use lattirust_poly::polynomials::{build_eq_x_r, VirtualPolynomial};
        use rand::thread_rng;
        use $crate::{
            arith::{r1cs::get_test_z_split, tests::get_test_ccs, Witness, CCCS},
            commitment::AjtaiCommitmentScheme,
            decomposition_parameters::DecompositionParams,
            nifs::linearization::{
                structs::LinearizationProver, utils::compute_u, LFLinearizationProver,
                LFLinearizationVerifier, LinearizationVerifier,
            },
            transcript::poseidon::PoseidonTranscript,
        };
        use $crate::{
            nifs::linearization::utils::prepare_lin_sumcheck_polynomial,
            //utils::mle::dense_vec_to_dense_mle,
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

        #[test]
        fn test_linearization_polynomial() {
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
                //M_z_mles.push(dense_vec_to_dense_mle(log_m, &mle));
                M_z_mles.push(DenseMultilinearExtension::from_slice(log_m, &mle));
            }

            let _ = g.add_mle_list(M_z_mles.clone().into_iter().map(|mle| Arc::new(mle)), c);
            let eq_b_r = build_eq_x_r(&beta).unwrap();
            let _ = g.mul_by_mle(eq_b_r, RqNTT::one());

            let polynomial =
                prepare_lin_sumcheck_polynomial(log_m, &[c], &M_z_mles, &[vec![0, 1, 2]], &beta)
                    .unwrap();

            for _ in 0..20 {
                let point: Vec<RqNTT> = (0..log_m).map(|_| RqNTT::rand(&mut rng)).collect();
                assert_eq!(
                    g.evaluate(&point).unwrap(),
                    polynomial.evaluate(&point).unwrap()
                )
            }
        }
        // Actual Tests
        #[test]
        fn test_linearization() {
            type R = RqNTT;
            type T = PoseidonTranscript<R, CS>;

            #[derive(Clone)]
            struct PP;
            impl DecompositionParams for PP {
                const B: u128 = $b;
                const L: usize = $l;
                const B_SMALL: usize = $b_small;
                const K: usize = $k;
            }

            const WIT_LEN: usize = 4; // 4 is the length of witness in this (Vitalik's) example
            const W: usize = WIT_LEN * PP::L; // the number of columns of the Ajtai matrix

            let ccs = get_test_ccs::<R>(W);
            let (_, x_ccs, w_ccs) = get_test_z_split::<R>(3);
            let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());

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
    generate_linearization_tests!(1024, 2, 2, 10);
}

#[cfg(test)]
mod tests_stark {

    use cyclotomic_rings::rings::StarkChallengeSet;
    use lattirust_ring::cyclotomic_ring::models::stark_prime::RqNTT;
    use num_bigint::BigUint;

    use crate::{
        arith::{r1cs::get_test_dummy_z_split, tests::get_test_dummy_ccs},
        utils::security_check::{check_ring_modulus_128_bits_security, check_witness_bound},
    };
    type CS = StarkChallengeSet;
    generate_linearization_tests!(1024, 2, 2, 10);

    #[test]
    fn test_dummy_linearization() {
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
    use lattirust_ring::cyclotomic_ring::models::goldilocks::RqNTT;
    type CS = GoldilocksChallengeSet;
    generate_linearization_tests!(1024, 2, 2, 10);
}

#[cfg(test)]
mod tests_frog {
    use cyclotomic_rings::rings::FrogChallengeSet;
    use lattirust_ring::cyclotomic_ring::models::frog_ring::RqNTT;
    type CS = FrogChallengeSet;
    generate_linearization_tests!(1024, 2, 2, 10);
}

#[cfg(test)]
mod tests_babybear {

    use cyclotomic_rings::rings::BabyBearChallengeSet;
    use lattirust_ring::cyclotomic_ring::models::babybear::RqNTT;
    type CS = BabyBearChallengeSet;

    generate_linearization_tests!(1024, 2, 2, 10);
}
