use crate::{
    arith::{r1cs::get_test_z_split, Witness, CCCS, CCS},
    ark_base::*,
    commitment::AjtaiCommitmentScheme,
    decomposition_parameters::{test_params::DP, DecompositionParams},
    nifs::linearization::{
        utils::{compute_u, prepare_lin_sumcheck_polynomial},
        LFLinearizationProver, LFLinearizationVerifier, LinearizationProof, LinearizationProver,
        LinearizationVerifier,
    },
    transcript::poseidon::PoseidonTranscript,
};

use crate::arith::tests::get_test_ccs;

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Compress, Validate};
use ark_std::io::Cursor;
use cyclotomic_rings::{challenge_set::LatticefoldChallengeSet, rings::SuitableRing};
use lattirust_poly::{
    mle::DenseMultilinearExtension,
    polynomials::{build_eq_x_r, RefCounter, VirtualPolynomial},
};
use lattirust_ring::OverField;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;

fn test_compute_ui<R: OverField>() {
    let mut mles = Vec::with_capacity(10);
    let mut rng = ark_std::test_rng();

    for _i in 0..10 {
        let evals: Vec<R> = (0..8).map(|_| R::rand(&mut rng)).collect();

        mles.push(DenseMultilinearExtension::from_evaluations_slice(3, &evals))
    }

    for b in 0..8_u8 {
        let us: Vec<R> = compute_u(
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
        //M_z_mles.push(dense_vec_to_dense_mle(log_m, &mle));
        M_z_mles.push(DenseMultilinearExtension::from_slice(log_m, &mle));
    }

    let _ = g.add_mle_list(M_z_mles.clone().into_iter().map(RefCounter::new), c);
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

fn generate_linearization_proof<RqNTT, CS>() -> (
    LinearizationProof<RqNTT>,
    PoseidonTranscript<RqNTT, CS>,
    CCCS<4, RqNTT>,
    CCS<RqNTT>,
)
where
    RqNTT: OverField + SuitableRing,
    CS: LatticefoldChallengeSet<RqNTT>,
{
    const WIT_LEN: usize = 4; // 4 is the length of witness in this (Vitalik's) example
    const W: usize = WIT_LEN * DP::L; // the number of columns of the Ajtai matrix

    let ccs = get_test_ccs::<RqNTT>(W);
    let (_, x_ccs, w_ccs) = get_test_z_split::<RqNTT>(3);
    let mut rng = ChaCha8Rng::seed_from_u64(0);
    let scheme = AjtaiCommitmentScheme::rand(&mut rng);

    let wit: Witness<RqNTT> = Witness::from_w_ccs::<DP>(w_ccs);
    let cm_i: CCCS<4, RqNTT> = CCCS {
        cm: wit.commit::<4, W, DP>(&scheme).unwrap(),
        x_ccs,
    };

    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();

    let res = LFLinearizationProver::<_, PoseidonTranscript<RqNTT, CS>>::prove(
        &cm_i,
        &wit,
        &mut transcript,
        &ccs,
    );

    let transcript = PoseidonTranscript::<RqNTT, CS>::default();

    let linearization_proof = res.unwrap().1;

    (linearization_proof, transcript, cm_i, ccs)
}

fn test_linearization<RqNTT, CS>()
where
    RqNTT: OverField + SuitableRing,
    CS: LatticefoldChallengeSet<RqNTT>,
{
    let (proof, mut transcript, cm_i, ccs) = generate_linearization_proof::<RqNTT, CS>();

    let res = LFLinearizationVerifier::<_, PoseidonTranscript<RqNTT, CS>>::verify(
        &cm_i,
        &proof,
        &mut transcript,
        &ccs,
    );

    res.expect("Failed to verify linearization proof");
}

fn test_linearization_proof_serialization<RqNTT, CS>()
where
    RqNTT: OverField + SuitableRing,
    CS: LatticefoldChallengeSet<RqNTT>,
{
    let proof = generate_linearization_proof::<RqNTT, CS>().0;

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

mod tests_stark {
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
    use lattirust_ring::cyclotomic_ring::models::stark_prime::RqNTT;
    use num_bigint::BigUint;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    #[test]
    fn test_compute_ui() {
        super::test_compute_ui::<RqNTT>();
    }

    #[test]
    fn test_linearization_polynomial() {
        super::test_linearization_polynomial::<RqNTT>();
    }

    #[test]
    fn test_linearization() {
        super::test_linearization::<RqNTT, StarkChallengeSet>();
    }
    #[test]
    fn test_linearization_proof_serialization() {
        super::test_linearization_proof_serialization::<RqNTT, StarkChallengeSet>();
    }

    #[test]
    fn test_dummy_linearization() {
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
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        let scheme = AjtaiCommitmentScheme::rand(&mut rng);

        let wit = Witness::from_w_ccs::<StarkDP>(w_ccs);

        // Make bound and securitty checks
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
            #[cfg(feature = "std")]
            println!(" Bound condition satisfied for 128 bits security");
        } else {
            #[cfg(feature = "std")]
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

mod tests_goldilocks {
    use cyclotomic_rings::rings::GoldilocksChallengeSet;
    use lattirust_ring::cyclotomic_ring::models::goldilocks::RqNTT;
    #[test]
    fn test_compute_ui() {
        super::test_compute_ui::<RqNTT>();
    }

    #[test]
    fn test_linearization_polynomial() {
        super::test_linearization_polynomial::<RqNTT>();
    }

    #[test]
    fn test_linearization() {
        super::test_linearization::<RqNTT, GoldilocksChallengeSet>();
    }
    #[test]
    fn test_linearization_proof_serialization() {
        super::test_linearization_proof_serialization::<RqNTT, GoldilocksChallengeSet>();
    }
}

mod tests_frog {
    use cyclotomic_rings::rings::FrogChallengeSet;
    use lattirust_ring::cyclotomic_ring::models::frog_ring::RqNTT;
    #[test]
    fn test_compute_ui() {
        super::test_compute_ui::<RqNTT>();
    }

    #[test]
    fn test_linearization_polynomial() {
        super::test_linearization_polynomial::<RqNTT>();
    }

    #[test]
    fn test_linearization() {
        super::test_linearization::<RqNTT, FrogChallengeSet>();
    }
    #[test]
    fn test_linearization_proof_serialization() {
        super::test_linearization_proof_serialization::<RqNTT, FrogChallengeSet>();
    }
}

mod tests_babybear {
    use cyclotomic_rings::rings::BabyBearChallengeSet;
    use lattirust_ring::cyclotomic_ring::models::babybear::RqNTT;
    #[test]
    fn test_compute_ui() {
        super::test_compute_ui::<RqNTT>();
    }

    #[test]
    fn test_linearization_polynomial() {
        super::test_linearization_polynomial::<RqNTT>();
    }

    #[test]
    fn test_linearization() {
        super::test_linearization::<RqNTT, BabyBearChallengeSet>();
    }

    #[test]
    fn test_linearization_proof_serialization() {
        super::test_linearization_proof_serialization::<RqNTT, BabyBearChallengeSet>();
    }
}
