use crate::{
    arith::{r1cs::get_test_z_split, tests::get_test_ccs, Witness, CCCS},
    commitment::AjtaiCommitmentScheme,
    decomposition_parameters::DecompositionParams,
    nifs::{
        decomposition::{
            DecompositionProver, DecompositionVerifier, LFDecompositionProver,
            LFDecompositionVerifier,
        },
        folding::{
            prepare_g1_i_mle_list, prepare_g2_i_mle_list, prepare_g3_i_mle_list, FoldingProver,
            FoldingVerifier, LFFoldingProver, LFFoldingVerifier,
        },
        linearization::{
            LFLinearizationProver, LFLinearizationVerifier, LinearizationProver,
            LinearizationVerifier,
        },
    },
    transcript::poseidon::PoseidonTranscript,
};
use ark_ff::{One, UniformRand};
use ark_std::iter::successors;
use ark_std::Zero;
use cyclotomic_rings::{StarkChallengeSet, StarkRingNTT};
use lattirust_poly::{
    mle::DenseMultilinearExtension,
    polynomials::{build_eq_x_r, eq_eval, VirtualPolynomial},
};
use rand::thread_rng;

// Boilerplate code to generate values needed for testing
type R = StarkRingNTT;
type CS = StarkChallengeSet;
type T = PoseidonTranscript<StarkRingNTT, CS>;

#[derive(Clone)]
struct PP;

impl DecompositionParams for PP {
    const B: u128 = 1_024;
    const L: usize = 1;
    const B_SMALL: usize = 2;
    const K: usize = 10;
}
#[test]
fn test_g_1() {
    let mut rng = thread_rng();
    let m = 8;
    let log_m = 3;

    let f_i: Vec<R> = (0..m).map(|_| R::rand(&mut rng)).collect();
    let r_i: Vec<R> = (0..log_m).map(|_| R::rand(&mut rng)).collect();

    let fi_mle = DenseMultilinearExtension::from_evaluations_vec(log_m, f_i.clone());
    let r_i_eq = build_eq_x_r(&r_i).unwrap();
    let mle_coeff = R::rand(&mut rng);
    fn evaluate(x: &[R], f_i: &DenseMultilinearExtension<R>, r_i: &[R], coeff: &R) -> R {
        eq_eval(r_i, x).unwrap() * f_i.evaluate(x).unwrap() * coeff
    }

    let mut g = VirtualPolynomial::new(log_m);

    let _ = prepare_g1_i_mle_list(&mut g, fi_mle.clone(), r_i_eq, mle_coeff);

    for _ in 0..20 {
        let point: Vec<RqNTT> = (0..log_m).map(|_| R::rand(&mut rng)).collect();
        assert_eq!(
            g.evaluate(&point).unwrap(),
            evaluate(&point, &fi_mle, &r_i, &mle_coeff)
        )
    }
}

#[test]
fn test_g_2() {
    let mut rng = thread_rng();
    let m = 8;
    let log_m = 3;
    let b = 8;
    let f_i: Vec<R> = (0..m).map(|_| R::rand(&mut rng)).collect();
    let beta: Vec<R> = (0..log_m).map(|_| R::rand(&mut rng)).collect();

    let fi_mle = DenseMultilinearExtension::from_evaluations_vec(log_m, f_i.clone());
    let beta_eq_x = build_eq_x_r(&beta).unwrap();
    let mu_i = R::rand(&mut rng);

    fn evaluate(x: &[R], f_i: &DenseMultilinearExtension<R>, b: usize, beta: &[R], mu_i: &R) -> R {
        let mut evaluation = R::one();

        for i in 1..b {
            let i_hat = R::from(i as u128);

            evaluation *= f_i.evaluate(x).unwrap() - i_hat;
            evaluation *= f_i.evaluate(x).unwrap() + i_hat;
        }
        evaluation *= f_i.evaluate(x).unwrap();
        evaluation *= eq_eval(beta, x).unwrap();
        evaluation * mu_i
    }

    let mut g = VirtualPolynomial::new(log_m);

    let _ = prepare_g2_i_mle_list(&mut g, fi_mle.clone(), b, mu_i, beta_eq_x);

    for _ in 0..20 {
        let point: Vec<RqNTT> = (0..log_m).map(|_| R::rand(&mut rng)).collect();
        assert_eq!(
            g.evaluate(&point).unwrap(),
            evaluate(&point, &fi_mle, b, &beta, &mu_i)
        )
    }
}

#[test]
fn test_g_3() {
    let mut rng = thread_rng();
    let m = 8;
    let log_m = 3;
    let t = 3;

    let mz_s: Vec<Vec<R>> = (0..t)
        .map(|_| (0..m).map(|_| R::rand(&mut rng)).collect())
        .collect();
    let r_i: Vec<R> = (0..log_m).map(|_| R::rand(&mut rng)).collect();

    let mz_mles: Vec<DenseMultilinearExtension<R>> = mz_s
        .into_iter()
        .map(|m_z| DenseMultilinearExtension::from_evaluations_vec(log_m, m_z))
        .collect();
    let r_i_eq = build_eq_x_r(&r_i).unwrap();
    let zeta_i = R::rand(&mut rng);
    fn evaluate(x: &[R], mz_mles: &[DenseMultilinearExtension<R>], r_i: &[R], zeta_i: &R) -> R {
        let mut evaluation = R::zero();

        for (zeta, M) in successors(Some(*zeta_i), |y| Some(*zeta_i * *y)).zip(mz_mles.iter()) {
            evaluation += zeta * M.evaluate(x).unwrap();
        }
        evaluation * eq_eval(x, r_i).unwrap()
    }

    let mut g = VirtualPolynomial::new(log_m);

    let _ = prepare_g3_i_mle_list(&mut g, &mz_mles, zeta_i, r_i_eq);

    for _ in 0..20 {
        let point: Vec<RqNTT> = (0..log_m).map(|_| R::rand(&mut rng)).collect();
        assert_eq!(
            g.evaluate(&point).unwrap(),
            evaluate(&point, &mz_mles, &r_i, &zeta_i)
        )
    }
}
#[test]
fn test_folding() {
    const WIT_LEN: usize = 4; // 4 is the length of witness in this (Vitalik's) example
    const W: usize = WIT_LEN * PP::L; // the number of columns of the Ajtai matrix

    let ccs = get_test_ccs::<R>(W);
    let (_, x_ccs, w_ccs) = get_test_z_split::<R>(3);
    let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());
    let wit: Witness<R> = Witness::from_w_ccs::<PP>(&w_ccs);
    let cm_i: CCCS<4, R> = CCCS {
        cm: wit.commit::<4, 4, PP>(&scheme).unwrap(),
        x_ccs,
    };

    let mut prover_transcript = PoseidonTranscript::<R, CS>::default();
    let mut verifier_transcript = PoseidonTranscript::<R, CS>::default();

    let (_, linearization_proof) =
        LFLinearizationProver::<_, T>::prove(&cm_i, &wit, &mut prover_transcript, &ccs).unwrap();

    let lcccs = LFLinearizationVerifier::<_, PoseidonTranscript<R, CS>>::verify(
        &cm_i,
        &linearization_proof,
        &mut verifier_transcript,
        &ccs,
    )
    .unwrap();

    let (_, vec_wit, decomposition_proof) = LFDecompositionProver::<_, T>::prove::<4, 4, PP>(
        &lcccs,
        &wit,
        &mut prover_transcript,
        &ccs,
        &scheme,
    )
    .unwrap();

    let vec_lcccs = LFDecompositionVerifier::<_, T>::verify::<4, PP>(
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
        LFFoldingProver::<_, T>::prove::<4, PP>(&lcccs, &wit_s, &mut prover_transcript, &ccs)
            .unwrap();

    let lcccs_verifier = LFFoldingVerifier::<_, T>::verify::<4, PP>(
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
    const W: usize = WIT_LEN * PP::L; // the number of columns of the Ajtai matrix

    let ccs = get_test_ccs::<R>(W);
    let (_, x_ccs, w_ccs) = get_test_z_split::<R>(3);
    let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());
    let wit: Witness<R> = Witness::from_w_ccs::<PP>(&w_ccs);
    let cm_i: CCCS<4, R> = CCCS {
        cm: wit.commit::<4, 4, PP>(&scheme).unwrap(),
        x_ccs,
    };

    let mut prover_transcript = PoseidonTranscript::<R, CS>::default();
    let mut verifier_transcript = PoseidonTranscript::<R, CS>::default();

    let (_, linearization_proof) =
        LFLinearizationProver::<_, T>::prove(&cm_i, &wit, &mut prover_transcript, &ccs).unwrap();

    let lcccs = LFLinearizationVerifier::<_, PoseidonTranscript<R, CS>>::verify(
        &cm_i,
        &linearization_proof,
        &mut verifier_transcript,
        &ccs,
    )
    .unwrap();

    let (_, mut vec_wit, decomposition_proof) = LFDecompositionProver::<_, T>::prove::<4, 4, PP>(
        &lcccs,
        &wit,
        &mut prover_transcript,
        &ccs,
        &scheme,
    )
    .unwrap();

    let vec_lcccs = LFDecompositionVerifier::<_, T>::verify::<4, PP>(
        &lcccs,
        &decomposition_proof,
        &mut verifier_transcript,
        &ccs,
    )
    .unwrap();

    vec_wit[0] = Witness::<R>::from_w_ccs::<PP>(&w_ccs);

    let res =
        LFFoldingProver::<_, T>::prove::<4, PP>(&vec_lcccs, &vec_wit, &mut prover_transcript, &ccs);

    assert!(res.is_err())
}

use lattirust_ring::cyclotomic_ring::models::stark_prime::RqNTT;
use num_bigint::BigUint;

use crate::{arith::r1cs::get_test_dummy_z_split, utils::security_check::check_witness_bound};
use crate::{
    arith::tests::get_test_dummy_ccs, utils::security_check::check_ring_modulus_128_bits_security,
};

#[test]
fn test_dummy_folding() {
    #[cfg(feature = "dhat-heap")]
    #[global_allocator]
    static ALLOC: dhat::Alloc = dhat::Alloc;

    type R = RqNTT;
    type CS = StarkChallengeSet;
    type T = PoseidonTranscript<R, CS>;

    #[derive(Clone)]
    struct PP;
    impl DecompositionParams for PP {
        const B: u128 = 3010936384;
        const L: usize = 8;
        const B_SMALL: usize = 38;
        const K: usize = 6;
    }

    const C: usize = 15;
    const X_LEN: usize = 1;
    const WIT_LEN: usize = 512;
    const W: usize = WIT_LEN * PP::L; // the number of columns of the Ajtai matrix
    let r1cs_rows_size = X_LEN + WIT_LEN + 1; // Let's have a square matrix

    #[cfg(feature = "dhat-heap")]
    let _profiler = dhat::Profiler::new_heap(); // Move a round to measure specific parts

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

    let mut prover_transcript = PoseidonTranscript::<R, CS>::default();

    let linearization_proof =
        LFLinearizationProver::<_, T>::prove(&cm_i, &wit, &mut prover_transcript, &ccs);

    let mut verifier_transcript = PoseidonTranscript::<R, CS>::default();

    let linearization_verification = LFLinearizationVerifier::<_, T>::verify(
        &cm_i,
        &linearization_proof
            .expect("Linearization proof generation error")
            .1,
        &mut verifier_transcript,
        &ccs,
    )
    .expect("Linearization Verification error");

    let lcccs = linearization_verification;

    let decomposition_prover = LFDecompositionProver::<_, T>::prove::<W, C, PP>(
        &lcccs,
        &wit,
        &mut prover_transcript,
        &ccs,
        &scheme,
    );

    let decomposition_proof = decomposition_prover.expect("Decomposition proof generation error");

    let decomposition_verification = LFDecompositionVerifier::<_, T>::verify::<C, PP>(
        &lcccs,
        &decomposition_proof.2,
        &mut verifier_transcript,
        &ccs,
    );

    let lcccs = decomposition_verification.expect("Decomposition Verification error");

    #[cfg(feature = "dhat-heap")]
    let _profiler = dhat::Profiler::new_heap();
    let (lcccs, wit_s) = {
        let mut lcccs = lcccs.clone();
        let mut lcccs_r = lcccs.clone();
        lcccs.append(&mut lcccs_r);

        let mut wit_s = decomposition_proof.1.clone();
        let mut wit_s_r = decomposition_proof.1;
        wit_s.append(&mut wit_s_r);

        (lcccs, wit_s)
    };
    let folding_prover =
        LFFoldingProver::<_, T>::prove::<C, PP>(&lcccs, &wit_s, &mut prover_transcript, &ccs);

    let folding_proof = folding_prover.expect("Folding proof generation error");

    let folding_verification = LFFoldingVerifier::<_, T>::verify::<C, PP>(
        &lcccs,
        &folding_proof.2,
        &mut verifier_transcript,
        &ccs,
    );

    folding_verification.expect("Folding Verification error");
}
