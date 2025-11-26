use ark_std::vec::Vec;
use latticefold_rings::rings::{BabyBearRingNTT, GoldilocksRingNTT, StarkRingNTT, SuitableRing};
use num_traits::One;
use rand::Rng;
#[cfg(feature = "parallel")]
use rayon::iter::{
    IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator, ParallelIterator,
};
use stark_rings_poly::mle::DenseMultilinearExtension;

use crate::{
    arith::{r1cs::get_test_z_split, tests::get_test_ccs, utils::mat_vec_mul, Witness, CCS, LCCCS},
    commitment::{AjtaiCommitmentScheme, Commitment},
    decomposition_parameters::{
        test_params::{BabyBearDP, GoldilocksDP, StarkDP},
        DecompositionParams,
    },
    nifs::{
        decomposition::{
            utils::{decompose_B_vec_into_k_vec, decompose_big_vec_into_k_vec_and_compose_back},
            DecompositionProver, DecompositionVerifier, LFDecompositionProver,
            LFDecompositionVerifier,
        },
        error::DecompositionError,
        linearization::utils::compute_u,
    },
    transcript::poseidon::PoseidonTranscript,
    utils::mle_helpers::{evaluate_mles, to_mles_err},
};

fn generate_decomposition_args<RqNTT, DP>(
    kappa: usize,
    n: usize,
    wit_len: usize,
) -> (
    LCCCS<RqNTT>,
    PoseidonTranscript<RqNTT>,
    PoseidonTranscript<RqNTT>,
    CCS<RqNTT>,
    Witness<RqNTT>,
    AjtaiCommitmentScheme<RqNTT>,
)
where
    RqNTT: SuitableRing,
    DP: DecompositionParams,
{
    let mut rng = ark_std::test_rng();
    let input: usize = rng.gen_range(1..101);
    let ccs = get_test_ccs(n, DP::L);
    let log_m = ccs.s;

    let scheme = AjtaiCommitmentScheme::rand(kappa, n, &mut rng);
    let (_, x_ccs, _) = get_test_z_split::<RqNTT>(input);

    let wit = Witness::rand::<_, DP>(&mut rng, wit_len);
    let mut z: Vec<RqNTT> = Vec::with_capacity(x_ccs.len() + wit_len + 1);

    z.extend_from_slice(&x_ccs);
    z.push(RqNTT::one());
    z.extend_from_slice(&wit.w_ccs);

    let cm: Commitment<RqNTT> = scheme.commit_ntt(&wit.f).unwrap();

    let r: Vec<RqNTT> = (0..log_m).map(|_| RqNTT::rand(&mut rng)).collect();
    let Mz_mles: Vec<DenseMultilinearExtension<RqNTT>> = ccs
        .M
        .iter()
        .map(|M| {
            DenseMultilinearExtension::from_evaluations_vec(log_m, mat_vec_mul(M, &z).unwrap())
        })
        .collect();

    let u = compute_u(&Mz_mles, &r).unwrap();

    let v = evaluate_mles::<RqNTT, &DenseMultilinearExtension<RqNTT>, _, DecompositionError>(
        &wit.f_hat, &r,
    )
    .unwrap();

    let lcccs = LCCCS {
        r,
        v,
        cm,
        u,
        x_w: x_ccs,
        h: RqNTT::one(),
    };

    (
        lcccs,
        PoseidonTranscript::<RqNTT>::default(),
        PoseidonTranscript::<RqNTT>::default(),
        ccs,
        wit,
        scheme,
    )
}

fn test_decomposition<RqNTT, DP>(kappa: usize, n: usize, wit_len: usize)
where
    RqNTT: SuitableRing,
    DP: DecompositionParams,
{
    let (lcccs, mut verifier_transcript, mut prover_transcript, ccs, wit, scheme) =
        generate_decomposition_args::<RqNTT, DP>(kappa, n, wit_len);

    let (_, _, _, decomposition_proof) =
        LFDecompositionProver::<_, PoseidonTranscript<RqNTT>>::prove::<DP>(
            &lcccs,
            &wit,
            &mut prover_transcript,
            &ccs,
            &scheme,
        )
        .unwrap();

    let res = LFDecompositionVerifier::<_, PoseidonTranscript<RqNTT>>::verify::<DP>(
        &lcccs,
        &decomposition_proof,
        &mut verifier_transcript,
        &ccs,
    );
    assert!(res.is_ok())
}

#[test]
fn test_decompose_witness() {
    type RqNTT = StarkRingNTT;
    type DP = StarkDP;
    const KAPPA: usize = 4;
    const WIT_LEN: usize = 4;
    const N: usize = WIT_LEN * DP::L;

    let (_, _, _, _, wit, _) = generate_decomposition_args::<RqNTT, DP>(KAPPA, N, WIT_LEN);

    let wit_vec =
        LFDecompositionProver::<_, PoseidonTranscript<RqNTT>>::decompose_witness::<DP>(&wit);

    // Compute expected result
    let f_s = decompose_B_vec_into_k_vec::<RqNTT, DP>(&wit.f_coeff);
    let expected_wit_vec: Vec<Witness<RqNTT>> = cfg_into_iter!(f_s)
        .map(Witness::from_f_coeff::<DP>)
        .collect();

    // Validate
    assert!(
        !wit_vec.is_empty(),
        "Decomposed witness vector should not be empty"
    );
    assert_eq!(
        wit_vec.len() * wit_vec[0].f_coeff.len(),
        wit.f_coeff.len() * DP::K,
        "Mismatch in decomposed witness vector length"
    );
    assert_eq!(
        wit_vec, expected_wit_vec,
        "Decomposed witness vector does not match expected evaluations"
    );
}

#[test]
fn test_compute_x_s() {
    type RqNTT = BabyBearRingNTT;
    type DP = BabyBearDP;
    const KAPPA: usize = 4;
    const WIT_LEN: usize = 4;
    const N: usize = WIT_LEN * DP::L;

    let (lcccs, _, _, _, _, _) = generate_decomposition_args::<RqNTT, DP>(KAPPA, N, WIT_LEN);
    let x_s = LFDecompositionProver::<_, PoseidonTranscript<RqNTT>>::compute_x_s::<DP>(
        lcccs.x_w.clone(),
        lcccs.h,
    );

    // Compute expected result
    let mut x_w_clone = lcccs.x_w.clone();
    x_w_clone.push(lcccs.h);
    let expected_x_s = decompose_big_vec_into_k_vec_and_compose_back::<RqNTT, DP>(x_w_clone);

    // Validate
    assert!(!x_s.is_empty(), "X_s vector should not be empty");
    assert_eq!(
        x_s.len(),
        lcccs.x_w.len() * DP::K,
        "Mismatch in X_s vector length"
    );
    assert_eq!(
        x_s, expected_x_s,
        "X_s vector does not match expected evaluations"
    );
}

#[test]
fn test_commit_witnesses() {
    type RqNTT = GoldilocksRingNTT;
    type DP = GoldilocksDP;
    const KAPPA: usize = 4;
    const WIT_LEN: usize = 4;
    const N: usize = WIT_LEN * DP::L;

    let (cm_i, _, _, _, wit, scheme) = generate_decomposition_args::<RqNTT, DP>(KAPPA, N, WIT_LEN);

    let wit_vec =
        LFDecompositionProver::<_, PoseidonTranscript<RqNTT>>::decompose_witness::<DP>(&wit);
    let y_s: Vec<Commitment<RqNTT>> =
        LFDecompositionProver::<_, PoseidonTranscript<RqNTT>>::commit_witnesses::<DP>(
            &wit_vec, &scheme, &cm_i,
        )
        .unwrap();

    // Compute expected result
    let expected_y_s: Vec<Commitment<RqNTT>> = cfg_iter!(wit_vec)
        .map(|wit| wit.commit::<DP>(&scheme))
        .collect::<Result<Vec<_>, _>>()
        .unwrap();

    // Validate
    assert!(!y_s.is_empty(), "Y_s vector should not be empty");
    assert_eq!(y_s.len(), wit_vec.len(), "Mismatch in Y_s vector length");
    assert_eq!(
        y_s, expected_y_s,
        "Y_s vector does not match expected evaluations"
    );
}

#[test]
fn test_compute_v_s() {
    type RqNTT = BabyBearRingNTT;
    type DP = BabyBearDP;
    const KAPPA: usize = 4;
    const WIT_LEN: usize = 4;
    const N: usize = WIT_LEN * DP::L;

    let (lcccs, _, _, _, wit, _) = generate_decomposition_args::<RqNTT, DP>(KAPPA, N, WIT_LEN);
    let wit_vec =
        LFDecompositionProver::<_, PoseidonTranscript<RqNTT>>::decompose_witness::<DP>(&wit);
    let v_s =
        LFDecompositionProver::<_, PoseidonTranscript<RqNTT>>::compute_v_s(&wit_vec, &lcccs.r)
            .unwrap();

    // Compute expected result
    let expected_v_s: Vec<Vec<RqNTT>> = cfg_iter!(wit_vec)
        .map(|wit| evaluate_mles::<RqNTT, _, _, DecompositionError>(&wit.f_hat, &lcccs.r).unwrap())
        .collect();

    // Validate
    assert!(!v_s.is_empty(), "V_s vector should not be empty");
    assert_eq!(v_s.len(), wit_vec.len(), "Mismatch in V_s vector length");
    assert_eq!(
        v_s, expected_v_s,
        "V_s vector does not match expected evaluations"
    );
}

#[test]
fn test_compute_u_s() {
    type RqNTT = GoldilocksRingNTT;
    type DP = GoldilocksDP;
    const KAPPA: usize = 4;
    const WIT_LEN: usize = 4;
    const N: usize = WIT_LEN * DP::L;

    let (lcccs, _, _, ccs, wit, _) = generate_decomposition_args::<RqNTT, DP>(KAPPA, N, WIT_LEN);
    let wit_vec =
        LFDecompositionProver::<_, PoseidonTranscript<RqNTT>>::decompose_witness::<DP>(&wit);
    let x_s = LFDecompositionProver::<_, PoseidonTranscript<RqNTT>>::compute_x_s::<DP>(
        lcccs.x_w.clone(),
        lcccs.h,
    );
    let mz_mles = LFDecompositionProver::<_, PoseidonTranscript<RqNTT>>::compute_mz_mles(
        &wit_vec, &ccs.M, &x_s, ccs.s,
    )
    .unwrap();

    let u_s =
        LFDecompositionProver::<_, PoseidonTranscript<RqNTT>>::compute_u_s(&mz_mles, &lcccs.r)
            .unwrap();

    // Compute expected result
    let expected_u_s: Vec<Vec<RqNTT>> = cfg_iter!(wit_vec)
        .enumerate()
        .map(|(i, wit)| {
            let z: Vec<_> = {
                let mut z = Vec::with_capacity(x_s[i].len() + wit.w_ccs.len());

                z.extend_from_slice(&x_s[i]);
                z.extend_from_slice(&wit.w_ccs);

                z
            };
            let mles = to_mles_err::<_, _, DecompositionError, _>(
                ccs.s,
                cfg_iter!(ccs.M).map(|M| mat_vec_mul(M, &z)),
            )
            .unwrap();

            evaluate_mles::<RqNTT, &DenseMultilinearExtension<_>, _, DecompositionError>(
                &mles, &lcccs.r,
            )
            .unwrap()
        })
        .collect();

    // Validate
    assert!(!u_s.is_empty(), "U_s vector should not be empty");
    assert_eq!(u_s.len(), wit_vec.len(), "Mismatch in U_s vector length");
    assert_eq!(
        u_s, expected_u_s,
        "U_s vector does not match expected evaluations"
    );
}

#[test]
fn test_test_decomposition() {
    type RqNTT = StarkRingNTT;
    type DP = StarkDP;
    const WIT_LEN: usize = 4;
    const KAPPA: usize = 4;
    const N: usize = WIT_LEN * DP::L;

    test_decomposition::<RqNTT, DP>(KAPPA, N, WIT_LEN);
}

#[test]
fn test_recompose_commitment() {
    type RqNTT = GoldilocksRingNTT;
    type DP = GoldilocksDP;
    type T = PoseidonTranscript<RqNTT>;
    type Verifier = LFDecompositionVerifier<RqNTT, T>;
    const WIT_LEN: usize = 4;
    const N: usize = WIT_LEN * DP::L;
    const KAPPA: usize = 4;

    let (lcccs, _, mut prover_transcript, ccs, wit, scheme) =
        generate_decomposition_args::<RqNTT, DP>(KAPPA, N, WIT_LEN);

    let (_, _, _, proof) = LFDecompositionProver::<_, T>::prove::<DP>(
        &lcccs,
        &wit,
        &mut prover_transcript,
        &ccs,
        &scheme,
    )
    .unwrap();

    let b_s = Verifier::calculate_b_s::<DP>();

    let should_equal_y0 =
        Verifier::recompose_commitment(&proof.y_s, &b_s).expect("Recomposing proof failed");

    assert_eq!(should_equal_y0, lcccs.cm);
}

#[test]
fn test_recompose_u() {
    type RqNTT = StarkRingNTT;
    type DP = StarkDP;
    type T = PoseidonTranscript<RqNTT>;
    type Verifier = LFDecompositionVerifier<RqNTT, T>;
    const WIT_LEN: usize = 4;
    const N: usize = WIT_LEN * DP::L;
    const KAPPA: usize = 4;

    let (lcccs, _, mut prover_transcript, ccs, wit, scheme) =
        generate_decomposition_args::<RqNTT, DP>(KAPPA, N, WIT_LEN);
    let (_, _, _, proof) = LFDecompositionProver::<_, T>::prove::<DP>(
        &lcccs,
        &wit,
        &mut prover_transcript,
        &ccs,
        &scheme,
    )
    .unwrap();

    let b_s = Verifier::calculate_b_s::<DP>();

    let u = Verifier::recompose(&proof.u_s, &b_s).expect("Recomposing proof u failed");

    assert_eq!(u, lcccs.u);
}

#[test]
fn test_recompose_v() {
    type RqNTT = BabyBearRingNTT;
    type DP = BabyBearDP;
    type T = PoseidonTranscript<RqNTT>;
    type Verifier = LFDecompositionVerifier<RqNTT, T>;
    const WIT_LEN: usize = 4;
    const N: usize = WIT_LEN * DP::L;
    const KAPPA: usize = 4;

    let (lcccs, _, mut prover_transcript, ccs, wit, scheme) =
        generate_decomposition_args::<RqNTT, DP>(KAPPA, N, WIT_LEN);

    let (_, _, _, proof) = LFDecompositionProver::<_, T>::prove::<DP>(
        &lcccs,
        &wit,
        &mut prover_transcript,
        &ccs,
        &scheme,
    )
    .unwrap();

    let b_s = Verifier::calculate_b_s::<DP>();

    let v = Verifier::recompose(&proof.v_s, &b_s).expect("Recomposing proof u failed");

    assert_eq!(v, lcccs.v);
}

#[test]
fn test_recompose_xw_and_h() {
    type RqNTT = GoldilocksRingNTT;
    type DP = GoldilocksDP;
    type T = PoseidonTranscript<RqNTT>;
    type Verifier = LFDecompositionVerifier<RqNTT, T>;
    const WIT_LEN: usize = 4;
    const N: usize = WIT_LEN * DP::L;
    const KAPPA: usize = 4;

    let (lcccs, _, mut prover_transcript, ccs, wit, scheme) =
        generate_decomposition_args::<RqNTT, DP>(KAPPA, N, WIT_LEN);

    let (_, _, _, proof) = LFDecompositionProver::<_, T>::prove::<DP>(
        &lcccs,
        &wit,
        &mut prover_transcript,
        &ccs,
        &scheme,
    )
    .unwrap();

    let b_s = Verifier::calculate_b_s::<DP>();

    let mut x_w = Verifier::recompose(&proof.x_s, &b_s).expect("Recomposing proof x_w failed");
    let h = x_w.pop().expect("x_w does not contain h");

    assert_eq!(x_w, lcccs.x_w);
    assert_eq!(h, lcccs.h);
}

#[test]
fn test_verify_full() {
    type RqNTT = StarkRingNTT;
    type DP = StarkDP;
    type T = PoseidonTranscript<RqNTT>;
    type Verifier = LFDecompositionVerifier<RqNTT, T>;
    const WIT_LEN: usize = 4;
    const N: usize = WIT_LEN * DP::L;
    const KAPPA: usize = 4;

    let (lcccs, mut verifier_transcript, mut prover_transcript, ccs, wit, scheme) =
        generate_decomposition_args::<RqNTT, DP>(KAPPA, N, WIT_LEN);

    let (_, _, _, proof) = LFDecompositionProver::<_, T>::prove::<DP>(
        &lcccs,
        &wit,
        &mut prover_transcript,
        &ccs,
        &scheme,
    )
    .unwrap();

    let _ = Verifier::verify::<DP>(&lcccs, &proof, &mut verifier_transcript, &ccs)
        .expect("Failed to verify decomposition proof");
}

#[test]
fn test_verify_invalid_proof() {
    type RqNTT = GoldilocksRingNTT;
    type DP = GoldilocksDP;
    type T = PoseidonTranscript<RqNTT>;
    type Verifier = LFDecompositionVerifier<RqNTT, T>;
    const WIT_LEN: usize = 4;
    const N: usize = WIT_LEN * DP::L;
    const KAPPA: usize = 4;

    let (lcccs, mut verifier_transcript, mut prover_transcript, ccs, wit, scheme) =
        generate_decomposition_args::<RqNTT, DP>(KAPPA, N, WIT_LEN);

    let (_, _, _, mut proof) = LFDecompositionProver::<_, T>::prove::<DP>(
        &lcccs,
        &wit,
        &mut prover_transcript,
        &ccs,
        &scheme,
    )
    .unwrap();

    // Make proof components have mismatched lengths
    if !proof.v_s.is_empty() {
        proof.v_s[0][0] += RqNTT::one();
    }

    let result = Verifier::verify::<DP>(&lcccs, &proof, &mut verifier_transcript, &ccs);
    assert!(
        result.is_err(),
        "Verification should fail with mismatched proof component lengths"
    );
}
