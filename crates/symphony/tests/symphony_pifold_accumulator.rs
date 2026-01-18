//! Accumulator tests for Π_fold.
//!
//! Tests multi-step accumulator workflow where instances are folded incrementally.

use ark_std::One;
use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
use std::sync::Arc;
use symphony::symphony_coins::derive_beta_chi;
use latticefold::transcript::Transcript;
use symphony::{
    pcs::{cmf_pcs, folding_pcs_l2},
    rp_rgchk::RPParams,
    symphony_pifold_batched::{verify_pi_fold_cp_poseidon_fs, PiFoldMatrices},
    symphony_pifold_streaming::{prove_pi_fold_poseidon_fs, PiFoldStreamingConfig},
    symphony_open::MultiAjtaiOpenVerifier,
};
use latticefold::commitment::AjtaiCommitmentScheme;
use stark_rings::{cyclotomic_ring::models::frog_ring::RqPoly as R, PolyRing, Ring};
use stark_rings_linalg::SparseMatrix;
use ark_ff::Field;

const MASTER_SEED: [u8; 32] = *b"SYMPHONY_AJTAI_SEED_V1_000000000";

fn fold_vec(beta0: R, beta1: R, a: &[R], b: &[R]) -> Vec<R> {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| beta0 * *x + beta1 * *y)
        .collect()
}

#[test]
fn test_pifold_accumulator_two_steps() {
    let n = 1 << 10;
    // NOTE: batchlin PCS (ℓ=2, r^3=n) currently requires log2(m*d) divisible by 3.
    // For Frog, d=16, so pick m=256 => m*d=4096 => log2=12.
    let m = 1 << 8;

    // Choose matrices so Π_had holds for any witness: M1=I, M2=0, M3=0.
    let mut m1 = SparseMatrix::<R>::identity(m);
    let mut m2 = SparseMatrix::<R>::identity(m);
    let mut m3 = SparseMatrix::<R>::identity(m);
    m1.pad_cols(n);
    m2.pad_cols(n);
    m3.pad_cols(n);
    for row in m2.coeffs.iter_mut() {
        row.clear();
    }
    for row in m3.coeffs.iter_mut() {
        row.clear();
    }
    let m1 = Arc::new(m1);
    let m2 = Arc::new(m2);
    let m3 = Arc::new(m3);

    let kappa_cm_f = 2usize;
    type BF = <<R as PolyRing>::BaseRing as Field>::BasePrimeField;

    let rg_params = RPParams {
        l_h: 64,
        // Keep m_J = (n/l_h)*lambda_pj <= m and m multiple of m_J for this small test.
        lambda_pj: 16,
        k_g: 4,
        d_prime: (R::dimension() as u128) / 2,
    };

    // CP commitment schemes for aux messages (needed for the canonical verifier API).
    let scheme_had =
        AjtaiCommitmentScheme::<R>::seeded(b"cfs_had_u", MASTER_SEED, 2, 3 * R::dimension());
    let scheme_mon =
        AjtaiCommitmentScheme::<R>::seeded(b"cfs_mon_b", MASTER_SEED, 2, rg_params.k_g);
    let scheme_g = AjtaiCommitmentScheme::<R>::seeded(b"cm_g", MASTER_SEED, 2, m * R::dimension());
    let open_cfs = MultiAjtaiOpenVerifier::<R>::new()
        .with_scheme("cfs_had_u", scheme_had.clone())
        .with_scheme("cfs_mon_b", scheme_mon.clone());

    // Accumulator starts at f_acc = 1.
    let mut f_acc = vec![R::one(); n];
    let pcs_params_f = {
        let flat_len = n * R::dimension();
        cmf_pcs::cmf_pcs_params_for_flat_len::<BF>(flat_len, kappa_cm_f).expect("cm_f pcs params")
    };
    let commit_cmf = |f: &[R]| -> Vec<R> {
        let flat: Vec<BF> = f
            .iter()
            .flat_map(|re| {
                re.coeffs()
                    .iter()
                    .map(|c| c.to_base_prime_field_elements().into_iter().next().expect("bf limb"))
            })
            .collect();
        let f_pcs = cmf_pcs::pad_flat_message(&pcs_params_f, &flat);
        let (t, _s) = folding_pcs_l2::commit(&pcs_params_f, &f_pcs).expect("cm_f pcs commit");
        cmf_pcs::pack_t_as_ring::<R>(&t)
    };
    let mut cm_acc = {
        commit_cmf(&f_acc)
    };

    // Step 1: fold in f1 = 0.
    let f1 = vec![R::ZERO; n];
    let cm1 = {
        let flat: Vec<BF> = f1
            .iter()
            .flat_map(|re| {
                re.coeffs()
                    .iter()
                    .map(|c| c.to_base_prime_field_elements().into_iter().next().expect("bf limb"))
            })
            .collect();
        let f_pcs = cmf_pcs::pad_flat_message(&pcs_params_f, &flat);
        let (t, _s) = folding_pcs_l2::commit(&pcs_params_f, &f_pcs).expect("cm_f pcs commit");
        cmf_pcs::pack_t_as_ring::<R>(&t)
    };

    // Build Ms: shared matrices repeated for each instance.
    let ms = vec![
        [m1.clone(), m2.clone(), m3.clone()],
        [m1.clone(), m2.clone(), m3.clone()],
    ];

    let cfg = PiFoldStreamingConfig::default();
    let pf1_out = prove_pi_fold_poseidon_fs::<R, PC>(
        ms.as_slice(),
        &[cm_acc.clone(), cm1.clone()],
        &[Arc::new(f_acc.clone()), Arc::new(f1.clone())],
        &[],
        Some(&scheme_had),
        Some(&scheme_mon),
        &scheme_g,
        rg_params.clone(),
        &cfg,
    )
    .unwrap();
    let pf1 = pf1_out.proof;

    // Recompute β for step 1 from the FS coin stream.
    let mut ts1 = symphony::public_coin_transcript::FixedTranscript::<R>::new_with_coins_and_events(
        pf1.coins.challenges.clone(),
        pf1.coins.bytes.clone(),
        pf1.coins.events.clone(),
    );
    ts1.absorb_field_element(&<R as stark_rings::PolyRing>::BaseRing::from(0x4c465053_50494250u128)); // "LFPS_PIBP"
    // absorb public inputs (empty)
    ts1.absorb_field_element(&<R as stark_rings::PolyRing>::BaseRing::from(0x4c465053_5055424cu128)); // "LFPS_PUBL"
    ts1.absorb_field_element(&<R as stark_rings::PolyRing>::BaseRing::from(0u128));
    let beta1_cts = derive_beta_chi::<R>(&mut ts1, 2);
    let beta10 = R::from(beta1_cts[0]);
    let beta11 = R::from(beta1_cts[1]);

    let attempt1 = verify_pi_fold_cp_poseidon_fs::<R, PC>(
        PiFoldMatrices::Shared([&*m1, &*m2, &*m3]),
        &[cm_acc.clone(), cm1.clone()],
        &pf1,
        &open_cfs,
        &pf1_out.cfs_had_u,
        &pf1_out.cfs_mon_b,
        &pf1_out.aux,
        &[],
    );
    let (folded1, _bat1) = attempt1.result.unwrap();

    // Update explicit accumulator witness and commitment.
    f_acc = fold_vec(beta10, beta11, &f_acc, &f1);
    cm_acc = commit_cmf(&f_acc);
    assert_eq!(folded1.c, cm_acc);

    // Step 2: fold in f2 = 1.
    let f2 = vec![R::one(); n];
    let cm2 = commit_cmf(&f2);

    let pf2_out = prove_pi_fold_poseidon_fs::<R, PC>(
        ms.as_slice(),
        &[cm_acc.clone(), cm2.clone()],
        &[Arc::new(f_acc.clone()), Arc::new(f2.clone())],
        &[],
        Some(&scheme_had),
        Some(&scheme_mon),
        &scheme_g,
        rg_params.clone(),
        &cfg,
    )
    .unwrap();
    let pf2 = pf2_out.proof;

    // Recompute β for step 2 from the FS coin stream.
    let mut ts2 = symphony::public_coin_transcript::FixedTranscript::<R>::new_with_coins_and_events(
        pf2.coins.challenges.clone(),
        pf2.coins.bytes.clone(),
        pf2.coins.events.clone(),
    );
    ts2.absorb_field_element(&<R as stark_rings::PolyRing>::BaseRing::from(0x4c465053_50494250u128)); // "LFPS_PIBP"
    // absorb public inputs (empty)
    ts2.absorb_field_element(&<R as stark_rings::PolyRing>::BaseRing::from(0x4c465053_5055424cu128)); // "LFPS_PUBL"
    ts2.absorb_field_element(&<R as stark_rings::PolyRing>::BaseRing::from(0u128));
    let beta2_cts = derive_beta_chi::<R>(&mut ts2, 2);
    let beta20 = R::from(beta2_cts[0]);
    let beta21 = R::from(beta2_cts[1]);

    let attempt2 = verify_pi_fold_cp_poseidon_fs::<R, PC>(
        PiFoldMatrices::Shared([&*m1, &*m2, &*m3]),
        &[cm_acc.clone(), cm2.clone()],
        &pf2,
        &open_cfs,
        &pf2_out.cfs_had_u,
        &pf2_out.cfs_mon_b,
        &pf2_out.aux,
        &[],
    );
    let (folded2, _bat2) = attempt2.result.unwrap();

    f_acc = fold_vec(beta20, beta21, &f_acc, &f2);
    cm_acc = commit_cmf(&f_acc);

    assert_eq!(folded2.c, cm_acc);
    assert!(!folded2.r.is_empty());
}
