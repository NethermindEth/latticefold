use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
use latticefold::commitment::AjtaiCommitmentScheme;
use symphony::{
    pcs::{cmf_pcs, folding_pcs_l2},
    rp_rgchk::RPParams,
    symphony_open::MultiAjtaiOpenVerifier,
    symphony_pifold_batched::{verify_pi_fold_cp_poseidon_fs, PiFoldMatrices},
    symphony_pifold_streaming::{prove_pi_fold_poseidon_fs, PiFoldStreamingConfig},
};
use ark_ff::Field;
use stark_rings::{cyclotomic_ring::models::frog_ring::RqPoly as R, PolyRing, Ring};
use stark_rings_linalg::SparseMatrix;
use std::sync::Arc;

const MASTER_SEED: [u8; 32] = *b"SYMPHONY_AJTAI_SEED_V1_000000000";

/// Round-trip test: the Poseidon-FS streaming prover output should be accepted by the standard
/// R_cp verifier path (which uses the batched verifier + CP commitment opening checks).
#[test]
fn test_streaming_pifold_fs_roundtrip_verifies() {
    let n = 1 << 10;
    // NOTE: batchlin PCS (ℓ=2) currently requires log2(m*d) divisible by 3.
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

    let rg_params = RPParams {
        l_h: 64,
        lambda_pj: 1,
        k_g: 3,
        d_prime: (R::dimension() as u128) - 2,
    };

    // Two instances to exercise batching/weights.
    let f0 = Arc::new(vec![R::ONE; n]);
    let f1 = Arc::new(vec![R::ZERO; n]);

    let kappa = 8;
    type BF = <<R as PolyRing>::BaseRing as Field>::BasePrimeField;
    let pcs_params_f = {
        let flat_len = n * R::dimension();
        cmf_pcs::cmf_pcs_params_for_flat_len::<BF>(flat_len, kappa).expect("cm_f pcs params")
    };
    let cm0 = {
        let flat: Vec<BF> = f0
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

    // CP commitment schemes for aux messages (opened in R_cp).
    let scheme_had =
        AjtaiCommitmentScheme::<R>::seeded(b"cfs_had_u", MASTER_SEED, kappa, 3 * R::dimension());
    let scheme_mon =
        AjtaiCommitmentScheme::<R>::seeded(b"cfs_mon_b", MASTER_SEED, kappa, rg_params.k_g);
    let scheme_g = AjtaiCommitmentScheme::<R>::seeded(b"cm_g", MASTER_SEED, kappa, m * R::dimension());

    let open = MultiAjtaiOpenVerifier::<R>::new()
        .with_scheme("cfs_had_u", scheme_had.clone())
        .with_scheme("cfs_mon_b", scheme_mon.clone());

    let cfg = PiFoldStreamingConfig::default();
    let ms = vec![
        [m1.clone(), m2.clone(), m3.clone()],
        [m1.clone(), m2.clone(), m3.clone()],
    ];
    let out = prove_pi_fold_poseidon_fs::<R, PC>(
        ms.as_slice(),
        &[cm0.clone(), cm1.clone()],
        &[f0.clone(), f1.clone()],
        &[],
        Some(&scheme_had),
        Some(&scheme_mon),
        &scheme_g,
        rg_params,
        &cfg,
    )
    .unwrap();

    // Canonical CP verifier path should accept.
    let attempt = verify_pi_fold_cp_poseidon_fs::<R, PC>(
        PiFoldMatrices::Shared([&*m1, &*m2, &*m3]),
        &[cm0, cm1],
        &out.proof,
        &open,
        &out.cfs_had_u,
        &out.cfs_mon_b,
        &out.aux,
        &[],
    );
    attempt.result.unwrap();
}

