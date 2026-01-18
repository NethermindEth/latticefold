use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
use latticefold::commitment::AjtaiCommitmentScheme;
use stark_rings::{cyclotomic_ring::models::frog_ring::RqPoly as R, PolyRing, Ring};
use stark_rings_linalg::SparseMatrix;
use ark_ff::Field;
use symphony::pcs::{cmf_pcs, folding_pcs_l2};
use symphony::rp_rgchk::RPParams;
use symphony::symphony_open::MultiAjtaiOpenVerifier;
use symphony::symphony_pifold_batched::{
    verify_pi_fold_cp_poseidon_fs, PiFoldMatrices,
};
use symphony::poseidon_trace::replay_poseidon_transcript_trace;
use symphony::symphony_pifold_streaming::{prove_pi_fold_poseidon_fs, PiFoldStreamingConfig};

const MASTER_SEED: [u8; 32] = *b"SYMPHONY_AJTAI_SEED_V1_000000000";

#[test]
fn test_pifold_streaming_hetero_m_roundtrip() {
    // Keep n >= m so SparseMatrix::pad_cols(n) doesn't need to "shrink" identity matrices.
    let n = 1 << 5; // 32 vars
    // NOTE: batchlin PCS (ℓ=2) currently requires log2(m*d) divisible by 3.
    // For Frog, d=16, so pick m=32 => m*d=512 => log2=9.
    let m = 1 << 5; // 32 rows per chunk/instance

    // Two different A-matrices; B=C=0 so A*f ∘ B*f - C*f == 0 holds for any witness.
    let mut a0 = SparseMatrix::<R>::identity(m);
    let mut b0 = SparseMatrix::<R>::identity(m);
    let mut c0 = SparseMatrix::<R>::identity(m);
    a0.pad_cols(n);
    b0.pad_cols(n);
    c0.pad_cols(n);
    for row in b0.coeffs.iter_mut() {
        row.clear();
    }
    for row in c0.coeffs.iter_mut() {
        row.clear();
    }

    let mut a1 = SparseMatrix::<R>::identity(m);
    let mut b1 = SparseMatrix::<R>::identity(m);
    let mut c1 = SparseMatrix::<R>::identity(m);
    a1.pad_cols(n);
    b1.pad_cols(n);
    c1.pad_cols(n);
    // Make a1 different from a0 by swapping two rows (still sparse and well-formed).
    a1.coeffs.swap(0, 1);
    for row in b1.coeffs.iter_mut() {
        row.clear();
    }
    for row in c1.coeffs.iter_mut() {
        row.clear();
    }

    let ms: Vec<[std::sync::Arc<SparseMatrix<R>>; 3]> = vec![
        [
            std::sync::Arc::new(a0.clone()),
            std::sync::Arc::new(b0.clone()),
            std::sync::Arc::new(c0.clone()),
        ],
        [
            std::sync::Arc::new(a1.clone()),
            std::sync::Arc::new(b1.clone()),
            std::sync::Arc::new(c1.clone()),
        ],
    ];

    let rg_params = RPParams {
        l_h: 4,
        lambda_pj: 1,
        k_g: 2,
        d_prime: (R::dimension() as u128) - 2,
    };

    let kappa_cm_f = 2usize;

    let f0 = vec![R::ONE; n];
    let f1 = (0..n).map(|i| if i % 2 == 0 { R::ONE } else { R::ZERO }).collect::<Vec<_>>();
    type BF = <<R as PolyRing>::BaseRing as Field>::BasePrimeField;
    let pcs_params_f = {
        let flat_len = n * R::dimension();
        cmf_pcs::cmf_pcs_params_for_flat_len::<BF>(flat_len, kappa_cm_f).expect("cm_f pcs params")
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

    let cms = vec![cm0.clone(), cm1.clone()];
    let witnesses = vec![std::sync::Arc::new(f0.clone()), std::sync::Arc::new(f1.clone())];
    let public_inputs: Vec<<R as PolyRing>::BaseRing> = vec![];

    // CP transcript-message commitment schemes (WE/DPP-facing path).
    let scheme_had = AjtaiCommitmentScheme::<R>::seeded(
        b"cfs_had_u",
        MASTER_SEED,
        2,
        3 * R::dimension(),
    );
    let scheme_mon = AjtaiCommitmentScheme::<R>::seeded(b"cfs_mon_b", MASTER_SEED, 2, rg_params.k_g);
    let scheme_g = AjtaiCommitmentScheme::<R>::seeded(b"cm_g", MASTER_SEED, 2, m * R::dimension());

    let cfg = PiFoldStreamingConfig::default();
    let out = prove_pi_fold_poseidon_fs::<R, PC>(
        ms.as_slice(),
        &cms,
        &witnesses,
        &public_inputs,
        Some(&scheme_had),
        Some(&scheme_mon),
        &scheme_g,
        rg_params.clone(),
        &cfg,
    )
    .expect("prove failed");

    let ms_refs: Vec<[&SparseMatrix<R>; 3]> = vec![
        [&a0, &b0, &c0],
        [&a1, &b1, &c1],
    ];

    // Canonical CP-style verification: verify using transcript-message commitments + openings to aux.
    let open_cfs = MultiAjtaiOpenVerifier::new()
        .with_scheme("cfs_had_u", scheme_had)
        .with_scheme("cfs_mon_b", scheme_mon);

    let attempt = verify_pi_fold_cp_poseidon_fs::<R, PC>(
        PiFoldMatrices::Hetero(ms_refs.as_slice()),
        &cms,
        &out.proof,
        &open_cfs,
        &out.cfs_had_u,
        &out.cfs_mon_b,
        &out.aux,
        &public_inputs,
    );
    attempt.result.expect("cp verify failed");
    let metrics = attempt.metrics;
    let trace = attempt.trace;

    assert_eq!(metrics.absorbed_elems as usize, trace.absorbed.len());
    assert_eq!(metrics.squeezed_field_elems as usize, trace.squeezed_field.len());
    assert_eq!(metrics.squeezed_bytes as usize, trace.squeezed_bytes.len());

    // Replay the trace to ensure it is a complete Poseidon transcript witness.
    let cfg = <PC as cyclotomic_rings::rings::GetPoseidonParams<<<R as PolyRing>::BaseRing as ark_ff::Field>::BasePrimeField>>::get_poseidon_config();
    let replay = replay_poseidon_transcript_trace(&cfg, &trace).expect("poseidon trace replay failed");
    // Sanity: permutation count should match simple upper-bound estimate order of magnitude.
    assert!(replay.permutes.len() > 10);
}

