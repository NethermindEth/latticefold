use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
use latticefold::commitment::AjtaiCommitmentScheme;
use symphony::{
    pcs::{cmf_pcs, folding_pcs_l2},
    rp_rgchk::RPParams,
    symphony_open::MultiAjtaiOpenVerifier,
    symphony_pifold_batched::{verify_pi_fold_cp_poseidon_fs, PiFoldMatrices},
    symphony_pifold_streaming::{prove_pi_fold_poseidon_fs, PiFoldStreamingConfig},
    symphony_we_relation::{FoldedOutput, ReducedRelation},
};
use stark_rings::{cyclotomic_ring::models::frog_ring::RqPoly as R, PolyRing, Ring, Zq};
use stark_rings_linalg::SparseMatrix;
use ark_ff::Field;

const MASTER_SEED: [u8; 32] = *b"SYMPHONY_AJTAI_SEED_V1_000000000";

struct RoCheckCEquals;

impl ReducedRelation<R> for RoCheckCEquals
where
    <R as PolyRing>::BaseRing: Zq,
{
    type Witness = Vec<R>;

    fn check(public: &FoldedOutput<R>, witness: &Self::Witness) -> Result<(), String> {
        if public.folded_inst.c != *witness {
            return Err("Ro: c mismatch".to_string());
        }
        Ok(())
    }
}

fn setup_two_instances_any_witness_holds() -> (
    [SparseMatrix<R>; 3],
    Vec<Vec<R>>,
    Vec<std::sync::Arc<Vec<R>>>,
    RPParams,
    AjtaiCommitmentScheme<R>,
    AjtaiCommitmentScheme<R>,
    AjtaiCommitmentScheme<R>,
) {
    let n = 1 << 10;
    // NOTE: batchlin PCS (ℓ=2) currently requires log2(m*d) divisible by 3.
    // For Frog, d=16, so pick m=256 => m*d=4096 => log_n=12.
    let m = 1 << 8;

    // Matrices so Π_had holds for any witness: M1=I, M2=0, M3=0.
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

    let f0 = std::sync::Arc::new(vec![R::ONE; n]);
    let f1 = std::sync::Arc::new(vec![R::ZERO; n]);
    let kappa_cm_f = 2usize;
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

    let rg_params = RPParams {
        l_h: 64,
        // Ensure m_J = (n/l_h)*lambda_pj <= m and m multiple of m_J.
        // Here blocks=(1024/64)=16, so lambda_pj=16 gives m_J=256 == m.
        lambda_pj: 16,
        k_g: 4,
        d_prime: (R::dimension() as u128) / 2,
    };

    let scheme_had = AjtaiCommitmentScheme::<R>::seeded(
        b"cfs_had_u",
        MASTER_SEED,
        2,
        3 * R::dimension(),
    );
    let scheme_mon = AjtaiCommitmentScheme::<R>::seeded(b"cfs_mon_b", MASTER_SEED, 2, rg_params.k_g);
    let scheme_g = AjtaiCommitmentScheme::<R>::seeded(b"cm_g", MASTER_SEED, 2, m * R::dimension());

    (
        [m1, m2, m3],
        vec![cm0, cm1],
        vec![f0, f1],
        rg_params,
        scheme_had,
        scheme_mon,
        scheme_g,
    )
}

#[test]
fn test_cp_verify_accepts_and_ro_check_accepts() {
    let (mats, cms, witnesses, rg_params, scheme_had, scheme_mon, scheme_g) =
        setup_two_instances_any_witness_holds();
    let [m1, m2, m3] = mats;

    let cfg = PiFoldStreamingConfig::default();
    let ms = vec![
        [std::sync::Arc::new(m1.clone()), std::sync::Arc::new(m2.clone()), std::sync::Arc::new(m3.clone())],
        [std::sync::Arc::new(m1.clone()), std::sync::Arc::new(m2.clone()), std::sync::Arc::new(m3.clone())],
    ];
    let out = prove_pi_fold_poseidon_fs::<R, PC>(
        ms.as_slice(),
        &cms,
        &witnesses,
        &[],
        Some(&scheme_had),
        Some(&scheme_mon),
        &scheme_g,
        rg_params,
        &cfg,
    )
    .unwrap();

    let open = MultiAjtaiOpenVerifier::<R>::new()
        .with_scheme("cfs_had_u", scheme_had)
        .with_scheme("cfs_mon_b", scheme_mon);

    let attempt = verify_pi_fold_cp_poseidon_fs::<R, PC>(
        PiFoldMatrices::Shared([&m1, &m2, &m3]),
        &cms,
        &out.proof,
        &open,
        &out.cfs_had_u,
        &out.cfs_mon_b,
        &out.aux,
        &[],
    )
    ;
    let (folded_inst, folded_bat) = attempt.result.unwrap();

    let folded = FoldedOutput {
        folded_inst: folded_inst.clone(),
        folded_bat,
    };
    RoCheckCEquals::check(&folded, &folded_inst.c).unwrap();
}

#[test]
fn test_ro_check_rejects_bad_witness() {
    let (mats, cms, witnesses, rg_params, scheme_had, scheme_mon, scheme_g) =
        setup_two_instances_any_witness_holds();
    let [m1, m2, m3] = mats;

    let out = prove_pi_fold_poseidon_fs::<R, PC>(
        &[
            [std::sync::Arc::new(m1.clone()), std::sync::Arc::new(m2.clone()), std::sync::Arc::new(m3.clone())],
            [std::sync::Arc::new(m1.clone()), std::sync::Arc::new(m2.clone()), std::sync::Arc::new(m3.clone())],
        ],
        &cms,
        &witnesses,
        &[],
        Some(&scheme_had),
        Some(&scheme_mon),
        &scheme_g,
        rg_params,
        &PiFoldStreamingConfig::default(),
    )
    .unwrap();

    let open = MultiAjtaiOpenVerifier::<R>::new()
        .with_scheme("cfs_had_u", scheme_had)
        .with_scheme("cfs_mon_b", scheme_mon);
    let attempt = verify_pi_fold_cp_poseidon_fs::<R, PC>(
        PiFoldMatrices::Shared([&m1, &m2, &m3]),
        &cms,
        &out.proof,
        &open,
        &out.cfs_had_u,
        &out.cfs_mon_b,
        &out.aux,
        &[],
    )
    ;
    let (folded_inst, folded_bat) = attempt.result.unwrap();

    let folded = FoldedOutput {
        folded_inst: folded_inst.clone(),
        folded_bat,
    };
    let bad = vec![R::ONE; folded_inst.c.len()];
    let err = RoCheckCEquals::check(&folded, &bad).unwrap_err();
    assert!(err.contains("Ro:"), "unexpected error: {err}");
}

#[test]
fn test_cp_verify_rejects_tampered_cfs_commitment() {
    let (mats, cms, witnesses, rg_params, scheme_had, scheme_mon, scheme_g) =
        setup_two_instances_any_witness_holds();
    let [m1, m2, m3] = mats;

    let out = prove_pi_fold_poseidon_fs::<R, PC>(
        &[
            [std::sync::Arc::new(m1.clone()), std::sync::Arc::new(m2.clone()), std::sync::Arc::new(m3.clone())],
            [std::sync::Arc::new(m1.clone()), std::sync::Arc::new(m2.clone()), std::sync::Arc::new(m3.clone())],
        ],
        &cms,
        &witnesses,
        &[],
        Some(&scheme_had),
        Some(&scheme_mon),
        &scheme_g,
        rg_params,
        &PiFoldStreamingConfig::default(),
    )
    .unwrap();

    let open = MultiAjtaiOpenVerifier::<R>::new()
        .with_scheme("cfs_had_u", scheme_had)
        .with_scheme("cfs_mon_b", scheme_mon);

    // Tamper with a CP commitment after computing it; opening verification must fail.
    let mut bad_cfs_had_u = out.cfs_had_u.clone();
    bad_cfs_had_u[0][0] += R::ONE;

    let err = verify_pi_fold_cp_poseidon_fs::<R, PC>(
        PiFoldMatrices::Shared([&m1, &m2, &m3]),
        &cms,
        &out.proof,
        &open,
        &bad_cfs_had_u,
        &out.cfs_mon_b,
        &out.aux,
        &[],
    )
    .result
    .unwrap_err();
    assert!(err.contains("AjtaiOpen") || err.contains("commitment"), "unexpected error: {err}");
}

