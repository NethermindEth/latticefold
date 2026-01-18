use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
use latticefold::commitment::AjtaiCommitmentScheme;
use stark_rings::{cyclotomic_ring::models::frog_ring::RqPoly as R, PolyRing};
use stark_rings_linalg::SparseMatrix;

use symphony::pcs::folding_pcs_l2::{BinMatrix, verify_folding_pcs_l2_with_c_matrices};
use symphony::rp_rgchk::RPParams;
use symphony::symphony_open::MultiAjtaiOpenVerifier;
use symphony::symphony_pifold_batched::{verify_pi_fold_cp_poseidon_fs, PiFoldMatrices};
use symphony::symphony_pifold_streaming::{prove_pi_fold_poseidon_fs, PiFoldStreamingConfig};
use symphony::transcript::PoseidonTraceOp;
use symphony::we_gate_arith::WeGateDr1csBuilder;
use symphony::pcs::batchlin_pcs::{batchlin_scalar_pcs_params, BATCHLIN_PCS_DOMAIN_SEP};
use symphony::pcs::cmf_pcs::CMF_PCS_DOMAIN_SEP;
use symphony::poseidon_trace::find_squeeze_bytes_idx_after_absorb_marker;
use symphony::pcs::{cmf_pcs, folding_pcs_l2};

use ark_ff::Field;

const MASTER_SEED: [u8; 32] = *b"SYMPHONY_AJTAI_SEED_V1_000000000";

fn bits_le_from_bytes(bytes: &[u8]) -> Vec<bool> {
    let mut out = Vec::with_capacity(bytes.len() * 8);
    for &b in bytes {
        for i in 0..8 {
            out.push(((b >> i) & 1) == 1);
        }
    }
    out
}

#[test]
fn we_gate_constraint_counts_real_pcs_in_gate() {
    // Keep this tiny: we just want a stable, deterministic constraint count smoke test.
    // NOTE: batchlin PCS (ℓ=2) currently requires log2(m*d) divisible by 3.
    // For Frog, d=16, so pick m=32 => m*d=512 => log_n=9.
    // Also keep n >= m so SparseMatrix::pad_cols(n) doesn't need to "shrink" identity matrices.
    let n = 1 << 5; // witness length
    let m = 1 << 5; // rows per instance

    // Two instances with sparse A, and B=C=0 so constraints hold for any witness.
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
    let ms_refs: Vec<[&SparseMatrix<R>; 3]> = vec![[&a0, &b0, &c0], [&a1, &b1, &c1]];

    let rg_params = RPParams {
        l_h: 4,
        lambda_pj: 1,
        k_g: 2,
        d_prime: (R::dimension() as u128) - 2,
    };

    let kappa_cm_f = 2usize;
    let f0 = vec![<R as stark_rings::Ring>::ONE; n];
    let f1 = (0..n)
        .map(|i| {
            if i % 2 == 0 {
                <R as stark_rings::Ring>::ONE
            } else {
                <R as stark_rings::Ring>::ZERO
            }
        })
        .collect::<Vec<_>>();
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
    let cms = vec![cm0, cm1];
    let witnesses = vec![std::sync::Arc::new(f0), std::sync::Arc::new(f1)];
    let public_inputs: Vec<<R as PolyRing>::BaseRing> = vec![];

    // CP transcript-message commitment schemes.
    // We build two copies so one can be moved into `open_cfs` and one can be used by the gate builder.
    let scheme_had_gate = AjtaiCommitmentScheme::<R>::seeded(
        b"cfs_had_u",
        MASTER_SEED,
        2,
        3 * R::dimension(),
    );
    let scheme_mon_gate =
        AjtaiCommitmentScheme::<R>::seeded(b"cfs_mon_b", MASTER_SEED, 2, rg_params.k_g);
    let scheme_had_verifier = AjtaiCommitmentScheme::<R>::seeded(
        b"cfs_had_u",
        MASTER_SEED,
        2,
        3 * R::dimension(),
    );
    let scheme_mon_verifier =
        AjtaiCommitmentScheme::<R>::seeded(b"cfs_mon_b", MASTER_SEED, 2, rg_params.k_g);
    let scheme_g = AjtaiCommitmentScheme::<R>::seeded(b"cm_g", MASTER_SEED, 2, m * R::dimension());

    let cfg = PiFoldStreamingConfig::default();
    let out = prove_pi_fold_poseidon_fs::<R, PC>(
        ms.as_slice(),
        &cms,
        &witnesses,
        &public_inputs,
        Some(&scheme_had_gate),
        Some(&scheme_mon_gate),
        &scheme_g,
        rg_params.clone(),
        &cfg,
    )
    .expect("prove failed");

    // Run the real verifier to obtain the Poseidon trace ops.
    let open_cfs = MultiAjtaiOpenVerifier::new()
        .with_scheme("cfs_had_u", scheme_had_verifier)
        .with_scheme("cfs_mon_b", scheme_mon_verifier);
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
    let trace = attempt.trace;

    let poseidon_cfg = <PC as cyclotomic_rings::rings::GetPoseidonParams<BF>>::get_poseidon_config();

    // ------------------------
    // R_cp (Poseidon + Π_fold-math + cfs-openings)
    // ------------------------
    let (rcp, rcp_asg) = WeGateDr1csBuilder::r_cp_poseidon_pifold_math_and_cfs_openings::<R>(
        &poseidon_cfg,
        &trace.ops,
        &cms,
        &out.proof,
        &scheme_had_gate,
        &scheme_mon_gate,
        &out.aux,
        &out.cfs_had_u,
        &out.cfs_mon_b,
    )
    .expect("build r_cp dr1cs failed");
    rcp.check(&rcp_asg).expect("r_cp unsat");

    // ------------------------
    // Full gate: R_cp × PCS (folding PCS ℓ=2, coins from Poseidon SqueezeBytes)
    // ------------------------
    let squeeze_bytes: Vec<Vec<u8>> = trace
        .ops
        .iter()
        .filter_map(|op| match op {
            PoseidonTraceOp::SqueezeBytes { out, .. } => Some(out.clone()),
            _ => None,
        })
        .collect();
    assert!(
        !squeeze_bytes.is_empty(),
        "expected verifier trace to contain at least one SqueezeBytes op"
    );
    let pcs_coin_squeeze_idx = find_squeeze_bytes_idx_after_absorb_marker(
        &trace.ops,
        BF::from(CMF_PCS_DOMAIN_SEP),
    )
    .expect("missing SqueezeBytes after CMF_PCS_DOMAIN_SEP");
    let c_bytes = &squeeze_bytes[pcs_coin_squeeze_idx];
    let bits = bits_le_from_bytes(c_bytes);

    // PCS#1 (cm_f PCS): MUST match the actual `cm_f` surface absorbed into the transcript.
    // Build a real opening proof for the first instance (cms[0] / witness f0) using the
    // transcript-derived C1/C2 coins (the `SqueezeBytes(N)` right after CMF_PCS_DOMAIN_SEP).
    let pcs_params = pcs_params_f.clone();

    let c1 = BinMatrix {
        rows: pcs_params.r * pcs_params.kappa,
        cols: pcs_params.kappa,
        data: (0..(pcs_params.r * pcs_params.kappa * pcs_params.kappa))
            .map(|i| if bits[i] { BF::ONE } else { BF::ZERO })
            .collect(),
    };
    let c2 = BinMatrix {
        rows: pcs_params.r * pcs_params.kappa,
        cols: pcs_params.kappa,
        data: (0..(pcs_params.r * pcs_params.kappa * pcs_params.kappa))
            .map(|i| {
                if bits[(pcs_params.r * pcs_params.kappa * pcs_params.kappa) + i] {
                    BF::ONE
                } else {
                    BF::ZERO
                }
            })
            .collect(),
    };

    let x0 = vec![BF::ONE; pcs_params.r];
    let x1 = vec![BF::ONE; pcs_params.r];
    let x2 = vec![BF::ONE; pcs_params.r];

    // Flatten witness f0 into BF and pad to PCS message length.
    let flat0: Vec<BF> = witnesses[0]
        .iter()
        .flat_map(|re| {
            re.coeffs()
                .iter()
                .map(|c| c.to_base_prime_field_elements().into_iter().next().expect("bf limb"))
        })
        .collect();
    let f_pcs0 = cmf_pcs::pad_flat_message(&pcs_params, &flat0);
    let (t_pcs, s0) = folding_pcs_l2::commit(&pcs_params, &f_pcs0).expect("cm_f pcs commit");
    assert_eq!(
        cms[0],
        cmf_pcs::pack_t_as_ring::<R>(&t_pcs),
        "cms[0] must equal pack_t_as_ring(commit(f0)) for cm_f PCS binding"
    );
    let (u_pcs, pcs_core) = folding_pcs_l2::open(
        &pcs_params,
        &f_pcs0,
        &s0,
        &x0,
        &x1,
        &x2,
        &c1,
        &c2,
    )
    .expect("cm_f pcs open");
    verify_folding_pcs_l2_with_c_matrices(&pcs_params, &t_pcs, &x0, &x1, &x2, &u_pcs, &pcs_core, &c1, &c2)
        .expect("native cm_f PCS sanity check failed");

    // Use the full trace ops here (production-shape).
    // PCS#2 (batchlin) must use the prover-produced commitment surface and transcript-bound coins.
    let log_n = ((out.proof.m * R::dimension()) as f64).log2() as usize;
    let pcs_params_g = batchlin_scalar_pcs_params::<BF>(log_n).expect("batchlin pcs params");
    let t_pcs_g = &out.proof.batchlin_pcs_t[0];
    let x0_g = &out.proof.batchlin_pcs_x0;
    let x1_g = &out.proof.batchlin_pcs_x1;
    let x2_g = &out.proof.batchlin_pcs_x2;
    let pcs_core_g = &out.proof.batchlin_pcs_core;

    let pcs_coin_squeeze_idx2 = find_squeeze_bytes_idx_after_absorb_marker(
        &trace.ops,
        BF::from(BATCHLIN_PCS_DOMAIN_SEP),
    )
    .expect("missing SqueezeBytes after BATCHLIN_PCS_DOMAIN_SEP");

    let (full, full_asg) = WeGateDr1csBuilder::poseidon_plus_pifold_plus_cfs_plus_pcs::<R>(
        &poseidon_cfg,
        &trace.ops,
        &cms,
        &out.proof,
        &scheme_had_gate,
        &scheme_mon_gate,
        &out.aux,
        &out.cfs_had_u,
        &out.cfs_mon_b,
        &pcs_params,
        &t_pcs,
        &x0,
        &x1,
        &x2,
        &pcs_core,
        pcs_coin_squeeze_idx,
        &pcs_params_g,
        t_pcs_g,
        x0_g,
        x1_g,
        x2_g,
        pcs_core_g,
        pcs_coin_squeeze_idx2,
    )
    .expect("build full gate dr1cs failed");
    full.check(&full_asg).expect("full gate unsat");

    // Print counts for humans (use `-- --nocapture`).
    eprintln!("==============================================================");
    eprintln!("WE gate constraint counts (real PCS-in-gate path)");
    eprintln!("  r_cp:   nvars={} constraints={}", rcp.nvars, rcp.constraints.len());
    eprintln!("  full:   nvars={} constraints={}", full.nvars, full.constraints.len());
    eprintln!("  delta:  nvars={} constraints={}", full.nvars.saturating_sub(rcp.nvars), full.constraints.len().saturating_sub(rcp.constraints.len()));
    eprintln!("==============================================================");

    // Keep a very weak invariant so the test fails if it accidentally becomes trivial.
    assert!(full.constraints.len() > rcp.constraints.len());
}

