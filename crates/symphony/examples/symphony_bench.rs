//! Micro-benchmark for Symphony Π_fold CP-style verification (Poseidon-FS).
//!
//! This is intentionally lightweight (no criterion) and meant for quick GM-style sanity checks:
//! measure how Π_fold proving/verifying scales with batch size ℓ under the current implementation,
//! including Ajtai opening checks for CP transcript-message commitments (`cfs_*`).

use std::sync::Arc;
use std::time::Instant;

fn main() {
    use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
    use latticefold::commitment::AjtaiCommitmentScheme;
    use symphony::{
        pcs::batchlin_pcs::{batchlin_scalar_pcs_params, BATCHLIN_PCS_DOMAIN_SEP},
        pcs::folding_pcs_l2::{BinMatrix, verify_folding_pcs_l2_with_c_matrices},
        pcs::{cmf_pcs, folding_pcs_l2},
        rp_rgchk::RPParams,
        symphony_open::MultiAjtaiOpenVerifier,
        symphony_pifold_batched::{verify_pi_fold_cp_poseidon_fs, PiFoldMatrices},
        symphony_pifold_streaming::{prove_pi_fold_poseidon_fs, PiFoldStreamingConfig},
        symphony_we_relation::{FoldedOutput, TrivialRo},
        transcript::PoseidonTraceOp,
        poseidon_trace::find_squeeze_bytes_idx_after_absorb_marker,
        we_gate_arith::WeGateDr1csBuilder,
    };
    use cyclotomic_rings::rings::GetPoseidonParams;
    use ark_ff::Field;
    use stark_rings::{cyclotomic_ring::models::frog_ring::RqPoly as R, PolyRing, Ring};
    use stark_rings_linalg::SparseMatrix;

    // Keep these at the test-friendly "toy" scale for interactive runs. Increase if you want.
    let n = 1 << 10;
    let m = 1 << 10;

    // Non-vacuous Π_had: M1=M2=M3=I.
    let mut m1 = SparseMatrix::<R>::identity(m);
    let mut m2 = SparseMatrix::<R>::identity(m);
    let mut m3 = SparseMatrix::<R>::identity(m);
    m1.pad_cols(n);
    m2.pad_cols(n);
    m3.pad_cols(n);
    let m1 = Arc::new(m1);
    let m2 = Arc::new(m2);
    let m3 = Arc::new(m3);

    let rg_params = RPParams {
        l_h: 64,
        lambda_pj: 32,
        k_g: 4,
        d_prime: (R::dimension() as u128) / 2,
    };

    // "Statement binding" placeholder (e.g. vk hash + claim digest).
    let public_inputs: Vec<<R as PolyRing>::BaseRing> = vec![
        <R as PolyRing>::BaseRing::from(123u128),
        <R as PolyRing>::BaseRing::from(456u128),
    ];

    // cm_f commitment surface: PCS-backed (packed into ring elements for existing Π_fold APIs).
    const MASTER_SEED: [u8; 32] = *b"SYMPHONY_AJTAI_SEED_V1_000000000";
    let kappa_cm_f = 8usize;

    // Try a few batch sizes ℓ.
    for &ell in &[2usize, 4, 8] {
        // Build witnesses: alternate 1 and 0 (both idempotent).
        let mut witnesses: Vec<Arc<Vec<R>>> = Vec::with_capacity(ell);
        let mut cm_f: Vec<Vec<R>> = Vec::with_capacity(ell);
        type BF = <<R as PolyRing>::BaseRing as Field>::BasePrimeField;
        // Save PCS#1 (cm_f PCS) material for instance 0 so the full WE gate can bind it.
        let mut pcs_params_f0: Option<symphony::pcs::folding_pcs_l2::FoldingPcsL2Params<BF>> = None;
        let mut f_pcs_f0: Option<Vec<BF>> = None;
        let mut s_pcs_f0: Option<[Vec<BF>; 3]> = None;
        let mut t_pcs_f0: Option<Vec<BF>> = None;
        for i in 0..ell {
            let f = if i % 2 == 0 {
                vec![R::ONE; n]
            } else {
                vec![R::ZERO; n]
            };
            let flat_witness: Vec<BF> = f
                .iter()
                .flat_map(|re| {
                    re.coeffs()
                        .iter()
                        .map(|c| c.to_base_prime_field_elements().into_iter().next().expect("bf limb"))
                })
                .collect();
            let pcs_params_f = cmf_pcs::cmf_pcs_params_for_flat_len::<BF>(
                flat_witness.len(),
                kappa_cm_f,
            )
            .expect("cm_f pcs params");
            let f_pcs_f = cmf_pcs::pad_flat_message(&pcs_params_f, &flat_witness);
            let (t_pcs_f, s_pcs_f) =
                folding_pcs_l2::commit(&pcs_params_f, &f_pcs_f).expect("cm_f pcs commit failed");
            let cm = cmf_pcs::pack_t_as_ring::<R>(&t_pcs_f);
            if i == 0 {
                pcs_params_f0 = Some(pcs_params_f.clone());
                f_pcs_f0 = Some(f_pcs_f.clone());
                s_pcs_f0 = Some(s_pcs_f);
                t_pcs_f0 = Some(t_pcs_f.clone());
            }
            witnesses.push(Arc::new(f));
            cm_f.push(cm);
        }

        // Build Ms: shared matrices repeated for each instance.
        let ms: Vec<[Arc<SparseMatrix<R>>; 3]> = (0..ell)
            .map(|_| [m1.clone(), m2.clone(), m3.clone()])
            .collect();

        // CP commitments for aux messages (cfs_*). Use a larger row count so we can see the
        // Ajtai-open verification cost.
        let scheme_had =
            AjtaiCommitmentScheme::<R>::seeded(b"cfs_had_u", MASTER_SEED, 32, 3 * R::dimension());
        let scheme_mon =
            AjtaiCommitmentScheme::<R>::seeded(b"cfs_mon_b", MASTER_SEED, 32, rg_params.k_g);
        let scheme_g =
            AjtaiCommitmentScheme::<R>::seeded(b"cm_g", MASTER_SEED, 32, m * R::dimension());

        let cfg = PiFoldStreamingConfig::default();
        let prove_start = Instant::now();
        let out = prove_pi_fold_poseidon_fs::<R, PC>(
            ms.as_slice(),
            &cm_f,
            &witnesses,
            &public_inputs,
            Some(&scheme_had),
            Some(&scheme_mon),
            &scheme_g,
            rg_params.clone(),
            &cfg,
        )
        .unwrap();
        let prove_time = prove_start.elapsed();

        let open = MultiAjtaiOpenVerifier::<R>::new()
            .with_scheme("cfs_had_u", scheme_had)
            .with_scheme("cfs_mon_b", scheme_mon);

        let verify_start = Instant::now();
        let attempt = verify_pi_fold_cp_poseidon_fs::<R, PC>(
            PiFoldMatrices::Shared([&*m1, &*m2, &*m3]),
            &cm_f,
            &out.proof,
            &open,
            &out.cfs_had_u,
            &out.cfs_mon_b,
            &out.aux,
            &public_inputs,
        )
        ;
        let (folded_inst, folded_bat) = attempt.result.unwrap();
        let verify_time = verify_start.elapsed();
        let trace = attempt.trace;

        // Sanity: `TrivialRo` exists just to show how R_WE would be checked; we don't benchmark it.
        let folded_out = FoldedOutput {
            folded_inst,
            folded_bat,
        };
        let _: Result<(), String> =
            <TrivialRo as symphony::symphony_we_relation::ReducedRelation<R>>::check(&folded_out, &());

        println!(
            "ℓ={ell:>2} | Π_fold prove: {:>8.3}s | R_cp verify: {:>8.3}s | proof bytes: {}",
            prove_time.as_secs_f64(),
            verify_time.as_secs_f64(),
            out.proof.coins.bytes.len()
        );

        // Gate size report: build r_cp and full (real PCS-in-gate).
        let poseidon_cfg = <PC as GetPoseidonParams<BF>>::get_poseidon_config();

        let (rcp, rcp_asg) = WeGateDr1csBuilder::r_cp_poseidon_pifold_math_and_cfs_openings::<R>(
            &poseidon_cfg,
            &trace.ops,
            &cm_f,
            &out.proof,
            &open.schemes["cfs_had_u"],
            &open.schemes["cfs_mon_b"],
            &out.aux,
            &out.cfs_had_u,
            &out.cfs_mon_b,
        )
        .expect("build r_cp dr1cs failed");
        rcp.check(&rcp_asg).expect("r_cp unsat");

        let squeeze_bytes: Vec<Vec<u8>> = trace
            .ops
            .iter()
            .filter_map(|op| match op {
                PoseidonTraceOp::SqueezeBytes { out, .. } => Some(out.clone()),
                _ => None,
            })
            .collect();
        if squeeze_bytes.is_empty() {
            println!("    gate dr1cs: r_cp(nvars={}, constraints={}) (no SqueezeBytes => skip full pcs)",
                rcp.nvars, rcp.constraints.len()
            );
            continue;
        }

        let pcs_coin_squeeze_idx = find_squeeze_bytes_idx_after_absorb_marker(
            &trace.ops,
            BF::from(cmf_pcs::CMF_PCS_DOMAIN_SEP),
        )
        .expect("missing SqueezeBytes after CMF_PCS_DOMAIN_SEP");
        let c_bytes = &squeeze_bytes[pcs_coin_squeeze_idx];
        let mut bits = Vec::with_capacity(c_bytes.len() * 8);
        for &b in c_bytes {
            for i in 0..8 {
                bits.push(((b >> i) & 1) == 1);
            }
        }
        // Build a real PCS#1 opening proof consistent with the absorbed cm_f surface.
        let pcs_params = pcs_params_f0.clone().expect("missing pcs_params_f0");
        let f_pcs = f_pcs_f0.clone().expect("missing f_pcs_f0");
        let s_pcs = s_pcs_f0.clone().expect("missing s_pcs_f0");
        let t_pcs = t_pcs_f0.clone().expect("missing t_pcs_f0");
        let r = pcs_params.r;
        let kappa = pcs_params.kappa;
        let c1 = BinMatrix {
            rows: r * kappa,
            cols: kappa,
            data: (0..(r * kappa * kappa))
                .map(|i| if bits[i] { <BF as Field>::ONE } else { <BF as Field>::ZERO })
                .collect(),
        };
        let c2 = BinMatrix {
            rows: r * kappa,
            cols: kappa,
            data: (0..(r * kappa * kappa))
                .map(|i| {
                    if bits[(r * kappa * kappa) + i] {
                        <BF as Field>::ONE
                    } else {
                        <BF as Field>::ZERO
                    }
                })
                .collect(),
        };
        let x0 = vec![<BF as Field>::ONE; r];
        let x1 = vec![<BF as Field>::ONE; r];
        let x2 = vec![<BF as Field>::ONE; r];
        let (u_pcs, pcs_core) =
            folding_pcs_l2::open(&pcs_params, &f_pcs, &s_pcs, &x0, &x1, &x2, &c1, &c2)
                .expect("cm_f pcs open failed");
        verify_folding_pcs_l2_with_c_matrices(&pcs_params, &t_pcs, &x0, &x1, &x2, &u_pcs, &pcs_core, &c1, &c2)
            .expect("native folding pcs sanity failed");

        let pcs_coin_squeeze_idx2 = find_squeeze_bytes_idx_after_absorb_marker(
            &trace.ops,
            BF::from(BATCHLIN_PCS_DOMAIN_SEP),
        )
        .expect("missing SqueezeBytes after BATCHLIN_PCS_DOMAIN_SEP");
        let log_n = ((out.proof.m * R::dimension()) as f64).log2() as usize;
        let pcs_params_g = batchlin_scalar_pcs_params::<BF>(log_n).expect("batchlin pcs params");
        let t_pcs_g = &out.proof.batchlin_pcs_t[0];
        let x0_g = &out.proof.batchlin_pcs_x0;
        let x1_g = &out.proof.batchlin_pcs_x1;
        let x2_g = &out.proof.batchlin_pcs_x2;
        let pcs_core_g = &out.proof.batchlin_pcs_core;

        let (full, full_asg) = WeGateDr1csBuilder::poseidon_plus_pifold_plus_cfs_plus_pcs::<R>(
            &poseidon_cfg,
            &trace.ops,
            &cm_f,
            &out.proof,
            &open.schemes["cfs_had_u"],
            &open.schemes["cfs_mon_b"],
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

        println!(
            "    gate dr1cs: r_cp(nvars={}, constraints={})  full(nvars={}, constraints={})  delta(nvars={}, constraints={})",
            rcp.nvars,
            rcp.constraints.len(),
            full.nvars,
            full.constraints.len(),
            full.nvars.saturating_sub(rcp.nvars),
            full.constraints.len().saturating_sub(rcp.constraints.len()),
        );
    }
}
