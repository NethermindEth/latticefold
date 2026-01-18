//! Chunked streaming proving for SP1 shrink verifier R1CS.
//!
//! This program does the following:
//! - open/build the SP1 chunk cache
//! - load all chunk matrices into memory
//! - run the canonical hetero-M streaming Π_fold prover (`prove_pi_fold_poseidon_fs`)
//!
//! Usage:
//!   SP1_R1CS=/path/to/shrink_verifier.r1cs \
//!     cargo run -p symphony --example symphony_sp1_oneproof --release
//!
//! To generate the R1CS file, run in the SP1 fork:
//!   OUT_R1CS=shrink_verifier.r1cs cargo run -p sp1-prover \
//!     --bin dump_shrink_verify_constraints --release

use std::sync::Arc;
use std::time::Instant;

use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
use cyclotomic_rings::rings::GetPoseidonParams;
use latticefold::commitment::AjtaiCommitmentScheme;
use ark_ff::Field;
use symphony::pcs::{cmf_pcs, folding_pcs_l2};
use symphony::pcs::dpp_folding_pcs_l2::folding_pcs_l2_params;
use symphony::pcs::folding_pcs_l2::{
    gadget_apply_digits, kron_ct_in_mul, kron_i_a_mul, kron_ikn_xt_mul, BinMatrix, DenseMatrix,
    FoldingPcsL2ProofCore,
    verify_folding_pcs_l2_with_c_matrices,
};
use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as R;
use stark_rings::PolyRing;
use stark_rings::Ring;
use symphony::rp_rgchk::RPParams;
use symphony::symphony_open::MultiAjtaiOpenVerifier;
use symphony::symphony_pifold_batched::{verify_pi_fold_cp_poseidon_fs, PiFoldMatrices};
use symphony::sp1_r1cs_loader::FieldFromU64;
use symphony::symphony_pifold_streaming::{
    prove_pi_fold_poseidon_fs, PiFoldStreamingConfig,
};
use symphony::symphony_sp1_r1cs::open_sp1_r1cs_chunk_cache;
use symphony::transcript::PoseidonTraceOp;
use symphony::poseidon_trace::find_squeeze_bytes_idx_after_absorb_marker;
use symphony::we_gate_arith::WeGateDr1csBuilder;

/// BabyBear field element for loading R1CS.
#[derive(Debug, Clone, Copy, Default)]
struct BabyBear(u64);

const BABYBEAR_P: u64 = 0x78000001; // 2013265921

impl FieldFromU64 for BabyBear {
    fn from_canonical_u64(val: u64) -> Self {
        BabyBear(val % BABYBEAR_P)
    }
    fn as_canonical_u64(&self) -> u64 {
        self.0
    }
}

fn main() {
    let r1cs_path = std::env::var("SP1_R1CS").expect("Set SP1_R1CS=/path/to/shrink.r1cs");
    let chunk_size: usize = std::env::var("CHUNK_SIZE")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(1 << 20); // 1M

    let l_h: usize = std::env::var("L_H")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(256);
    let lambda_pj: usize = std::env::var("LAMBDA_PJ")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(1);

    let mut cfg = PiFoldStreamingConfig::default();
    cfg.profile = std::env::var("SYMPHONY_PROFILE").ok().as_deref() == Some("1");

    println!("=========================================================");
    println!("Symphony SP1 One-Proof");
    println!("=========================================================");
    println!("  CHUNK_SIZE={chunk_size}  L_H={l_h}  LAMBDA_PJ={lambda_pj}");
    println!(
        "  rayon_threads={}",
        rayon::current_num_threads()
    );

    let t_load = Instant::now();
    let pad_cols_to_multiple_of = l_h;
    let cache = open_sp1_r1cs_chunk_cache::<R, BabyBear>(&r1cs_path, chunk_size, pad_cols_to_multiple_of)
        .expect("Failed to open/build chunk cache");
    println!("  cache open: {:?}", t_load.elapsed());

    let num_chunks = cache.num_chunks;
    let ncols = cache.ncols;
    println!("  chunks={num_chunks} ncols={ncols}");

    // Witness (constant-coeff embedded)
    let mut witness: Vec<R> = vec![R::ZERO; ncols];
    witness[0] = R::ONE;
    let witness = Arc::new(witness);

    // Params
    let rg_params = RPParams {
        l_h,
        lambda_pj,
        k_g: 3,
        d_prime: (R::dimension() as u128) - 2,
    };

    // Commit (`cm_f` is PCS-backed; packed into ring elements for existing Π_fold APIs).
    const MASTER_SEED: [u8; 32] = *b"SYMPHONY_AJTAI_SEED_V1_000000000";
    let kappa_cm_f = 8usize;
    type BF = <<R as PolyRing>::BaseRing as Field>::BasePrimeField;
    let flat_witness: Vec<BF> = witness
        .iter()
        .flat_map(|re| {
            re.coeffs()
                .iter()
                .map(|c| c.to_base_prime_field_elements().into_iter().next().expect("bf limb"))
        })
        .collect();
    let pcs_params_f =
        cmf_pcs::cmf_pcs_params_for_flat_len::<BF>(flat_witness.len(), kappa_cm_f)
            .expect("cm_f pcs params");
    let f_pcs_f = cmf_pcs::pad_flat_message(&pcs_params_f, &flat_witness);
    let (t_pcs_f, _s_pcs_f) =
        folding_pcs_l2::commit(&pcs_params_f, &f_pcs_f).expect("cm_f pcs commit failed");
    let cm_main = cmf_pcs::pack_t_as_ring::<R>(&t_pcs_f);
    let scheme_had = Arc::new(AjtaiCommitmentScheme::<R>::seeded(
        b"cfs_had_u",
        MASTER_SEED,
        8,
        3 * R::dimension(),
    ));
    let scheme_mon = Arc::new(AjtaiCommitmentScheme::<R>::seeded(
        b"cfs_mon_b",
        MASTER_SEED,
        8,
        rg_params.k_g,
    ));

    // Public statement binding (same placeholder as other bench/example)
    let public_inputs: Vec<<R as PolyRing>::BaseRing> = vec![
        <R as PolyRing>::BaseRing::from(1u128),
    ];

    // Load all mats
    let t_mats = Instant::now();
    let mut all_mats: Vec<[Arc<stark_rings_linalg::SparseMatrix<R>>; 3]> = Vec::with_capacity(num_chunks);
    for i in 0..num_chunks {
        let [m1, m2, m3] = cache.read_chunk(i).expect("read_chunk failed");
        all_mats.push([Arc::new(m1), Arc::new(m2), Arc::new(m3)]);
    }
    println!("  load all mats: {:?}", t_mats.elapsed());

    // Commitment scheme for the monomial vectors g^(i): length = m*d, kappa matches other schemes.
    let m = all_mats[0][0].nrows;
    let scheme_g = AjtaiCommitmentScheme::<R>::seeded(b"cm_g", MASTER_SEED, 8, m * R::dimension());

    let cms_all: Vec<Vec<R>> = vec![cm_main; num_chunks];
    // Clone the Arc (refcount bump) so we can still use `witness` later for optional verification.
    let witnesses_all: Vec<Arc<Vec<R>>> = vec![witness.clone(); num_chunks];

    let t_prove = Instant::now();
    let out = prove_pi_fold_poseidon_fs::<R, PC>(
        all_mats.as_slice(),
        &cms_all,
        &witnesses_all,
        &public_inputs,
        Some(scheme_had.as_ref()),
        Some(scheme_mon.as_ref()),
        &scheme_g,
        rg_params,
        &cfg,
    )
    .expect("prove failed");
    println!(
        "  prove total: {:?} (proof_bytes={})",
        t_prove.elapsed(),
        out.proof.coins.bytes.len()
    );

    // Optional: verify with transcript metrics (useful to estimate DPP-friendly cost).
    //
    // This runs the CP/aux verifier path:
    // - Poseidon-FS for challenges/bytes (records metrics + trace)
    // - verifies CP commitments `cfs_*` open to `aux`
    // - runs the core algebraic checks using `aux` (does NOT require full witness openings)
    //
    // Enable with:
    //   VERIFY=1 cargo run -p symphony --example symphony_sp1_oneproof --release
    if std::env::var("VERIFY").ok().as_deref() == Some("1") {
        let ms_ref: Vec<[&stark_rings_linalg::SparseMatrix<R>; 3]> = all_mats
            .iter()
            .map(|ms| [&*ms[0], &*ms[1], &*ms[2]])
            .collect();

        // Opening verifier for CP transcript-message commitments.
        let open_cfs = MultiAjtaiOpenVerifier::<R>::new()
            .with_scheme("cfs_had_u", (*scheme_had).clone())
            .with_scheme("cfs_mon_b", (*scheme_mon).clone());

        let t_vfy = Instant::now();
        let attempt = verify_pi_fold_cp_poseidon_fs::<R, PC>(
            PiFoldMatrices::Hetero(ms_ref.as_slice()),
            &cms_all,
            &out.proof,
            &open_cfs,
            &out.cfs_had_u,
            &out.cfs_mon_b,
            &out.aux,
            &public_inputs,
        );
        let res = attempt.result;
        let metrics = attempt.metrics;
        let trace = attempt.trace;
        println!("  verify (cp/aux): {:?}", t_vfy.elapsed());
        if let Err(e) = &res {
            println!("  verify (cp/aux) failed (expected with dummy witness): {e}");
        }

        // Same estimator we use in other logs: rate=20 field elems; 160 bytes per perm block.
        let perms_absorb = (metrics.absorbed_elems + 19) / 20;
        let perms_squeeze_field = (metrics.squeezed_field_elems + 19) / 20;
        let perms_squeeze_bytes = (metrics.squeezed_bytes + 159) / 160;
        println!(
            "  transcript metrics: absorbed_elems={} squeezed_field_elems={} squeezed_bytes={}",
            metrics.absorbed_elems, metrics.squeezed_field_elems, metrics.squeezed_bytes
        );
        println!(
            "  est poseidon perms: ceil(absorb/20)={} + ceil(sq_field/20)={} + ceil(bytes/160)={} => {}",
            perms_absorb,
            perms_squeeze_field,
            perms_squeeze_bytes,
            perms_absorb + perms_squeeze_field + perms_squeeze_bytes
        );

        // Extra: time the cm_f PCS recompute against the witness (linear work in witness size).
        let t_cm = Instant::now();
        let cm_re = {
            let flat_witness: Vec<BF> = witness
                .iter()
                .flat_map(|re| {
                    re.coeffs()
                        .iter()
                        .map(|c| c.to_base_prime_field_elements().into_iter().next().expect("bf limb"))
                })
                .collect();
            let f_pcs = cmf_pcs::pad_flat_message(&pcs_params_f, &flat_witness);
            let (t, _s) = folding_pcs_l2::commit(&pcs_params_f, &f_pcs).expect("cm_f pcs commit");
            cmf_pcs::pack_t_as_ring::<R>(&t)
        };
        assert_eq!(cm_re, cms_all[0], "cm_f binding mismatch");
        println!("  cm_f pcs recompute: {:?}", t_cm.elapsed());

        // ---------------------------------------------------------------------
        // dR1CS constraint counts for the WE gate (R_cp and full with PCS).
        //
        // This uses the *real verifier transcript trace* to build the dR1CS instance(s),
        // so the counts reflect the actual coin schedule / #rounds for this run's params.
        // ---------------------------------------------------------------------
        type BF = <<R as PolyRing>::BaseRing as Field>::BasePrimeField;
        let poseidon_cfg = <PC as GetPoseidonParams<BF>>::get_poseidon_config();

        let (rcp, rcp_asg) = WeGateDr1csBuilder::r_cp_poseidon_pifold_math_and_cfs_openings::<R>(
            &poseidon_cfg,
            &trace.ops,
            &cms_all,
            &out.proof,
            scheme_had.as_ref(),
            scheme_mon.as_ref(),
            &out.aux,
            &out.cfs_had_u,
            &out.cfs_mon_b,
        )
        .expect("build r_cp dr1cs failed");
        rcp.check(&rcp_asg).expect("r_cp unsat");

        // Build a tiny FoldingPCS(ℓ=2) instance just to measure incremental PCS arithmetization cost.
        // (This is not yet the full production PCS-over-SP1-witness integration.)
        let squeeze_bytes: Vec<Vec<u8>> = trace
            .ops
            .iter()
            .filter_map(|op| match op {
                PoseidonTraceOp::SqueezeBytes { out, .. } => Some(out.clone()),
                _ => None,
            })
            .collect();
        if squeeze_bytes.is_empty() {
            println!("  gate dr1cs: no SqueezeBytes in trace; skipping full-gate (pcs) count");
            return;
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
        let r = 1usize;
        let kappa = 2usize;
        let pcs_n = 4usize;
        // Must satisfy delta^alpha >= modulus (enforced by folding_pcs_l2 exactness guard).
        let delta = 1u64 << 32;
        let alpha = 2usize;
        let beta0 = 1u64 << 10;
        let beta1 = 2 * beta0;
        let beta2 = 2 * beta1;
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
        let mut a_data = vec![<BF as Field>::ZERO; pcs_n * (r * pcs_n * alpha)];
        for i in 0..pcs_n {
            a_data[i * (r * pcs_n * alpha) + i] = <BF as Field>::ONE;
        }
        let a = DenseMatrix::new(pcs_n, r * pcs_n * alpha, a_data);
        let pcs_params = folding_pcs_l2_params(r, kappa, pcs_n, delta, alpha, beta0, beta1, beta2, a);
        let x0 = vec![<BF as Field>::ONE; r];
        let x1 = vec![<BF as Field>::ONE; r];
        let x2 = vec![<BF as Field>::ONE; r];
        let y0 = vec![<BF as Field>::ONE; pcs_params.y0_len()];
        let y1 = kron_ct_in_mul(&c1, pcs_n, &y0);
        let y2 = kron_ct_in_mul(&c2, pcs_n, &y1);
        let t_pcs = kron_i_a_mul(&pcs_params.a, pcs_params.kappa, pcs_params.r * pcs_params.n * pcs_params.alpha, &y0);
        let mut delta_pows = Vec::with_capacity(alpha);
        let mut acc = <BF as Field>::ONE;
        let delta_f = <BF as From<u64>>::from(delta);
        for _ in 0..alpha {
            delta_pows.push(acc);
            acc *= delta_f;
        }
        let v0 = gadget_apply_digits(&delta_pows, r * kappa * pcs_n, &y0);
        let v1 = gadget_apply_digits(&delta_pows, r * kappa * pcs_n, &y1);
        let v2 = gadget_apply_digits(&delta_pows, r * kappa * pcs_n, &y2);
        let u_pcs = kron_ikn_xt_mul(&x2, kappa, pcs_n, &v0);
        let pcs_core = FoldingPcsL2ProofCore { y0, v0, y1, v1, y2, v2 };
        verify_folding_pcs_l2_with_c_matrices(&pcs_params, &t_pcs, &x0, &x1, &x2, &u_pcs, &pcs_core, &c1, &c2)
            .expect("native folding pcs sanity failed");

        let (full, full_asg) = WeGateDr1csBuilder::poseidon_plus_pifold_plus_cfs_plus_pcs::<R>(
            &poseidon_cfg,
            &trace.ops,
            &cms_all,
            &out.proof,
            scheme_had.as_ref(),
            scheme_mon.as_ref(),
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
            &pcs_params,
            &t_pcs,
            &x0,
            &x1,
            &x2,
            &pcs_core,
            pcs_coin_squeeze_idx.saturating_add(1),
        )
        .expect("build full gate dr1cs failed");
        full.check(&full_asg).expect("full gate unsat");

        println!(
            "  gate dr1cs: r_cp(nvars={}, constraints={})  full(nvars={}, constraints={})  delta(nvars={}, constraints={})",
            rcp.nvars,
            rcp.constraints.len(),
            full.nvars,
            full.constraints.len(),
            full.nvars.saturating_sub(rcp.nvars),
            full.constraints.len().saturating_sub(rcp.constraints.len()),
        );
    }
}

