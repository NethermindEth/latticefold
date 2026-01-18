#![allow(non_local_definitions)]

use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
use latticefold::commitment::AjtaiCommitmentScheme;
use stark_rings::{cyclotomic_ring::models::frog_ring::RqPoly as R, PolyRing};
use stark_rings_linalg::SparseMatrix;

use symphony::dpp_poseidon::poseidon_sponge_dr1cs_from_trace;
use symphony::we_gate_arith::WeGateDr1csBuilder;
use symphony::poseidon_trace::find_squeeze_bytes_idx_after_absorb_marker;
use symphony::poseidon_trace::replay_poseidon_transcript_trace;
use symphony::rp_rgchk::RPParams;
use symphony::symphony_open::MultiAjtaiOpenVerifier;
use symphony::symphony_pifold_batched::{verify_pi_fold_cp_poseidon_fs, PiFoldMatrices};
use symphony::symphony_pifold_streaming::{prove_pi_fold_poseidon_fs, PiFoldStreamingConfig};
use symphony::transcript::PoseidonTraceOp;
use symphony::pcs::cmf_pcs::CMF_PCS_DOMAIN_SEP;
use symphony::pcs::folding_pcs_l2::{
    BinMatrix, verify_folding_pcs_l2_with_c_matrices,
};
use symphony::pcs::{cmf_pcs, folding_pcs_l2};
use dpp::{
    dr1cs_flpcp::{Dr1csInstanceSparse as DppDr1csInstanceSparse, RsDr1csFlpcpSparse, RsDr1csNpFlpcpSparse},
    embedding::EmbeddingParams,
    packing::{BoundedFlpcpSparse, PackedDppParams},
    sparse::SparseVec,
    pipeline::build_rev2_dpp_sparse_boolean,
};
use ark_ff::{BigInteger, Field, Fp, Fp256, MontBackend, MontConfig, PrimeField};
use rand::{rngs::StdRng, SeedableRng};

const MASTER_SEED: [u8; 32] = *b"SYMPHONY_AJTAI_SEED_V1_000000000";

/// End-to-end harness: prove Π_fold → verify (WE/CP path) → extract Poseidon trace →
/// arithmetize that trace into sparse dR1CS → check satisfiable.
///
/// This used to time out when the arithmetization builder densified matrices. With the
/// sparse frontend, it should be fast enough to run normally.
#[test]
fn test_poseidon_trace_arithmetization_roundtrip_real_verifier() {
    // NOTE: batchlin PCS (ℓ=2, r^3=n) currently requires log2(m*d) divisible by 3.
    // For Frog, d=16, so pick m=32 => m*d=512 => log2=9.
    let n = 1 << 5; // 32 vars
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
    let ms_refs: Vec<[&SparseMatrix<R>; 3]> = vec![[&a0, &b0, &c0], [&a1, &b1, &c1]];

    let rg_params = RPParams {
        l_h: 4,
        lambda_pj: 1,
        k_g: 2,
        d_prime: (R::dimension() as u128) - 2,
    };

    // Choose `kappa_commit` in base-field elements (not ring elements) to match the packed `cm_f`
    // surface length used in the transcript: `kappa_commit = cm_f.len() * d`.
    // Since `cm_f` is packed into ring elements of dimension d, using `kappa_commit=d` keeps
    // the `cm_f` surface unambiguous (no trailing zero coefficient padding).
    let kappa_cm_f = R::dimension();
    let f0 = vec![<R as stark_rings::Ring>::ONE; n];
    let f1 = (0..n)
        .map(|i| if i % 2 == 0 { <R as stark_rings::Ring>::ONE } else { <R as stark_rings::Ring>::ZERO })
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
    let scheme_mon =
        AjtaiCommitmentScheme::<R>::seeded(b"cfs_mon_b", MASTER_SEED, 2, rg_params.k_g);
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

    // WE/CP-facing verification: verify using CP transcript-message commitments + openings to aux.
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

    // Replay check: the trace must be a valid Poseidon transcript witness.
    let poseidon_cfg =
        <PC as cyclotomic_rings::rings::GetPoseidonParams<BF>>::get_poseidon_config();
    let _replay = replay_poseidon_transcript_trace(&poseidon_cfg, &trace)
        .expect("poseidon trace replay failed");

    // Arithmetize the trace into sparse dR1CS and check satisfiable.
    let (inst, assignment, _replay2, _byte_wit) =
        poseidon_sponge_dr1cs_from_trace::<BF>(&poseidon_cfg, &trace.ops)
            .expect("build poseidon-sponge dr1cs failed");
    inst.check(&assignment).expect("dR1CS not satisfied");
}

// 256-bit prime field for the large-field embedding path.
// secp256k1 prime: 2^256 - 2^32 - 977
#[derive(MontConfig)]
#[modulus = "115792089237316195423570985008687907852837564279074904382605163141518161494337"]
#[generator = "7"]
pub struct Secp256k1Config;
type FLarge = Fp256<MontBackend<Secp256k1Config, 4>>;

// 521-bit Mersenne prime field for “no-booleanization” full-scale packing tests.
// p = 2^521 - 1
#[derive(MontConfig)]
#[modulus = "6864797660130609714981900799081393217269435300143305409394463459185543183397656052122559640661454554977296311391480858037121987999716643812574028291115057151"]
#[generator = "2"]
pub struct P521Config;
type FPack521 = Fp<MontBackend<P521Config, 9>, 9>;

/// Full sparse “Embed → Pack” check on a real verifier-derived Poseidon trace.
///
/// This exercises:
/// - Symphony verifier trace extraction
/// - Poseidon arithmetization into sparse dR1CS
/// - sparse RS dR1CS FLPCP
/// - sparse bounded embedding
/// - sparse query packing into 1-query DPP
/// - `sample_query` + `verify_with_query`
#[test]
fn test_poseidon_trace_sparse_dpp_end_to_end_accepts() {
    // Reuse the existing harness to produce a real Poseidon transcript trace.
    // (We keep parameters tiny so this runs in CI.)
    let n = 1 << 5;
    let m = 1 << 5;

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
        .map(|i| if i % 2 == 0 { <R as stark_rings::Ring>::ONE } else { <R as stark_rings::Ring>::ZERO })
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

    let scheme_had = AjtaiCommitmentScheme::<R>::seeded(
        b"cfs_had_u",
        MASTER_SEED,
        2,
        3 * R::dimension(),
    );
    let scheme_mon =
        AjtaiCommitmentScheme::<R>::seeded(b"cfs_mon_b", MASTER_SEED, 2, rg_params.k_g);
    let scheme_g = AjtaiCommitmentScheme::<R>::seeded(b"cm_g", MASTER_SEED, 2, m * R::dimension());

    let out = prove_pi_fold_poseidon_fs::<R, PC>(
        ms.as_slice(),
        &cms,
        &witnesses,
        &public_inputs,
        Some(&scheme_had),
        Some(&scheme_mon),
        &scheme_g,
        rg_params.clone(),
        &PiFoldStreamingConfig::default(),
    )
    .expect("prove failed");

    let open_cfs = MultiAjtaiOpenVerifier::new()
        .with_scheme("cfs_had_u", scheme_had.clone())
        .with_scheme("cfs_mon_b", scheme_mon.clone());
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

    // IMPORTANT: for this DPP integration test we only take a short prefix of the transcript ops,
    // to keep the resulting dR1CS small (and thus keep the RS FLPCP proof length `m=2k` small).
    //
    // The correctness property we want to exercise here is the end-to-end plumbing:
    //   (real verifier trace prefix) -> arithmetize -> RS-FLPCP -> embed -> pack -> verify_with_query.
    let take_ops = core::cmp::min(25usize, trace.ops.len());
    let ops = &trace.ops[..take_ops];

    // Arithmetize the trace prefix into sparse dR1CS (over BF).
    let poseidon_cfg =
        <PC as cyclotomic_rings::rings::GetPoseidonParams<BF>>::get_poseidon_config();

    // Poseidon-only arithmetization (used below for the RS-FLPCP plumbing part of this test).
    let (inst, assignment, _replay2, _byte_wit) =
        poseidon_sponge_dr1cs_from_trace::<BF>(&poseidon_cfg, ops)
            .expect("build poseidon-sponge dr1cs failed");
    inst.check(&assignment).expect("poseidon dR1CS not satisfied");

    // Poseidon + Ajtai-open(cfs_*) arithmetization (this is the WE gate “binding” part).
    let (merged, merged_asg) = WeGateDr1csBuilder::poseidon_plus_cfs_openings::<R>(
        &poseidon_cfg,
        ops,
        &scheme_had,
        &scheme_mon,
        &out.aux,
        &out.cfs_had_u,
        &out.cfs_mon_b,
    )
    .expect("build merged (poseidon+ajtai-open) dr1cs failed");
    merged.check(&merged_asg).expect("merged (poseidon+ajtai) dr1cs unsat");

    // Negative test: tamper with a CP commitment (but keep the same opened message).
    // The resulting merged system must become unsatisfiable.
    if !cms.is_empty() {
        let i = 0usize;
        // Rebuild using a tampered commitment for cfs_had_u[0] while keeping the same opened aux,
        // and ensure the merged system becomes UNSAT.
        let mut bad_cfs_had_u = out.cfs_had_u.clone();
        bad_cfs_had_u[i][0] += <R as stark_rings::Ring>::ONE;
        let (bad_merged, bad_asg) = WeGateDr1csBuilder::poseidon_plus_cfs_openings::<R>(
            &poseidon_cfg,
            ops,
            &scheme_had,
            &scheme_mon,
            &out.aux,
            &bad_cfs_had_u,
            &out.cfs_mon_b,
        )
        .expect("build bad merged dr1cs failed");
        assert!(bad_merged.check(&bad_asg).is_err(), "tampered Ajtai-open should be UNSAT");
    }

    // Production-shape binding (current architecture): use PCS-in-gate instead of AjtaiOpen(cm_f, f).
    //
    // For this integration test we build a **real cm_f PCS opening** for the first instance,
    // and take its `C1/C2` coins from the *real* verifier transcript:
    //   `CMF_PCS_DOMAIN_SEP` -> `SqueezeBytes(n_bytes_cmf)`.
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
    // Little-endian bits per byte (same convention as `bytes_to_bitstream_le_dr1cs`).
    let mut bits = Vec::with_capacity(c_bytes.len() * 8);
    for &b in c_bytes {
        for i in 0..8 {
            bits.push(((b >> i) & 1) == 1);
        }
    }

    // Recompute the cm_f PCS commitment surface for the first instance and generate a valid opening.
    let flat0: Vec<BF> = witnesses[0]
        .iter()
        .flat_map(|re| {
            re.coeffs()
                .iter()
                .map(|c| c.to_base_prime_field_elements().into_iter().next().expect("bf limb"))
        })
        .collect();
    let f_pcs0 = cmf_pcs::pad_flat_message(&pcs_params_f, &flat0);
    let (t_pcs, s_pcs) = folding_pcs_l2::commit(&pcs_params_f, &f_pcs0).expect("cm_f pcs commit");

    let r = pcs_params_f.r;
    let kappa = pcs_params_f.kappa;

    let need_bits = 2usize * r * kappa * kappa;
    assert!(
        bits.len() >= need_bits,
        "expected {} bits for cm_f PCS C1/C2, got {}",
        need_bits,
        bits.len()
    );
    let c1 = BinMatrix {
        rows: r * kappa,
        cols: kappa,
        data: (0..(r * kappa * kappa))
            .map(|i| if bits[i] { BF::ONE } else { BF::ZERO })
            .collect(),
    };
    let c2 = BinMatrix {
        rows: r * kappa,
        cols: kappa,
        data: (0..(r * kappa * kappa))
            .map(|i| if bits[(r * kappa * kappa) + i] { BF::ONE } else { BF::ZERO })
            .collect(),
    };

    let x0 = vec![BF::ONE; r];
    let x1 = vec![BF::ONE; r];
    let x2 = vec![BF::ONE; r];
    let (u_pcs, pcs_core) =
        folding_pcs_l2::open(&pcs_params_f, &f_pcs0, &s_pcs, &x0, &x1, &x2, &c1, &c2)
            .expect("cm_f pcs open");
    verify_folding_pcs_l2_with_c_matrices(
        &pcs_params_f,
        &t_pcs,
        &x0,
        &x1,
        &x2,
        &u_pcs,
        &pcs_core,
        &c1,
        &c2,
    )
    .expect("native cm_f pcs sanity check failed");

    // PCS#2 (batchlin) uses a real prover-produced commitment surface + proof core.
    use symphony::pcs::batchlin_pcs::{batchlin_scalar_pcs_params, BATCHLIN_PCS_DOMAIN_SEP};
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
    .expect("missing SqueezeBytes after BATCHLIN_PCS_DOMAIN_SEP (PCS#2 coin source)");

    let (merged3, merged3_asg) = WeGateDr1csBuilder::poseidon_plus_pifold_plus_cfs_plus_pcs::<R>(
        &poseidon_cfg,
        &trace.ops,
        &cms,
        &out.proof,
        &scheme_had,
        &scheme_mon,
        &out.aux,
        &out.cfs_had_u,
        &out.cfs_mon_b,
        &pcs_params_f,
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
    .expect("build merged (poseidon+pifold-math+cfs+pcs) failed");
    if let Err(e) = merged3.check(&merged3_asg) {
        // Helpful debug on failures: print the failing constraint and its evaluated values.
        let idx = e
            .split_whitespace()
            .find_map(|tok| tok.parse::<usize>().ok())
            .unwrap_or(0);
        let row = &merged3.constraints[idx];
        let eval_lc = |lc: &Vec<(BF, usize)>| -> BF {
            lc.iter().fold(BF::ZERO, |acc, (c, i)| acc + (*c) * merged3_asg[*i])
        };
        let a = eval_lc(&row.a);
        let b = eval_lc(&row.b);
        let c = eval_lc(&row.c);
        let a_terms = row
            .a
            .iter()
            .take(8)
            .map(|(coef, vi)| format!("({coef:?} * asg[{vi}]={:?})", merged3_asg[*vi]))
            .collect::<Vec<_>>()
            .join(" + ");
        let c_terms = row
            .c
            .iter()
            .take(8)
            .map(|(coef, vi)| format!("({coef:?} * asg[{vi}]={:?})", merged3_asg[*vi]))
            .collect::<Vec<_>>()
            .join(" + ");
        panic!(
            "merged (poseidon+pifold-math+cfs+pcs) dr1cs unsat: {e}\n  idx={idx}\n  a={a:?}\n  b={b:?}\n  c={c:?}\n  a*b={ab:?}\n  row.a.len={la} row.b.len={lb} row.c.len={lc}\n  row.a(head)={a_terms}\n  row.c(head)={c_terms}",
            ab = a * b,
            la = row.a.len(),
            lb = row.b.len(),
            lc = row.c.len(),
        );
    }

    // Convert Symphony sparse dR1CS to DPP sparse dR1CS rows.
    //
    // NOTE: The RS-FLPCP prover is now (asymptotically) **~O(k log k)** via consecutive-point
    // extrapolation using convolution (FFT when the field supports a radix-2 domain, otherwise
    // CRT+NTT fallback). We still keep `k` small here to keep this integration test fast in CI.
    let k_full = inst.constraints.len();
    let k = core::cmp::min(64usize, k_full);
    let mut a_rows = Vec::with_capacity(k);
    let mut b_rows = Vec::with_capacity(k);
    let mut c_rows = Vec::with_capacity(k);
    for row in inst.constraints.iter().take(k) {
        a_rows.push(SparseVec::new(row.a.clone()));
        b_rows.push(SparseVec::new(row.b.clone()));
        c_rows.push(SparseVec::new(row.c.clone()));
    }
    // Compress variable space to only variables actually used by these k constraints.
    // This keeps the NP witness length small so the embedding+packing modulus condition is feasible.
    let mut used = std::collections::BTreeSet::<usize>::new();
    for row in a_rows.iter().chain(b_rows.iter()).chain(c_rows.iter()) {
        for (_, idx) in row.terms.iter() {
            used.insert(*idx);
        }
    }
    // Ensure constant slot 0 is included.
    used.insert(0);
    let used_list = used.into_iter().collect::<Vec<_>>();
    let mut map = std::collections::BTreeMap::<usize, usize>::new();
    for (new, old) in used_list.iter().enumerate() {
        map.insert(*old, new);
    }

    let assignment_reduced = used_list.iter().map(|old| assignment[*old]).collect::<Vec<_>>();
    let remap = |row: &SparseVec<BF>| -> SparseVec<BF> {
        SparseVec::new(
            row.terms
                .iter()
                .map(|(c, idx)| (*c, *map.get(idx).expect("idx in map")))
                .collect::<Vec<_>>(),
        )
    };
    let a_rows = a_rows.iter().map(remap).collect::<Vec<_>>();
    let b_rows = b_rows.iter().map(remap).collect::<Vec<_>>();
    let c_rows = c_rows.iter().map(remap).collect::<Vec<_>>();

    let dr1cs = DppDr1csInstanceSparse::<BF> { n: assignment_reduced.len(), a: a_rows, b: b_rows, c: c_rows };

    // NP-style RS FLPCP over BF: take x = [] (public), witness z_w = full assignment.
    //
    // This matches the WE setting: statement is public, and the witness is private.
    let ell = 2 * k; // minimal, good enough for this test harness
    let flpcp_small = RsDr1csNpFlpcpSparse::<BF>::new(dr1cs, 0, ell);
    let x_small: Vec<BF> = vec![];
    let pi_field = flpcp_small.prove(&x_small, &assignment_reduced);
    assert_eq!(pi_field.len(), assignment_reduced.len() + 2 * k);

    // For this small integration test we keep the Rev2 “Booleanize → Embed → Pack” path,
    // so that the bounded decoding condition holds over a 256-bit packed field.
    let boolized = dpp::BooleanProofFlpcpSparse::<BF, _>::new(flpcp_small.clone());
    let pi_bits = boolized.encode_proof_bits(&pi_field);

    let dpp = build_rev2_dpp_sparse_boolean::<BF, FLarge, _>(
        flpcp_small,
        EmbeddingParams { gamma: 2, assume_boolean_proof: true, k_prime: 0 },
        PackedDppParams { ell: 2 },
    );

    // Lift x and π into FLarge.
    let x_large: Vec<FLarge> = vec![];
    let pi_large = pi_bits
        .iter()
        .map(|wi| FLarge::from_le_bytes_mod_order(&wi.into_bigint().to_bytes_le()))
        .collect::<Vec<_>>();

    // Sample a packed query and verify with the provided proof (arm/decap style).
    let mut rng = StdRng::seed_from_u64(12345);
    let q = dpp.sample_query(&mut rng, &x_large).expect("sample_query");
    let ok = dpp.verify_with_query(&x_large, &pi_large, &q).expect("verify_with_query");
    assert!(ok);
}

/// Full-scale RS-FLPCP prover check on the full Poseidon trace arithmetization.
///
/// This flips the “take_ops / k caps” and is meant to be run in release:
/// `cargo test -p symphony --release -- --ignored test_poseidon_trace_rs_flpcp_full_trace_honest_accepts`
#[test]
#[ignore]
fn test_poseidon_trace_rs_flpcp_full_trace_honest_accepts() {
    // Produce a real verifier trace (same harness as the other tests).
    let n = 1 << 5;
    let m = 1 << 5;

    let mut a0 = SparseMatrix::<R>::identity(m);
    let mut b0 = SparseMatrix::<R>::identity(m);
    let mut c0 = SparseMatrix::<R>::identity(m);
    a0.pad_cols(n);
    b0.pad_cols(n);
    c0.pad_cols(n);
    for row in b0.coeffs.iter_mut() { row.clear(); }
    for row in c0.coeffs.iter_mut() { row.clear(); }

    let mut a1 = SparseMatrix::<R>::identity(m);
    let mut b1 = SparseMatrix::<R>::identity(m);
    let mut c1 = SparseMatrix::<R>::identity(m);
    a1.pad_cols(n);
    b1.pad_cols(n);
    c1.pad_cols(n);
    a1.coeffs.swap(0, 1);
    for row in b1.coeffs.iter_mut() { row.clear(); }
    for row in c1.coeffs.iter_mut() { row.clear(); }

    let ms: Vec<[std::sync::Arc<SparseMatrix<R>>; 3]> = vec![
        [std::sync::Arc::new(a0.clone()), std::sync::Arc::new(b0.clone()), std::sync::Arc::new(c0.clone())],
        [std::sync::Arc::new(a1.clone()), std::sync::Arc::new(b1.clone()), std::sync::Arc::new(c1.clone())],
    ];
    let ms_refs: Vec<[&SparseMatrix<R>; 3]> = vec![[&a0, &b0, &c0], [&a1, &b1, &c1]];

    let rg_params = RPParams { l_h: 4, lambda_pj: 1, k_g: 2, d_prime: (R::dimension() as u128) - 2 };

    let kappa_cm_f = 2usize;
    let f0 = vec![<R as stark_rings::Ring>::ONE; n];
    let f1 = (0..n)
        .map(|i| if i % 2 == 0 { <R as stark_rings::Ring>::ONE } else { <R as stark_rings::Ring>::ZERO })
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

    let scheme_had = AjtaiCommitmentScheme::<R>::seeded(b"cfs_had_u", MASTER_SEED, 2, 3 * R::dimension());
    let scheme_mon = AjtaiCommitmentScheme::<R>::seeded(b"cfs_mon_b", MASTER_SEED, 2, rg_params.k_g);
    let scheme_g = AjtaiCommitmentScheme::<R>::seeded(b"cm_g", MASTER_SEED, 2, m * R::dimension());

    let out = prove_pi_fold_poseidon_fs::<R, PC>(
        ms.as_slice(),
        &cms,
        &witnesses,
        &public_inputs,
        Some(&scheme_had),
        Some(&scheme_mon),
        &scheme_g,
        rg_params.clone(),
        &PiFoldStreamingConfig::default(),
    )
    .expect("prove failed");

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
    let trace = attempt.trace;

    // Arithmetize the FULL trace into sparse dR1CS.
    let poseidon_cfg = <PC as cyclotomic_rings::rings::GetPoseidonParams<BF>>::get_poseidon_config();
    let (inst, assignment, _replay2, _byte_wit) =
        poseidon_sponge_dr1cs_from_trace::<BF>(&poseidon_cfg, &trace.ops).expect("build dr1cs failed");
    inst.check(&assignment).expect("poseidon dR1CS not satisfied");

    // Bridge to DPP sparse dR1CS.
    let k = inst.constraints.len();
    let nvars = inst.nvars;
    let mut a_rows = Vec::with_capacity(k);
    let mut b_rows = Vec::with_capacity(k);
    let mut c_rows = Vec::with_capacity(k);
    for row in &inst.constraints {
        a_rows.push(SparseVec::new(row.a.clone()));
        b_rows.push(SparseVec::new(row.b.clone()));
        c_rows.push(SparseVec::new(row.c.clone()));
    }
    let dr1cs = DppDr1csInstanceSparse::<BF> { n: nvars, a: a_rows, b: b_rows, c: c_rows };

    // Full-size RS FLPCP proof (π = w of length 2k).
    let ell = 2 * k;
    let flpcp = RsDr1csFlpcpSparse::<BF>::new(dr1cs, ell);
    let w = flpcp.prove(&assignment);
    assert_eq!(w.len(), 2 * k);

    // Spot-check a few random verifier queries accept on the honest proof.
    let mut rng = StdRng::seed_from_u64(999);
    for _ in 0..5 {
        let (qs, pred) = flpcp.sample_queries_and_predicate_sparse(&mut rng, &assignment).unwrap();
        assert_eq!(qs.len(), 3);
        let v = [assignment.as_slice(), w.as_slice()].concat();
        let a1 = qs[0].dot(&v);
        let a2 = qs[1].dot(&v);
        let a3 = qs[2].dot(&v);
        assert!(pred.check(&[a1, a2, a3]));
    }
}

/// Full-trace/full-k DPP pipeline (release, ignored), **without Boolean proof expansion**.
///
/// This is the production-relevant direction for large instances: Booleanization would multiply
/// proof length by ~64× and is not currently practical at the 500k-constraint scale.
///
/// Run with:
/// `cargo test -p symphony --release -- --ignored test_poseidon_trace_full_dpp_end_to_end_no_boolean_full_trace`
#[test]
#[ignore]
fn test_poseidon_trace_full_dpp_end_to_end_no_boolean_full_trace() {
    // Same real-trace harness.
    let n = 1 << 5;
    let m = 1 << 5;

    let mut a0 = SparseMatrix::<R>::identity(m);
    let mut b0 = SparseMatrix::<R>::identity(m);
    let mut c0 = SparseMatrix::<R>::identity(m);
    a0.pad_cols(n);
    b0.pad_cols(n);
    c0.pad_cols(n);
    for row in b0.coeffs.iter_mut() { row.clear(); }
    for row in c0.coeffs.iter_mut() { row.clear(); }

    let mut a1 = SparseMatrix::<R>::identity(m);
    let mut b1 = SparseMatrix::<R>::identity(m);
    let mut c1 = SparseMatrix::<R>::identity(m);
    a1.pad_cols(n);
    b1.pad_cols(n);
    c1.pad_cols(n);
    a1.coeffs.swap(0, 1);
    for row in b1.coeffs.iter_mut() { row.clear(); }
    for row in c1.coeffs.iter_mut() { row.clear(); }

    let ms: Vec<[std::sync::Arc<SparseMatrix<R>>; 3]> = vec![
        [std::sync::Arc::new(a0.clone()), std::sync::Arc::new(b0.clone()), std::sync::Arc::new(c0.clone())],
        [std::sync::Arc::new(a1.clone()), std::sync::Arc::new(b1.clone()), std::sync::Arc::new(c1.clone())],
    ];
    let ms_refs: Vec<[&SparseMatrix<R>; 3]> = vec![[&a0, &b0, &c0], [&a1, &b1, &c1]];

    let rg_params = RPParams { l_h: 4, lambda_pj: 1, k_g: 2, d_prime: (R::dimension() as u128) - 2 };

    let kappa_cm_f = 2usize;
    let f0 = vec![<R as stark_rings::Ring>::ONE; n];
    let f1 = (0..n)
        .map(|i| if i % 2 == 0 { <R as stark_rings::Ring>::ONE } else { <R as stark_rings::Ring>::ZERO })
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

    let scheme_had = AjtaiCommitmentScheme::<R>::seeded(b"cfs_had_u", MASTER_SEED, 2, 3 * R::dimension());
    let scheme_mon = AjtaiCommitmentScheme::<R>::seeded(b"cfs_mon_b", MASTER_SEED, 2, rg_params.k_g);
    let scheme_g = AjtaiCommitmentScheme::<R>::seeded(b"cm_g", MASTER_SEED, 2, m * R::dimension());

    let out = prove_pi_fold_poseidon_fs::<R, PC>(
        ms.as_slice(),
        &cms,
        &witnesses,
        &public_inputs,
        Some(&scheme_had),
        Some(&scheme_mon),
        &scheme_g,
        rg_params.clone(),
        &PiFoldStreamingConfig::default(),
    )
    .expect("prove failed");

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
    let trace = attempt.trace;

    // Arithmetize FULL trace to sparse dR1CS.
    let poseidon_cfg = <PC as cyclotomic_rings::rings::GetPoseidonParams<BF>>::get_poseidon_config();
    let (inst, assignment, _replay2, _byte_wit) =
        poseidon_sponge_dr1cs_from_trace::<BF>(&poseidon_cfg, &trace.ops).expect("build dr1cs failed");
    inst.check(&assignment).expect("poseidon dR1CS not satisfied");

    // Bridge to DPP sparse dR1CS.
    let k = inst.constraints.len();
    let nvars = inst.nvars;
    let mut a_rows = Vec::with_capacity(k);
    let mut b_rows = Vec::with_capacity(k);
    let mut c_rows = Vec::with_capacity(k);
    for row in &inst.constraints {
        a_rows.push(SparseVec::new(row.a.clone()));
        b_rows.push(SparseVec::new(row.b.clone()));
        c_rows.push(SparseVec::new(row.c.clone()));
    }
    let dr1cs = DppDr1csInstanceSparse::<BF> { n: nvars, a: a_rows, b: b_rows, c: c_rows };

    // NP-style RS FLPCP: statement x=[], witness is full assignment.
    let ell = 2 * k;
    let flpcp = RsDr1csNpFlpcpSparse::<BF>::new(dr1cs, 0, ell);
    let pi_field = flpcp.prove(&[], &assignment);

    // Embed+pack (NO Booleanization): use a large packed field to satisfy Claim 5.22 modulus condition.
    let dpp = dpp::pipeline::build_rev2_dpp_sparse::<BF, FPack521, _>(
        flpcp,
        EmbeddingParams { gamma: 2, assume_boolean_proof: false, k_prime: 0 },
        PackedDppParams { ell: 2 },
    );

    let x_large: Vec<FPack521> = vec![];
    let pi_large = pi_field
        .iter()
        .map(|wi| FPack521::from_le_bytes_mod_order(&wi.into_bigint().to_bytes_le()))
        .collect::<Vec<_>>();

    let mut rng = StdRng::seed_from_u64(424242);
    let q = dpp.sample_query(&mut rng, &x_large).expect("sample_query");
    let ok = dpp.verify_with_query(&x_large, &pi_large, &q).expect("verify_with_query");
    assert!(ok);
}

