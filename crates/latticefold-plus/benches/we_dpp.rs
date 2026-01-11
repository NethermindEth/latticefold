//! WE-gate + DPP integration bench (research).
//!
//! Current scope:
//! - Build a WE sparse dR1CS for verifying one `CmProof` (commitment transform / Π_cm)
//! - Convert it into the prototype dpp::dr1cs_flpcp pipeline and run verification
//!
//! This is not yet the full LF+ WE gate (DecompProof still TODO).

#![allow(non_snake_case)]
#![allow(non_local_definitions)]

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
use cyclotomic_rings::rings::GetPoseidonParams;

use ark_ff::{BigInteger, Field, Fp384, MontBackend, MontConfig, PrimeField};
use rand::{rngs::StdRng, RngCore, SeedableRng};

use latticefold_plus::cm::Cm;
use latticefold_plus::rgchk::{DecompParameters, Rg, RgInstance};
use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as R;
use stark_rings::PolyRing;
use stark_rings_linalg::{Matrix, SparseMatrix};

use latticefold_plus::recording_transcript::TracePoseidonTranscript;
use latticefold_plus::we_gate_arith::{build_we_dr1cs_for_cm_proof_debug, WeCmBuildDebug};
use latticefold_plus::we_statement::{
    digest32_to_field, we_statement_hash_lf_plus, WeParams, LFP_WE_GATE_DIGEST_V1,
};

use latticefold::transcript::Transcript;
use dpp::BoundedFlpcpSparse;
use dpp::packing::{
    centered_bigint_to_field, field_to_centered_bigint, sample_packing_weights, FlpcpPredicate,
    PackedDppQuerySparse,
};
use sha2::{Digest, Sha256};

// -----------------------------------------------------------------------------
// Big field for Rev2 embedding (p' large enough for packing).
// -----------------------------------------------------------------------------

#[derive(MontConfig)]
// NIST P-384 prime (as used by Symphony’s Rev2 embedding bench).
#[modulus = "39402006196394479212279040100143613805079739270465446667948293404245721771496870329047266088258938001861606973112319"]
#[generator = "2"]
pub struct Secp384r1Config;
type FBig = Fp384<MontBackend<Secp384r1Config, 6>>;

fn lift_to_big<Fs: PrimeField>(x: Fs) -> FBig {
    FBig::from_le_bytes_mod_order(&x.into_bigint().to_bytes_le())
}

fn bench_we_dpp(c: &mut Criterion) {
    // Keep defaults small-ish so local runs work; override on server by editing this file for now.
    // Toy params, but must still satisfy decomposition constraints:
    // - `RgInstance::from_f` uses `split(..., padding_size = ell)` to gadget-decompose commitment entries.
    // - If `ell` is too small, `balanced_decomposition` can panic (needs enough digits to represent a typical field element).
    // Use a conservative `ell=32` as in other benches.
    //
    // We also keep `f=0` so the *cf(f)* decomposition with `k=1` is safe.
    let k = 1usize;
    let kappa = 1usize;
    let ell = 32usize;
    let b = 2u128;
    // Ensure `n >= tau_unpadded_len` for `split`:
    // tau_unpadded_len = kappa * (k*d) * ell * d.
    let d = R::dimension();
    let tau_unpadded_len = kappa * (k * d) * ell * d;
    let n = tau_unpadded_len.next_power_of_two();
    let nvars = ark_std::log2(n) as usize;

    let dparams = DecompParameters { b, k, l: ell };
    let mut rng = ark_std::test_rng();

    // Single-instance Cm setup.
    let f = vec![R::from(<R as PolyRing>::BaseRing::ZERO); n];
    let A = Matrix::<R>::rand(&mut rng, kappa, n);
    let inst = RgInstance::from_f(f, &A, &dparams);
    let rg = Rg {
        nvars,
        instances: vec![inst],
        dparams: dparams.clone(),
    };
    let cm = Cm { rg };
    // NOTE: `setchk` expects a witness matrix `M` whose entries are unit monomials.
    // Use plain identity (all-ones on diagonal) to stay in the monomial set.
    let M: Vec<SparseMatrix<R>> = vec![SparseMatrix::identity(n)];

    // Prover-side Cm proof.
    let mut ts = latticefold_plus::transcript::PoseidonTranscript::empty::<PC>();
    // Model SP1: one public input digest (statement-defined) absorbed into the transcript *before* proving.
    // (In production this comes from SP1 public inputs.)
    type FSmall = <<R as PolyRing>::BaseRing as ark_ff::Field>::BasePrimeField;
    // Use a "random-looking" in-field digest (so we don't accidentally rely on small constants).
    let sp1_public_input_digest: FSmall = {
        let d: [u8; 32] = Sha256::digest(b"LFP_SP1_PUBLIC_INPUT_DIGEST_V1").into();
        digest32_to_field::<FSmall>(d)
    };
    ts.absorb_field_element(&sp1_public_input_digest);
    let (_com, proof) = cm.prove(&M, &mut ts);

    // Record verifier transcript ops.
    let mut rec = TracePoseidonTranscript::<R>::empty::<PC>();
    rec.absorb_field_element(&sp1_public_input_digest);
    proof.verify(&M, &mut rec).expect("cm proof verify");
    let trace = rec.trace().clone();

    // Statement params prefix (placeholder values; we only bind layout in this bench).
    let params = WeParams {
        nvars_setchk: nvars as u64,
        degree_setchk: 3,
        nvars_cm: nvars as u64,
        degree_cm: 2,
        kappa: kappa as u64,
        ring_dim_d: R::dimension() as u64,
        k: k as u64,
        l: ell as u64,
        mlen: M.len() as u64,
    };

    let poseidon_cfg = PC::get_poseidon_config();

    let mut group = c.benchmark_group("we_dpp");
    group.sample_size(10);

    group.bench_function(BenchmarkId::new("build_we_dr1cs_cm_proof", n), |bch| {
        bch.iter(|| {
            let (out, dbg) =
            build_we_dr1cs_for_cm_proof_debug::<R>(
                &poseidon_cfg,
                &trace,
                &params,
                &[sp1_public_input_digest],
                &proof,
                M.len(),
            )
                    .expect("build_we_dr1cs_for_cm_proof_debug");
            if let Err(e) = out.inst.check(&out.assignment) {
                let msg = explain_failed_constraint(&out, &dbg, &e);
                panic!("dr1cs satisfied: {e}\n{msg}");
            }
        })
    });

    group.bench_function(BenchmarkId::new("dpp_verify_cm_proof", n), |bch| {
        // Build once outside the timed loop.
        let (out, dbg) = build_we_dr1cs_for_cm_proof_debug::<R>(
            &poseidon_cfg,
            &trace,
            &params,
            &[sp1_public_input_digest],
            &proof,
            M.len(),
        )
        .expect("build_we_dr1cs_for_cm_proof_debug");
        if let Err(e) = out.inst.check(&out.assignment) {
            let msg = explain_failed_constraint(&out, &dbg, &e);
            panic!("dr1cs satisfied: {e}\n{msg}");
        }

        // Convert sparse dR1CS -> sparse dR1CS instance for the prototype RS FLPCP.
        let inst_sparse = dpp::dr1cs_flpcp::Dr1csInstanceSparse::<FSmall> {
            n: out.inst.nvars,
            a: out
                .inst
                .constraints
                .iter()
                .map(|row| dpp::SparseVec::new(row.a.clone()))
                .collect(),
            b: out
                .inst
                .constraints
                .iter()
                .map(|row| dpp::SparseVec::new(row.b.clone()))
                .collect(),
            c: out
                .inst
                .constraints
                .iter()
                .map(|row| dpp::SparseVec::new(row.c.clone()))
                .collect(),
        };
        let k_rows = inst_sparse.k();
        let ell = 2 * k_rows;
        // IMPORTANT (WE/DPP path):
        // Use the NP-style FLPCP (statement+ witness), but expose the WE statement prefix
        // as public input `x` (length = out.public_len).
        let l_public = out.public_len;
        let flpcp = dpp::dr1cs_flpcp::RsDr1csNpFlpcpSparse::<FSmall>::new(inst_sparse, l_public, ell);

        let x_small = out.assignment[..l_public].to_vec();
        let z_w_small = out.assignment[l_public..].to_vec();
        let (_pi_field_small, cw) = flpcp.prove_with_codewords(&x_small, &z_w_small);

        // Rev2 pipeline (Booleanize -> Embed -> Pack) into a large field.
        //
        // Use the same builder as Symphony to match bounds/packing behavior exactly.
        let dppv = dpp::pipeline::build_rev2_dpp_sparse_boolean_auto::<FSmall, FBig, _>(
            flpcp,
            dpp::EmbeddingParams {
                gamma: 2,
                assume_boolean_proof: true,
                k_prime: 0,
            },
        )
        .expect("build_rev2_dpp_sparse_boolean_auto");

        // Proof-agnostic arming model: derive query coins from a statement digest (no per-proof artifacts).
        // (In production, `vk_hash` and `r1cs_digest` are provided by SP1, and `gate_digest` is a fixed per-gate constant.)
        let vk_hash = [1u8; 32];
        let r1cs_digest = [2u8; 32];
        // Gate digest: production model is a precomputed constant per WE gate version.
        // (Do NOT hash over 10^8+ nonzeros at runtime.)
        let gate_digest: [u8; 32] = LFP_WE_GATE_DIGEST_V1;
        // In SP1, "public inputs" for statement arming are just the SP1 public I/O digest(s).
        let public_inputs_small = vec![sp1_public_input_digest];
        let stmt_digest = we_statement_hash_lf_plus::<R>(vk_hash, r1cs_digest, gate_digest, &public_inputs_small);

        const ARMER_SEED: [u8; 32] = *b"LFP_ARMER_SEED_V1_00000000000000";
        let lock_j: u64 = 0;
        let coin_seed: [u8; 32] = {
            let mut h = Sha256::new();
            h.update(b"LFP_LOCK_COIN_V1");
            h.update(&ARMER_SEED);
            h.update(&stmt_digest);
            h.update(&lock_j.to_le_bytes());
            h.finalize().into()
        };
        let mut rng = StdRng::from_seed(coin_seed);

        // do NOT expand packed query vectors via `sample_query()` (O(k) work).
        // Sample coins (idx, λ) + packing weights, answer the 3 RS queries in coin form by indexing
        // the cached codewords, then pack/decode via `verify_packed_answer`.
        let b = dppv.flpcp.bounds_b();
        let w = sample_packing_weights::<FBig>(&mut rng, dppv.params.ell, &b)
            .expect("sample_packing_weights");
        let pred = FlpcpPredicate::MulEqModP {
            p_small: num_bigint::BigInt::from_bytes_le(
                num_bigint::Sign::Plus,
                &FSmall::MODULUS.to_bytes_le(),
            ),
        };
        let ell_rs = 2 * k_rows;
        let idx = (rng.next_u64() as usize) % ell_rs;
        let lambda_small = FSmall::from(rng.next_u64());

        let (a_small, b_small, c_small) = if idx < k_rows {
            let a = cw.y_a[idx];
            let b0 = cw.y_b[idx];
            let wv = cw.w[idx];
            let cx_minus = cw.y_c[idx] - wv;
            let c = wv + lambda_small * cx_minus;
            (a, b0, c)
        } else {
            let j = idx - k_rows;
            let a = cw.y_a_tail[j];
            let b0 = cw.y_b_tail[j];
            let wv = cw.w[idx];
            // Tail-half: C-part unused; answer is w(α)=a*b.
            let c = wv;
            (a, b0, c)
        };

        let ans_field: [FBig; 3] = [
            lift_to_big::<FSmall>(a_small),
            lift_to_big::<FSmall>(b_small),
            lift_to_big::<FSmall>(c_small),
        ];
        let mut a_int = num_bigint::BigInt::from(0);
        for (wi, ai) in w.iter().zip(ans_field.iter()) {
            let ai_int = field_to_centered_bigint::<FBig>(ai);
            a_int += wi * ai_int;
        }
        let a = centered_bigint_to_field::<FBig>(&a_int);

        let q_meta = PackedDppQuerySparse::<FBig> { q: dpp::SparseVec::default(), w, b, pred };
        bch.iter(|| {
            let ok = dppv.verify_packed_answer(&a, &q_meta).expect("verify_packed_answer");
            assert!(ok);
        })
    });

    group.finish();
}

fn parse_failed_constraint_idx(msg: &str) -> Option<usize> {
    // expected "constraint {i} failed"
    let msg = msg.trim();
    let msg = msg.strip_prefix("constraint ")?;
    let msg = msg.strip_suffix(" failed")?;
    msg.parse::<usize>().ok()
}

fn explain_failed_constraint(
    out: &latticefold_plus::we_gate_arith::WeDr1csOutput<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
    dbg: &WeCmBuildDebug,
    err: &str,
) -> String {
    let Some(i) = parse_failed_constraint_idx(err) else {
        return "[we_dpp] could not parse failed constraint index".to_string();
    };
    let mut acc = 0usize;
    let names = [
        "poseidon",
        "params",
        "setchk_verify",
        "dcom_absorb",
        "cm_short_bytes",
        "cm_field_chals",
        "cm_verify",
    ];
    for (part_idx, &cnt) in dbg.part_constraints.iter().enumerate() {
        if i < acc + cnt {
            let name = names.get(part_idx).copied().unwrap_or("unknown");
            let mut msg = format!(
                "[we_dpp] failed constraint {i} is in PART {part_idx} ({name}), start={acc}, len={cnt}"
            );
            if part_idx == 6 && !dbg.cm_phase_marks.is_empty() {
                let local = i - acc;
                let mut phase = "unknown";
                for (j, &m) in dbg.cm_phase_marks.iter().enumerate() {
                    if local < m {
                        phase = dbg.cm_phase_names.get(j).map(|s| s.as_str()).unwrap_or("unknown");
                        break;
                    }
                }
                if phase == "unknown" {
                    if let Some(last) = dbg.cm_phase_names.last() {
                        phase = last;
                    }
                }
                msg.push_str(&format!("\n[we_dpp] cm_verify local_idx={local}, phase≈{phase}"));
            }
            return msg;
        }
        acc += cnt;
    }
    // Glue constraints
    let glue_idx = i.saturating_sub(dbg.base_constraints);
    if glue_idx < dbg.glue.len() {
        let (pa, xa, pb, xb) = dbg.glue[glue_idx];
        // Compute merged-space indices to show witness mismatch.
        let mut offsets = Vec::with_capacity(dbg.part_nvars.len());
        let mut off = 0usize;
        for &nv in &dbg.part_nvars {
            offsets.push(off);
            off += nv - 1;
        }
        let ga = if xa == 0 { 0 } else { xa + offsets[pa] };
        let gb = if xb == 0 { 0 } else { xb + offsets[pb] };
        let va = out.assignment[ga];
        let vb = out.assignment[gb];
        return format!(
            "[we_dpp] failed constraint {i} is GLUE #{glue_idx}: (part {pa}, var {xa}) == (part {pb}, var {xb})\n\
             merged idxs: {ga} vs {gb}\n\
             values: {va:?} vs {vb:?}"
        );
    } else {
        format!("[we_dpp] failed constraint {i} is after all parts+glue??")
    }
}

criterion_group!(benches, bench_we_dpp);
criterion_main!(benches);

