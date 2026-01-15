//! LF+ one-proof harness for SP1 shrink verifier R1LF (production path).
//!
//! This produces a real `PlusProof<R, ComR1CSProof<R>>` so the **existing** LF+ WE/DPP gate
//! (`build_we_dr1cs_for_plus_proof`) can arithmetize and verify it unchanged.
//!
//! Implementation strategy (Salsa/Symphony-style):
//! - load `.r1lf` chunk cache and materialize A/B/C into in-memory sparse matrices (const-coeff)
//! - load SP1 witness and embed into `R` as constant-coeff ring elements
//! - run `PlusProver` to produce a `PlusProof`
//! - record the verifier transcript trace and sanity-check that the WE gate dR1CS is satisfied
//!
//! Usage:
//!   SP1_R1LF=/path/to/shrink_verifier.r1lf \
//!   SP1_WITNESS=/path/to/shrink_verifier.witness.u64le \
//!     cargo run -p latticefold-plus --example lf_plus_sp1_oneproof --features we_gate --release

#![cfg(feature = "we_gate")]
#![allow(non_local_definitions)]

use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
use cyclotomic_rings::rings::GetPoseidonParams;
use ark_ff::{BigInteger, Fp384, MontBackend, MontConfig, PrimeField};
use latticefold::commitment::AjtaiCommitmentScheme;
use latticefold::transcript::Transcript;
use latticefold_plus::lin::LinearizedVerify;
use latticefold_plus::utils::maybe_print_rss;
use latticefold_plus::we_statement::{we_statement_hash_lf_plus, LFP_WE_GATE_DIGEST_V1};
use cyclotomic_rings::rings::FrogRingPoly as R;
use stark_rings::PolyRing;
use stark_rings_linalg::SparseMatrix;
use std::sync::Arc;
use std::time::Instant;
use sha2::{Digest, Sha256};
use rand::{rngs::StdRng, RngCore, SeedableRng};

use dpp::BoundedFlpcpSparse;
use dpp::dr1cs_flpcp::{Dr1csInstanceSparse as DppInst, RsDr1csNpFlpcpSparse};
use dpp::packing::{
    centered_bigint_to_field, field_to_centered_bigint, sample_packing_weights, FlpcpPredicate,
    PackedDppQuerySparse,
};
use dpp::pipeline::build_rev2_dpp_sparse_boolean_auto;
use dpp::sparse::SparseVec;

// Big field for Rev2 embedding (same as `test_large_trace` / `benches/we_dpp.rs`).
#[derive(MontConfig)]
#[modulus = "39402006196394479212279040100143613805079739270465446667948293404245721771496870329047266088258938001861606973112319"]
#[generator = "2"]
pub struct Secp384r1Config;
type FBig = Fp384<MontBackend<Secp384r1Config, 6>>;

type F = <R as PolyRing>::BaseRing;

fn lift_to_big<Fs: PrimeField>(x: Fs) -> FBig {
    FBig::from_le_bytes_mod_order(&x.into_bigint().to_bytes_le())
}

#[inline]
fn babybear_u64_to_centered_host(x: u64, p_bb: u64) -> F {
    debug_assert!(p_bb > 1);
    let half = p_bb / 2;
    if x > half {
        let neg = p_bb - x;
        -F::from(neg)
    } else {
        F::from(x)
    }
}

fn main() {
    let r1lf_path = std::env::var("SP1_R1LF").expect("Set SP1_R1LF=/path/to/shrink.r1lf");
    let witness_path =
        std::env::var("SP1_WITNESS").expect("Set SP1_WITNESS=/path/to/shrink_verifier.witness.u64le");
    let chunk_size: usize = std::env::var("CHUNK_SIZE")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(1 << 20);
    let pad_cols_to_multiple_of: usize = std::env::var("PAD_COLS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(256);

    println!("=========================================================");
    println!("LF+ SP1 One-Proof (R1LF -> full PlusProof -> WE gate check)");
    println!("=========================================================");
    println!("  CHUNK_SIZE={chunk_size} PAD_COLS={pad_cols_to_multiple_of}");

    let t0 = Instant::now();
    let cache =
        latticefold_plus::sp1_r1lf::open_sp1_r1lf_chunk_cache::<R>(&r1lf_path, chunk_size, pad_cols_to_multiple_of)
            .expect("open_sp1_r1lf_chunk_cache");
    println!("  cache open: {:?}", t0.elapsed());
    maybe_print_rss("after cache open");
    println!("  chunks={} ncols={}", cache.num_chunks, cache.ncols);
    println!(
        "  stats: num_vars={} num_constraints={} num_public={} p_bb={} total_nonzeros={}",
        cache.stats.num_vars,
        cache.stats.num_constraints,
        cache.stats.num_public,
        cache.stats.p_bb,
        cache.stats.total_nonzeros
    );
    println!("  digest={:02x?}...", &cache.stats.digest[..8]);

    // Materialize full (A,B,C) as SparseMatrix<F> (constant-coeff) by concatenating chunk rows.
    let t_mats = Instant::now();
    let total_rows = cache.num_chunks * chunk_size;
    let mut a_rows: Vec<Vec<(F, usize)>> = Vec::with_capacity(total_rows);
    let mut b_rows: Vec<Vec<(F, usize)>> = Vec::with_capacity(total_rows);
    let mut c_rows: Vec<Vec<(F, usize)>> = Vec::with_capacity(total_rows);
    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        type Rows = Vec<Vec<(F, usize)>>;

        // NOTE: `into_par_iter()` over a range is an IndexedParallelIterator, so `collect::<Vec<_>>()`
        // preserves order by `chunk_idx`. This keeps row ordering deterministic.
        let chunks: Vec<(Rows, Rows, Rows)> = (0..cache.num_chunks)
            .into_par_iter()
            .map(|chunk_idx| {
                let [a, b, c] = cache.read_chunk(chunk_idx).expect("read_chunk");
                debug_assert_eq!(a.nrows, chunk_size);

                let conv = |m: stark_rings_linalg::SparseMatrix<F>| m.coeffs;

                (conv(a), conv(b), conv(c))
            })
            .collect();

        for (ar, br, cr) in chunks {
            a_rows.extend(ar);
            b_rows.extend(br);
            c_rows.extend(cr);
        }
    }

    #[cfg(not(feature = "parallel"))]
    {
        for chunk_idx in 0..cache.num_chunks {
            let [a, b, c] = cache.read_chunk(chunk_idx).expect("read_chunk");
            debug_assert_eq!(a.nrows, chunk_size);
            for row in a.coeffs {
                a_rows.push(row);
            }
            for row in b.coeffs {
                b_rows.push(row);
            }
            for row in c.coeffs {
                c_rows.push(row);
            }
        }
    }
    let m_a = SparseMatrix::<F> { nrows: total_rows, ncols: cache.ncols, coeffs: a_rows };
    let m_b = SparseMatrix::<F> { nrows: total_rows, ncols: cache.ncols, coeffs: b_rows };
    let m_c = SparseMatrix::<F> { nrows: total_rows, ncols: cache.ncols, coeffs: c_rows };
    println!(
        "  build full mats: {:?} (nrows={} ncols={})",
        t_mats.elapsed(),
        total_rows,
        cache.ncols
    );
    maybe_print_rss("after build full mats (A,B,C)");

    let (w_u64, base_len, aux_len) =
        latticefold_plus::sp1_witness_io::load_sp1_witness_any(&witness_path, cache.stats.num_vars)
            .expect("load witness");
    println!("  loaded witness: base={} aux={} full={}", base_len, aux_len, w_u64.len());
    assert!(!w_u64.is_empty() && w_u64[0] == 1, "witness must have w[0]=1");
    maybe_print_rss("after load witness u64");

    // Map u64 witness -> Frog base field scalars once, matching SP1 lift semantics:
    // - **all vars (including aux)**: centered embedding mod p_bb
    let p_bb = cache.stats.p_bb;
    let t_w = Instant::now();
    let w_host: Arc<Vec<F>> = Arc::new(
        w_u64
            .iter()
            .copied()
            .enumerate()
            .map(|(i, x)| {
                if x >= p_bb {
                    panic!("witness word out of [0,p_bb) range at idx={i}: x={x} p_bb={p_bb}");
                }
                babybear_u64_to_centered_host(x, p_bb)
            })
            .collect(),
    );
    println!("  map witness u64->F: {:?}", t_w.elapsed());
    maybe_print_rss("after map witness u64->F");

    // Keep witness as **base scalars** (const-coeff embedding).
    //
    // IMPORTANT: do NOT pad to `ncols`. We treat missing columns as implicit zeros throughout
    // the prover, while still committing / sampling challenges over the full `ncols` domain.
    let t_f0 = Instant::now();
    let mut f0 = (*w_host).clone();
    f0.truncate(cache.stats.num_vars);
    let f0: Arc<Vec<F>> = Arc::new(f0);
    println!("  build f0 (base scalars, padded): {:?}", t_f0.elapsed());
    maybe_print_rss("after build f0 padded");

    // Build `ComR1CS` instance and run the full LF+ prover to produce a `PlusProof`.
    let t_setup = Instant::now();
    let r1cs = latticefold::arith::r1cs::R1CS::<F> { l: 0, A: m_a, B: m_b, C: m_c };
    maybe_print_rss("after build r1cs struct");

    // Deterministic Ajtai commitment scheme (system parameter). Keep kappa=1 for now.
    let kappa: usize = 1;
    const AJTAI_SEED: [u8; 32] = *b"LFP_SP1_AJTAI_SEED_V1_0000000000";
    let ajtai = AjtaiCommitmentScheme::<R>::seeded(b"lf_plus_ajtai", AJTAI_SEED, kappa, cache.ncols);
    maybe_print_rss("after init Ajtai scheme");

    let cr1cs = latticefold_plus::r1cs::ComR1CSBase::<R>::from_f0_seeded_base(r1cs, f0, 0, &ajtai);
    maybe_print_rss("after ComR1CS::from_f0_seeded");
    let m0 = cr1cs.x.matrices_arc_base();
    maybe_print_rss("after matrices_arc");

    // LF+ parameters: boundedness base b=2^16,k=2, and a conservative decomp base B for Π_decomp.
    let we_params =
        latticefold_plus::sp1_r1lf::sp1_default_we_params_for_r1lf_cache::<R>(&cache, kappa as u64, m0.len() as u64)
            .expect("sp1_default_we_params_for_r1lf_cache");
    let dparams = latticefold_plus::rgchk::DecompParameters {
        b: (we_params.decomp_b as u128),
        k: (we_params.k as usize),
        l: (we_params.l as usize),
    };
    let lin_params = latticefold_plus::lin::LinParameters { kappa, decomp: dparams };
    // Π_decomp splits each base-field element into exactly 2 digits in base `B` (see `decomp.rs`),
    // so we must pick `B` large enough that **every** value we decompose fits in 2 digits.
    //
    // The balanced decomposition uses signed representatives in [-q/2, q/2], so it's sufficient to
    // choose B > sqrt(q). We round up to a power-of-two (even) for fast divisions.
    fn isqrt_u128(x: u128) -> u128 {
        // Integer floor sqrt via binary search (no floats).
        let mut lo: u128 = 0;
        let mut hi: u128 = 1u128 << 64; // sqrt(u128::MAX) < 2^64
        while lo + 1 < hi {
            let mid = (lo + hi) >> 1;
            if mid.saturating_mul(mid) <= x {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        lo
    }
    fn next_pow2_u128(x: u128) -> u128 {
        if x <= 1 {
            return 1;
        }
        let p = 128 - (x - 1).leading_zeros() as u32;
        1u128 << p
    }
    let q_u128: u128 = <F as ark_ff::PrimeField>::MODULUS.0[0] as u128;
    let b_min = isqrt_u128(q_u128) + 1;
    let b_decomp: u128 = next_pow2_u128(b_min);
    let pparams = latticefold_plus::plus::PlusParameters { lin: lin_params, B: b_decomp };

    // Public statement binding: use the SP1 r1lf digest bits as public inputs (boolean field elems).
    type BFSmall = <<R as PolyRing>::BaseRing as ark_ff::Field>::BasePrimeField;
    let public_inputs: Vec<BFSmall> = {
        let d: [u8; 32] = cache.stats.digest;
        latticefold_plus::we_statement::digest32_to_bits_field::<BFSmall>(d)
    };

    // Proof-agnostic arming statement digest (binds vk, r1cs, gate version, **params**, and public inputs).
    // This is what an honest armer/decapper should use to derive lock coins.
    let vk_hash = [0u8; 32]; // TODO: wire SP1 program/vk digest here (statement-defined)
    let r1cs_digest = cache.stats.digest; // SP1 R1LF instance digest (statement-defined)
    let stmt_digest = we_statement_hash_lf_plus::<R>(vk_hash, r1cs_digest, LFP_WE_GATE_DIGEST_V1, &we_params, &public_inputs);
    println!("  stmt_digest={:02x?}...", &stmt_digest[..8]);

    // Demonstrate how an honest armer derives lock/query coins from the statement digest.
    // (This is *outside* the LF+ transcript, so it does not affect prover/verifier behavior.)
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
    println!("  lock_coin_seed={:02x?}... (j={lock_j})", &coin_seed[..8]);

    let mut prover = latticefold_plus::plus::PlusProverSparseBase::init_seeded_base(
        ajtai.clone(),
        m0.clone(),
        1,
        pparams.clone(),
        latticefold_plus::transcript::PoseidonTranscript::empty::<PC>(),
    );
    for b in &public_inputs {
        prover.transcript.absorb_field_element(b);
    }
    println!("  setup full LF+: {:?}", t_setup.elapsed());
    maybe_print_rss("after setup full LF+");

    let t_prove = Instant::now();
    let proof = prover.prove_sparse_base(std::slice::from_ref(&cr1cs));
    println!("  PlusProverSparseBase::prove_sparse_base: {:?}", t_prove.elapsed());
    maybe_print_rss("after prove_sparse_base");

    // Record verifier trace and ensure the existing WE gate arithmetization is satisfied.
    let poseidon_cfg = PC::get_poseidon_config();
    let mut rec = latticefold_plus::recording_transcript::TracePoseidonTranscript::<R>::empty::<PC>();
    for b in &public_inputs {
        rec.absorb_field_element(b);
    }
    let t_verify_record = Instant::now();
    for lp in &proof.lproof {
        lp.verify(&mut rec);
    }
    proof.cmproof
        .verify_with_mlen(m0.len(), &mut rec)
        .expect("cm proof verify");
    println!("  PlusVerifier::verify(record trace): {:?}", t_verify_record.elapsed());
    maybe_print_rss("after verify(record)");
    let trace = rec.trace().clone();

    let t_we = Instant::now();
    let out = latticefold_plus::we_gate_arith::build_we_dr1cs_for_plus_proof::<R>(
        &poseidon_cfg,
        &trace,
        &we_params,
        &public_inputs,
        &proof,
        m0.len(),
        b_decomp,
    )
    .expect("build_we_dr1cs_for_plus_proof");
    println!("  WE gate build_dr1cs: {:?}", t_we.elapsed());
    maybe_print_rss("after WE build_dr1cs");

    let t_sat = Instant::now();
    out.inst.check(&out.assignment).expect("we gate dr1cs satisfied");
    println!("  WE gate dr1cs sat check: {:?}", t_sat.elapsed());
    maybe_print_rss("after WE sat check");

    // -------------------------------------------------------------------------
    // Armer-side randomness: match `test_large_trace` coin sampling.
    // -------------------------------------------------------------------------
    let t_arm = Instant::now();
    let (inst, assignment, public_len) = (out.inst, out.assignment, out.public_len);
    let n = inst.nvars;
    let mut a = Vec::with_capacity(inst.constraints.len());
    let mut b = Vec::with_capacity(inst.constraints.len());
    let mut c = Vec::with_capacity(inst.constraints.len());
    for mut row in inst.constraints {
        a.push(SparseVec::new(std::mem::take(&mut row.a)));
        b.push(SparseVec::new(std::mem::take(&mut row.b)));
        c.push(SparseVec::new(std::mem::take(&mut row.c)));
    }
    let inst_sparse = DppInst::<BFSmall> { n, a, b, c };
    let k_rows = inst_sparse.k();
    let ell_rs = 2 * k_rows;
    let flpcp = RsDr1csNpFlpcpSparse::<BFSmall>::new(inst_sparse, public_len, ell_rs);

    // Produce codewords once so we can answer the RS-FLPCP query in coin form
    // (matching `test_large_trace`).
    let x_small = assignment[..public_len].to_vec();
    let z_w_small = assignment[public_len..].to_vec();
    let (_pi_field, cw) = flpcp.prove_with_codewords(&x_small, &z_w_small);

    let dppv = build_rev2_dpp_sparse_boolean_auto::<BFSmall, FBig, _>(
        flpcp,
        dpp::EmbeddingParams {
            gamma: 2,
            assume_boolean_proof: true,
            k_prime: 0,
        },
    )
    .expect("build dpp");

    // Sample coins exactly like `test_large_trace`:
    // - packing weights first
    // - then idx and lambda
    let mut rng = StdRng::from_seed(coin_seed);
    let bnds = dppv.flpcp.bounds_b();
    let w = sample_packing_weights::<FBig>(&mut rng, dppv.params.ell, &bnds)
        .expect("sample_packing_weights");
    let idx = (rng.next_u64() as usize) % ell_rs;
    let lambda_small = BFSmall::from(rng.next_u64());

    // Coin-form RS-FLPCP answers (exactly like `test_large_trace`).
    let (a_small, b_small, c_small) = if idx < k_rows {
        let a0 = cw.y_a[idx];
        let b0 = cw.y_b[idx];
        let wv = cw.w[idx];
        let cx_minus = cw.y_c[idx] - wv;
        let c0 = wv + lambda_small * cx_minus;
        (a0, b0, c0)
    } else {
        let j = idx - k_rows;
        let a0 = cw.y_a_tail[j];
        let b0 = cw.y_b_tail[j];
        let wv = cw.w[idx];
        // Tail-half: the C-part is unused in q3; answer is just w(α)=a*b.
        (a0, b0, wv)
    };

    // Pack into one big-field element and verify the accepting predicate.
    let ans_field: [FBig; 3] = [lift_to_big(a_small), lift_to_big(b_small), lift_to_big(c_small)];
    let mut a_int = num_bigint::BigInt::from(0);
    for (wi, ai) in w.iter().zip(ans_field.iter()) {
        let ai_int = field_to_centered_bigint::<FBig>(ai);
        a_int += wi * ai_int;
    }
    let a_big = centered_bigint_to_field::<FBig>(&a_int);
    let pred = FlpcpPredicate::MulEqModP {
        p_small: num_bigint::BigInt::from_bytes_le(
            num_bigint::Sign::Plus,
            &BFSmall::MODULUS.to_bytes_le(),
        ),
    };
    let q_meta = PackedDppQuerySparse::<FBig> { q: dpp::sparse::SparseVec::default(), w, b: bnds, pred };
    let ok = dppv.verify_packed_answer(&a_big, &q_meta).expect("verify_packed_answer");

    println!(
        "  armer/decap: idx={idx} ell_rs={ell_rs} lambda_small={:?} ok={ok} (arm_time={:?})",
        lambda_small,
        t_arm.elapsed()
    );

    // Non-transcript (local) consistency check for Π_decomp.
    // This does not affect the recorded verifier trace; WE gate enforces Π_decomp separately.
    let t_decomp_local = Instant::now();
    proof
        .dproof
        .verify(&proof.linb2x.cm_g, &proof.linb2x.vo, b_decomp);
    println!("  Π_decomp local verify (non-trace): {:?}", t_decomp_local.elapsed());

    println!("  OK: WE gate DR1CS satisfied (existing gate, unchanged)");
}

