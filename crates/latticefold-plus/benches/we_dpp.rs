//! WE-gate + DPP integration bench (research).
//!
//! Current scope:
//! - Build a WE sparse dR1CS for verifying one **full** LF+ `PlusProof`
//!   (Π_lin + Π_cm + Π_decomp, i.e. the full LF+ verifier)
//! - Convert it into the prototype dpp::dr1cs_flpcp pipeline and run verification
//!
//! This is intended to be the “apples-to-apples” WE/DPP benchmark surface for LF+:
//! arithmetize the verifier trace, then run the Rev2 (Booleanize → Embed → Pack) DPP pipeline.

#![allow(non_snake_case)]
#![allow(non_local_definitions)]

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
use cyclotomic_rings::rings::GetPoseidonParams;

use ark_ff::{BigInteger, Field, Fp384, MontBackend, MontConfig, PrimeField};
use rand::{rngs::StdRng, RngCore, SeedableRng};

use latticefold_plus::lin::LinearizedVerify;
use latticefold_plus::lin::LinParameters;
use latticefold_plus::plus::{PlusParameters, PlusProver};
use latticefold_plus::r1cs::{r1cs_decomposed_square, ComR1CS};
use latticefold_plus::rgchk::DecompParameters;
use latticefold_plus::utils::estimate_bound;
use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as R;
use stark_rings::PolyRing;
use stark_rings_linalg::{Matrix, SparseMatrix};
use std::sync::Arc;

use latticefold_plus::recording_transcript::TracePoseidonTranscript;
use latticefold_plus::we_gate_arith::build_we_dr1cs_for_plus_proof;
use latticefold_plus::we_statement::{
    digest32_to_bits_field, we_statement_hash_lf_plus, WeParams, LFP_WE_GATE_DIGEST_V1,
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
    // Toy-ish params, but must still satisfy decomposition constraints. Use a conservative `ell=32`.
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

    // Model SP1: one public input digest (statement-defined) absorbed into the transcript *before* proving.
    // (In production this comes from SP1 public inputs.)
    type FSmall = <<R as PolyRing>::BaseRing as ark_ff::Field>::BasePrimeField;
    // Use a "random-looking" in-field digest (so we don't accidentally rely on small constants).
    let sp1_public_input_digest_bits: Vec<FSmall> = {
        let d: [u8; 32] = Sha256::digest(b"LFP_SP1_PUBLIC_INPUT_DIGEST_V1").into();
        digest32_to_bits_field::<FSmall>(d)
    };

    // Ajtai matrix + monomial witness matrices (identity keeps `setchk` happy).
    //
    // IMPORTANT: `build_we_dr1cs_for_plus_proof` enforces a prefix binding constraint:
    // `cm_f[j] == public_inputs[j]` for `j < min(kappa, 8)`.
    //
    // For this bench we keep `kappa=1`, so we must ensure `cm_f[0]` is a constant-coeff ring
    // element equal to `public_inputs[0]` (a digest bit in {0,1}). We do that by making the
    // first Ajtai row expose `f[0]` directly: `A[0,0]=1`, all other entries 0.
    let A: Matrix<R> = {
        let mut rows = vec![vec![<R as stark_rings::Ring>::ZERO; n]; kappa];
        rows[0][0] = <R as stark_rings::Ring>::ONE;
        rows.into()
    };
    let M: Vec<Arc<SparseMatrix<R>>> = vec![Arc::new(SparseMatrix::identity(n))];

    // Minimal Π_lin component so the transcript prefix is exercised.
    //
    // We use the same “decomposed square” harness as `test_large_trace` in `we_gate_arith.rs`.
    // This is a bench harness, not a tight bound.
    let sop = R::dimension() * 128;
    let B_bound = estimate_bound(sop, 1, R::dimension(), k) + 1;
    let m = n / k;
    let mut z: Vec<R> = (0..m)
        .map(|_| R::from((rng.next_u64() & 1) as u128))
        .collect();
    // Satisfy the enforced prefix binding for this toy bench: ensure `f[0] == public_inputs[0]`.
    // With `k=1`, the gadget decomposition leaves `f[0] == z[0]`.
    let bit0 = if sp1_public_input_digest_bits[0] == <FSmall as Field>::ZERO {
        0u128
    } else {
        1u128
    };
    z[0] = R::from(bit0);
    let r1cs0 = r1cs_decomposed_square(
        latticefold::arith::r1cs::R1CS::<R> {
            l: 1,
            A: SparseMatrix::identity(m),
            B: SparseMatrix::identity(m),
            C: SparseMatrix::identity(m),
        },
        n,
        B_bound,
        k,
    );
    let cr1cs = ComR1CS::new(r1cs0, z, 1, B_bound, k, &A);
    let lin_params = LinParameters { kappa, decomp: dparams.clone() };
    let pparams = PlusParameters { lin: lin_params, B: B_bound };

    // Prover-side full Plus proof.
    let transcript = latticefold_plus::transcript::PoseidonTranscript::empty::<PC>();
    let mut prover = PlusProver::init(A.clone(), M.clone(), 1, pparams.clone(), transcript);
    for b in &sp1_public_input_digest_bits {
        prover.transcript.absorb_field_element(b);
    }
    let proof = prover.prove(&[cr1cs]);

    // Record verifier transcript ops.
    let mut rec = TracePoseidonTranscript::<R>::empty::<PC>();
    for b in &sp1_public_input_digest_bits {
        rec.absorb_field_element(b);
    }
    // Mirror `PlusVerifier::verify` to record the full verifier trace.
    for lp in &proof.lproof {
        lp.verify(&mut rec);
    }
    proof.cmproof.verify(&M, &mut rec).expect("cm proof verify");
    proof
        .dproof
        .verify(&proof.linb2x.cm_g, &proof.linb2x.vo, B_bound);
    let trace = rec.trace().clone();

    // Statement params prefix (placeholder values; we only bind layout in this bench).
    let params = WeParams {
        nvars_setchk: nvars as u64,
        degree_setchk: 3,
        nvars_cm: nvars as u64,
        degree_cm: 2,
        kappa: kappa as u64,
        ring_dim_d: R::dimension() as u64,
        // Must match the decomposition parameters actually used by the LF+ proof in this bench.
        decomp_b: b as u64,
        k: k as u64,
        l: ell as u64,
        mlen: M.len() as u64,
    };

    let poseidon_cfg = PC::get_poseidon_config();

    let mut group = c.benchmark_group("we_dpp");
    group.sample_size(10);

    group.bench_function(BenchmarkId::new("build_we_dr1cs_plus_proof", n), |bch| {
        bch.iter(|| {
            let out = build_we_dr1cs_for_plus_proof::<R>(
                &poseidon_cfg,
                &trace,
                &params,
                &sp1_public_input_digest_bits,
                &proof,
                M.len(),
                B_bound,
            )
            .expect("build_we_dr1cs_for_plus_proof");
            out.inst.check(&out.assignment).expect("dr1cs satisfied");
        })
    });

    group.bench_function(BenchmarkId::new("dpp_verify_plus_proof", n), |bch| {
        // Build once outside the timed loop.
        let out = build_we_dr1cs_for_plus_proof::<R>(
            &poseidon_cfg,
            &trace,
            &params,
            &sp1_public_input_digest_bits,
            &proof,
            M.len(),
            B_bound,
        )
        .expect("build_we_dr1cs_for_plus_proof");
        out.inst.check(&out.assignment).expect("dr1cs satisfied");

        // Convert sparse dR1CS -> sparse dR1CS instance for the prototype RS FLPCP.
        //
        // IMPORTANT: avoid cloning multi-million sparse rows. Consume `out.inst.constraints` and
        // move `(a,b,c)` out of each row.
        let (inst, assignment, public_len) = (out.inst, out.assignment, out.public_len);
        let n = inst.nvars;
        let mut a = Vec::with_capacity(inst.constraints.len());
        let mut b = Vec::with_capacity(inst.constraints.len());
        let mut c = Vec::with_capacity(inst.constraints.len());
        for mut row in inst.constraints {
            a.push(dpp::SparseVec::new(std::mem::take(&mut row.a)));
            b.push(dpp::SparseVec::new(std::mem::take(&mut row.b)));
            c.push(dpp::SparseVec::new(std::mem::take(&mut row.c)));
        }
        let inst_sparse = dpp::dr1cs_flpcp::Dr1csInstanceSparse::<FSmall> { n, a, b, c };
        let k_rows = inst_sparse.k();
        let ell = 2 * k_rows;
        // IMPORTANT (WE/DPP path):
        // Use the NP-style FLPCP (statement+ witness), but expose the WE statement prefix
        // as public input `x` (length = out.public_len).
        let l_public = public_len;
        let flpcp = dpp::dr1cs_flpcp::RsDr1csNpFlpcpSparse::<FSmall>::new(inst_sparse, l_public, ell);

        let x_small = assignment[..l_public].to_vec();
        let z_w_small = assignment[l_public..].to_vec();
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
        let stmt_digest =
            we_statement_hash_lf_plus::<R>(vk_hash, r1cs_digest, gate_digest, &params, &sp1_public_input_digest_bits);

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

criterion_group!(benches, bench_we_dpp);
criterion_main!(benches);

