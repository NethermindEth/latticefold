//! LF+ one-proof harness for SP1 shrink verifier R1LF (research).
//!
//! This is the intended entrypoint for provers:
//! - load SP1’s `.r1lf` + `{path}.chunks` cache
//! - load SP1 witness (base+aux)
//! - enforce boundedness (base 2^16, k=2) on the full witness (incl. lift `t_i`)
//! - produce + verify a cryptographic streaming-sumcheck proof that each R1LF chunk satisfies
//!   the lifted R1CS relation: (A·w) * (B·w) == (C·w)
//!
//! Usage:
//!   SP1_R1LF=/path/to/shrink_verifier.r1lf \
//!   SP1_WITNESS=/path/to/shrink_verifier.witness.u64le \
//!     cargo run -p latticefold-plus --example lf_plus_sp1_oneproof --features we_gate --release

#![cfg(feature = "we_gate")]

use ark_ff::{BigInteger, PrimeField};
use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
use latticefold::transcript::Transcript;
use latticefold::utils::sumcheck::utils::eq_eval;
use latticefold::utils::sumcheck::MLSumcheck;
use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as R;
use stark_rings::PolyRing;
use stark_rings::Ring;
use std::sync::Arc;
use std::time::Instant;

type F = <R as PolyRing>::BaseRing;

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

fn modulus_u128<P: PrimeField>() -> u128 {
    let le = P::MODULUS.to_bytes_le();
    let mut out = 0u128;
    for (i, b) in le.iter().enumerate().take(16) {
        out |= (*b as u128) << (8 * i);
    }
    out
}

#[inline]
fn centered_i128_from_canonical_u64(x: u64, q: u128) -> i128 {
    let xu = x as u128;
    let half = q >> 1;
    if xu > half {
        (xu as i128) - (q as i128)
    } else {
        xu as i128
    }
}

#[inline]
fn check_base2_16_k2(centered: i128) -> bool {
    const B: i128 = 1i128 << 16;
    const HALF: i128 = 1i128 << 15;
    let mut x = centered;
    for _ in 0..2 {
        let mut d = x.rem_euclid(B);
        if d >= HALF {
            d -= B;
        }
        if d < i16::MIN as i128 || d > i16::MAX as i128 {
            return false;
        }
        x = (x - d) / B;
    }
    x == 0
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
    println!("LF+ SP1 One-Proof (R1LF prove+verify, streaming sumcheck per chunk)");
    println!("=========================================================");
    println!("  CHUNK_SIZE={chunk_size} PAD_COLS={pad_cols_to_multiple_of}");

    let t0 = Instant::now();
    let cache =
        latticefold_plus::sp1_r1lf::open_sp1_r1lf_chunk_cache::<R>(&r1lf_path, chunk_size, pad_cols_to_multiple_of)
            .expect("open_sp1_r1lf_chunk_cache");
    println!("  cache open: {:?}", t0.elapsed());
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

    let (w_u64, base_len, aux_len) =
        latticefold_plus::sp1_witness_io::load_sp1_witness_any(&witness_path, cache.stats.num_vars)
            .expect("load witness");
    println!("  loaded witness: base={} aux={} full={}", base_len, aux_len, w_u64.len());
    assert!(!w_u64.is_empty() && w_u64[0] == 1, "witness must have w[0]=1");

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

    // Boundedness (base 2^16, k=2) over the full witness (incl. aux t_i).
    let t_b = Instant::now();
    let q = modulus_u128::<F>();
    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        let bad = w_host.par_iter().enumerate().find_any(|(_i, fx)| {
            let bi = fx.into_bigint();
            let limbs = bi.as_ref();
            let x0 = limbs.get(0).copied().unwrap_or(0);
            let centered = centered_i128_from_canonical_u64(x0, q);
            !check_base2_16_k2(centered)
        });
        if let Some((i, fx)) = bad {
            panic!("boundedness failed at idx={i}: fx={fx:?}");
        }
    }
    #[cfg(not(feature = "parallel"))]
    {
        for (i, fx) in w_host.iter().enumerate() {
            let bi = fx.into_bigint();
            let limbs = bi.as_ref();
            let x0 = limbs.get(0).copied().unwrap_or(0);
            let centered = centered_i128_from_canonical_u64(x0, q);
            if !check_base2_16_k2(centered) {
                panic!("boundedness failed at idx={i}: fx={fx:?}");
            }
        }
    }
    println!("  boundedness (base2^16,k=2) OK: {:?}", t_b.elapsed());

    // Prove+verify per chunk: sumcheck over rows proving eq(row,r0)*(Aw*Bw - Cw) sums to 0.
    //
    // This is memory-friendly (one chunk in RAM at a time), and cryptographically binds the claim
    // via the Fiat–Shamir transcript.
    let t_pv = Instant::now();
    for chunk_idx in 0..cache.num_chunks {
        let t_chunk = Instant::now();
        let [a, b, c] = cache.read_chunk(chunk_idx).expect("read_chunk");
        let nrows = a.nrows;
        assert!(nrows.is_power_of_two(), "chunk nrows must be power-of-two");
        let nvars = ark_std::log2(nrows) as usize;

        let a = Arc::new(a);
        let b = Arc::new(b);
        let c = Arc::new(c);

        // Prover transcript.
        let mut ts = latticefold_plus::transcript::PoseidonTranscript::<R>::empty::<PC>();
        ts.absorb_field_element(&F::from(0x4c46502b_53503152u128)); // "LFP+SP1R"
        ts.absorb_field_element(&F::from(chunk_idx as u128));
        // Bind to the statement digest (cheap field reduction).
        ts.absorb_field_element(&F::from_le_bytes_mod_order(&cache.stats.digest));

        let r0 = ts.get_challenges(nvars);
        let one = F::from(1u64);
        let one_minus_r0 = r0.iter().copied().map(|x| one - x).collect::<Vec<_>>();

        let mles = vec![
            latticefold_plus::streaming_sumcheck::StreamingMleEnum::<R>::EqBase {
                scale: one,
                r: r0.clone(),
                one_minus_r: one_minus_r0,
            },
            latticefold_plus::streaming_sumcheck::StreamingMleEnum::<R>::SparseMatVecConstCoeffBase {
                matrix: a.clone(),
                witness0: w_host.clone(),
                num_vars: nvars,
            },
            latticefold_plus::streaming_sumcheck::StreamingMleEnum::<R>::SparseMatVecConstCoeffBase {
                matrix: b.clone(),
                witness0: w_host.clone(),
                num_vars: nvars,
            },
            latticefold_plus::streaming_sumcheck::StreamingMleEnum::<R>::SparseMatVecConstCoeffBase {
                matrix: c.clone(),
                witness0: w_host.clone(),
                num_vars: nvars,
            },
        ];

        let comb_fn = |vals: &[R]| -> R { vals[0] * (vals[1] * vals[2] - vals[3]) };

        let (proof, _rand, final_vals) =
            latticefold_plus::streaming_sumcheck::StreamingSumcheck::prove_as_subprotocol(
                &mut ts,
                mles,
                nvars,
                3,
                comb_fn,
            );

        let va = final_vals[1];
        let vb = final_vals[2];
        let vc = final_vals[3];
        ts.absorb_slice(&[va, vb, vc]);

        // Verifier transcript (must match).
        let mut tv = latticefold_plus::transcript::PoseidonTranscript::<R>::empty::<PC>();
        tv.absorb_field_element(&F::from(0x4c46502b_53503152u128));
        tv.absorb_field_element(&F::from(chunk_idx as u128));
        tv.absorb_field_element(&F::from_le_bytes_mod_order(&cache.stats.digest));
        let r0_v: Vec<R> = tv
            .get_challenges(nvars)
            .into_iter()
            .map(R::from)
            .collect();

        let sub = MLSumcheck::verify_as_subprotocol(&mut tv, nvars, 3, R::ZERO, &proof)
            .expect("sumcheck verify");
        let ro: Vec<R> = sub.point.into_iter().map(R::from).collect();
        let s = sub.expected_evaluation;

        tv.absorb_slice(&[va, vb, vc]);

        let eq = eq_eval(&r0_v, &ro).expect("eq_eval");
        assert_eq!(eq * (va * vb - vc), s, "chunk {chunk_idx} failed");

        if chunk_idx == 0 || chunk_idx + 1 == cache.num_chunks {
            println!(
                "  chunk {}/{} OK: {:?} (nrows={}, nvars={})",
                chunk_idx + 1,
                cache.num_chunks,
                t_chunk.elapsed(),
                nrows,
                nvars
            );
        }
    }
    println!("  total prove+verify chunks: {:?}", t_pv.elapsed());
}

