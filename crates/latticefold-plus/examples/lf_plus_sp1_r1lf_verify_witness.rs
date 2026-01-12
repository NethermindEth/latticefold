//! Verify that an SP1 `.r1lf` lifted R1CS is satisfied by a provided witness (research).
//!
//! This is **not** a full LF+ `PlusProver/PlusVerifier` run (that requires materializing full
//! matrices and commitments). Instead, it performs the two correctness checks we need before
//! attempting folding:
//! - **Boundedness**: each witness value is representable in balanced base \(2^16\) with `k=2`.
//! - **Lifted R1CS satisfaction**: for each constraint row, (A·w) * (B·w) == (C·w) in Frog base field.
//!
//! Usage:
//!   SP1_R1LF=/tmp/shrink_verifier.r1lf \
//!   SP1_WITNESS=/tmp/shrink_verifier.witness.u64le \
//!     cargo run -p latticefold-plus --example lf_plus_sp1_r1lf_verify_witness --features we_gate --release
//!
//! Optional env:
//! - CHECK_ALL=1         (default: 0) verify all chunks; otherwise only chunk 0.
//! - CHECK_BOUNDS=1      (default: 1) verify base2^16,k=2 boundedness on witness.
//! - FAIL_FAST=1         (default: 1) stop at first failing constraint; otherwise count failures.
//! - CHUNK_SIZE=...      (default: 1<<20) must match cache build.
//! - PAD_COLS=...        (default: 256) must match cache build.

#![cfg(feature = "we_gate")]

use ark_ff::{BigInteger, Field, PrimeField};
use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as R;
use stark_rings::PolyRing;

type F = <R as PolyRing>::BaseRing;

#[inline]
fn coeff0(c: &R) -> F {
    c.coeffs()[0]
}

#[inline]
fn row_dot_base(row: &[(R, usize)], w: &[F]) -> F {
    let mut acc = F::ZERO;
    for (c, j) in row {
        let wj = w.get(*j).copied().unwrap_or(F::ZERO);
        acc += coeff0(c) * wj;
    }
    acc
}

fn modulus_u128<P: PrimeField>() -> u128 {
    let le = P::MODULUS.to_bytes_le();
    let mut out = 0u128;
    for (i, b) in le.iter().enumerate().take(16) {
        out |= (*b as u128) << (8 * i);
    }
    out
}

/// Centered embedding: map x in [0,q) to a signed integer in roughly (-q/2, q/2].
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

/// Check representability in balanced base 2^16 with k=2 digits (i16 range).
#[inline]
fn check_base2_16_k2(centered: i128) -> bool {
    const B: i128 = 1i128 << 16;
    const HALF: i128 = 1i128 << 15;
    let mut x = centered;
    for _ in 0..2 {
        let mut d = x.rem_euclid(B); // 0..B-1
        if d >= HALF {
            d -= B; // now in [-HALF, HALF)
        }
        // Digits must fit i16 in the DigitMatrix representation.
        if d < i16::MIN as i128 || d > i16::MAX as i128 {
            return false;
        }
        x = (x - d) / B;
    }
    x == 0
}

fn main() {
    let path = std::env::var("SP1_R1LF").expect("Set SP1_R1LF=/path/to/file.r1lf");
    let witness_path = std::env::var("SP1_WITNESS").expect("Set SP1_WITNESS=/path/to/witness.u64le");
    let chunk_size: usize = std::env::var("CHUNK_SIZE")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(1 << 20);
    let pad_cols_to_multiple_of: usize = std::env::var("PAD_COLS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(256);
    let check_all: bool = std::env::var("CHECK_ALL").ok().as_deref() == Some("1");
    let check_bounds: bool = std::env::var("CHECK_BOUNDS")
        .ok()
        .as_deref()
        .map(|s| s != "0")
        .unwrap_or(true);
    let fail_fast: bool = std::env::var("FAIL_FAST")
        .ok()
        .as_deref()
        .map(|s| s != "0")
        .unwrap_or(true);

    println!("=========================================================");
    println!("LF+ SP1 R1LF Verify Witness (boundedness + R1CS check)");
    println!("=========================================================");
    println!("  CHUNK_SIZE={chunk_size} PAD_COLS={pad_cols_to_multiple_of} CHECK_ALL={check_all} CHECK_BOUNDS={check_bounds} FAIL_FAST={fail_fast}");

    let t0 = std::time::Instant::now();
    let cache =
        latticefold_plus::sp1_r1lf::open_sp1_r1lf_chunk_cache::<R>(&path, chunk_size, pad_cols_to_multiple_of)
            .expect("open R1LF chunk cache");
    println!("  open/build cache: {:?}", t0.elapsed());
    println!(
        "  header: num_vars={} num_constraints={} num_public={} p_bb={} chunks={} ncols={}",
        cache.stats.num_vars, cache.stats.num_constraints, cache.stats.num_public, cache.stats.p_bb, cache.num_chunks, cache.ncols
    );
    println!("  digest={:02x?}...", &cache.stats.digest[..8]);

    let (w_u64, base_len, aux_len) =
        latticefold_plus::sp1_witness_io::load_sp1_witness_any(&witness_path, cache.stats.num_vars)
            .expect("load witness");
    println!("  loaded witness: base={} aux={} full={}", base_len, aux_len, w_u64.len());
    assert!(!w_u64.is_empty() && w_u64[0] == 1, "witness must have w[0]=1");

    // Map u64 witness -> Frog base field scalars once.
    let t_w = std::time::Instant::now();
    let w: Vec<F> = w_u64.into_iter().map(F::from).collect();
    println!("  map witness u64->F: {:?}", t_w.elapsed());

    if check_bounds {
        let q = modulus_u128::<F>();
        let t_b = std::time::Instant::now();

        #[cfg(feature = "parallel")]
        {
            use rayon::prelude::*;
            let bad = w
                .par_iter()
                .enumerate()
                .find_any(|(_i, fx)| {
                    // Recover canonical u64 via bigint limbs (fast path for 64-bit fields).
                    let bi = fx.into_bigint();
                    let x = bi.as_ref()[0];
                    let centered = centered_i128_from_canonical_u64(x, q);
                    !check_base2_16_k2(centered)
                });
            if let Some((i, fx)) = bad {
                let bi = fx.into_bigint();
                let x = bi.as_ref()[0];
                let centered = centered_i128_from_canonical_u64(x, q);
                panic!("boundedness failed at idx={i}: x={x} centered={centered}");
            }
        }
        #[cfg(not(feature = "parallel"))]
        {
            for (i, fx) in w.iter().enumerate() {
                let bi = fx.into_bigint();
                let x = bi.as_ref()[0];
                let centered = centered_i128_from_canonical_u64(x, q);
                if !check_base2_16_k2(centered) {
                    panic!("boundedness failed at idx={i}: x={x} centered={centered}");
                }
            }
        }

        println!("  boundedness (base2^16,k=2) OK: {:?}", t_b.elapsed());
    }

    let num_chunks = if check_all { cache.num_chunks } else { 1 };
    let mut total_rows_checked: u64 = 0;
    let mut total_failures: u64 = 0;
    let t_check = std::time::Instant::now();

    for chunk_idx in 0..num_chunks {
        let t_c = std::time::Instant::now();
        let [a, b, c] = cache.read_chunk(chunk_idx).expect("read_chunk");
        println!(
            "  chunk {chunk_idx}/{num_chunks}: read {:?} (nrows={})",
            t_c.elapsed(),
            a.nrows
        );

        let nrows = a.nrows;
        total_rows_checked += nrows as u64;

        #[cfg(feature = "parallel")]
        {
            use rayon::prelude::*;
            if fail_fast {
                let bad = (0..nrows).into_par_iter().find_any(|&i| {
                    let av = row_dot_base(&a.coeffs[i], &w);
                    let bv = row_dot_base(&b.coeffs[i], &w);
                    let cv = row_dot_base(&c.coeffs[i], &w);
                    av * bv != cv
                });
                if let Some(i) = bad {
                    let global = (chunk_idx as u64) * (chunk_size as u64) + (i as u64);
                    let av = row_dot_base(&a.coeffs[i], &w);
                    let bv = row_dot_base(&b.coeffs[i], &w);
                    let cv = row_dot_base(&c.coeffs[i], &w);
                    panic!("R1CS check failed at global_row={global} (chunk={chunk_idx} row={i}): (A·w)*(B·w) != (C·w)");
                }
            } else {
                let fails = (0..nrows)
                    .into_par_iter()
                    .map(|i| {
                        let av = row_dot_base(&a.coeffs[i], &w);
                        let bv = row_dot_base(&b.coeffs[i], &w);
                        let cv = row_dot_base(&c.coeffs[i], &w);
                        if av * bv != cv { 1u64 } else { 0u64 }
                    })
                    .sum::<u64>();
                total_failures += fails;
                if fails > 0 {
                    println!("    chunk {chunk_idx}: failures={fails}");
                }
            }
        }

        #[cfg(not(feature = "parallel"))]
        {
            for i in 0..nrows {
                let av = row_dot_base(&a.coeffs[i], &w);
                let bv = row_dot_base(&b.coeffs[i], &w);
                let cv = row_dot_base(&c.coeffs[i], &w);
                if av * bv != cv {
                    let global = (chunk_idx as u64) * (chunk_size as u64) + (i as u64);
                    if fail_fast {
                        panic!("R1CS check failed at global_row={global} (chunk={chunk_idx} row={i})");
                    } else {
                        total_failures += 1;
                    }
                }
            }
        }
    }

    println!(
        "  R1CS check done: rows_checked={} failures={} time={:?}",
        total_rows_checked,
        total_failures,
        t_check.elapsed()
    );
    if total_failures == 0 {
        println!("  OK");
    } else {
        panic!("R1CS check failed (failures={total_failures})");
    }
}

