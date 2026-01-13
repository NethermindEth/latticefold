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

/// SP1 lift semantics: interpret a BabyBear residue `x in [0,p_bb)` as a centered integer in
/// (-(p_bb-1)/2 .. (p_bb-1)/2], then embed into the host field `F`.
///
/// IMPORTANT: this is **not** the same as `F::from(x)` because the host modulus is not `p_bb`.
#[inline]
fn babybear_u64_to_centered_host(x: u64, p_bb: u64) -> F {
    debug_assert!(p_bb > 1);
    let half = p_bb / 2;
    if x > half {
        // centered negative: x - p_bb
        let neg = p_bb - x;
        -F::from(neg)
    } else {
        F::from(x)
    }
}

#[inline]
fn bb_centered_i128(x: u64, p_bb: u64) -> i128 {
    let x = x as i128;
    let p = p_bb as i128;
    let half = p / 2;
    if x > half { x - p } else { x }
}

#[inline]
fn row_dot_base(row: &[(F, usize)], w: &[F]) -> F {
    let mut acc = F::ZERO;
    for (c, j) in row {
        let wj = w.get(*j).copied().unwrap_or(F::ZERO);
        acc += (*c) * wj;
    }
    acc
}

fn read_row_terms_from_chunk_cache(
    cache_path: &str,
    chunk_idx: usize,
    row_idx: usize,
    which_matrix: usize, // 0=A,1=B,2=C
) -> Result<(usize, Vec<(u32, i64)>), String> {
    use std::io::{Read, Seek, SeekFrom};

    let mut f = std::fs::File::open(cache_path).map_err(|e| format!("open {cache_path}: {e}"))?;
    let mut buf4 = [0u8; 4];
    let mut buf8 = [0u8; 8];

    // Header: magic(4) version(4) digest(32) + 8 u64 fields (p_bb, num_vars, num_constraints, num_public,
    // nnz, chunk_size, ncols, num_chunks) => 104 bytes total.
    let mut hdr = [0u8; 4];
    f.read_exact(&mut hdr).map_err(|e| format!("read magic: {e}"))?;
    if &hdr != b"LFC1" {
        return Err("bad chunk cache magic (expected LFC1)".to_string());
    }
    f.read_exact(&mut buf4).map_err(|e| format!("read version: {e}"))?;
    let version = u32::from_le_bytes(buf4);
    if version != 1 {
        return Err(format!("bad chunk cache version {version} (expected 1)"));
    }
    let mut digest = [0u8; 32];
    f.read_exact(&mut digest).map_err(|e| format!("read digest: {e}"))?;
    // Skip p_bb..ncols (7 u64s)
    for _ in 0..7 {
        f.read_exact(&mut buf8).map_err(|e| format!("read hdr u64: {e}"))?;
    }
    // num_chunks
    f.read_exact(&mut buf8).map_err(|e| format!("read num_chunks: {e}"))?;
    let num_chunks = u64::from_le_bytes(buf8) as usize;
    if chunk_idx >= num_chunks {
        return Err(format!("chunk_idx out of range: {chunk_idx} >= {num_chunks}"));
    }

    let offsets_start = 104u64;
    f.seek(SeekFrom::Start(offsets_start + (chunk_idx as u64) * 8))
        .map_err(|e| format!("seek offsets: {e}"))?;
    f.read_exact(&mut buf8).map_err(|e| format!("read chunk offset: {e}"))?;
    let chunk_off = u64::from_le_bytes(buf8);
    f.seek(SeekFrom::Start(chunk_off))
        .map_err(|e| format!("seek chunk: {e}"))?;
    f.read_exact(&mut buf8).map_err(|e| format!("read nrows: {e}"))?;
    let nrows = u64::from_le_bytes(buf8) as usize;
    if row_idx >= nrows {
        return Err(format!("row_idx out of range: {row_idx} >= {nrows}"));
    }

    fn skip_row(
        f: &mut std::fs::File,
        buf4: &mut [u8; 4],
    ) -> Result<(), String> {
        use std::io::{Read, Seek, SeekFrom};
        f.read_exact(buf4).map_err(|e| format!("read num_terms: {e}"))?;
        let nt = u32::from_le_bytes(*buf4) as u64;
        // each term: col(u32)+coeff(i64) = 12 bytes
        f.seek(SeekFrom::Current((nt * 12) as i64))
            .map_err(|e| format!("skip terms: {e}"))?;
        Ok(())
    }

    fn read_row(
        f: &mut std::fs::File,
        buf4: &mut [u8; 4],
        buf8: &mut [u8; 8],
    ) -> Result<Vec<(u32, i64)>, String> {
        use std::io::Read;
        f.read_exact(buf4).map_err(|e| format!("read num_terms: {e}"))?;
        let nt = u32::from_le_bytes(*buf4) as usize;
        let mut out = Vec::with_capacity(nt);
        for _ in 0..nt {
            f.read_exact(buf4).map_err(|e| format!("read col: {e}"))?;
            let col = u32::from_le_bytes(*buf4);
            f.read_exact(buf8).map_err(|e| format!("read coeff: {e}"))?;
            let coeff = i64::from_le_bytes(*buf8);
            if coeff != 0 {
                out.push((col, coeff));
            }
        }
        Ok(out)
    }

    for m in 0..3 {
        for r in 0..nrows {
            if m == which_matrix && r == row_idx {
                let row = read_row(&mut f, &mut buf4, &mut buf8)?;
                return Ok((nrows, row));
            } else {
                skip_row(&mut f, &mut buf4)?;
            }
        }
    }
    Err("unreachable: failed to locate row".to_string())
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

    // Infer `base_vars` from the first (smallest) column index that appears with coeff ±p_bb in C.
    // This is useful for sanity/telemetry, but **SP1 evaluates all witness slots using centered
    // BabyBear integers**, including aux vars (see `eval_row_i128` in `sp1/.../r1cs/lf.rs`).
    let p_bb = cache.stats.p_bb;
    let t_inf = std::time::Instant::now();
    let mut base_vars = cache.stats.num_vars; // fallback if no aux terms exist
    'outer: for chunk_idx in 0..cache.num_chunks {
        let [_a, _b, c] = cache
            .read_chunk(chunk_idx)
            .expect("read_chunk for infer_base_vars");
        for row in &c.coeffs {
            for (coeff, col_idx) in row {
                // Coefficients are stored as base-ring scalars in the chunk cache.
                if *coeff == F::from(p_bb) || *coeff == -F::from(p_bb) {
                    base_vars = base_vars.min(*col_idx);
                }
            }
        }
        if base_vars != cache.stats.num_vars {
            break 'outer;
        }
    }
    println!(
        "  inferred base_vars={} (infer {:?})",
        base_vars,
        t_inf.elapsed()
    );

    // Map u64 witness -> Frog base field scalars once, matching SP1 lift semantics:
    // - **all vars (including aux)**: centered embedding mod p_bb
    let t_w = std::time::Instant::now();
    let w_host: Vec<F> = w_u64
        .iter()
        .copied()
        .enumerate()
        .map(|(i, x)| {
            if x >= p_bb {
                panic!("witness word out of [0,p_bb) range at idx={i}: x={x} p_bb={p_bb}");
            }
            let _ = i; // keep index for error context above
            babybear_u64_to_centered_host(x, p_bb)
        })
        .collect();
    println!("  map witness u64->F: {:?}", t_w.elapsed());

    if check_bounds {
        let q = modulus_u128::<F>();
        let t_b = std::time::Instant::now();

        #[cfg(feature = "parallel")]
        {
            use rayon::prelude::*;
            let bad = w_host
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
            for (i, fx) in w_host.iter().enumerate() {
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
                // Deterministic fail-fast: pick the *smallest* failing row index.
                // (Rayon's `find_any` is intentionally nondeterministic.)
                let bad = (0..nrows)
                    .into_par_iter()
                    .filter(|&i| {
                    let av = row_dot_base(&a.coeffs[i], &w_host);
                    let bv = row_dot_base(&b.coeffs[i], &w_host);
                    let cv = row_dot_base(&c.coeffs[i], &w_host);
                    av * bv != cv
                    })
                    .reduce_with(|x, y| x.min(y));
                if let Some(i) = bad {
                    let global = (chunk_idx as u64) * (chunk_size as u64) + (i as u64);
                    let av = row_dot_base(&a.coeffs[i], &w_host);
                    let bv = row_dot_base(&b.coeffs[i], &w_host);
                    let cv = row_dot_base(&c.coeffs[i], &w_host);

                    // Deep debug: re-read raw i64 coefficients from the `{path}.chunks` cache and
                    // evaluate in SP1 integer semantics (centered mod p_bb).
                    let cache_path = format!("{path}.chunks");
                    let (nrows0, a_i64) =
                        read_row_terms_from_chunk_cache(&cache_path, chunk_idx, i, 0).expect("read A row i64");
                    let (_nrows1, b_i64) =
                        read_row_terms_from_chunk_cache(&cache_path, chunk_idx, i, 1).expect("read B row i64");
                    let (_nrows2, c_i64) =
                        read_row_terms_from_chunk_cache(&cache_path, chunk_idx, i, 2).expect("read C row i64");
                    eprintln!("  [debug] chunk_nrows={nrows0} row_in_chunk={i} global_row={global}");
                    eprintln!("  [debug] A terms={} B terms={} C terms={}", a_i64.len(), b_i64.len(), c_i64.len());
                    for (col, coeff) in &a_i64 {
                        let wu = w_u64[*col as usize];
                        eprintln!(
                            "  [debug] A term: coeff={} col={} w_u64={} w_centered={}",
                            coeff,
                            col,
                            wu,
                            bb_centered_i128(wu, p_bb)
                        );
                    }
                    for (col, coeff) in &b_i64 {
                        let wu = w_u64[*col as usize];
                        eprintln!(
                            "  [debug] B term: coeff={} col={} w_u64={} w_centered={}",
                            coeff,
                            col,
                            wu,
                            bb_centered_i128(wu, p_bb)
                        );
                    }
                    // Print any ±p_bb term in C.
                    for (col, coeff) in &c_i64 {
                        if *coeff == p_bb as i64 || *coeff == -(p_bb as i64) {
                            let wu = w_u64[*col as usize];
                            eprintln!(
                                "  [debug] C has p term: coeff={} col={} w_u64={} w_centered={}",
                                coeff,
                                col,
                                wu,
                                bb_centered_i128(wu, p_bb)
                            );
                        }
                    }
                    let eval_int = |terms: &[(u32, i64)]| -> i128 {
                        terms
                            .iter()
                            .map(|(col, coeff)| (*coeff as i128) * bb_centered_i128(w_u64[*col as usize], p_bb))
                            .sum()
                    };
                    let a_int = eval_int(&a_i64);
                    let b_int = eval_int(&b_i64);
                    let c_int = eval_int(&c_i64);
                    let diff = a_int * b_int - c_int;
                    eprintln!("  [debug] int: a={} b={} c={} a*b-c={}", a_int, b_int, c_int, diff);
                    eprintln!("  [debug] int mod p_bb: {}", diff.rem_euclid(p_bb as i128));
                    // Also show host-field lhs/rhs limbs (u64 canonical) for quick inspection.
                    let av_u = av.into_bigint().as_ref()[0];
                    let bv_u = bv.into_bigint().as_ref()[0];
                    let cv_u = cv.into_bigint().as_ref()[0];
                    eprintln!("  [debug] host(F): av={} bv={} cv={}", av_u, bv_u, cv_u);

                    panic!("R1CS check failed at global_row={global} (chunk={chunk_idx} row={i}): (A·w)*(B·w) != (C·w)");
                }
            } else {
                let fails = (0..nrows)
                    .into_par_iter()
                    .map(|i| {
                        let av = row_dot_base(&a.coeffs[i], &w_host);
                        let bv = row_dot_base(&b.coeffs[i], &w_host);
                        let cv = row_dot_base(&c.coeffs[i], &w_host);
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
                let av = row_dot_base(&a.coeffs[i], &w_host);
                let bv = row_dot_base(&b.coeffs[i], &w_host);
                let cv = row_dot_base(&c.coeffs[i], &w_host);
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

