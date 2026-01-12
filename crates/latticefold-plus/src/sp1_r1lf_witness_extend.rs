//! Witness extension for SP1 R1LF lift vars (q/carry).
//!
//! Given:
//! - an R1LF `.chunks` cache (so we can stream rows cheaply),
//! - and a base witness for the **original** variables (indices `[0..orig_num_vars)`),
//! this computes the auxiliary lift variables (indices `[orig_num_vars..num_vars)`) such that each
//! lifted constraint row satisfies:
//!
//!   (A·w) * (B·w) = (C·w) + p_bb * t_i
//!
//! where `t_i` is the unique aux var identified by a `±p_bb` term in the C-row.
//!
//! Notes:
//! - This assumes the SP1 R1LF lift allocates aux vars by appending indices after the original
//!   witness (true for the current SP1 lift pass).
//! - We compute `t_i` as an **integer** quotient when possible:
//!     t_i = (A·w * B·w - C_noaux·w) / p_bb
//!   and require divisibility.
//! - We use i128 arithmetic for safety (row dots are sparse).

#![cfg(feature = "we_gate")]

use crate::sp1_r1lf::R1LfChunkCache;

#[derive(Debug, Clone)]
pub struct ExtendWitnessStats {
    pub num_constraints: usize,
    pub num_aux_vars: usize,
    pub num_assigned: usize,
    pub num_linear: usize,
    pub num_mul: usize,
}

/// Compute `orig_num_vars = num_vars - num_aux_vars`.
///
/// This matches the current SP1 lift which appends exactly one aux var per lifted constraint.
pub fn orig_num_vars_from_counts(num_vars: usize, num_aux_vars: usize) -> Result<usize, String> {
    num_vars
        .checked_sub(num_aux_vars)
        .ok_or_else(|| "num_aux_vars > num_vars".to_string())
}

/// Extend the witness by computing aux vars from the R1LF constraints.
///
/// - `w_base` must have length `orig_num_vars`
/// - Returns `w_full` of length `cache.stats.num_vars` where aux slots are filled.
///
/// This function is intentionally strict: it errors if an aux var is assigned inconsistently or if
/// divisibility by `p_bb` fails for any lifted row.
pub fn extend_witness_from_r1lf_chunks(
    cache: &R1LfChunkCache<stark_rings::cyclotomic_ring::models::frog_ring::RqPoly>,
    w_base: &[u64],
    num_aux_vars: usize,
) -> Result<(Vec<u64>, ExtendWitnessStats), String> {
    use stark_rings::cyclotomic_ring::models::frog_ring::Fq;
    use stark_rings::Zq;

    let num_vars = cache.stats.num_vars;
    let orig_num_vars = orig_num_vars_from_counts(num_vars, num_aux_vars)?;
    if w_base.len() != orig_num_vars {
        return Err(format!(
            "w_base length mismatch: expected orig_num_vars={} got {}",
            orig_num_vars,
            w_base.len()
        ));
    }
    if w_base.is_empty() || w_base[0] != 1 {
        return Err("w_base[0] must be 1 (R1CS constant ONE slot)".to_string());
    }

    let p_bb = cache.stats.p_bb as i128;
    let p_bb_u64 = cache.stats.p_bb;

    // Output witness: copy base, append zeros for aux.
    let mut w_full = vec![0u64; num_vars];
    w_full[..orig_num_vars].copy_from_slice(w_base);

    // Track assigned aux vars (only for aux range).
    let mut assigned = vec![false; num_aux_vars];

    // Parse the chunk cache directly for speed (raw i64 coeffs).
    use std::io::{BufReader, Read, Seek, SeekFrom};
    const IO_BUFFER_SIZE: usize = 256 * 1024 * 1024;
    let file = std::fs::File::open(&cache.cache_path).map_err(|e| format!("{e}"))?;
    let mut r = BufReader::with_capacity(IO_BUFFER_SIZE, file);

    let mut buf4 = [0u8; 4];
    let mut buf8 = [0u8; 8];

    let mut stats = ExtendWitnessStats {
        num_constraints: cache.stats.num_constraints,
        num_aux_vars,
        num_assigned: 0,
        num_linear: 0,
        num_mul: 0,
    };

    // Helper to read one row and compute dot product with current witness.
    // Returns (dot, is_one_row) where is_one_row means the row is exactly "1" (num_terms=1, (col=0,coeff=1)).
    fn read_row_dot(
        r: &mut BufReader<std::fs::File>,
        buf4: &mut [u8; 4],
        buf8: &mut [u8; 8],
        witness: &[u64],
    ) -> Result<(i128, bool), String> {
        r.read_exact(buf4).map_err(|e| format!("{e}"))?;
        let num_terms = u32::from_le_bytes(*buf4) as usize;
        let mut is_one = num_terms == 1;
        let mut acc: i128 = 0;
        for _ in 0..num_terms {
            r.read_exact(buf4).map_err(|e| format!("{e}"))?;
            let col_idx = u32::from_le_bytes(*buf4) as usize;
            r.read_exact(buf8).map_err(|e| format!("{e}"))?;
            let coeff = i64::from_le_bytes(*buf8) as i128;
            if is_one {
                is_one = col_idx == 0 && coeff == 1;
            }
            let wv = *witness
                .get(col_idx)
                .ok_or_else(|| format!("witness index out of range: {col_idx}"))? as i128;
            acc += coeff * wv;
        }
        Ok((acc, is_one))
    }

    for chunk_idx in 0..cache.num_chunks {
        // Seek to chunk start.
        r.seek(SeekFrom::Start(cache.chunk_offsets[chunk_idx]))
            .map_err(|e| format!("{e}"))?;
        r.read_exact(&mut buf8).map_err(|e| format!("{e}"))?;
        let nrows = u64::from_le_bytes(buf8) as usize;

        // Read A dots and "is_one" flags
        let mut a_dot = vec![0i128; nrows];
        let mut a_is_one = vec![false; nrows];
        for i in 0..nrows {
            let (d, one) = read_row_dot(&mut r, &mut buf4, &mut buf8, &w_full)?;
            a_dot[i] = d;
            a_is_one[i] = one;
        }

        // Read B dots and flags
        let mut b_dot = vec![0i128; nrows];
        let mut b_is_one = vec![false; nrows];
        for i in 0..nrows {
            let (d, one) = read_row_dot(&mut r, &mut buf4, &mut buf8, &w_full)?;
            b_dot[i] = d;
            b_is_one[i] = one;
        }

        // Read C rows: compute c_noaux dot and capture aux var index + sign.
        for row_idx in 0..nrows {
            r.read_exact(&mut buf4).map_err(|e| format!("{e}"))?;
            let num_terms = u32::from_le_bytes(buf4) as usize;
            let mut c_acc: i128 = 0;
            let mut aux_col: Option<usize> = None;
            let mut aux_sign: i128 = 0;
            for _ in 0..num_terms {
                r.read_exact(&mut buf4).map_err(|e| format!("{e}"))?;
                let col_idx = u32::from_le_bytes(buf4) as usize;
                r.read_exact(&mut buf8).map_err(|e| format!("{e}"))?;
                let coeff = i64::from_le_bytes(buf8) as i128;

                if coeff == p_bb || coeff == -p_bb {
                    aux_col = Some(col_idx);
                    aux_sign = if coeff == p_bb { 1 } else { -1 };
                    continue; // exclude aux term from c_acc
                }

                let wv = *w_full
                    .get(col_idx)
                    .ok_or_else(|| format!("witness index out of range: {col_idx}"))? as i128;
                c_acc += coeff * wv;
            }

            // Only lifted rows have an aux var.
            let Some(v_idx) = aux_col else { continue };

            if v_idx < orig_num_vars {
                return Err(format!("aux var index {v_idx} is < orig_num_vars {orig_num_vars}"));
            }
            if v_idx >= num_vars {
                return Err(format!("aux var index {v_idx} out of range num_vars {num_vars}"));
            }
            let aux_slot = v_idx - orig_num_vars;

            let is_linear = a_is_one[row_idx] || b_is_one[row_idx];
            if is_linear {
                stats.num_linear += 1;
            } else {
                stats.num_mul += 1;
            }

            // Compute t as integer quotient:
            //   aux_sign * p * t = a*b - c_noaux
            // so t = (a*b - c_noaux) / (aux_sign * p)
            let num = a_dot[row_idx] * b_dot[row_idx] - c_acc;
            let denom = aux_sign * p_bb;
            if denom == 0 {
                return Err("denom=0 for aux var".to_string());
            }
            if num % denom != 0 {
                return Err(format!(
                    "non-divisible lifted row: (a*b-c) mod (p_bb) != 0 (row_in_chunk={}, chunk_idx={})",
                    row_idx, chunk_idx
                ));
            }
            let t_int = num / denom;
            // Map back into field element representative in [0,p_bb) if possible.
            if t_int < 0 || t_int >= (p_bb_u64 as i128) {
                // For now, require canonical quotient/carry in [0,p_bb). This matches the intended lift semantics.
                return Err(format!("t out of canonical range [0,p_bb): t={t_int}"));
            }
            let t_u64 = t_int as u64;

            if assigned[aux_slot] {
                if w_full[v_idx] != t_u64 {
                    return Err(format!("aux var assigned inconsistently at idx={v_idx}"));
                }
            } else {
                assigned[aux_slot] = true;
                w_full[v_idx] = t_u64;
                stats.num_assigned += 1;
            }
        }
    }

    // Ensure all aux vars got assigned.
    if stats.num_assigned != num_aux_vars {
        return Err(format!(
            "not all aux vars assigned: assigned={} expected={}",
            stats.num_assigned, num_aux_vars
        ));
    }

    // Optional: quick boundedness sanity for b=2^16,k=2 (centered magnitude < 2^31).
    // This is NOT a full LF+ rgchk proof; just a cheap preflight check.
    let bound = 1i128 << 31;
    for &x in &w_full {
        let fq = Fq::from(x);
        let mag = fq.center().to_u64().unwrap_or(u64::MAX) as i128;
        if mag >= bound {
            return Err("witness entry exceeds 2^31 bound; would not fit k=2,b=2^16".to_string());
        }
    }

    Ok((w_full, stats))
}

