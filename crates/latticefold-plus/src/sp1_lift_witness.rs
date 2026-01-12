//! Helpers to compute (or at least locate) the auxiliary `q_i` / carry vars introduced by the
//! SP1 BabyBear→LF lift (R1LF).
//!
//! Why this exists:
//! - The lifted relation adds one fresh witness variable per lifted constraint via a `(+p_bb) * v`
//!   term in the C-row.
//! - **Soundness requires boundedness checks to cover these fresh vars too** (otherwise the `+p_bb*v`
//!   term is vacuous mod q).
//!
//! This module is intentionally minimal: it provides a streaming scan over a `.r1lf.chunks` cache to
//! identify which columns are aux vars and which constraints are "linear" vs "true-mul" (matching the
//! SP1 lift heuristic).

#![cfg(feature = "we_gate")]

use crate::sp1_r1lf::R1LfChunkCache;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// Classification of an aux var introduced by the lift.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AuxVarKind {
    /// Added for a linear constraint (A==1 || B==1 in SP1 lift).
    Carry,
    /// Added for a true multiplication constraint.
    Quotient,
}

/// Summary of aux variables implied by an R1LF relation.
#[derive(Debug, Clone, Default)]
pub struct AuxVarSummary {
    /// Total number of constraints scanned.
    pub num_constraints: usize,
    /// Distinct aux vars found.
    pub num_aux_vars: usize,
    /// Count of aux vars classified as `Carry`.
    pub num_carry: usize,
    /// Count of aux vars classified as `Quotient`.
    pub num_quotient: usize,
}

/// Scan a `.chunks` cache and extract which witness columns are aux vars, classified by kind.
///
/// Returns:
/// - `kinds`: vector of length `cache.stats.num_vars` where `Some(kind)` marks an aux var column.
/// - `summary`: counts.
///
/// Notes:
/// - This does **not** compute the aux witness *values* (only locates/classifies the columns).
/// - This is enough to wire boundedness checks over the correct committed witness surface.
pub fn scan_aux_vars_from_r1lf_chunks<R>(
    cache: &R1LfChunkCache<R>,
) -> std::io::Result<(Vec<Option<AuxVarKind>>, AuxVarSummary)> {
    let profile = std::env::var("LF_PLUS_PROFILE").ok().as_deref() == Some("1");
    let p_bb_i64: i64 = cache.stats.p_bb as i64;

    // Pull the minimal scan inputs out of `cache` so parallel code doesn't capture `R`.
    let cache_path = cache.cache_path.clone();
    let chunk_offsets = cache.chunk_offsets.clone();
    let num_chunks = cache.num_chunks;

    // We need raw i64 coefficients; use the chunk-cache file directly.
    //
    // Cache file format (per chunk):
    //   nrows(u64), then 3 matrices, each is nrows rows of (num_terms(u32), [col(u32), coeff(i64)]...).
    //
    // Performance note:
    // - DO NOT allocate per-row vectors; for SP1 scale (nrows~1<<20) this causes huge allocator overhead.
    // - We only need to know if A-row or B-row is exactly `1` (i.e., [(col=0, coeff=1)]).
    // - So we store `a_is_one[row]` and `b_is_one[row]` as booleans, then scan C for the ±p_bb term.
    use std::io::{BufReader, Read, Seek, SeekFrom};
    const IO_BUFFER_SIZE: usize = 256 * 1024 * 1024;

    #[inline]
    fn scan_one_chunk(
        cache_path: &str,
        chunk_offset: u64,
        p_bb_i64: i64,
    ) -> std::io::Result<Vec<(usize, AuxVarKind)>> {
        let file = std::fs::File::open(cache_path)?;
        let mut r = BufReader::with_capacity(IO_BUFFER_SIZE, file);
        r.seek(SeekFrom::Start(chunk_offset))?;

        let mut buf4 = [0u8; 4];
        let mut buf8 = [0u8; 8];

        r.read_exact(&mut buf8)?;
        let nrows = u64::from_le_bytes(buf8) as usize;

        let mut a_is_one = vec![false; nrows];
        let mut b_is_one = vec![false; nrows];

        // A rows
        for row_idx in 0..nrows {
            r.read_exact(&mut buf4)?;
            let num_terms = u32::from_le_bytes(buf4) as usize;
            let mut is_one = num_terms == 1;
            for _ in 0..num_terms {
                r.read_exact(&mut buf4)?;
                let col_idx = u32::from_le_bytes(buf4) as usize;
                r.read_exact(&mut buf8)?;
                let coeff = i64::from_le_bytes(buf8);
                if is_one {
                    is_one = col_idx == 0 && coeff == 1;
                }
            }
            a_is_one[row_idx] = is_one;
        }

        // B rows
        for row_idx in 0..nrows {
            r.read_exact(&mut buf4)?;
            let num_terms = u32::from_le_bytes(buf4) as usize;
            let mut is_one = num_terms == 1;
            for _ in 0..num_terms {
                r.read_exact(&mut buf4)?;
                let col_idx = u32::from_le_bytes(buf4) as usize;
                r.read_exact(&mut buf8)?;
                let coeff = i64::from_le_bytes(buf8);
                if is_one {
                    is_one = col_idx == 0 && coeff == 1;
                }
            }
            b_is_one[row_idx] = is_one;
        }

        let mut out: Vec<(usize, AuxVarKind)> = Vec::new();
        // Now C rows: identify the aux var per row.
        for row_idx in 0..nrows {
            r.read_exact(&mut buf4)?;
            let num_terms = u32::from_le_bytes(buf4) as usize;
            let mut aux_col: Option<usize> = None;
            for _ in 0..num_terms {
                r.read_exact(&mut buf4)?;
                let col_idx = u32::from_le_bytes(buf4) as usize;
                r.read_exact(&mut buf8)?;
                let coeff = i64::from_le_bytes(buf8);
                if coeff == p_bb_i64 || coeff == -p_bb_i64 {
                    aux_col = Some(col_idx);
                }
            }

            if let Some(v) = aux_col {
                let is_linear = a_is_one[row_idx] || b_is_one[row_idx];
                let kind = if is_linear { AuxVarKind::Carry } else { AuxVarKind::Quotient };
                out.push((v, kind));
            }
        }

        Ok(out)
    }

    // 0=unset, 1=carry, 2=quotient
    let kinds_atomic: Vec<std::sync::atomic::AtomicU8> =
        (0..cache.stats.num_vars).map(|_| std::sync::atomic::AtomicU8::new(0)).collect();

    let (num_aux_vars, num_carry, num_quotient) = {
        #[cfg(feature = "parallel")]
        {
            if profile {
                eprintln!(
                    "[LF+ scan_aux] parallel scan (rayon_threads={})",
                    rayon::current_num_threads()
                );
            }
            let per_chunk = (0..num_chunks)
                .into_par_iter()
                .map(|chunk_idx| scan_one_chunk(&cache_path, chunk_offsets[chunk_idx], p_bb_i64))
                .collect::<Result<Vec<_>, _>>()?;
            for entries in per_chunk {
                for (v, kind) in entries {
                    if v >= kinds_atomic.len() {
                        continue;
                    }
                    let new = match kind {
                        AuxVarKind::Carry => 1u8,
                        AuxVarKind::Quotient => 2u8,
                    };
                    // Each aux var is unique per constraint in the lift, so this should be a one-time set.
                    // If it repeats, it should be consistent; we accept "already set".
                    let prev = kinds_atomic[v].load(std::sync::atomic::Ordering::Relaxed);
                    if prev == 0 {
                        let _ = kinds_atomic[v].compare_exchange(
                            0,
                            new,
                            std::sync::atomic::Ordering::Relaxed,
                            std::sync::atomic::Ordering::Relaxed,
                        );
                    }
                }
            }
        }

        #[cfg(not(feature = "parallel"))]
        {
            // Sequential scan when rayon is disabled.
            if profile {
                eprintln!("[LF+ scan_aux] sequential scan (feature=parallel is disabled)");
            }
            for chunk_idx in 0..num_chunks {
                if profile {
                    eprintln!(
                        "[LF+ scan_aux] chunk {}/{} ({:.1}%)",
                        chunk_idx + 1,
                        num_chunks,
                        100.0 * (chunk_idx as f64 + 1.0) / (num_chunks as f64)
                    );
                }
                let entries = scan_one_chunk(&cache_path, chunk_offsets[chunk_idx], p_bb_i64)?;
                for (v, kind) in entries {
                    if v >= kinds_atomic.len() {
                        continue;
                    }
                    let new = match kind {
                        AuxVarKind::Carry => 1u8,
                        AuxVarKind::Quotient => 2u8,
                    };
                    let prev = kinds_atomic[v].load(std::sync::atomic::Ordering::Relaxed);
                    if prev == 0 {
                        let _ = kinds_atomic[v].compare_exchange(
                            0,
                            new,
                            std::sync::atomic::Ordering::Relaxed,
                            std::sync::atomic::Ordering::Relaxed,
                        );
                    }
                }
            }
        }

        // Reduce atomics into counts.
        let mut num_aux = 0usize;
        let mut ncarry = 0usize;
        let mut nquot = 0usize;
        for a in &kinds_atomic {
            match a.load(std::sync::atomic::Ordering::Relaxed) {
                1 => {
                    num_aux += 1;
                    ncarry += 1;
                }
                2 => {
                    num_aux += 1;
                    nquot += 1;
                }
                _ => {}
            }
        }
        (num_aux, ncarry, nquot)
    };

    let kinds: Vec<Option<AuxVarKind>> = kinds_atomic
        .into_iter()
        .map(|a| match a.load(std::sync::atomic::Ordering::Relaxed) {
            1 => Some(AuxVarKind::Carry),
            2 => Some(AuxVarKind::Quotient),
            _ => None,
        })
        .collect();

    let summary = AuxVarSummary {
        num_constraints: cache.stats.num_constraints,
        num_aux_vars,
        num_carry,
        num_quotient,
    };

    Ok((kinds, summary))
}

