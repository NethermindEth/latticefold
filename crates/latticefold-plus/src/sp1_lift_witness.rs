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

/// Lightweight "raw" sparse matrix row term used for scanning.
///
/// (coeff_i64, col_idx)
type RowTermsI64 = Vec<(i64, usize)>;

#[inline]
fn row_is_one(row: &RowTermsI64) -> bool {
    row.len() == 1 && row[0].0 == 1 && row[0].1 == 0
}

#[inline]
fn find_p_bb_term(row: &RowTermsI64, p_bb: i64) -> Option<usize> {
    // The lift adds (+p_bb)*v to C, and coefficients are centered for all original terms,
    // so encountering ±p_bb is an unambiguous marker for the aux var.
    row.iter()
        .find_map(|(coeff, col)| (*coeff == p_bb || *coeff == -p_bb).then_some(*col))
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
    let p_bb_i64: i64 = cache.stats.p_bb as i64;
    let mut kinds: Vec<Option<AuxVarKind>> = vec![None; cache.stats.num_vars];
    let mut summary = AuxVarSummary {
        num_constraints: cache.stats.num_constraints,
        ..Default::default()
    };

    // We need raw i64 coefficients; use the chunk-cache file directly.
    //
    // Cache file format (per chunk):
    //   nrows(u64), then 3 matrices, each is nrows rows of (num_terms(u32), [col(u32), coeff(i64)]...).
    //
    // We only need to scan rows, so we parse per-row without building ring elements.
    use std::io::{BufReader, Read, Seek, SeekFrom};
    const IO_BUFFER_SIZE: usize = 256 * 1024 * 1024;

    let file = std::fs::File::open(&cache.cache_path)?;
    let mut r = BufReader::with_capacity(IO_BUFFER_SIZE, file);

    // Jump to offset table end for chunk parsing: we already have offsets, so just seek per chunk.
    let mut buf4 = [0u8; 4];
    let mut buf8 = [0u8; 8];

    for chunk_idx in 0..cache.num_chunks {
        r.seek(SeekFrom::Start(cache.chunk_offsets[chunk_idx]))?;
        r.read_exact(&mut buf8)?;
        let nrows = u64::from_le_bytes(buf8) as usize;

        // Scan A,B,C row-by-row for this chunk.
        let mut a_rows: Vec<RowTermsI64> = Vec::with_capacity(nrows);
        let mut b_rows: Vec<RowTermsI64> = Vec::with_capacity(nrows);

        for rows_out in [&mut a_rows, &mut b_rows] {
            for _ in 0..nrows {
                r.read_exact(&mut buf4)?;
                let num_terms = u32::from_le_bytes(buf4) as usize;
                let mut row: RowTermsI64 = Vec::with_capacity(num_terms);
                for _ in 0..num_terms {
                    r.read_exact(&mut buf4)?;
                    let col_idx = u32::from_le_bytes(buf4) as usize;
                    r.read_exact(&mut buf8)?;
                    let coeff = i64::from_le_bytes(buf8);
                    if coeff != 0 {
                        row.push((coeff, col_idx));
                    }
                }
                rows_out.push(row);
            }
        }

        // Now C rows: identify the aux var per row.
        for row_idx in 0..nrows {
            r.read_exact(&mut buf4)?;
            let num_terms = u32::from_le_bytes(buf4) as usize;
            let mut crow: RowTermsI64 = Vec::with_capacity(num_terms);
            for _ in 0..num_terms {
                r.read_exact(&mut buf4)?;
                let col_idx = u32::from_le_bytes(buf4) as usize;
                r.read_exact(&mut buf8)?;
                let coeff = i64::from_le_bytes(buf8);
                if coeff != 0 {
                    crow.push((coeff, col_idx));
                }
            }

            if let Some(v) = find_p_bb_term(&crow, p_bb_i64) {
                let is_linear = row_is_one(&a_rows[row_idx]) || row_is_one(&b_rows[row_idx]);
                let kind = if is_linear { AuxVarKind::Carry } else { AuxVarKind::Quotient };
                if v < kinds.len() {
                    if kinds[v].is_none() {
                        summary.num_aux_vars += 1;
                        match kind {
                            AuxVarKind::Carry => summary.num_carry += 1,
                            AuxVarKind::Quotient => summary.num_quotient += 1,
                        }
                    }
                    kinds[v] = Some(kind);
                }
            }
        }
    }

    Ok((kinds, summary))
}

