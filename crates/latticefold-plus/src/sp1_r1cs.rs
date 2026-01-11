//! SP1 shrink verifier R1CS loader helpers (research).
//!
//! This module intentionally reuses Symphony’s chunk-cache loader so we can do apples-to-apples
//! experiments without duplicating the file format parsing logic.
//!
//! Enable with `--features sp1_import`.

#![cfg(feature = "sp1_import")]

use symphony::sp1_r1cs_loader::FieldFromU64;
use symphony::symphony_sp1_r1cs::open_sp1_r1cs_chunk_cache;
use ark_ff::PrimeField;
use stark_rings::{OverField, PolyRing, Zq};

/// Open/build the SP1 R1CS chunk cache, returning Symphony’s cache object.
///
/// The cache can stream chunks (A,B,C) as `stark_rings_linalg::SparseMatrix<R>` over the provided ring `R`.
pub fn open_sp1_cache<R, F>(
    r1cs_path: &str,
    chunk_size: usize,
    pad_cols_to_multiple_of: usize,
) -> Result<symphony::symphony_sp1_r1cs::ChunkCache<R>, String>
where
    R: OverField + PolyRing,
    // Matches Symphony loader expectations (chunk cache stores u64 coeffs embedded into BaseRing).
    R::BaseRing: Zq + PrimeField + From<u64> + Send + Sync,
    F: FieldFromU64 + Clone + Send + Sync,
{
    open_sp1_r1cs_chunk_cache::<R, F>(r1cs_path, chunk_size, pad_cols_to_multiple_of)
        .map_err(|e| format!("{e}"))
}

