//! SP1 witness file I/O helpers.
//!
//! Convention:
//! - `WITNESS_PATH` is a **single-file full witness** of length `num_vars` from the R1LF header
//!   (u64-le words).

#![cfg(feature = "we_gate")]

use std::io::Read;

pub fn read_u64le_vec(path: &str) -> std::io::Result<Vec<u64>> {
    let mut f = std::fs::File::open(path)?;
    let mut bytes = Vec::new();
    f.read_to_end(&mut bytes)?;
    if bytes.len() % 8 != 0 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "witness file length not a multiple of 8",
        ));
    }
    Ok(bytes
        .chunks_exact(8)
        .map(|c| u64::from_le_bytes(c.try_into().unwrap()))
        .collect())
}

/// Load a witness in either of these formats:
/// - **Single file only**: `witness_path` is the full witness of length `num_vars`.
///
/// Returns `(full_witness, base_len, aux_len)`.
pub fn load_sp1_witness_any(
    witness_path: &str,
    num_vars: usize,
) -> Result<(Vec<u64>, usize, usize), String> {
    let base = read_u64le_vec(witness_path).map_err(|e| format!("read {witness_path}: {e}"))?;
    if base.len() == num_vars {
        // Single-file full witness.
        if base.is_empty() || base[0] != 1 {
            return Err("witness must have w[0]=1 (constant ONE slot)".to_string());
        }
        return Ok((base, num_vars, 0));
    }

    Err(format!(
        "witness length mismatch: got len={} expected num_vars={}. \
Provide the single full witness file emitted by the SP1 exporter.",
        base.len(),
        num_vars
    ))
}

