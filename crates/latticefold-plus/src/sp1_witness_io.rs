//! SP1 witness file I/O helpers.
//!
//! Convention:
//! - If `WITNESS_PATH` is provided, we require a sibling `WITNESS_PATH.aux`.
//! - `WITNESS_PATH` contains the **base** witness for original vars (u64-le words).
//! - `WITNESS_PATH.aux` contains the **aux** witness values (q/carry) appended by the lift.
//! - Full witness is `base || aux` and must have length `num_vars` from the R1LF header.

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
/// - **Single file**: `witness_path` is the full witness of length `num_vars`.
/// - **Split files**: `witness_path` is the base witness and `witness_path.aux` is the aux tail,
///   and `base_len + aux_len == num_vars`.
///
/// Returns `(full_witness, base_len, aux_len)`.
pub fn load_sp1_witness_any(
    witness_path: &str,
    num_vars: usize,
) -> Result<(Vec<u64>, usize, usize), String> {
    let aux_path = format!("{witness_path}.aux");

    let base = read_u64le_vec(witness_path).map_err(|e| format!("read {witness_path}: {e}"))?;
    if base.len() == num_vars {
        // Single-file full witness.
        if base.is_empty() || base[0] != 1 {
            return Err("witness must have w[0]=1 (constant ONE slot)".to_string());
        }
        return Ok((base, num_vars, 0));
    }

    // Split-file witness (base + aux).
    if !std::path::Path::new(&aux_path).exists() {
        return Err(format!(
            "witness file is not full-length (len={}, expected={}), so we require a sibling aux file: `{aux_path}`",
            base.len(),
            num_vars
        ));
    }
    let aux = read_u64le_vec(&aux_path).map_err(|e| format!("read {aux_path}: {e}"))?;

    if base.len() + aux.len() != num_vars {
        return Err(format!(
            "witness length mismatch: base_len={} aux_len={} total={} expected_num_vars={}",
            base.len(),
            aux.len(),
            base.len() + aux.len(),
            num_vars
        ));
    }
    if base.is_empty() || base[0] != 1 {
        return Err("base witness must have w[0]=1 (constant ONE slot)".to_string());
    }

    let mut full = Vec::with_capacity(num_vars);
    full.extend_from_slice(&base);
    full.extend_from_slice(&aux);
    Ok((full, base.len(), aux.len()))
}

