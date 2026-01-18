//! SP1 witness file I/O helpers.
//!
//! Convention:
//! - `WITNESS_PATH` is a **bundle** that includes public inputs and the R1LF digest.

#![cfg(feature = "we_gate")]

use std::io::Read;

const SP1_WITNESS_BUNDLE_MAGIC: &[u8; 4] = b"SP1W";
const SP1_WITNESS_BUNDLE_VERSION: u32 = 1;
const SP1_WITNESS_BUNDLE_HEADER_LEN: usize = 112;

#[derive(Debug, Clone)]
pub struct Sp1WitnessBundle {
    pub witness: Vec<u64>,
    pub base_len: usize,
    pub aux_len: usize,
    pub public_inputs: ([u8; 32], [u8; 32]),
    pub r1lf_digest: [u8; 32],
}

/// Load a witness in either of these formats:
/// - **Single file only**: `witness_path` is the full witness of length `num_vars`.
///
/// Returns `(full_witness, base_len, aux_len)`.
pub fn load_sp1_witness_any(
    witness_path: &str,
    num_vars: usize,
) -> Result<Sp1WitnessBundle, String> {
    let mut bytes = Vec::new();
    std::fs::File::open(witness_path)
        .and_then(|mut f| f.read_to_end(&mut bytes))
        .map_err(|e| format!("read {witness_path}: {e}"))?;

    if bytes.len() >= SP1_WITNESS_BUNDLE_HEADER_LEN
        && bytes[0..4] == *SP1_WITNESS_BUNDLE_MAGIC
    {
        let version = u32::from_le_bytes(bytes[4..8].try_into().unwrap());
        if version != SP1_WITNESS_BUNDLE_VERSION {
            return Err(format!(
                "unsupported witness bundle version {version} (expected {SP1_WITNESS_BUNDLE_VERSION})"
            ));
        }
        let mut r1lf_digest = [0u8; 32];
        r1lf_digest.copy_from_slice(&bytes[8..40]);
        let bundle_num_vars = u64::from_le_bytes(bytes[40..48].try_into().unwrap()) as usize;
        let mut vk_hash = [0u8; 32];
        vk_hash.copy_from_slice(&bytes[48..80]);
        let mut committed_values_digest = [0u8; 32];
        committed_values_digest.copy_from_slice(&bytes[80..112]);
        let witness_len = bundle_num_vars;
        if witness_len != num_vars {
            return Err(format!(
                "witness bundle num_vars mismatch: file has {}, expected {}",
                witness_len, num_vars
            ));
        }
        let expected_bytes = SP1_WITNESS_BUNDLE_HEADER_LEN + witness_len * 8;
        if bytes.len() != expected_bytes {
            return Err(format!(
                "witness bundle byte length mismatch: got {} expected {}",
                bytes.len(),
                expected_bytes
            ));
        }
        let mut witness = Vec::with_capacity(witness_len);
        for chunk in bytes[SP1_WITNESS_BUNDLE_HEADER_LEN..].chunks_exact(8) {
            witness.push(u64::from_le_bytes(chunk.try_into().unwrap()));
        }
        if witness.is_empty() || witness[0] != 1 {
            return Err("witness must have w[0]=1 (constant ONE slot)".to_string());
        }
        return Ok(Sp1WitnessBundle {
            witness,
            base_len: num_vars,
            aux_len: 0,
            public_inputs: (vk_hash, committed_values_digest),
            r1lf_digest,
        });
    }

    Err(format!(
        "unrecognized witness format: expected SP1 witness bundle. \
Set SP1_WITNESS when exporting from SP1."
    ))
}

