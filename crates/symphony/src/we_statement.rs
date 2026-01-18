//! Helpers to bind WE locks to a Symphony statement.
//!
//! In Architecture‑T (statement-only arming), the armer derives per-lock queries from a
//! **statement hash** that binds:
//! - the SP1 verifier key / program id (vk hash),
//! - the SP1 R1CS digest,
//! - the WE gate digest (exact gate relation version),
//! - and the public inputs.
//!
//! This module provides a canonical `sha256` hash for the **public** portion of `R_WE`.

use ark_ff::{BigInteger, Field, PrimeField};
use sha2::{Digest, Sha256};
use stark_rings::OverField;

/// SHA256 hash of the public statement for `R_WE` (hetero‑M).
///
/// This intentionally does **not** include the folding proof (`PiFoldBatchedProof`) nor auxiliary
/// witness messages (`aux`) nor reduced witness (`ro_witness`): those are witnesses to the relation.
///
/// Callers should include **only statement-defined values**:
/// - `vk_hash`: a 32-byte digest identifying the SP1 shrink verifier circuit/version.
/// - `r1cs_digest`: a 32-byte digest identifying the SP1 R1CS being folded (from SP1 R1CS header).
/// - `gate_digest`: a 32-byte digest identifying the **exact** WE gate relation `R_WE`
///   (i.e., the dR1CS constraint system / verifier circuit shape the DPP targets).
/// - `public_inputs`: the public inputs absorbed by the Poseidon transcript.
///
/// NOTE: This intentionally does NOT include prover-chosen proof artifacts (e.g. `cm_f`, `cfs_*`),
/// because locks must be armable without observing a particular proving run. Those objects are
/// checked in-gate as part of the witness to `R_WE`, but are not statement-defining here.
pub fn we_statement_hash_hetero_m<R: OverField>(
    vk_hash: [u8; 32],
    r1cs_digest: [u8; 32],
    gate_digest: [u8; 32],
    public_inputs: &[R::BaseRing],
) -> [u8; 32]
where
    R::BaseRing: Field,
{
    let mut h = Sha256::new();
    h.update(b"SYMPHONY_WE_STATEMENT_V1");
    h.update(&vk_hash);
    h.update(&r1cs_digest);
    h.update(&gate_digest);

    // Public inputs (base field elements).
    h.update(&(public_inputs.len() as u64).to_le_bytes());
    for x in public_inputs {
        absorb_field_elem::<R>(x, &mut h);
    }

    h.finalize().into()
}

fn absorb_field_elem<R: OverField>(x: &R::BaseRing, h: &mut Sha256)
where
    R::BaseRing: Field,
{
    for fp in x.to_base_prime_field_elements() {
        h.update(fp.into_bigint().to_bytes_le());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use stark_rings::{cyclotomic_ring::models::frog_ring::RqPoly as R, PolyRing};

    #[test]
    fn test_statement_hash_deterministic() {
        let vk = [1u8; 32];
        let r1cs = [2u8; 32];
        let gate = [3u8; 32];
        let public_inputs: Vec<<R as PolyRing>::BaseRing> = vec![<R as PolyRing>::BaseRing::from(5u128)];

        let a = we_statement_hash_hetero_m::<R>(vk, r1cs, gate, &public_inputs);
        let b = we_statement_hash_hetero_m::<R>(vk, r1cs, gate, &public_inputs);
        assert_eq!(a, b);
    }

    #[test]
    fn test_statement_hash_changes_on_public_inputs() {
        let vk = [1u8; 32];
        let r1cs = [2u8; 32];
        let gate = [3u8; 32];

        let p0: Vec<<R as PolyRing>::BaseRing> = vec![<R as PolyRing>::BaseRing::from(5u128)];
        let p1: Vec<<R as PolyRing>::BaseRing> = vec![<R as PolyRing>::BaseRing::from(6u128)];
        let a = we_statement_hash_hetero_m::<R>(vk, r1cs, gate, &p0);
        let b = we_statement_hash_hetero_m::<R>(vk, r1cs, gate, &p1);
        assert_ne!(a, b);
    }

    #[test]
    fn test_statement_hash_changes_on_gate_digest() {
        let vk = [1u8; 32];
        let r1cs = [2u8; 32];
        let gate0 = [3u8; 32];
        let gate1 = [4u8; 32];
        let public_inputs: Vec<<R as PolyRing>::BaseRing> = vec![<R as PolyRing>::BaseRing::from(5u128)];

        let a = we_statement_hash_hetero_m::<R>(vk, r1cs, gate0, &public_inputs);
        let b = we_statement_hash_hetero_m::<R>(vk, r1cs, gate1, &public_inputs);
        assert_ne!(a, b);
    }
}

