//! Explicit accepting-set enumeration for packed DPP queries (only for small domains).
//!
//! This module exists to support **arm-before-proof** locks in the GPT‑PRO interface:
//! the lock layer wants `(q, A)` where `A` is an explicit (small) accepting set such that:
//!   accept ⇔ <q, (x||π)> ∈ A.
//!
//! Our current packed verifier is predicate-style:
//!   accept ⇔ pred(Decode(<q,(x||π)>; w,b)) = 1.
//!
//! In general the induced accepting set is huge. Here we provide a *guarded* enumerator that
//! only works when the decoded-answer domain is small enough to brute-force.

use ark_ff::{BigInteger, PrimeField};
use num_bigint::BigInt;
use num_traits::{One, Signed, ToPrimitive, Zero};
use thiserror::Error;

use crate::packing::{centered_bigint_to_field, FlpcpPredicate, PackedDppQuery, PackedDppQuerySparse};

#[derive(Debug, Error)]
pub enum AcceptingSetError {
    #[error("invalid parameters")]
    InvalidParams,
    #[error("decoded-answer domain too large to enumerate (estimated {estimated} > limit {limit})")]
    DomainTooLarge { estimated: u128, limit: u128 },
}

/// Estimate the number of centered integer tuples `ans` satisfying `|ans_i| <= b_i-1`.
fn estimate_domain_size(b: &[BigInt]) -> Option<u128> {
    let mut acc: u128 = 1;
    for bi in b {
        // range size = 2*(b_i-1)+1
        let bound = bi - BigInt::one();
        if bound.is_negative() {
            return None;
        }
        let r = (bound * BigInt::from(2u64)) + BigInt::one();
        let r_u = r.to_u128()?;
        acc = acc.checked_mul(r_u)?;
    }
    Some(acc)
}

/// Enumerate the induced accepting set
///   A = { <ans, w> : ans in [-b_i+1 .. b_i-1], pred(ans)=true } ⊂ F
/// for a dense packed query.
///
/// This is **only feasible** when the decoded-answer domain is tiny.
pub fn accepting_set_for_packed_query<F: PrimeField>(
    query: &PackedDppQuery<F>,
    limit: u128,
) -> Result<Vec<F>, AcceptingSetError> {
    let k = query.w.len();
    if k == 0 || query.b.len() != k {
        return Err(AcceptingSetError::InvalidParams);
    }
    let est = estimate_domain_size(&query.b).ok_or(AcceptingSetError::InvalidParams)?;
    if est > limit {
        return Err(AcceptingSetError::DomainTooLarge { estimated: est, limit });
    }

    // Enumerate centered integer answers within bounds.
    let mut cur: Vec<BigInt> = vec![BigInt::zero(); k];
    let mut out: Vec<F> = Vec::new();

    fn rec<F: PrimeField>(
        i: usize,
        cur: &mut [BigInt],
        w: &[BigInt],
        b: &[BigInt],
        pred: &FlpcpPredicate<F>,
        out: &mut Vec<F>,
    ) {
        let k = b.len();
        if i == k {
            // Convert to field for predicate check.
            let ans_field: Vec<F> = cur.iter().map(|z| centered_bigint_to_field::<F>(z)).collect();
            if pred.check(&ans_field) {
                // Compute packed answer a = <cur, w> (integer), then map into field.
                let mut a_int = BigInt::zero();
                for j in 0..k {
                    a_int += &cur[j] * &w[j];
                }
                out.push(centered_bigint_to_field::<F>(&a_int));
            }
            return;
        }

        let bound = &b[i] - BigInt::one();
        // iterate t ∈ [-bound .. bound]
        let mut t = -bound.clone();
        while t <= bound {
            cur[i] = t.clone();
            rec::<F>(i + 1, cur, w, b, pred, out);
            t += BigInt::one();
        }
    }

    rec::<F>(0, &mut cur, &query.w, &query.b, &query.pred, &mut out);

    // Dedup and return in canonical order (by field integer rep).
    //
    // Note: this is only for small sets; O(|A| log |A|) is fine.
    let p = BigInt::from_bytes_le(num_bigint::Sign::Plus, &F::MODULUS.to_bytes_le());
    out.sort_by(|a, b| {
        let aa = BigInt::from_bytes_le(num_bigint::Sign::Plus, &a.into_bigint().to_bytes_le()) % &p;
        let bb = BigInt::from_bytes_le(num_bigint::Sign::Plus, &b.into_bigint().to_bytes_le()) % &p;
        aa.cmp(&bb)
    });
    out.dedup();
    Ok(out)
}

/// Sparse packed query variant.
pub fn accepting_set_for_packed_query_sparse<F: PrimeField>(
    query: &PackedDppQuerySparse<F>,
    limit: u128,
) -> Result<Vec<F>, AcceptingSetError> {
    // Same logic; only metadata differs.
    let dense = PackedDppQuery::<F> {
        q: vec![], // unused
        w: query.w.clone(),
        b: query.b.clone(),
        pred: query.pred.clone(),
    };
    accepting_set_for_packed_query::<F>(&dense, limit)
}

