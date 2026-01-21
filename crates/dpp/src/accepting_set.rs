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

    // Fast path for k=3 and small numeric domains.
    //
    // The generic enumerator below is BigInt-heavy (nested loops doing BigInt increments and
    // conversions) and becomes slow even around ~1e5–1e6 points. For tiny-field demo tests we
    // want this to run quickly, so when:
    // - k == 3,
    // - bounds and weights fit in i64/i128,
    // - and the predicate is MulEqModP over a small modulus,
    // we enumerate in i64/i128 and only convert accepted tuples into field elements.
    if query.w.len() == 3 {
        if let FlpcpPredicate::MulEqModP { p_small } = &query.pred {
            if let (Some(p_small_u128), Some(p_u128)) = (p_small.to_u128(), modulus_u128::<F>()) {
                if p_small_u128 > 1 && p_u128 > 1 {
                    // Precompute p (mod p_small) for handling negative centered reps.
                    let p_mod_ps = (p_u128 % p_small_u128) as u128;

                    // Convert bounds b_i-1 into i64.
                    let mut bound: [i64; 3] = [0; 3];
                    for i in 0..3 {
                        let bi = &query.b[i];
                        let b1 = (bi - BigInt::one()).to_i64();
                        if b1.is_none() {
                            bound = [0; 3];
                            break;
                        }
                        bound[i] = b1.unwrap();
                    }
                    if bound[0] > 0 && bound[1] > 0 && bound[2] > 0 {
                        // Convert weights w_i into i128 (safe for toy params).
                        let mut w: [i128; 3] = [0; 3];
                        for i in 0..3 {
                            let wi = query.w[i].to_i128();
                            if wi.is_none() {
                                w = [0; 3];
                                break;
                            }
                            w[i] = wi.unwrap();
                        }
                        if w[0] != 0 && w[1] != 0 && w[2] != 0 {
                            #[inline]
                            fn mod_small(z: i64, p_mod_ps: u128, ps: u128) -> u128 {
                                if z >= 0 {
                                    (z as u128) % ps
                                } else {
                                    let t = ((-z) as u128) % ps;
                                    if t == 0 {
                                        0
                                    } else {
                                        (p_mod_ps + ps - t) % ps
                                    }
                                }
                            }

                            #[inline]
                            fn mod_p_u128(z: i128, p: u128) -> u128 {
                                // z mod p in [0,p)
                                let p_i = p as i128;
                                let mut r = z % p_i;
                                if r < 0 {
                                    r += p_i;
                                }
                                r as u128
                            }

                            // Store residues in u128 for fast sort/dedup, then map to F at end.
                            let mut out_u128: Vec<u128> = Vec::new();
                            let b0 = bound[0];
                            let b1 = bound[1];
                            let b2 = bound[2];
                            for t0 in -b0..=b0 {
                                let a0 = mod_small(t0, p_mod_ps, p_small_u128);
                                for t1 in -b1..=b1 {
                                    let a1 = mod_small(t1, p_mod_ps, p_small_u128);
                                    let a01 = (a0 * a1) % p_small_u128;
                                    for t2 in -b2..=b2 {
                                        let a2 = mod_small(t2, p_mod_ps, p_small_u128);
                                        if (a01 + p_small_u128 - a2) % p_small_u128 == 0 {
                                            // Accept: compute packed a = <t, w> as integer, then reduce mod p.
                                            let a_int = (w[0] * (t0 as i128))
                                                + (w[1] * (t1 as i128))
                                                + (w[2] * (t2 as i128));
                                            let a_mod = mod_p_u128(a_int, p_u128);
                                            out_u128.push(a_mod);
                                        }
                                    }
                                }
                            }

                            out_u128.sort_unstable();
                            out_u128.dedup();
                            let out = out_u128
                                .into_iter()
                                .map(|a| F::from_le_bytes_mod_order(&a.to_le_bytes()))
                                .collect::<Vec<_>>();
                            return Ok(out);
                        }
                    }
                }
            }
        }
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

fn modulus_u128<F: PrimeField>() -> Option<u128> {
    // Parse modulus as u128 if it fits (true for Goldilocks / Fp64).
    let p = BigInt::from_bytes_le(num_bigint::Sign::Plus, &F::MODULUS.to_bytes_le());
    p.to_u128()
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

