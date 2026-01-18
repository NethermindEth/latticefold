//! Boolean proof encoding wrapper (Claim 5.3).
//!
//! The Rev2 paper introduces a “binary decomposition optimization” (see Section 1.3),
//! formalized as Claim 5.3: any FLPCP can be transformed to have **Boolean proofs**
//! by encoding each proof field element in binary.
//!
//! This tightens naive integer bounds: without Booleanization, integer representatives can involve
//! products of field elements (≈ (p-1)^2); with Boolean proofs, contributions become linear in p.

use ark_ff::{BigInteger, PrimeField};
use num_bigint::BigInt;
use num_traits::One;
use rand::RngCore;
use rayon::prelude::*;

use crate::packing::{BoundedFlpcp, BoundedFlpcpSparse, FlpcpPredicate};
use crate::sparse::SparseVec;

/// Wrapper that replaces an m-element proof over F_p with an (m * bit_len)-bit proof in {0,1}.
#[derive(Clone, Debug)]
pub struct BooleanProofFlpcp<F: PrimeField, V: BoundedFlpcp<F>> {
    pub inner: V,
    pub bit_len: usize,
    _pd: core::marker::PhantomData<F>,
}

impl<F: PrimeField, V: BoundedFlpcp<F>> BooleanProofFlpcp<F, V> {
    pub fn new(inner: V) -> Self {
        // Use modulus bit size as the fixed-width binary encoding length.
        let bit_len = F::MODULUS_BIT_SIZE as usize;
        Self { inner, bit_len, _pd: core::marker::PhantomData }
    }

    /// Compute the Boolean proof length m' = m * bit_len.
    pub fn m_bits(&self) -> usize {
        self.inner.m() * self.bit_len
    }

    /// Encode a proof vector into little-endian bits (length m*bit_len).
    pub fn encode_proof_bits(&self, pi: &[F]) -> Vec<F> {
        assert_eq!(pi.len(), self.inner.m());
        let mut out = Vec::with_capacity(self.m_bits());
        for x in pi {
            let bits = x.into_bigint().to_bits_le();
            // Truncate/pad to bit_len.
            for i in 0..self.bit_len {
                let b = bits.get(i).copied().unwrap_or(false);
                out.push(if b { F::ONE } else { F::ZERO });
            }
        }
        out
    }
}

impl<F: PrimeField, V: BoundedFlpcp<F>> BoundedFlpcp<F> for BooleanProofFlpcp<F, V> {
    fn n(&self) -> usize {
        self.inner.n()
    }

    fn m(&self) -> usize {
        self.m_bits()
    }

    fn k(&self) -> usize {
        self.inner.k()
    }

    fn bounds_b(&self) -> Vec<BigInt> {
        // Boolean proof improves boundedness only for the proof part: bits are in {0,1}.
        //
        // For a generic query vector `q` with centered reps bounded by (p-1)/2:
        // - x-part contribution per coordinate is bounded by ((p-1)/2)^2
        // - proof-bit contribution per coordinate is bounded by (p-1)/2
        //
        // So a safe bound is:
        //   b = n * ((p-1)/2)^2 + m_bits * ((p-1)/2)
        let p = BigInt::from_bytes_le(num_bigint::Sign::Plus, &F::MODULUS.to_bytes_le());
        let half = (&p - BigInt::one()) / BigInt::from(2u64);
        let n = BigInt::from(self.n() as u64);
        let m_bits = BigInt::from(self.m() as u64);
        let b = n * (&half * &half) + m_bits * half;
        vec![b; self.k()]
    }

    fn sample_queries_and_predicate(
        &self,
        rng: &mut dyn RngCore,
        x: &[F],
    ) -> Result<(Vec<Vec<F>>, FlpcpPredicate<F>), String> {
        // Sample queries for the original (x || pi_field) space.
        let (q_rows, pred) = self.inner.sample_queries_and_predicate(rng, x)?;
        let n = self.inner.n();
        let m = self.inner.m();
        let bit_len = self.bit_len;
        let mut out_rows = Vec::with_capacity(q_rows.len());

        // For each query row, expand proof-part coefficients using powers of 2:
        // q_pi[j] * pi[j]  ==  Σ_t (q_pi[j] * 2^t) * bit_{j,t}.
        let mut pow2 = Vec::with_capacity(bit_len);
        let mut cur = F::ONE;
        for _ in 0..bit_len {
            pow2.push(cur);
            cur = cur + cur;
        }

        for row in q_rows {
            if row.len() != n + m {
                return Err("BooleanProofFlpcp: inner query length mismatch".to_string());
            }
            let mut new_row = Vec::with_capacity(n + m * bit_len);
            // x-part unchanged
            new_row.extend_from_slice(&row[..n]);
            // proof-part expanded
            for j in 0..m {
                let coeff = row[n + j];
                for t in 0..bit_len {
                    new_row.push(coeff * pow2[t]);
                }
            }
            out_rows.push(new_row);
        }

        Ok((out_rows, pred))
    }
}

/// Sparse-query version of the Boolean proof wrapper (Claim 5.3).
#[derive(Clone, Debug)]
pub struct BooleanProofFlpcpSparse<F: PrimeField, V: BoundedFlpcpSparse<F>> {
    pub inner: V,
    pub bit_len: usize,
    _pd: core::marker::PhantomData<F>,
}

impl<F: PrimeField, V: BoundedFlpcpSparse<F>> BooleanProofFlpcpSparse<F, V> {
    pub fn new(inner: V) -> Self {
        let bit_len = F::MODULUS_BIT_SIZE as usize;
        Self { inner, bit_len, _pd: core::marker::PhantomData }
    }

    pub fn m_bits(&self) -> usize {
        self.inner.m() * self.bit_len
    }

    pub fn encode_proof_bits(&self, pi: &[F]) -> Vec<F> {
        assert_eq!(pi.len(), self.inner.m());
        let mut out = Vec::with_capacity(self.m_bits());
        for x in pi {
            let bits = x.into_bigint().to_bits_le();
            for i in 0..self.bit_len {
                let b = bits.get(i).copied().unwrap_or(false);
                out.push(if b { F::ONE } else { F::ZERO });
            }
        }
        out
    }

    /// Encode proof into a **bitpacked** byte array (little-endian bits within each proof element).
    ///
    /// This is a memory-critical optimization for large instances: the unpacked `Vec<F>` form uses
    /// one full field element per bit, which explodes RAM in benchmarks.
    pub fn encode_proof_bits_packed(&self, pi: &[F]) -> Vec<u8> {
        assert_eq!(pi.len(), self.inner.m());
        // Fast path: bit_len is a multiple of 8 (true for BF=64-bit, FLarge=384-bit, etc).
        // Then bitpacking is just a fixed-width little-endian byte encoding per element.
        if self.bit_len % 8 == 0 {
            let bytes_per = self.bit_len / 8;
            let mut out = vec![0u8; pi.len() * bytes_per];
            out.par_chunks_mut(bytes_per)
                .zip(pi.par_iter())
                .for_each(|(chunk, x)| {
                    // Ensure fixed-width output.
                    chunk.fill(0);
                    // Avoid per-element heap allocation from `to_bytes_le()` by writing limbs directly.
                    let xb = x.into_bigint();
                    let limbs = xb.as_ref(); // little-endian u64 limbs
                    let mut off = 0usize;
                    for limb in limbs {
                        if off >= bytes_per {
                            break;
                        }
                        let limb_bytes = limb.to_le_bytes();
                        let take = core::cmp::min(bytes_per - off, 8);
                        chunk[off..off + take].copy_from_slice(&limb_bytes[..take]);
                        off += take;
                    }
                });
            return out;
        }

        // Generic slow path (non-byte-aligned bit_len): keep for completeness.
        let nbits = self.m_bits();
        let nbytes = (nbits + 7) / 8;
        let mut out = vec![0u8; nbytes];
        let mut bit_pos = 0usize;
        for x in pi {
            let bytes = x.into_bigint().to_bytes_le();
            for i in 0..self.bit_len {
                let b = (bytes.get(i / 8).copied().unwrap_or(0) >> (i % 8)) & 1;
                if b == 1 {
                    out[bit_pos / 8] |= 1u8 << (bit_pos % 8);
                }
                bit_pos += 1;
            }
        }
        debug_assert_eq!(bit_pos, nbits);
        out
    }
}

impl<F: PrimeField, V: BoundedFlpcpSparse<F>> BoundedFlpcpSparse<F> for BooleanProofFlpcpSparse<F, V> {
    fn n(&self) -> usize {
        self.inner.n()
    }

    fn m(&self) -> usize {
        self.m_bits()
    }

    fn k(&self) -> usize {
        self.inner.k()
    }

    fn bounds_b(&self) -> Vec<BigInt> {
        // Same bound logic as dense wrapper.
        let p = BigInt::from_bytes_le(num_bigint::Sign::Plus, &F::MODULUS.to_bytes_le());
        let half = (&p - BigInt::one()) / BigInt::from(2u64);
        let n = BigInt::from(self.n() as u64);
        let m_bits = BigInt::from(self.m() as u64);
        let b = n * (&half * &half) + m_bits * half;
        vec![b; self.k()]
    }

    fn sample_queries_and_predicate_sparse(
        &self,
        rng: &mut dyn RngCore,
        x: &[F],
    ) -> Result<(Vec<SparseVec<F>>, FlpcpPredicate<F>), String> {
        let (q_rows, pred) = self.inner.sample_queries_and_predicate_sparse(rng, x)?;
        let n = self.inner.n();
        let m = self.inner.m();
        let bit_len = self.bit_len;

        let mut pow2 = Vec::with_capacity(bit_len);
        let mut cur = F::ONE;
        for _ in 0..bit_len {
            pow2.push(cur);
            cur = cur + cur;
        }

        let mut out_rows = Vec::with_capacity(q_rows.len());
        for row in q_rows {
            let mut new_terms: Vec<(F, usize)> = Vec::new();
            for (coeff, idx) in row.terms {
                if idx < n {
                    new_terms.push((coeff, idx));
                } else {
                    let j = idx - n;
                    if j >= m {
                        return Err("BooleanProofFlpcpSparse: inner index out of range".to_string());
                    }
                    for t in 0..bit_len {
                        let c = coeff * pow2[t];
                        if !c.is_zero() {
                            new_terms.push((c, n + j * bit_len + t));
                        }
                    }
                }
            }
            out_rows.push(SparseVec::new(new_terms));
        }

        Ok((out_rows, pred))
    }
}

