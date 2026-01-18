//! FLPCP for deterministic R1CS (dR1CS) — prototype (Section 4.1).
//!
//! The paper constructs efficient 3-query FLPCPs for dR1CS using multiplication codes.
//! This module will implement the RS-based instantiation (Corollary 4.8) and expose it
//! through the `BoundedFlpcp` interface so it can be packed into a DPP (Section 5.2).

use ark_ff::{BigInteger, Field, FftField, PrimeField};
use num_bigint::BigInt;
use num_traits::One;
use rand::RngCore;

use rayon::prelude::*;
use rayon::join;

use crate::packing::{BoundedFlpcp, BoundedFlpcpSparse, FlpcpPredicate};
use crate::rs::{barycentric_weights_consecutive, extrapolate_consecutive_next_block, lagrange_coeffs_at};
use crate::sparse::SparseVec;

/// Dense dR1CS instance: check (A x) ⊙ (B x) == (C x).
#[derive(Clone, Debug)]
pub struct Dr1csInstance<F: Field> {
    pub a: Vec<Vec<F>>, // k x n
    pub b: Vec<Vec<F>>, // k x n
    pub c: Vec<Vec<F>>, // k x n
}

impl<F: Field> Dr1csInstance<F> {
    pub fn k(&self) -> usize {
        self.a.len()
    }
    pub fn n(&self) -> usize {
        self.a.first().map(|r| r.len()).unwrap_or(0)
    }
}

/// RS-based 3-query FLPCP for dR1CS (Theorem 4.6 / Corollary 4.8).
///
/// Proof length is `m = 2k` field elements (the systematic prefix of the square code).
#[derive(Clone, Debug)]
pub struct RsDr1csFlpcp<F: PrimeField + FftField> {
    pub inst: Dr1csInstance<F>,
    /// Codeword length ℓ (must satisfy ℓ >= 2k and ℓ <= |F|).
    pub ell: usize,
    /// Evaluation points α_0..α_{ℓ-1}.
    pub points: Vec<F>,
    ws_k: Vec<F>,
    ws_2k: Vec<F>,
}

impl<F: PrimeField + FftField> RsDr1csFlpcp<F> {
    pub fn new(inst: Dr1csInstance<F>, ell: usize) -> Self {
        let k = inst.k();
        assert!(k > 0);
        assert!(ell >= 2 * k);
        // Simple deterministic points: 1,2,...,ell
        let points = (0..ell).map(|i| F::from((i as u64) + 1)).collect::<Vec<_>>();
        // O(k) weights for consecutive points.
        let ws_k = barycentric_weights_consecutive::<F>(k, 1);
        let ws_2k = barycentric_weights_consecutive::<F>(2 * k, 1);
        Self { inst, ell, points, ws_k, ws_2k }
    }

    /// Prover for the RS multiplication-code FLPCP:
    /// computes `w[i] = E(Ax)[i] * E(Bx)[i]` for i in [0..2k).
    pub fn prove(&self, x: &[F]) -> Vec<F> {
        let k = self.inst.k();
        assert_eq!(x.len(), self.inst.n());

        // Compute yA = A x, yB = B x (length k). Independent: run concurrently.
        let (y_a, y_b) = join(|| mat_vec(&self.inst.a, x), || mat_vec(&self.inst.b, x));

        // Fast systematic RS extrapolation on consecutive points:
        // - we already have f(1..k) as y_*
        // - compute f(k+1..2k) in O(k log k) via one convolution
        //
        // Memory note:
        // Extrapolation uses large NTT buffers when k is huge; running multiple extrapolations
        // concurrently can blow up RSS. For large sizes, compute tails sequentially.
        // Two extrapolations is usually an acceptable memory multiplier; keep them parallel.
        let (y_a_tail, y_b_tail) = join(
            || extrapolate_consecutive_next_block::<F>(&y_a),
            || extrapolate_consecutive_next_block::<F>(&y_b),
        );
        debug_assert_eq!(y_a_tail.len(), k);
        debug_assert_eq!(y_b_tail.len(), k);

        let mut w = Vec::with_capacity(2 * k);
        for i in 0..k {
            w.push(y_a[i] * y_b[i]);
        }
        for i in 0..k {
            w.push(y_a_tail[i] * y_b_tail[i]);
        }
        w
    }
}

impl<F: PrimeField + FftField> BoundedFlpcp<F> for RsDr1csFlpcp<F> {
    fn n(&self) -> usize {
        self.inst.n()
    }

    fn m(&self) -> usize {
        2 * self.inst.k()
    }

    fn k(&self) -> usize {
        3
    }

    fn bounds_b(&self) -> Vec<BigInt> {
        // Unbounded FLPCP: return the trivial bound b_i = (n+m) * ((p-1)/2)^2 as in Def 5.2 note.
        // This makes the packing modulus condition fail unless the field is enormous; the bounded
        // embedding (Theorem 5.6) is the intended way to get small b.
        let len = self.n() + self.m();
        let p = BigInt::from_bytes_le(num_bigint::Sign::Plus, &F::MODULUS.to_bytes_le());
        let half = (&p - BigInt::one()) / BigInt::from(2u64);
        let b = BigInt::from(len as u64) * &half * &half;
        vec![b.clone(), b.clone(), b]
    }

    fn sample_queries_and_predicate(
        &self,
        rng: &mut dyn RngCore,
        x: &[F],
    ) -> Result<(Vec<Vec<F>>, FlpcpPredicate<F>), String> {
        let k = self.inst.k();
        if x.len() != self.inst.n() || k == 0 {
            return Err("bad input".to_string());
        }
        // Prover's claimed proof π is w (length 2k).
        // Query schedule (Theorem 4.6):
        // Pick random i in [ℓ], random λ in F; query:
        //  a = E(Ax)[i]
        //  b = E(Bx)[i]
        //  c = E*(w)[i] + λ * E*(Cx - w≤k || 0)[i]
        // Accept if a*b == c.

        let idx = (rng.next_u64() as usize) % self.ell;
        let lambda = F::from(rng.next_u64());
        let alpha = self.points[idx];

        // Lagrange coeffs over first k points (for E(Ax), E(Bx), E(Cx)).
        let lam_k = lagrange_coeffs_at(&self.points[..k], &self.ws_k, alpha);
        // Lagrange coeffs over first 2k points (for E*(w), and for E*(Cx-w||0)).
        let lam_2k = lagrange_coeffs_at(&self.points[..2 * k], &self.ws_2k, alpha);

        // Build query vectors over v = (x || w) of length n + 2k.
        // q_a: coefficients on x to compute dot(lam_k, A x)
        let q_a = lin_combo_rows(&self.inst.a, &lam_k);
        let q_b = lin_combo_rows(&self.inst.b, &lam_k);
        let _q_cx = lin_combo_rows(&self.inst.c, &lam_k); // gives dot(lam_k, Cx)

        let n = self.inst.n();
        let mut q1 = vec![F::zero(); n + 2 * k];
        let mut q2 = vec![F::zero(); n + 2 * k];
        let mut q3 = vec![F::zero(); n + 2 * k];

        // q1 = (q_a || 0)
        q1[..n].copy_from_slice(&q_a);
        // q2 = (q_b || 0)
        q2[..n].copy_from_slice(&q_b);

        // q3 = ( ? || ? )
        // Part 1: E*(w)[i] is linear in w with coeffs lam_2k.
        // Part 2: E*(Cx - w≤k || 0)[i] = dot(lam_2k[..k], (Cx - w_prefix)) since last k entries are 0.
        // => equals dot(lam_2k[..k], Cx) - dot(lam_2k[..k], w_prefix).
        // Combine: c = dot(lam_2k, w) + λ*( dot(lam_2k[..k], Cx) - dot(lam_2k[..k], w_prefix) )
        // = λ * dot(lam_2k[..k], Cx)  + dot(lam_2k, w) - λ * dot(lam_2k[..k], w_prefix)
        // So w coeffs: for j<k: lam_2k[j] - λ*lam_2k[j]; for j>=k: lam_2k[j].
        // x coeffs: λ * (lin_combo_rows(C, lam_2k[..k])).
        // Compute lin combo using lam_2k[..k] over C rows.
        let q_cx2 = lin_combo_rows(&self.inst.c, &lam_2k[..k]);

        for j in 0..n {
            q3[j] = lambda * q_cx2[j];
        }
        for j in 0..2 * k {
            let coeff = if j < k {
                lam_2k[j] - (lambda * lam_2k[j])
            } else {
                lam_2k[j]
            };
            q3[n + j] = coeff;
        }

        Ok((vec![q1, q2, q3], FlpcpPredicate::MulEq))
    }
}

/// Sparse dR1CS instance: check (A x) ⊙ (B x) == (C x) with sparse rows.
#[derive(Clone, Debug)]
pub struct Dr1csInstanceSparse<F: PrimeField> {
    pub n: usize,
    pub a: Vec<SparseVec<F>>, // k rows
    pub b: Vec<SparseVec<F>>, // k rows
    pub c: Vec<SparseVec<F>>, // k rows
}

impl<F: PrimeField> Dr1csInstanceSparse<F> {
    pub fn k(&self) -> usize {
        self.a.len()
    }
}

/// RS-based 3-query FLPCP for sparse dR1CS rows.
///
/// Same protocol as `RsDr1csFlpcp`, but query vectors are emitted as `SparseVec` over v=(x||w).
#[derive(Clone, Debug)]
pub struct RsDr1csFlpcpSparse<F: PrimeField + FftField> {
    pub inst: Dr1csInstanceSparse<F>,
    pub ell: usize,
    pub points: Vec<F>,
    ws_k: Vec<F>,
    ws_2k: Vec<F>,
}

impl<F: PrimeField + FftField> RsDr1csFlpcpSparse<F> {
    pub fn new(inst: Dr1csInstanceSparse<F>, ell: usize) -> Self {
        let k = inst.k();
        assert!(k > 0);
        assert!(ell >= 2 * k);
        let points = (0..ell).map(|i| F::from((i as u64) + 1)).collect::<Vec<_>>();
        let ws_k = barycentric_weights_consecutive::<F>(k, 1);
        let ws_2k = barycentric_weights_consecutive::<F>(2 * k, 1);
        Self { inst, ell, points, ws_k, ws_2k }
    }

    pub fn prove(&self, x: &[F]) -> Vec<F> {
        let k = self.inst.k();
        assert_eq!(x.len(), self.inst.n);
        // A and B mat-vecs are independent: run concurrently.
        let (y_a, y_b) = join(
            || mat_vec_sparse(&self.inst.a, x),
            || mat_vec_sparse(&self.inst.b, x),
        );
        // Extrapolations are independent and often dominate: run concurrently.
        let (y_a_tail, y_b_tail) = join(
            || extrapolate_consecutive_next_block::<F>(&y_a),
            || extrapolate_consecutive_next_block::<F>(&y_b),
        );
        let mut w = Vec::with_capacity(2 * k);
        for i in 0..k {
            w.push(y_a[i] * y_b[i]);
        }
        for i in 0..k {
            w.push(y_a_tail[i] * y_b_tail[i]);
        }
        w
    }
}

impl<F: PrimeField + FftField> BoundedFlpcpSparse<F> for RsDr1csFlpcpSparse<F> {
    fn n(&self) -> usize {
        self.inst.n
    }

    fn m(&self) -> usize {
        2 * self.inst.k()
    }

    fn k(&self) -> usize {
        3
    }

    fn bounds_b(&self) -> Vec<BigInt> {
        // Same unbounded bound as dense prototype.
        let len = self.n() + self.m();
        let p = BigInt::from_bytes_le(num_bigint::Sign::Plus, &F::MODULUS.to_bytes_le());
        let half = (&p - BigInt::one()) / BigInt::from(2u64);
        let b = BigInt::from(len as u64) * &half * &half;
        vec![b.clone(), b.clone(), b]
    }

    fn sample_queries_and_predicate_sparse(
        &self,
        rng: &mut dyn RngCore,
        x: &[F],
    ) -> Result<(Vec<SparseVec<F>>, FlpcpPredicate<F>), String> {
        let k = self.inst.k();
        let n = self.inst.n;
        if x.len() != n || k == 0 {
            return Err("bad input".to_string());
        }

        let idx = (rng.next_u64() as usize) % self.ell;
        let lambda = F::from(rng.next_u64());
        let alpha = self.points[idx];

        let lam_k = lagrange_coeffs_at(&self.points[..k], &self.ws_k, alpha);
        let lam_2k = lagrange_coeffs_at(&self.points[..2 * k], &self.ws_2k, alpha);

        let q_a = lin_combo_rows_sparse(&self.inst.a, n, &lam_k);
        let q_b = lin_combo_rows_sparse(&self.inst.b, n, &lam_k);
        let q_cx2 = lin_combo_rows_sparse(&self.inst.c, n, &lam_2k[..k]);

        // q1 = (q_a || 0), q2 = (q_b || 0)
        let q1 = q_a;
        let q2 = q_b;

        // q3 has x-part and w-part (shifted by n).
        let mut q3_terms: Vec<(F, usize)> = Vec::new();
        for (c, idx) in q_cx2.terms.iter() {
            let cc = lambda * *c;
            if !cc.is_zero() {
                q3_terms.push((cc, *idx));
            }
        }
        for j in 0..2 * k {
            let coeff = if j < k {
                lam_2k[j] - (lambda * lam_2k[j])
            } else {
                lam_2k[j]
            };
            if !coeff.is_zero() {
                q3_terms.push((coeff, n + j));
            }
        }

        Ok((vec![q1, q2, SparseVec::new(q3_terms)], FlpcpPredicate::MulEq))
    }
}

/// RS-based 3-query FLPCP for NP-style dR1CS, where the witness vector is part of the proof.
///
/// - Public input `x` has length `l` (can be 0).
/// - Private witness `z_w` has length `n_total - l`.
/// - The FLPCP proof is `(z_w || w)`, where `w` is the systematic square-code prefix (length 2k).
///
/// This matches the WE/DPP use-case: the statement is public, and the witness is private.
#[derive(Clone, Debug)]
pub struct RsDr1csNpFlpcpSparse<F: PrimeField + FftField> {
    pub inst: Dr1csInstanceSparse<F>,
    /// Number of public variables in `z` (prefix length).
    pub l: usize,
    pub ell: usize,
    pub points: Vec<F>,
    ws_k: Vec<F>,
    ws_2k: Vec<F>,
}

/// Auxiliary RS codeword data for coin-form query evaluation.
///
/// This is NOT part of the published proof; it is prover-side cached data useful for answering
/// many locks (queries) efficiently after proving.
#[derive(Clone, Debug)]
pub struct RsNpCodewords<F: PrimeField> {
    pub y_a: Vec<F>,      // values at points 1..k
    pub y_a_tail: Vec<F>, // values at points k+1..2k
    pub y_b: Vec<F>,
    pub y_b_tail: Vec<F>,
    pub y_c: Vec<F>,
    pub w: Vec<F>, // length 2k: w[i]=y_a_full[i]*y_b_full[i]
}

impl<F: PrimeField + FftField> RsDr1csNpFlpcpSparse<F> {
    pub fn new(inst: Dr1csInstanceSparse<F>, l: usize, ell: usize) -> Self {
        let k = inst.k();
        assert!(k > 0);
        assert!(l <= inst.n);
        assert!(ell >= 2 * k);
        let points = (0..ell).map(|i| F::from((i as u64) + 1)).collect::<Vec<_>>();
        let ws_k = barycentric_weights_consecutive::<F>(k, 1);
        let ws_2k = barycentric_weights_consecutive::<F>(2 * k, 1);
        Self { inst, l, ell, points, ws_k, ws_2k }
    }

    /// Prover: given public `x` and private witness `z_w`, output π = (z_w || w).
    pub fn prove(&self, x: &[F], z_w: &[F]) -> Vec<F> {
        assert_eq!(x.len(), self.l);
        assert_eq!(z_w.len(), self.inst.n - self.l);
        let k = self.inst.k();
        // IMPORTANT: avoid materializing z = (x || z_w) for large instances (multi-million entries).
        // This copy can dominate runtime and look “single-threaded”.
        // A and B mat-vecs are independent: run concurrently.
        let (y_a, y_b) = join(
            || mat_vec_sparse_np(&self.inst.a, x, z_w, self.l),
            || mat_vec_sparse_np(&self.inst.b, x, z_w, self.l),
        );
        // Extrapolations are independent and often dominate: run concurrently.
        let (y_a_tail, y_b_tail) = join(
            || extrapolate_consecutive_next_block::<F>(&y_a),
            || extrapolate_consecutive_next_block::<F>(&y_b),
        );
        let mut w = Vec::with_capacity(2 * k);
        for i in 0..k {
            w.push(y_a[i] * y_b[i]);
        }
        for i in 0..k {
            w.push(y_a_tail[i] * y_b_tail[i]);
        }
        // Output π = (z_w || w) without an extra concat allocation.
        let mut pi = Vec::with_capacity(z_w.len() + w.len());
        pi.extend_from_slice(z_w);
        pi.extend_from_slice(&w);
        pi
    }

    /// Prover + cache: output π = (z_w || w) and also return the RS codeword values needed to
    /// answer verifier coins `(idx, λ)` in **O(1)** time (by indexing), assuming `ell == 2k`.
    pub fn prove_with_codewords(&self, x: &[F], z_w: &[F]) -> (Vec<F>, RsNpCodewords<F>) {
        assert_eq!(x.len(), self.l);
        assert_eq!(z_w.len(), self.inst.n - self.l);
        let k = self.inst.k();
        // Compute y_a, y_b, y_c concurrently.
        let (y_a, (y_b, y_c)) = join(
            || mat_vec_sparse_np(&self.inst.a, x, z_w, self.l),
            || join(
                || mat_vec_sparse_np(&self.inst.b, x, z_w, self.l),
                || mat_vec_sparse_np(&self.inst.c, x, z_w, self.l),
            ),
        );
        // Compute tails.
        //
        // Memory note:
        // For huge k, each extrapolation allocates large NTT buffers. Running 3 extrapolations
        // concurrently can easily exceed machine memory. For large sizes, compute tails sequentially.
        // Memory note:
        // Running 3 huge extrapolations concurrently can explode RSS. Keep bounded parallelism:
        // compute 2 in parallel, then the third.
        let (y_a_tail, y_b_tail) = join(
            || extrapolate_consecutive_next_block::<F>(&y_a),
            || extrapolate_consecutive_next_block::<F>(&y_b),
        );
        // NOTE: we intentionally do **not** compute `y_c_tail`.
        //
        // For the RS-FLPCP query used by the packed lock check, the tail-half query does not use the
        // `C`-part (see coin-form logic in `we_gate_arith::tests::test_large_trace`), so computing
        // `y_c_tail` is wasted work and can be extremely expensive at SP1 scale.

        // Build w (length 2k).
        let mut w = Vec::with_capacity(2 * k);
        for i in 0..k {
            w.push(y_a[i] * y_b[i]);
        }
        for i in 0..k {
            w.push(y_a_tail[i] * y_b_tail[i]);
        }

        // Proof π = (z_w || w).
        let mut pi = Vec::with_capacity(z_w.len() + w.len());
        pi.extend_from_slice(z_w);
        pi.extend_from_slice(&w);

        (
            pi,
            RsNpCodewords {
                y_a,
                y_a_tail,
                y_b,
                y_b_tail,
                y_c,
                w,
            },
        )
    }

    fn map_z_index_to_v(&self, idx: usize) -> (bool, usize) {
        // Returns (is_public, mapped_index)
        if idx < self.l {
            (true, idx)
        } else {
            (false, idx - self.l)
        }
    }
}

impl<F: PrimeField + FftField> BoundedFlpcpSparse<F> for RsDr1csNpFlpcpSparse<F> {
    fn n(&self) -> usize {
        self.l
    }

    fn m(&self) -> usize {
        // proof = z_w (n-l) || w (2k)
        (self.inst.n - self.l) + 2 * self.inst.k()
    }

    fn k(&self) -> usize {
        3
    }

    fn bounds_b(&self) -> Vec<BigInt> {
        // Unbounded FLPCP: same trivial bound.
        let len = self.n() + self.m();
        let p = BigInt::from_bytes_le(num_bigint::Sign::Plus, &F::MODULUS.to_bytes_le());
        let half = (&p - BigInt::one()) / BigInt::from(2u64);
        let b = BigInt::from(len as u64) * &half * &half;
        vec![b.clone(), b.clone(), b]
    }

    fn sample_queries_and_predicate_sparse(
        &self,
        rng: &mut dyn RngCore,
        x: &[F],
    ) -> Result<(Vec<SparseVec<F>>, FlpcpPredicate<F>), String> {
        if x.len() != self.l || self.inst.k() == 0 {
            return Err("bad input".to_string());
        }
        let k = self.inst.k();
        let idx = (rng.next_u64() as usize) % self.ell;
        let lambda = F::from(rng.next_u64());
        let alpha = self.points[idx];

        let lam_k = lagrange_coeffs_at(&self.points[..k], &self.ws_k, alpha);
        let lam_2k = lagrange_coeffs_at(&self.points[..2 * k], &self.ws_2k, alpha);

        // Build sparse linear combos over z indices.
        let q_a_z = lin_combo_rows_sparse(&self.inst.a, self.inst.n, &lam_k);
        let q_b_z = lin_combo_rows_sparse(&self.inst.b, self.inst.n, &lam_k);
        let q_cx2_z = lin_combo_rows_sparse(&self.inst.c, self.inst.n, &lam_2k[..k]);

        // Map z-indices into v=(x || z_w || w) indices.
        let mut q1_terms: Vec<(F, usize)> = Vec::new();
        for (c, idx) in q_a_z.terms {
            let (is_pub, j) = self.map_z_index_to_v(idx);
            let v_idx = if is_pub { j } else { self.l + j };
            q1_terms.push((c, v_idx));
        }
        let mut q2_terms: Vec<(F, usize)> = Vec::new();
        for (c, idx) in q_b_z.terms {
            let (is_pub, j) = self.map_z_index_to_v(idx);
            let v_idx = if is_pub { j } else { self.l + j };
            q2_terms.push((c, v_idx));
        }

        // q3: x/witness part + w part (shift by l + (n-l) = n_total).
        let z_w_len = self.inst.n - self.l;
        let base = self.l + z_w_len;
        let mut q3_terms: Vec<(F, usize)> = Vec::new();
        for (c, idx) in q_cx2_z.terms.iter() {
            let cc = lambda * *c;
            if cc.is_zero() {
                continue;
            }
            let (is_pub, j) = self.map_z_index_to_v(*idx);
            let v_idx = if is_pub { j } else { self.l + j };
            q3_terms.push((cc, v_idx));
        }
        // w part indices are after z_w in the proof: position base + j
        for j in 0..2 * k {
            let coeff = if j < k {
                lam_2k[j] - (lambda * lam_2k[j])
            } else {
                lam_2k[j]
            };
            if !coeff.is_zero() {
                q3_terms.push((coeff, base + j));
            }
        }

        Ok((
            vec![
                SparseVec::new(q1_terms),
                SparseVec::new(q2_terms),
                SparseVec::new(q3_terms),
            ],
            FlpcpPredicate::MulEq,
        ))
    }
}

fn mat_vec<F: Field>(m: &[Vec<F>], x: &[F]) -> Vec<F> {
    // Hot path for large dR1CS instances: each row dot is independent.
    if m.len() >= 256 {
        m.par_iter()
            .map(|row| row.iter().zip(x.iter()).fold(F::ZERO, |acc, (a, b)| acc + (*a * *b)))
            .collect()
    } else {
        m.iter()
            .map(|row| row.iter().zip(x.iter()).fold(F::ZERO, |acc, (a, b)| acc + (*a * *b)))
            .collect()
    }
}

fn mat_vec_sparse<F: PrimeField>(m: &[SparseVec<F>], x: &[F]) -> Vec<F> {
    if m.len() >= 256 {
        m.par_iter().map(|row| row.dot(x)).collect()
    } else {
        m.iter().map(|row| row.dot(x)).collect()
    }
}

fn mat_vec_sparse_np<F: PrimeField>(m: &[SparseVec<F>], x: &[F], z_w: &[F], l: usize) -> Vec<F> {
    debug_assert_eq!(x.len(), l);
    debug_assert_eq!(l + z_w.len(), m.first().map(|_| l + z_w.len()).unwrap_or(l + z_w.len()));
    if m.len() >= 256 {
        m.par_iter()
            .map(|row| {
                row.terms.iter().fold(F::ZERO, |acc, (c, idx)| {
                    let v = if *idx < l { x[*idx] } else { z_w[*idx - l] };
                    acc + (*c * v)
                })
            })
            .collect()
    } else {
        m.iter()
            .map(|row| {
                row.terms.iter().fold(F::ZERO, |acc, (c, idx)| {
                    let v = if *idx < l { x[*idx] } else { z_w[*idx - l] };
                    acc + (*c * v)
                })
            })
            .collect()
    }
}

fn lin_combo_rows<F: Field>(m: &[Vec<F>], coeffs: &[F]) -> Vec<F> {
    let n = m.first().map(|r| r.len()).unwrap_or(0);
    let mut out = vec![F::ZERO; n];
    // This is a matrix-transpose times vector. Parallelize over columns (independent outputs).
    if n >= 1024 && m.len() >= 64 {
        out.par_iter_mut().enumerate().for_each(|(j, out_j)| {
            let mut acc = F::ZERO;
            for (row, c) in m.iter().zip(coeffs.iter()) {
                acc += *c * row[j];
            }
            *out_j = acc;
        });
        return out;
    }
    for (row, c) in m.iter().zip(coeffs.iter()) {
        for j in 0..n {
            out[j] += *c * row[j];
        }
    }
    out
}

fn lin_combo_rows_sparse<F: PrimeField>(m: &[SparseVec<F>], n: usize, coeffs: &[F]) -> SparseVec<F> {
    debug_assert_eq!(m.len(), coeffs.len());
    let mut acc: std::collections::BTreeMap<usize, F> = std::collections::BTreeMap::new();
    for (row, c) in m.iter().zip(coeffs.iter()) {
        if c.is_zero() {
            continue;
        }
        for (aij, idx) in row.terms.iter() {
            debug_assert!(*idx < n);
            let entry = acc.entry(*idx).or_insert(F::ZERO);
            *entry += *c * *aij;
        }
    }
    let terms = acc
        .into_iter()
        .filter_map(|(idx, v)| if v.is_zero() { None } else { Some((v, idx)) })
        .collect::<Vec<_>>();
    SparseVec::new(terms)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::{Fp64, MontBackend, MontConfig};
    use rand_chacha::ChaCha20Rng;
    use rand::SeedableRng;
    use crate::packing::{DppFromBoundedFlpcp, DppFromBoundedFlpcpSparse, PackingError, PackedDppParams};

    #[derive(MontConfig)]
    #[modulus = "10007"]
    #[generator = "5"]
    pub struct F10007Config;
    type F = Fp64<MontBackend<F10007Config, 1>>;

    #[test]
    fn test_rs_dr1cs_flpcp_honest_accepts() {
        // Tiny instance: k=2, n=2; choose A,B,C so that relation holds for x=[1,2].
        // Let A = I, B = I, C = I => require x_i^2 == x_i, which doesn't hold, so instead:
        // Set A=I, B=0, C=0 so relation holds trivially.
        let a = vec![vec![F::ONE, F::ZERO], vec![F::ZERO, F::ONE]];
        let b = vec![vec![F::ZERO, F::ZERO], vec![F::ZERO, F::ZERO]];
        let c = vec![vec![F::ZERO, F::ZERO], vec![F::ZERO, F::ZERO]];
        let inst = Dr1csInstance { a, b, c };
        let flpcp = RsDr1csFlpcp::new(inst, 8);
        let x = vec![F::from(3u64), F::from(5u64)];
        let w = flpcp.prove(&x);

        // Sample verifier queries and check predicate holds for honest answers.
        let mut rng = ChaCha20Rng::seed_from_u64(7);
        let (qs, pred) = flpcp.sample_queries_and_predicate(&mut rng, &x).unwrap();
        let v = [x.clone(), w.clone()].concat();
        let ans = qs
            .iter()
            .map(|q| q.iter().zip(v.iter()).fold(F::ZERO, |acc, (a, b)| acc + (*a * *b)))
            .collect::<Vec<_>>();
        assert!(pred.check(&ans));
    }

    #[test]
    fn test_packing_reports_modulus_too_small_for_unbounded_flpcp() {
        // Same as above, but try to pack directly: should usually fail modulus bound check because b is huge.
        let a = vec![vec![F::ONE, F::ZERO], vec![F::ZERO, F::ONE]];
        let b = vec![vec![F::ZERO, F::ZERO], vec![F::ZERO, F::ZERO]];
        let c = vec![vec![F::ZERO, F::ZERO], vec![F::ZERO, F::ZERO]];
        let inst = Dr1csInstance { a, b, c };
        let flpcp = RsDr1csFlpcp::new(inst, 8);
        let dpp = DppFromBoundedFlpcp::<F, _>::new(flpcp.clone(), crate::packing::PackedDppParams { ell: 32 });

        let x = vec![F::from(3u64), F::from(5u64)];
        let pi = flpcp.prove(&x);
        let mut rng = ChaCha20Rng::seed_from_u64(8);
        let res = dpp.verify(&mut rng, &x, &pi);
        assert!(matches!(res, Err(crate::packing::PackingError::ModulusTooSmall)));
    }

    #[test]
    fn test_rs_dr1cs_flpcp_sparse_honest_accepts() {
        let n = 2usize;
        let a = vec![
            SparseVec::new(vec![(F::ONE, 0)]),
            SparseVec::new(vec![(F::ONE, 1)]),
        ];
        let b = vec![SparseVec::default(), SparseVec::default()];
        let c = vec![SparseVec::default(), SparseVec::default()];
        let inst = Dr1csInstanceSparse { n, a, b, c };
        let flpcp = RsDr1csFlpcpSparse::new(inst, 8);

        let x = vec![F::from(3u64), F::from(5u64)];
        let w = flpcp.prove(&x);

        let mut rng = ChaCha20Rng::seed_from_u64(9);
        let (qs, pred) = flpcp.sample_queries_and_predicate_sparse(&mut rng, &x).unwrap();
        let v = [x.clone(), w.clone()].concat();
        let ans = qs.iter().map(|q| q.dot(&v)).collect::<Vec<_>>();
        assert!(pred.check(&ans));
    }

    #[test]
    fn test_sparse_packing_reports_modulus_too_small_for_unbounded_flpcp() {
        let n = 2usize;
        let a = vec![
            SparseVec::new(vec![(F::ONE, 0)]),
            SparseVec::new(vec![(F::ONE, 1)]),
        ];
        let b = vec![SparseVec::default(), SparseVec::default()];
        let c = vec![SparseVec::default(), SparseVec::default()];
        let inst = Dr1csInstanceSparse { n, a, b, c };
        let flpcp = RsDr1csFlpcpSparse::new(inst, 8);
        let dpp = DppFromBoundedFlpcpSparse::<F, _>::new(flpcp.clone(), PackedDppParams { ell: 32 });

        let x = vec![F::from(3u64), F::from(5u64)];
        let mut rng = ChaCha20Rng::seed_from_u64(10);
        let q = dpp.sample_query(&mut rng, &x);
        assert!(matches!(q, Err(PackingError::ModulusTooSmall)));
    }
}


