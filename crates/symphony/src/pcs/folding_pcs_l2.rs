//! Folding PCS base protocol for ℓ=2 (paper 2024-281, Figure 5).
//!
//! We keep this as a standalone, **paper-faithful** algebraic verifier (and helper ops)
//! because it is the exact object we later arithmetize into dR1CS/DPP.
//!
//! Important modeling choices in Symphony:
//! - We represent vectors in a **row-block** layout compatible with `C^T ⊗ I_n` operations.
//! - Verifier challenges `C1,C2` are treated as Fiat–Shamir coins (Poseidon-derived in the WE gate),
//!   so the canonical proof object should NOT include them.

use ark_ff::{BigInteger, PrimeField};
use ark_std::vec::Vec;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use rayon::prelude::*;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BinMatrix<F: PrimeField> {
    pub rows: usize,
    pub cols: usize,
    /// Row-major entries, each must be 0/1 in `F`.
    pub data: Vec<F>,
}

impl<F: PrimeField> BinMatrix<F> {
    pub fn get(&self, row: usize, col: usize) -> F {
        self.data[row * self.cols + col]
    }
}

/// Dense matrix in row-major form.
#[derive(Clone, Debug)]
pub struct DenseMatrix<F: PrimeField> {
    pub rows: usize,
    pub cols: usize,
    pub data: Vec<F>,
}

impl<F: PrimeField> DenseMatrix<F> {
    pub fn new(rows: usize, cols: usize, data: Vec<F>) -> Self {
        assert_eq!(rows * cols, data.len());
        Self { rows, cols, data }
    }

    #[inline]
    pub fn get(&self, r: usize, c: usize) -> F {
        self.data[r * self.cols + c]
    }

    pub fn mul_vec_par(&self, x: &[F]) -> Vec<F> {
        assert_eq!(x.len(), self.cols);
        (0..self.rows)
            .into_par_iter()
            .map(|r| {
                let mut acc = F::ZERO;
                let base = r * self.cols;
                for c in 0..self.cols {
                    acc += self.data[base + c] * x[c];
                }
                acc
            })
            .collect()
    }
}

#[derive(Clone, Debug)]
pub struct FoldingPcsL2Params<F: PrimeField> {
    pub r: usize,
    pub kappa: usize,
    pub n: usize,
    /// Base-δ gadget: each packed element is represented by `alpha` digits.
    pub delta: u64,
    pub alpha: usize,
    /// Norm bounds for the three short vectors y0,y1,y2 (paper: β0,β1,β2).
    ///
    /// Interpreted as a **per-coordinate signed bound**: each coefficient must be in [-β, β].
    pub beta0: u64,
    pub beta1: u64,
    pub beta2: u64,
    /// Matrix A ∈ F^{n × (r*n*alpha)} (paper: A ∈ Z_q^{n×r n log q}).
    pub a: DenseMatrix<F>,
}

impl<F: PrimeField> FoldingPcsL2Params<F> {
    pub fn y0_len(&self) -> usize { self.kappa * self.r * self.n * self.alpha }
    pub fn y1_len(&self) -> usize { self.kappa * self.r * self.n * self.alpha }
    pub fn y2_len(&self) -> usize { self.kappa * self.r * self.n * self.alpha }

    pub fn f_len(&self) -> usize { self.r * self.r * self.r * self.kappa * self.n }
    pub fn t_len(&self) -> usize { self.kappa * self.n }

    pub fn v0_len(&self) -> usize { self.r * self.kappa * self.n }
    pub fn v1_len(&self) -> usize { self.r * self.kappa * self.n }
    pub fn v2_len(&self) -> usize { self.r * self.kappa * self.n }

    #[inline]
    fn delta_pows(&self) -> Vec<F> {
        let mut pows = Vec::with_capacity(self.alpha);
        let mut acc = F::ONE;
        let delta = F::from(self.delta);
        for _ in 0..self.alpha {
            pows.push(acc);
            acc *= delta;
        }
        pows
    }
}

/// Canonical proof object *without* verifier coins.
#[derive(Clone, Debug)]
pub struct FoldingPcsL2ProofCore<F: PrimeField> {
    pub y0: Vec<F>,
    pub v0: Vec<F>,
    pub y1: Vec<F>,
    pub v1: Vec<F>,
    pub y2: Vec<F>,
    pub v2: Vec<F>,
}

/// Multiply by (C^T ⊗ I_block).
pub fn kron_ct_in_mul<F: PrimeField>(c: &BinMatrix<F>, block_len: usize, v: &[F]) -> Vec<F> {
    assert_eq!(v.len(), c.rows * block_len);
    let mut out = vec![F::ZERO; c.cols * block_len];
    out.par_chunks_mut(block_len)
        .enumerate()
        .for_each(|(col, out_block)| {
            for row in 0..c.rows {
                let bit = c.get(row, col);
                if bit != F::ZERO {
                    let in_block = &v[row * block_len..(row + 1) * block_len];
                    for (o, x) in out_block.iter_mut().zip(in_block.iter()) {
                        *o += bit * *x;
                    }
                }
            }
        });
    out
}

/// Multiply by \((I_{\kappa n} \otimes x^T)\) on a vector in **row-block** layout.
pub fn kron_ikn_xt_mul<F: PrimeField>(x: &[F], kappa: usize, n: usize, v: &[F]) -> Vec<F> {
    let r = x.len();
    assert!(r > 0);
    assert_eq!(v.len(), r * kappa * n);
    (0..kappa)
        .into_par_iter()
        .flat_map_iter(|k| {
            (0..n).map(move |j| {
                let mut acc = F::ZERO;
                for a in 0..r {
                    acc += x[a] * v[((k * r + a) * n) + j];
                }
                acc
            })
        })
        .collect()
}

/// Apply gadget matrix G_{len} to digit vector y_digits, producing packed length `len`.
pub fn gadget_apply_digits<F: PrimeField>(delta_pows: &[F], len: usize, y_digits: &[F]) -> Vec<F> {
    let alpha = delta_pows.len();
    assert_eq!(y_digits.len(), len * alpha);
    (0..len)
        .into_par_iter()
        .map(|i| {
            let mut acc = F::ZERO;
            let base = i * alpha;
            for j in 0..alpha {
                acc += y_digits[base + j] * delta_pows[j];
            }
            acc
        })
        .collect()
}

/// Multiply by (I_κ ⊗ A) where A: n × (r*n*alpha), and input y is κ blocks of length r*n*alpha.
pub fn kron_i_a_mul<F: PrimeField>(a: &DenseMatrix<F>, kappa: usize, in_block_len: usize, y: &[F]) -> Vec<F> {
    assert_eq!(a.cols, in_block_len);
    assert_eq!(y.len(), kappa * in_block_len);
    (0..kappa)
        .into_par_iter()
        .flat_map_iter(|i| {
            let block = &y[i * in_block_len..(i + 1) * in_block_len];
            a.mul_vec_par(block)
        })
        .collect()
}

fn field_to_centered_mag_u64<Fp: PrimeField>(x: Fp) -> u64 {
    let v = BigUint::from_bytes_le(&x.into_bigint().to_bytes_le());
    let p = BigUint::from_bytes_le(&Fp::MODULUS.to_bytes_le());
    let half = &p >> 1;
    if v > half { (&p - v).to_u64().unwrap_or(u64::MAX) } else { v.to_u64().unwrap_or(u64::MAX) }
}

fn check_vec_bound<Fp: PrimeField>(v: &[Fp], beta: u64) -> Result<(), String> {
    for (i, &x) in v.iter().enumerate() {
        let mag = field_to_centered_mag_u64::<Fp>(x);
        if mag > beta {
            return Err(format!("verify: norm bound violated at idx {i}: |x|={mag} > {beta}"));
        }
    }
    Ok(())
}

/// Plain verifier for ℓ=2 folding PCS base protocol, with explicit verifier coins `C1,C2`.
pub fn verify_folding_pcs_l2_with_c_matrices<F: PrimeField>(
    p: &FoldingPcsL2Params<F>,
    t: &[F],
    x0: &[F],
    x1: &[F],
    x2: &[F],
    u: &[F],
    proof: &FoldingPcsL2ProofCore<F>,
    c1: &BinMatrix<F>,
    c2: &BinMatrix<F>,
) -> Result<(), String> {
    if t.len() != p.t_len() || u.len() != p.t_len() {
        return Err("dim: t/u".to_string());
    }
    if x0.len() != p.r || x1.len() != p.r || x2.len() != p.r {
        return Err("dim: x".to_string());
    }
    if proof.y0.len() != p.y0_len()
        || proof.y1.len() != p.y1_len()
        || proof.y2.len() != p.y2_len()
        || proof.v0.len() != p.v0_len()
        || proof.v1.len() != p.v1_len()
        || proof.v2.len() != p.v2_len()
    {
        return Err("dim: proof".to_string());
    }
    if c1.rows != p.r * p.kappa || c1.cols != p.kappa {
        return Err("dim: c1".to_string());
    }
    if c2.rows != p.r * p.kappa || c2.cols != p.kappa {
        return Err("dim: c2".to_string());
    }

    // Norm bounds for y0,y1,y2 (per-coordinate signed bounds).
    check_vec_bound::<F>(&proof.y0, p.beta0)?;
    check_vec_bound::<F>(&proof.y1, p.beta1)?;
    check_vec_bound::<F>(&proof.y2, p.beta2)?;

    let delta_pows = p.delta_pows();

    // (I_k ⊗ A) y0 = t
    let lhs0 = kron_i_a_mul(&p.a, p.kappa, p.r * p.n * p.alpha, &proof.y0);
    if lhs0 != t {
        return Err("verify: lhs0".to_string());
    }

    // (I_k ⊗ A) y1 = (C1^T ⊗ I_n) G y0
    let gy0 = gadget_apply_digits(&delta_pows, p.r * p.kappa * p.n, &proof.y0);
    let rhs1 = kron_ct_in_mul(c1, p.n, &gy0);
    let lhs1 = kron_i_a_mul(&p.a, p.kappa, p.r * p.n * p.alpha, &proof.y1);
    if lhs1 != rhs1 {
        return Err("verify: lhs1".to_string());
    }

    // (I_k ⊗ A) y2 = (C2^T ⊗ I_n) G y1
    let gy1 = gadget_apply_digits(&delta_pows, p.r * p.kappa * p.n, &proof.y1);
    let rhs2 = kron_ct_in_mul(c2, p.n, &gy1);
    let lhs2 = kron_i_a_mul(&p.a, p.kappa, p.r * p.n * p.alpha, &proof.y2);
    if lhs2 != rhs2 {
        return Err("verify: lhs2".to_string());
    }

    // G y2 = v2
    let gy2 = gadget_apply_digits(&delta_pows, p.r * p.kappa * p.n, &proof.y2);
    if gy2 != proof.v2 {
        return Err("verify: gy2".to_string());
    }

    // (I_{κn} ⊗ x2^T) v0 = u
    let u_re = kron_ikn_xt_mul(x2, p.kappa, p.n, &proof.v0);
    if u_re != u {
        return Err("verify: u".to_string());
    }

    // (I_{κn} ⊗ x1^T) v1 = (C1^T ⊗ I_n) v0
    let lhs_v1 = kron_ikn_xt_mul(x1, p.kappa, p.n, &proof.v1);
    let rhs_v1 = kron_ct_in_mul(c1, p.n, &proof.v0);
    if lhs_v1 != rhs_v1 {
        return Err("verify: v1".to_string());
    }

    // (I_{κn} ⊗ x0^T) v2 = (C2^T ⊗ I_n) v1
    let lhs_v2 = kron_ikn_xt_mul(x0, p.kappa, p.n, &proof.v2);
    let rhs_v2 = kron_ct_in_mul(c2, p.n, &proof.v1);
    if lhs_v2 != rhs_v2 {
        return Err("verify: v2".to_string());
    }

    Ok(())
}

// NOTE: We intentionally do NOT provide a wrapper that takes prover-supplied `C1/C2`.
// In this setting, `C1/C2` must be derived from the transcript and provided explicitly.

// ============================================================================
// PROVER FUNCTIONS
// ============================================================================

#[derive(Clone, Debug)]
struct GadgetExactParams {
    delta: u64,
    alpha: usize,
    delta_big: BigUint,
    modulus: BigUint,
    half_modulus: BigUint,
    delta_pow_alpha: BigUint,
}

fn gadget_exact_params<F: PrimeField>(delta: u64, alpha: usize) -> Result<GadgetExactParams, String> {
    if delta < 2 {
        return Err("gadget params: require delta>=2".to_string());
    }
    if alpha == 0 {
        return Err("gadget params: require alpha>=1".to_string());
    }
    let modulus = BigUint::from_bytes_le(&F::MODULUS.to_bytes_le());
    let half_modulus = &modulus >> 1;
    let delta_big = BigUint::from(delta);
    let delta_pow_alpha = delta_big.pow(alpha as u32);
    if delta_pow_alpha < modulus {
        return Err(format!(
            "gadget params: delta^alpha too small (delta={delta}, alpha={alpha}); need delta^alpha >= modulus"
        ));
    }
    Ok(GadgetExactParams {
        delta,
        alpha,
        delta_big,
        modulus,
        half_modulus,
        delta_pow_alpha,
    })
}

fn gadget_decompose_checked_with_params<F: PrimeField>(
    x: &[F],
    gp: &GadgetExactParams,
) -> Result<Vec<F>, String> {
    let len = x.len();
    let alpha = gp.alpha;
    let mut y = vec![F::ZERO; len * alpha];

    y.par_chunks_mut(alpha)
        .enumerate()
        .try_for_each(|(i, digits)| -> Result<(), String> {
            let val_bytes = x[i].into_bigint().to_bytes_le();
            let mut val = BigUint::from_bytes_le(&val_bytes);

            let is_negative = val > gp.half_modulus;
            if is_negative {
                val = &gp.modulus - &val;
            }

            for j in 0..alpha {
                let digit = (&val % &gp.delta_big).to_u64().unwrap_or(0);
                val /= &gp.delta_big;
                if is_negative && digit > 0 {
                    digits[j] = -F::from(digit);
                } else {
                    digits[j] = F::from(digit);
                }
            }

            // Critical soundness guard: reject truncation.
            if val != BigUint::ZERO {
                return Err(format!(
                    "gadget_decompose_checked: truncation at coord i={i} (delta={}, alpha={}, delta^alpha_bits={}, remaining_bits={})",
                    gp.delta,
                    gp.alpha,
                    gp.delta_pow_alpha.bits(),
                    val.bits(),
                ));
            }
            Ok(())
        })?;

    Ok(y)
}

fn hash_level<F: PrimeField>(
    p: &FoldingPcsL2Params<F>,
    gp: &GadgetExactParams,
    cur: &[F],
    groups: usize,
) -> Result<(Vec<F>, Vec<F>), String> {
    // Input is `groups` blocks of length r*n.
    let in_block = p.r * p.n;
    if cur.len() != groups * in_block {
        return Err("commit: hash_level dim mismatch".to_string());
    }

    // For each group, compute s = G^{-1}(x) and then y = A*s.
    let per_group = (0..groups)
        .into_par_iter()
        .map(|gidx| {
            let x = &cur[gidx * in_block..(gidx + 1) * in_block];
            let s = gadget_decompose_checked_with_params::<F>(x, gp)?;
            let y = p.a.mul_vec_par(&s);
            Ok::<_, String>((y, s))
        })
        .collect::<Result<Vec<_>, String>>()?;

    let mut next = Vec::with_capacity(groups * p.n);
    let mut s_all = Vec::with_capacity(cur.len() * p.alpha);
    for (y, s) in per_group {
        next.extend_from_slice(&y);
        s_all.extend_from_slice(&s);
    }
    Ok((next, s_all))
}

/// Commit to a vector `f` using the ℓ=2 folding PCS (paper 2024-281, Fig. 5).
///
/// Returns `(t, s)` where:
/// - `t` is the public commitment (length κ*n)
/// - `s = [s0, s1, s2]` are the nested gadget openings (digits), used by `open`.
pub fn commit<F: PrimeField>(
    p: &FoldingPcsL2Params<F>,
    f: &[F],
) -> Result<(Vec<F>, [Vec<F>; 3]), String> {
    if f.len() != p.f_len() {
        return Err(format!("commit: f length {} != f_len {}", f.len(), p.f_len()));
    }

    let r = p.r;
    let kappa = p.kappa;
    // Precompute and validate gadget params once (avoid per-block BigUint pow/modulus work).
    let gp = gadget_exact_params::<F>(p.delta, p.alpha)?;

    // Level-2: f -> t2
    let (t2, s2) = hash_level(p, &gp, f, r * r * kappa)?;
    // Level-1: t2 -> t1
    let (t1, s1) = hash_level(p, &gp, &t2, r * kappa)?;
    // Level-0: t1 -> t
    let (t, s0) = hash_level(p, &gp, &t1, kappa)?;

    Ok((t, [s0, s1, s2]))
}

/// Open a PCS commitment at evaluation point `(x0, x1, x2)` with folding challenges `C1, C2`.
///
/// Returns `(u, proof)` where:
/// - `u` is the claimed evaluation (length κ*n)
/// - `proof` is the opening proof
pub fn open<F: PrimeField>(
    p: &FoldingPcsL2Params<F>,
    f: &[F],
    s: &[Vec<F>; 3],
    x0: &[F],
    x1: &[F],
    x2: &[F],
    c1: &BinMatrix<F>,
    c2: &BinMatrix<F>,
) -> Result<(Vec<F>, FoldingPcsL2ProofCore<F>), String> {
    if f.len() != p.f_len() {
        return Err("open: f length mismatch".to_string());
    }
    if x0.len() != p.r || x1.len() != p.r || x2.len() != p.r {
        return Err("open: x dimension mismatch".to_string());
    }

    // Openings must be present (these are the nested gadget digits produced by `commit`).
    let s0 = &s[0];
    let s1 = &s[1];
    let s2 = &s[2];
    if s0.len() != p.y0_len() {
        return Err("open: s0 length mismatch".to_string());
    }
    if s1.len() != (p.r * p.r * p.kappa * p.n * p.alpha) {
        return Err("open: s1 length mismatch".to_string());
    }
    if s2.len() != (p.f_len() * p.alpha) {
        return Err("open: s2 length mismatch".to_string());
    }

    let r = p.r;
    let n = p.n;
    let kappa = p.kappa;
    let alpha = p.alpha;

    // Prover messages y0, v0.
    let y0 = s0.clone();
    let tmp0 = kron_ikn_xt_mul(x0, r * r * kappa, n, f); // len r^2 κ n
    let v0 = kron_ikn_xt_mul(x1, r * kappa, n, &tmp0); // len r κ n
    let u = kron_ikn_xt_mul(x2, kappa, n, &v0); // len κ n

    // y1, v1.
    let y1 = kron_ct_in_mul(c1, r * n * alpha, s1);
    let t_c1_f = kron_ct_in_mul(c1, r * r * n, f);
    let v1 = kron_ikn_xt_mul(x0, r * kappa, n, &t_c1_f);

    // y2, v2.
    let t_c1_s2 = kron_ct_in_mul(c1, r * r * n * alpha, s2);
    let y2 = kron_ct_in_mul(c2, r * n * alpha, &t_c1_s2);
    let v2 = kron_ct_in_mul(c2, r * n, &t_c1_f);

    let proof = FoldingPcsL2ProofCore {
        y0,
        v0,
        y1,
        v1,
        y2,
        v2,
    };

    Ok((u, proof))
}

/// Compute MLE evaluation at point using tensor-product structure.
///
/// For a vector `f` of length `r^3`, interprets it as an MLE and evaluates at `(x0, x1, x2)`.
/// This is: `⟨x0 ⊗ x1 ⊗ x2, f⟩ = Σ_{a,b,c} x0[a] * x1[b] * x2[c] * f[a*r^2 + b*r + c]`
pub fn eval_mle_tensor<F: PrimeField>(f: &[F], x0: &[F], x1: &[F], x2: &[F]) -> F {
    let r = x0.len();
    assert_eq!(x1.len(), r);
    assert_eq!(x2.len(), r);
    assert_eq!(f.len(), r * r * r);

    (0..r)
        .into_par_iter()
        .map(|a| {
            let mut sum_bc = F::ZERO;
            for b in 0..r {
                for c in 0..r {
                    let idx = a * r * r + b * r + c;
                    sum_bc += x1[b] * x2[c] * f[idx];
                }
            }
            x0[a] * sum_bc
        })
        .sum()
}
