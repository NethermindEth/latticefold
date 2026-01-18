//! dR1CS arithmetization of the ℓ=2 folding PCS base verifier (paper 2024-281, Fig. 5).
//!
//! DPP-friendly immediately:
//! - C1, C2 are witness-provided bits (constrained boolean).
//! - No Poseidon-bit extraction / transcript binding yet.

use ark_ff::{BigInteger, PrimeField};
use ark_std::vec::Vec;
use num_bigint::BigUint;
use num_traits::ToPrimitive;

use crate::dpp_sumcheck::Dr1csBuilder;

use super::folding_pcs_l2::{BinMatrix, DenseMatrix, FoldingPcsL2Params, FoldingPcsL2ProofCore};

fn constrain_binary<F: PrimeField>(b: &mut Dr1csBuilder<F>, bit: usize) {
    let one_minus = b.new_var(F::ONE - b.assignment[bit]);
    b.add_constraint(
        vec![(F::ONE, b.one()), (-F::ONE, bit)],
        vec![(F::ONE, b.one())],
        vec![(F::ONE, one_minus)],
    );
    b.add_constraint(
        vec![(F::ONE, bit)],
        vec![(F::ONE, one_minus)],
        vec![(F::ZERO, b.one())],
    );
}

fn constraint_mul_eq_zero<F: PrimeField>(b: &mut Dr1csBuilder<F>, x: usize, y: usize) {
    // x * y = 0
    b.add_constraint(vec![(F::ONE, x)], vec![(F::ONE, y)], vec![(F::ZERO, b.one())]);
}

fn constrain_bit_reconstruction<F: PrimeField>(b: &mut Dr1csBuilder<F>, bits: &[usize]) -> usize {
    // abs = Σ 2^i * bits[i]
    let mut power = F::ONE;
    let mut sum = F::ZERO;
    for &bit in bits {
        sum += power * b.assignment[bit];
        power = power.double();
    }
    let abs = b.new_var(sum);
    let mut lc: Vec<(F, usize)> = Vec::with_capacity(bits.len() + 1);
    lc.push((F::ONE, abs));
    let mut power2 = F::ONE;
    for &bit in bits {
        lc.push((-power2, bit));
        power2 = power2.double();
    }
    b.enforce_lc_times_one_eq_const(lc);
    abs
}

fn constrain_abs_leq_const<F: PrimeField>(b: &mut Dr1csBuilder<F>, bits_le: &[usize], beta: u64) -> Result<(), String> {
    // Enforce value(bits_le) <= beta using prefix-equality trick.
    //
    // Let bits be little-endian; we scan MSB->LSB.
    let num_bits = bits_le.len();
    if num_bits == 0 {
        return Ok(());
    }
    let mut prefix_eq = b.one(); // 1 if all higher bits equal beta so far
    for i in (0..num_bits).rev() {
        let bit = bits_le[i];
        let beta_bit = ((beta >> i) & 1) as u64;
        if beta_bit == 0 {
            // If prefix_eq == 1, then bit must be 0: prefix_eq * bit = 0.
            constraint_mul_eq_zero(b, prefix_eq, bit);
            // Update prefix_eq' = prefix_eq * (1 - bit)
            let one_minus = b.new_var(F::ONE - b.assignment[bit]);
            b.add_constraint(
                vec![(F::ONE, b.one()), (-F::ONE, bit)],
                vec![(F::ONE, b.one())],
                vec![(F::ONE, one_minus)],
            );
            let next = b.new_var(b.assignment[prefix_eq] * b.assignment[one_minus]);
            b.enforce_mul(prefix_eq, one_minus, next);
            prefix_eq = next;
        } else {
            // beta_bit == 1: update prefix_eq' = prefix_eq * bit
            let next = b.new_var(b.assignment[prefix_eq] * b.assignment[bit]);
            b.enforce_mul(prefix_eq, bit, next);
            prefix_eq = next;
        }
    }
    Ok(())
}

fn field_to_centered_sign_mag<F: PrimeField>(val: F) -> (bool, u64) {
    // Convert val ∈ F to centered integer in [-p/2, p/2], return (sign, magnitude).
    // This is only used for witness generation (bits/sign vars).
    let v = BigUint::from_bytes_le(&val.into_bigint().to_bytes_le());
    let p = BigUint::from_bytes_le(&F::MODULUS.to_bytes_le());
    let half = &p >> 1;
    if v > half {
        // negative
        let mag = (&p - v).to_u64().unwrap_or(u64::MAX);
        (true, mag)
    } else {
        let mag = v.to_u64().unwrap_or(u64::MAX);
        (false, mag)
    }
}

fn constrain_signed_in_range<F: PrimeField>(
    b: &mut Dr1csBuilder<F>,
    coeff_var: usize,
    beta: u64,
) -> Result<(), String> {
    // Signed decomposition: coeff = (1 - 2*sign) * abs, with abs in [0..beta].
    //
    // - Allocate sign bit (boolean)
    // - Allocate abs bits (boolean) for abs in [0..2^k)
    // - Enforce abs <= beta (exact comparator)
    // - Enforce coeff = abs - 2*(sign*abs)
    let (sign, mag) = field_to_centered_sign_mag::<F>(b.assignment[coeff_var]);
    let sign_bit = b.new_var(if sign { F::ONE } else { F::ZERO });
    constrain_binary(b, sign_bit);

    // k = ceil(log2(beta+1)) with k>=1. (If beta=0, k=1.)
    let num_bits = if beta == 0 { 1 } else { 64usize - beta.leading_zeros() as usize };

    let mut bits = Vec::with_capacity(num_bits);
    let mut tmp = mag;
    for _ in 0..num_bits {
        let bit = b.new_var(if (tmp & 1) == 1 { F::ONE } else { F::ZERO });
        constrain_binary(b, bit);
        bits.push(bit);
        tmp >>= 1;
    }

    let abs = constrain_bit_reconstruction(b, &bits);
    constrain_abs_leq_const(b, &bits, beta)?;

    // t = sign * abs
    let t_val = b.assignment[sign_bit] * b.assignment[abs];
    let t = b.new_var(t_val);
    b.enforce_mul(sign_bit, abs, t);

    // coeff = abs - 2*t
    let two_t_val = t_val.double();
    let two_t = b.new_var(two_t_val);
    b.add_constraint(
        vec![(F::ONE.double(), t)],
        vec![(F::ONE, b.one())],
        vec![(F::ONE, two_t)],
    );
    b.enforce_lc_times_one_eq_const(vec![(F::ONE, abs), (-F::ONE, two_t), (-F::ONE, coeff_var)]);

    Ok(())
}

fn constrain_vec_signed_bounds<F: PrimeField>(
    b: &mut Dr1csBuilder<F>,
    vars: &[usize],
    beta: u64,
) -> Result<(), String> {
    for &v in vars {
        constrain_signed_in_range::<F>(b, v, beta)?;
    }
    Ok(())
}

fn alloc_vec<F: PrimeField>(b: &mut Dr1csBuilder<F>, v: &[F]) -> Vec<usize> {
    v.iter().map(|&x| b.new_var(x)).collect()
}

fn enforce_vec_eq<F: PrimeField>(b: &mut Dr1csBuilder<F>, lhs: &[usize], rhs: &[usize]) -> Result<(), String> {
    if lhs.len() != rhs.len() {
        return Err("len mismatch".to_string());
    }
    for (&a, &c) in lhs.iter().zip(rhs.iter()) {
        b.enforce_lc_times_one_eq_const(vec![(F::ONE, a), (-F::ONE, c)]);
    }
    Ok(())
}


/// out = (I_kappa ⊗ A) * y, where y is kappa blocks of length in_block_len.
fn kron_i_a_mul_dr1cs<F: PrimeField>(
    b: &mut Dr1csBuilder<F>,
    a: &DenseMatrix<F>,
    kappa: usize,
    in_block_len: usize,
    y: &[usize],
) -> Result<Vec<usize>, String> {
    if y.len() != kappa * in_block_len {
        return Err("kron_i_a_mul: y len mismatch".to_string());
    }
    if a.cols != in_block_len {
        return Err("kron_i_a_mul: A cols mismatch".to_string());
    }
    let mut out = Vec::with_capacity(kappa * a.rows);
    for blk in 0..kappa {
        let yblk = &y[blk * in_block_len..(blk + 1) * in_block_len];
        for r in 0..a.rows {
            let mut acc = F::ZERO;
            let mut lc: Vec<(F, usize)> = Vec::with_capacity(in_block_len + 1);
            for c in 0..in_block_len {
                let coeff = a.get(r, c);
                if coeff != F::ZERO {
                    acc += coeff * b.assignment[yblk[c]];
                    lc.push((coeff, yblk[c]));
                }
            }
            let v = b.new_var(acc);
            lc.push((-F::ONE, v));
            b.enforce_lc_times_one_eq_const(lc);
            out.push(v);
        }
    }
    Ok(out)
}

/// out[i] = sum_{j=0..alpha-1} y_digits[i*alpha+j] * delta^j.
fn gadget_apply_digits_dr1cs<F: PrimeField>(
    b: &mut Dr1csBuilder<F>,
    delta: u64,
    alpha: usize,
    y_digits: &[usize],
) -> Result<Vec<usize>, String> {
    if y_digits.len() % alpha != 0 {
        return Err("gadget_apply_digits: len % alpha != 0".to_string());
    }
    let len = y_digits.len() / alpha;
    let mut pows = Vec::with_capacity(alpha);
    let mut accp = F::ONE;
    let d = F::from(delta);
    for _ in 0..alpha {
        pows.push(accp);
        accp *= d;
    }

    let mut out = Vec::with_capacity(len);
    for i in 0..len {
        let mut acc = F::ZERO;
        let mut lc: Vec<(F, usize)> = Vec::with_capacity(alpha + 1);
        for j in 0..alpha {
            let vj = y_digits[i * alpha + j];
            acc += pows[j] * b.assignment[vj];
            lc.push((pows[j], vj));
        }
        let v = b.new_var(acc);
        lc.push((-F::ONE, v));
        b.enforce_lc_times_one_eq_const(lc);
        out.push(v);
    }
    Ok(out)
}

/// out = (C^T ⊗ I_block) * v
///
/// C entries are witness bit variables.
fn kron_ct_in_mul_dr1cs<F: PrimeField>(
    b: &mut Dr1csBuilder<F>,
    c_rows: usize,
    c_cols: usize,
    c_bits: &[usize], // length rows*cols
    block_len: usize,
    v: &[usize], // length rows*block_len
) -> Result<Vec<usize>, String> {
    if c_bits.len() != c_rows * c_cols {
        return Err("kron_ct_in_mul: C len mismatch".to_string());
    }
    if v.len() != c_rows * block_len {
        return Err("kron_ct_in_mul: v len mismatch".to_string());
    }

    let mut out = Vec::with_capacity(c_cols * block_len);
    for col in 0..c_cols {
        for i in 0..block_len {
            let mut acc = F::ZERO;
            let mut sum_terms: Vec<(F, usize)> = Vec::new();
            for row in 0..c_rows {
                let bit = c_bits[row * c_cols + col];
                let val = v[row * block_len + i];
                let prod_val = b.assignment[bit] * b.assignment[val];
                let prod = b.new_var(prod_val);
                b.enforce_mul(bit, val, prod);
                acc += prod_val;
                sum_terms.push((F::ONE, prod));
            }
            let outv = b.new_var(acc);
            sum_terms.push((-F::ONE, outv));
            b.enforce_lc_times_one_eq_const(sum_terms);
            out.push(outv);
        }
    }
    Ok(out)
}

/// out = (I_blocks ⊗ x^T) v, where x are public constants.
fn kron_ikn_xt_mul_dr1cs<F: PrimeField>(
    b: &mut Dr1csBuilder<F>,
    x: &[F],      // length r
    kappa: usize,
    n: usize,
    v: &[usize],  // length r*kappa*n in row-block layout
) -> Result<Vec<usize>, String> {
    let r = x.len();
    if r == 0 {
        return Err("kron_ikn_xt_mul: r=0".to_string());
    }
    if v.len() != r * kappa * n {
        return Err("kron_ikn_xt_mul: dim".to_string());
    }

    // This is purely linear in the witness (x are public constants), so no mul constraints needed.
    let mut out = Vec::with_capacity(kappa * n);
    for k in 0..kappa {
        for j in 0..n {
            let mut acc = F::ZERO;
            let mut lc: Vec<(F, usize)> = Vec::with_capacity(r + 1);
            for a in 0..r {
                let idx = v[((k * r + a) * n) + j];
                if x[a] != F::ZERO {
                    acc += x[a] * b.assignment[idx];
                    lc.push((x[a], idx));
                }
            }
            let outv = b.new_var(acc);
            lc.push((-F::ONE, outv));
            b.enforce_lc_times_one_eq_const(lc);
            out.push(outv);
        }
    }
    Ok(out)
}

#[derive(Clone, Debug)]
pub struct FoldingPcsL2Wiring {
    pub c1_bits: Vec<usize>,
    pub c2_bits: Vec<usize>,
    /// Witness variables for the public commitment vector `t` (so callers can glue it to Poseidon Absorb).
    pub t_vars: Vec<usize>,
    /// Witness variables for the computed evaluation `u_re = (I_{κn} ⊗ x2^T) v0`.
    pub u_re_vars: Vec<usize>,
}

/// Build a row-major bit-matrix (flattened) from a bitstream.
///
/// This is intended for transcript integration: once Poseidon byte/bit extraction is wired,
/// you can feed the resulting bitstream here to obtain `C1/C2` bits without the prover
/// choosing them.
pub fn bin_matrix_bits_from_stream<F: PrimeField>(
    bits: &[usize],
    rows: usize,
    cols: usize,
) -> Result<Vec<usize>, String> {
    let need = rows * cols;
    if bits.len() < need {
        return Err(format!("bitstream too short: need {need}, have {}", bits.len()));
    }
    Ok(bits[..need].to_vec())
}

fn byte_to_bits_le_dr1cs<F: PrimeField>(b: &mut Dr1csBuilder<F>, byte_var: usize) -> Vec<usize> {
    // Allocate 8 bits, constrain boolean, and constrain byte = Σ 2^i * bit_i.
    let byte_val_u64 = b.assignment[byte_var].into_bigint().as_ref().get(0).copied().unwrap_or(0);
    let mut bits = Vec::with_capacity(8);
    for i in 0..8 {
        let bi = (byte_val_u64 >> i) & 1;
        let v = b.new_var(if bi == 1 { F::ONE } else { F::ZERO });
        constrain_binary(b, v);
        bits.push(v);
    }
    let recon = constrain_bit_reconstruction(b, &bits);
    // byte - recon = 0
    b.enforce_lc_times_one_eq_const(vec![(F::ONE, byte_var), (-F::ONE, recon)]);
    bits
}

/// Expand byte variables into a little-endian bitstream (LSB-first per byte).
pub fn bytes_to_bitstream_le_dr1cs<F: PrimeField>(
    b: &mut Dr1csBuilder<F>,
    byte_vars: &[usize],
) -> Vec<usize> {
    let mut out = Vec::with_capacity(byte_vars.len() * 8);
    for &bv in byte_vars {
        out.extend(byte_to_bits_le_dr1cs::<F>(b, bv));
    }
    out
}

/// Verify folding PCS ℓ=2 in dR1CS, with `C1/C2` supplied as bit-variables.
///
/// This is the interface needed for WE/DPP integration: `C1/C2` should ultimately
/// be derived from the transcript (Poseidon), not prover-chosen.
pub fn folding_pcs_l2_verify_dr1cs_with_c_bits<F: PrimeField>(
    b: &mut Dr1csBuilder<F>,
    p: &FoldingPcsL2Params<F>,
    t: &[F],
    x0: &[F],
    x1: &[F],
    x2: &[F],
    proof: &FoldingPcsL2ProofCore<F>,
    c1_bits: &[usize], // length (r*kappa)*kappa
    c2_bits: &[usize], // length (r*kappa)*kappa
) -> Result<FoldingPcsL2Wiring, String> {
    let y0 = alloc_vec(b, &proof.y0);
    let y1 = alloc_vec(b, &proof.y1);
    let y2 = alloc_vec(b, &proof.y2);
    let v0 = alloc_vec(b, &proof.v0);
    let v1 = alloc_vec(b, &proof.v1);
    let v2 = alloc_vec(b, &proof.v2);
    let t_vars = alloc_vec(b, t);

    if c1_bits.len() != (p.r * p.kappa) * p.kappa {
        return Err("C1 bits length mismatch".to_string());
    }
    if c2_bits.len() != (p.r * p.kappa) * p.kappa {
        return Err("C2 bits length mismatch".to_string());
    }
    for &bit in c1_bits.iter().chain(c2_bits.iter()) {
        constrain_binary(b, bit);
    }

    // Norm bounds (per-coordinate signed bounds) for y0,y1,y2.
    constrain_vec_signed_bounds::<F>(b, &y0, p.beta0)?;
    constrain_vec_signed_bounds::<F>(b, &y1, p.beta1)?;
    constrain_vec_signed_bounds::<F>(b, &y2, p.beta2)?;

    let lhs0 = kron_i_a_mul_dr1cs(b, &p.a, p.kappa, p.r * p.n * p.alpha, &y0)?;
    enforce_vec_eq(b, &lhs0, &t_vars)?;

    // NOTE: In Fig. 5, v0 is *not* required to equal G(y0).
    // It is a prover message that is later tied to u / v1 / v2 via the eval-chain checks.
    let gy0 = gadget_apply_digits_dr1cs(b, p.delta, p.alpha, &y0)?;
    let rhs1 = kron_ct_in_mul_dr1cs(b, p.r * p.kappa, p.kappa, &c1_bits, p.n, &gy0)?;
    let lhs1 = kron_i_a_mul_dr1cs(b, &p.a, p.kappa, p.r * p.n * p.alpha, &y1)?;
    enforce_vec_eq(b, &lhs1, &rhs1)?;

    // NOTE: In Fig. 5, v1 is *not* required to equal G(y1).
    let gy1 = gadget_apply_digits_dr1cs(b, p.delta, p.alpha, &y1)?;
    let rhs2 = kron_ct_in_mul_dr1cs(b, p.r * p.kappa, p.kappa, &c2_bits, p.n, &gy1)?;
    let lhs2 = kron_i_a_mul_dr1cs(b, &p.a, p.kappa, p.r * p.n * p.alpha, &y2)?;
    enforce_vec_eq(b, &lhs2, &rhs2)?;

    let gy2 = gadget_apply_digits_dr1cs(b, p.delta, p.alpha, &y2)?;
    enforce_vec_eq(b, &gy2, &v2)?;

    let u_re = kron_ikn_xt_mul_dr1cs(b, x2, p.kappa, p.n, &v0)?;

    let lhs_v1 = kron_ikn_xt_mul_dr1cs(b, x1, p.kappa, p.n, &v1)?;
    let rhs_v1 = kron_ct_in_mul_dr1cs(b, p.r * p.kappa, p.kappa, &c1_bits, p.n, &v0)?;
    enforce_vec_eq(b, &lhs_v1, &rhs_v1)?;

    let lhs_v2 = kron_ikn_xt_mul_dr1cs(b, x0, p.kappa, p.n, &v2)?;
    let rhs_v2 = kron_ct_in_mul_dr1cs(b, p.r * p.kappa, p.kappa, &c2_bits, p.n, &v1)?;
    enforce_vec_eq(b, &lhs_v2, &rhs_v2)?;

    Ok(FoldingPcsL2Wiring {
        c1_bits: c1_bits.to_vec(),
        c2_bits: c2_bits.to_vec(),
        t_vars,
        u_re_vars: u_re,
    })
}

/// Convenience wrapper: derive `C1/C2` from a provided byte stream (little-endian bits),
/// then verify.
///
/// This is meant for the upcoming Poseidon integration, where the transcript exposes
/// a stream of bytes and we interpret them as C-matrix bits.
pub fn folding_pcs_l2_verify_dr1cs_with_c_bytes<F: PrimeField>(
    b: &mut Dr1csBuilder<F>,
    p: &FoldingPcsL2Params<F>,
    t: &[F],
    x0: &[F],
    x1: &[F],
    x2: &[F],
    proof: &FoldingPcsL2ProofCore<F>,
    c_bytes: &[usize],
) -> Result<FoldingPcsL2Wiring, String> {
    let bits = bytes_to_bitstream_le_dr1cs::<F>(b, c_bytes);
    let need = (p.r * p.kappa) * p.kappa;
    let c1_bits = bin_matrix_bits_from_stream::<F>(&bits, p.r * p.kappa, p.kappa)?;
    let c2_bits = bin_matrix_bits_from_stream::<F>(&bits[need..], p.r * p.kappa, p.kappa)?;
    folding_pcs_l2_verify_dr1cs_with_c_bits(b, p, t, x0, x1, x2, proof, &c1_bits, &c2_bits)
}

// NOTE: no wrapper that takes prover-supplied `C1/C2` — coins must come from transcript.

pub fn folding_pcs_l2_params<F: PrimeField>(
    r: usize,
    kappa: usize,
    n: usize,
    delta: u64,
    alpha: usize,
    beta0: u64,
    beta1: u64,
    beta2: u64,
    a: DenseMatrix<F>,
) -> FoldingPcsL2Params<F> {
    FoldingPcsL2Params { r, kappa, n, delta, alpha, beta0, beta1, beta2, a }
}

pub fn bin_matrix<F: PrimeField>(rows: usize, cols: usize, data: Vec<F>) -> BinMatrix<F> {
    BinMatrix { rows, cols, data }
}

