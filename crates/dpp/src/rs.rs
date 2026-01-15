//! Systematic Reed–Solomon utilities (prototype).
//!
//! We implement enough RS machinery to support the dR1CS FLPCP construction in Section 4.1:
//! - fixed evaluation points α_0..α_{ℓ-1}
//! - barycentric weights for fast evaluation from the first `k` (or `2k`) “systematic” values
//! - Lagrange coefficients at a given evaluation point

use ark_ff::{BigInteger, Field, FftField, PrimeField};
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};

use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::{Mutex, OnceLock};

/// Compute barycentric weights for distinct points `xs`.
///
/// w_i = 1 / ∏_{j≠i} (x_i - x_j)
pub fn barycentric_weights<F: Field>(xs: &[F]) -> Vec<F> {
    let n = xs.len();
    let mut w = vec![F::ZERO; n];
    for i in 0..n {
        let mut denom = F::ONE;
        for j in 0..n {
            if i == j {
                continue;
            }
            denom *= xs[i] - xs[j];
        }
        w[i] = denom.inverse().expect("distinct points => invertible denom");
    }
    w
}

/// Batch-invert nonzero elements in-place using 1 field inversion.
fn batch_inverse_in_place<F: Field>(v: &mut [F]) {
    let n = v.len();
    if n == 0 {
        return;
    }
    // prefix[i] = ∏_{j< i} v[j]
    let mut prefix = Vec::with_capacity(n);
    let mut acc = F::ONE;
    for x in v.iter() {
        prefix.push(acc);
        acc *= *x;
    }
    let mut inv_acc = acc.inverse().expect("batch inversion: nonzero product");
    for i in (0..n).rev() {
        let xi = v[i];
        v[i] = inv_acc * prefix[i];
        inv_acc *= xi;
    }
}

/// Compute barycentric weights for consecutive integer points.
///
/// Points are `xs[i] = start + i` interpreted in the field.
///
/// For consecutive points, we have the closed form:
/// \[
///   \prod_{j\neq i} (x_i - x_j) = i!\,(-1)^{n-1-i}\,(n-1-i)!
/// \]
/// so
/// \[
///   w_i = (-1)^{n-1-i} / (i!\,(n-1-i)!)
/// \]
///
/// This computes all weights in **O(n)** time using factorial and inverse factorial tables.
pub fn barycentric_weights_consecutive<F: Field>(n: usize, start: u64) -> Vec<F> {
    if n == 0 {
        return Vec::new();
    }
    // factorials up to n-1 in the field
    let mut fact = vec![F::ONE; n];
    for i in 1..n {
        fact[i] = fact[i - 1] * F::from(i as u64);
    }
    // inv factorials
    let mut inv_fact = vec![F::ONE; n];
    inv_fact[n - 1] = fact[n - 1].inverse().expect("field nonzero factorial invertible");
    for i in (1..n).rev() {
        inv_fact[i - 1] = inv_fact[i] * F::from(i as u64);
    }
    // weights
    let mut w = vec![F::ZERO; n];
    for i in 0..n {
        let mut wi = inv_fact[i] * inv_fact[n - 1 - i];
        // (-1)^{n-1-i}
        if ((n - 1 - i) & 1) == 1 {
            wi = -wi;
        }
        w[i] = wi;
    }
    // Touch `start` so callers can’t accidentally ignore it when passing other sequences.
    // (The weights for consecutive points are translation-invariant.)
    let _ = start;
    w
}

/// Evaluate the unique degree < `xs.len()` polynomial interpolating `(xs[i], ys[i])` at `x`.
///
/// Uses barycentric interpolation.
pub fn barycentric_eval<F: Field>(xs: &[F], ws: &[F], ys: &[F], x: F) -> F {
    debug_assert_eq!(xs.len(), ws.len());
    debug_assert_eq!(xs.len(), ys.len());
    // If x is one of the xs, return the corresponding y.
    for (xi, yi) in xs.iter().zip(ys.iter()) {
        if *xi == x {
            return *yi;
        }
    }
    let mut num = F::ZERO;
    let mut den = F::ZERO;
    for i in 0..xs.len() {
        let inv = (x - xs[i]).inverse().expect("x != xs[i] so invertible");
        let t = ws[i] * inv;
        num += t * ys[i];
        den += t;
    }
    num * den.inverse().expect("den != 0 for distinct points")
}

/// Compute Lagrange coefficients `λ_i(x)` for the basis over points `xs` evaluated at `x`.
///
/// Returns vector `lambda` such that for any `ys`, `f(x) = Σ_i lambda[i] * ys[i]`.
pub fn lagrange_coeffs_at<F: Field>(xs: &[F], ws: &[F], x: F) -> Vec<F> {
    debug_assert_eq!(xs.len(), ws.len());
    // If x == xs[i], lambda is unit vector.
    for (i, xi) in xs.iter().enumerate() {
        if *xi == x {
            let mut out = vec![F::ZERO; xs.len()];
            out[i] = F::ONE;
            return out;
        }
    }
    // Compute inv(x - x_i) via batch inversion (1 inversion total).
    let mut diffs = xs.iter().map(|xi| x - *xi).collect::<Vec<_>>();
    batch_inverse_in_place(&mut diffs);

    let mut denom = F::ZERO;
    let mut ts = Vec::with_capacity(xs.len());
    for i in 0..xs.len() {
        let t = ws[i] * diffs[i];
        ts.push(t);
        denom += t;
    }
    let denom_inv = denom.inverse().expect("den != 0");
    ts.into_iter().map(|t| t * denom_inv).collect()
}

/// Convolution (polynomial multiplication).
///
/// Uses radix-2 FFT when the field supports a large enough 2-adic domain; otherwise falls back
/// to the naive O(n^2) algorithm (only intended for small test sizes).
fn convolution<F: FftField + PrimeField>(a: &[F], b: &[F]) -> Vec<F> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let out_len = a.len() + b.len() - 1;
    let size = out_len.next_power_of_two();
    let domain = Radix2EvaluationDomain::<F>::new(size);

    if let Some(domain) = domain {
        let mut fa = vec![F::ZERO; size];
        fa[..a.len()].copy_from_slice(a);
        let mut fb = vec![F::ZERO; size];
        fb[..b.len()].copy_from_slice(b);

        domain.fft_in_place(&mut fa);
        domain.fft_in_place(&mut fb);
        for i in 0..size {
            fa[i] *= fb[i];
        }
        domain.ifft_in_place(&mut fa);
        fa.truncate(out_len);
        fa
    } else {
        // No large 2-adic domain in this field (e.g. Frog prime has v2(p-1)=3).
        //
        // For large instances, we use CRT+NTT convolution over a few NTT-friendly primes
        // and reconstruct the exact integer coefficient within a known bound, then reduce mod p.
        //
        // This keeps the RS extrapolation scalable even when `F` itself is not FFT-friendly.
        if out_len > (1 << 16) {
            convolution_crt_ntt::<F>(a, b, out_len, size)
        } else {
            // Small-size fallback.
            let mut out = vec![F::ZERO; out_len];
            for (i, ai) in a.iter().enumerate() {
                for (j, bj) in b.iter().enumerate() {
                    out[i + j] += *ai * *bj;
                }
            }
            out
        }
    }
}

// ---------------- CRT+NTT convolution for 64-bit prime fields ----------------

#[inline]
fn mul_mod_u32(a: u32, b: u32, m: u32) -> u32 {
    ((a as u64 * b as u64) % (m as u64)) as u32
}

#[inline]
fn add_mod_u32(a: u32, b: u32, m: u32) -> u32 {
    let s = a as u64 + b as u64;
    (if s >= m as u64 { s - m as u64 } else { s }) as u32
}

#[inline]
fn sub_mod_u32(a: u32, b: u32, m: u32) -> u32 {
    (if a >= b { a - b } else { a + m - b }) as u32
}

fn pow_mod_u32(mut a: u32, mut e: u32, m: u32) -> u32 {
    let mut r: u32 = 1;
    while e > 0 {
        if (e & 1) == 1 {
            r = mul_mod_u32(r, a, m);
        }
        a = mul_mod_u32(a, a, m);
        e >>= 1;
    }
    r
}

fn inv_mod_u32(a: u32, m: u32) -> u32 {
    // m is prime
    pow_mod_u32(a, m - 2, m)
}

fn bit_reverse_permute(a: &mut [u32]) {
    let n = a.len();
    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            a.swap(i, j);
        }
    }
}

fn ntt(a: &mut [u32], invert: bool, modulus: u32, primitive_root: u32) {
    let n = a.len();
    bit_reverse_permute(a);

    // Compute n-th primitive root: g^{(p-1)/n}
    let root = pow_mod_u32(primitive_root, (modulus - 1) / (n as u32), modulus);
    let root_inv = inv_mod_u32(root, modulus);

    // Stage-parallel NTT:
    // For a fixed stage length `len`, the butterflies operate independently on each
    // contiguous block of `len` coefficients. So we can parallelize over blocks via
    // `par_chunks_mut(len)` safely.
    //
    // Thresholds:
    // - for small sizes, rayon overhead outweighs benefit
    // - for large sizes (our RS extrapolation), this unlocks multi-core scaling
    const PAR_NTT_MIN_N: usize = 1 << 18; // 262k
    const PAR_NTT_MIN_LEN: usize = 1 << 10; // 1024
    let use_par = n >= PAR_NTT_MIN_N && rayon::current_num_threads() > 1;

    // Reuse twiddle buffer to avoid per-stage allocations.
    let mut twiddles: Vec<u32> = Vec::new();

    let mut len = 2usize;
    while len <= n {
        let wlen = if invert {
            pow_mod_u32(root_inv, (n / len) as u32, modulus)
        } else {
            pow_mod_u32(root, (n / len) as u32, modulus)
        };
        // Precompute twiddles for this stage: tw[j] = wlen^j.
        twiddles.resize(len / 2, 0u32);
        twiddles[0] = 1u32;
        for j in 1..(len / 2) {
            twiddles[j] = mul_mod_u32(twiddles[j - 1], wlen, modulus);
        }

        if use_par && len >= PAR_NTT_MIN_LEN {
            a.par_chunks_mut(len).for_each(|chunk| {
                debug_assert_eq!(chunk.len(), len);
                for j in 0..(len / 2) {
                    let u = chunk[j];
                    let v = mul_mod_u32(chunk[j + len / 2], twiddles[j], modulus);
                    chunk[j] = add_mod_u32(u, v, modulus);
                    chunk[j + len / 2] = sub_mod_u32(u, v, modulus);
                }
            });
        } else {
            for i in (0..n).step_by(len) {
                for j in 0..(len / 2) {
                    let u = a[i + j];
                    let v = mul_mod_u32(a[i + j + len / 2], twiddles[j], modulus);
                    a[i + j] = add_mod_u32(u, v, modulus);
                    a[i + j + len / 2] = sub_mod_u32(u, v, modulus);
                }
            }
        }
        len <<= 1;
    }

    if invert {
        let n_inv = inv_mod_u32(n as u32, modulus);
        if use_par {
            a.par_iter_mut().for_each(|x| {
                *x = mul_mod_u32(*x, n_inv, modulus);
            });
        } else {
            for x in a.iter_mut() {
                *x = mul_mod_u32(*x, n_inv, modulus);
            }
        }
    }
}

fn convolution_ntt_mod(a: &[u32], b: &[u32], out_len: usize, size: usize, modulus: u32, primitive_root: u32) -> Vec<u32> {
    let mut fa = vec![0u32; size];
    let mut fb = vec![0u32; size];
    fa[..a.len()].copy_from_slice(a);
    fb[..b.len()].copy_from_slice(b);

    ntt(&mut fa, false, modulus, primitive_root);
    ntt(&mut fb, false, modulus, primitive_root);
    for i in 0..size {
        fa[i] = mul_mod_u32(fa[i], fb[i], modulus);
    }
    ntt(&mut fa, true, modulus, primitive_root);
    fa.truncate(out_len);
    fa
}

/// CRT+NTT convolution fallback.
///
/// Uses 5 NTT-friendly 32-bit primes whose product is ~155 bits, enough to reconstruct
/// coefficients bounded by `O(n * p^2)` for our RS extrapolation use-case at ~1e6 scale.
fn convolution_crt_ntt<F: PrimeField>(a: &[F], b: &[F], out_len: usize, size: usize) -> Vec<F> {
    // IMPORTANT:
    // RS extrapolation wants a convolution length `size` which can be as large as 2^24 and beyond.
    // Hardcoding a small set of “popular” NTT primes is brittle: some support only up to 2^22/2^23.
    //
    // Instead, we *deterministically synthesize* a set of 32-bit primes p ≡ 1 (mod size) at runtime,
    // find a primitive root for each, and use them for CRT+NTT convolution.
    //
    // This keeps the fallback working as we change parameters, and avoids coupling correctness to
    // an ad-hoc modulus list.
    let (mods, roots) = ntt_moduli_for_size(size, /*target_count=*/ 8);
    let k = mods.len();
    if k < 4 {
        // For very large power-of-two sizes (e.g. 2^28), there may be too few 32-bit primes p ≡ 1 (mod size).
        // Fall back to a 64-bit CRT+NTT path which has many more valid primes and needs fewer moduli
        // to reach the required reconstruction bitwidth.
        return convolution_crt_ntt_u64::<F>(a, b, out_len, size);
    }

    // Extract modulus of F as u64 (Frog prime fits in u64).
    let p_bytes = F::MODULUS.to_bytes_le();
    let mut p_u64: u64 = 0;
    for (i, byte) in p_bytes.iter().enumerate().take(8) {
        p_u64 |= (*byte as u64) << (8 * i);
    }
    assert!(p_u64 > 1, "expected small (<=64-bit) prime modulus");

    // Map field elements to u64 in [0,p).
    let to_u64 = |x: &F| -> u64 {
        let xb = x.into_bigint().to_bytes_le();
        let mut v: u64 = 0;
        for (i, byte) in xb.iter().enumerate().take(8) {
            v |= (*byte as u64) << (8 * i);
        }
        v % p_u64
    };

    let a_u = a.par_iter().map(to_u64).collect::<Vec<_>>();
    let b_u = b.par_iter().map(to_u64).collect::<Vec<_>>();

    // For each modulus, reduce inputs and convolve.
    //
    // IMPORTANT (memory):
    // Each convolution allocates O(size) NTT buffers; parallelizing across moduli can easily
    // explode RSS for large sizes (e.g. 2^27, 2^28). We keep NTT *internally* parallel and
    // run moduli sequentially to cap memory.
    let residues_vec: Vec<Vec<u32>> = (0..k)
        .into_iter()
        .map(|i| {
            let modulus = mods[i];
            let root = roots[i];
            debug_assert!(
                ((modulus - 1) as usize) % size == 0,
                "internal: synthesized NTT modulus does not support size {size}"
            );
            let aa = a_u.iter().map(|&v| (v % (modulus as u64)) as u32).collect::<Vec<_>>();
            let bb = b_u.iter().map(|&v| (v % (modulus as u64)) as u32).collect::<Vec<_>>();
            convolution_ntt_mod(&aa, &bb, out_len, size, modulus, root)
        })
        .collect();

    let residues = residues_vec;

    // Precompute prefix products modulo each modulus and modulo p.
    let mut prefix_mod = vec![vec![0u32; k]; k];
    let mut inv_prefix = vec![0u32; k];
    for i in 0..k {
        let mi = mods[i];
        let mut prod: u64 = 1;
        for j in 0..i {
            prefix_mod[i][j] = (prod % (mi as u64)) as u32;
            prod = (prod * mods[j] as u64) % (mi as u64);
        }
        let prod_i = (prod % (mi as u64)) as u32;
        inv_prefix[i] = inv_mod_u32(prod_i, mi);
    }
    let mut prefix_p = vec![0u64; k];
    prefix_p[0] = 1;
    for i in 1..k {
        prefix_p[i] =
            (((prefix_p[i - 1] as u128) * (mods[i - 1] as u128)) % (p_u64 as u128)) as u64;
    }

    // Garner reconstruction per coefficient, then reduce mod p and map back to F.
    let out: Vec<F> = (0..out_len)
        .into_par_iter()
        .map(|t| {
            // mixed radix digits c[i] in [0, MODS[i])
            let mut c = vec![0u32; k];
            c[0] = residues[0][t];
            for i in 1..k {
                let mi = mods[i];
                let mut acc = residues[i][t];
                for j in 0..i {
                    let term = mul_mod_u32(c[j], prefix_mod[i][j], mi);
                    acc = sub_mod_u32(acc, term, mi);
                }
                c[i] = mul_mod_u32(acc, inv_prefix[i], mi);
            }

            // x mod p
            let mut x_mod_p: u64 = 0;
            for i in 0..k {
                let add = (((c[i] as u128) * (prefix_p[i] as u128)) % (p_u64 as u128)) as u64;
                x_mod_p = (((x_mod_p as u128) + (add as u128)) % (p_u64 as u128)) as u64;
            }
            F::from_le_bytes_mod_order(&x_mod_p.to_le_bytes())
        })
        .collect();

    out
}

// ---------------- CRT+NTT convolution for 64-bit NTT primes (fallback for huge sizes) ----------------

#[inline]
fn mul_mod_u64(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % (m as u128)) as u64
}

#[inline]
fn add_mod_u64(a: u64, b: u64, m: u64) -> u64 {
    let s = a as u128 + b as u128;
    (if s >= m as u128 { s - m as u128 } else { s }) as u64
}

#[inline]
fn sub_mod_u64(a: u64, b: u64, m: u64) -> u64 {
    if a >= b { a - b } else { a + m - b }
}

fn pow_mod_u64_wide(mut a: u64, mut e: u64, m: u64) -> u64 {
    let mut r: u64 = 1;
    while e > 0 {
        if (e & 1) == 1 {
            r = mul_mod_u64(r, a, m);
        }
        a = mul_mod_u64(a, a, m);
        e >>= 1;
    }
    r
}

fn inv_mod_u64(a: u64, m: u64) -> u64 {
    // m is prime
    pow_mod_u64_wide(a, m - 2, m)
}

fn bit_reverse_permute_u64(a: &mut [u64]) {
    let n = a.len();
    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            a.swap(i, j);
        }
    }
}

/// Find an n-th primitive root of unity in F_p for power-of-two n dividing (p-1).
fn root_of_unity_u64(p: u64, n: usize) -> u64 {
    debug_assert!(n.is_power_of_two());
    debug_assert!(((p - 1) as usize) % n == 0);
    let exp = (p - 1) / (n as u64);
    let half = n >> 1;
    for g in 2u64.. {
        let w = pow_mod_u64_wide(g, exp, p);
        if w == 1 {
            continue;
        }
        // Since n is a power-of-two, w has order n iff w^(n/2) != 1.
        if pow_mod_u64_wide(w, half as u64, p) != 1 {
            return w;
        }
    }
    unreachable!("failed to find root of unity for p={p} n={n}");
}

fn ntt_u64(a: &mut [u64], invert: bool, modulus: u64, root_n: u64) {
    let n = a.len();
    bit_reverse_permute_u64(a);

    let root = root_n;
    let root_inv = inv_mod_u64(root, modulus);

    const PAR_NTT_MIN_N: usize = 1 << 18;
    const PAR_NTT_MIN_LEN: usize = 1 << 10;
    let use_par = n >= PAR_NTT_MIN_N && rayon::current_num_threads() > 1;

    let mut twiddles: Vec<u64> = Vec::new();

    let mut len = 2usize;
    while len <= n {
        let wlen = if invert {
            pow_mod_u64_wide(root_inv, (n / len) as u64, modulus)
        } else {
            pow_mod_u64_wide(root, (n / len) as u64, modulus)
        };
        twiddles.resize(len / 2, 0u64);
        twiddles[0] = 1;
        for j in 1..(len / 2) {
            twiddles[j] = mul_mod_u64(twiddles[j - 1], wlen, modulus);
        }

        if use_par && len >= PAR_NTT_MIN_LEN {
            a.par_chunks_mut(len).for_each(|chunk| {
                for j in 0..(len / 2) {
                    let u = chunk[j];
                    let v = mul_mod_u64(chunk[j + len / 2], twiddles[j], modulus);
                    chunk[j] = add_mod_u64(u, v, modulus);
                    chunk[j + len / 2] = sub_mod_u64(u, v, modulus);
                }
            });
        } else {
            for i in (0..n).step_by(len) {
                for j in 0..(len / 2) {
                    let u = a[i + j];
                    let v = mul_mod_u64(a[i + j + len / 2], twiddles[j], modulus);
                    a[i + j] = add_mod_u64(u, v, modulus);
                    a[i + j + len / 2] = sub_mod_u64(u, v, modulus);
                }
            }
        }
        len <<= 1;
    }

    if invert {
        let n_inv = inv_mod_u64(n as u64, modulus);
        if use_par {
            a.par_iter_mut().for_each(|x| *x = mul_mod_u64(*x, n_inv, modulus));
        } else {
            for x in a.iter_mut() {
                *x = mul_mod_u64(*x, n_inv, modulus);
            }
        }
    }
}

fn convolution_ntt_mod_u64(a: &[u64], b: &[u64], out_len: usize, size: usize, modulus: u64, root_n: u64) -> Vec<u64> {
    let mut fa = vec![0u64; size];
    let mut fb = vec![0u64; size];
    fa[..a.len()].copy_from_slice(a);
    fb[..b.len()].copy_from_slice(b);

    ntt_u64(&mut fa, false, modulus, root_n);
    ntt_u64(&mut fb, false, modulus, root_n);
    for i in 0..size {
        fa[i] = mul_mod_u64(fa[i], fb[i], modulus);
    }
    ntt_u64(&mut fa, true, modulus, root_n);
    fa.truncate(out_len);
    fa
}

fn is_prime_u64(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    if n % 2 == 0 {
        return n == 2;
    }
    // Deterministic Miller-Rabin for u64.
    // See https://miller-rabin.appspot.com/ : these bases are sufficient for < 2^64.
    const BASES: [u64; 7] = [2, 325, 9375, 28178, 450775, 9780504, 1795265022];
    let d = n - 1;
    let s = d.trailing_zeros();
    let d_odd = d >> s;
    'outer: for &a in &BASES {
        let a = a % n;
        if a == 0 {
            continue;
        }
        let mut x = pow_mod_u64_wide(a, d_odd, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mul_mod_u64(x, x, n);
            if x == n - 1 {
                continue 'outer;
            }
        }
        return false;
    }
    true
}

fn ntt_moduli_for_size_u64(size: usize, target_count: usize) -> Vec<u64> {
    assert!(size.is_power_of_two());
    // Search primes of form p = m*size + 1 in the 64-bit range.
    // Start at a multiplier that makes p ~ 2^63 so we get large primes (better CRT bitwidth).
    let min_p: u128 = 1u128 << 63;
    let mut m: u128 = (min_p.saturating_sub(1) / (size as u128)) + 1;
    if m == 0 {
        m = 1;
    }
    let mut out = Vec::new();
    while out.len() < target_count {
        let p = m.saturating_mul(size as u128).saturating_add(1);
        if p >= (u64::MAX as u128) {
            break;
        }
        let p64 = p as u64;
        if is_prime_u64(p64) {
            out.push(p64);
        }
        m += 1;
    }
    out
}

/// Deterministically synthesize `target_count` 64-bit primes p such that `p ≡ 1 (mod size)`,
/// along with a size-th root of unity for each prime.
///
/// This is used as a fallback when there are too few 32-bit NTT primes for large `size` (e.g. 2^28).
fn ntt_moduli_roots_for_size_u64(size: usize, target_count: usize) -> (Vec<u64>, Vec<u64>) {
    assert!(size.is_power_of_two(), "NTT size must be power-of-two");
    assert!(size >= 2, "NTT size too small");

    // Cache per (size,target_count) so RS extrapolation doesn't repeatedly search primes/roots.
    static CACHE: OnceLock<Mutex<HashMap<(usize, usize), (Vec<u64>, Vec<u64>)>>> = OnceLock::new();
    let cache = CACHE.get_or_init(|| Mutex::new(HashMap::new()));
    if let Some(v) = cache
        .lock()
        .ok()
        .and_then(|m| m.get(&(size, target_count)).cloned())
    {
        return v;
    }

    let mods = ntt_moduli_for_size_u64(size, target_count);
    let roots_n: Vec<u64> = mods.iter().map(|&m| root_of_unity_u64(m, size)).collect();
    let out = (mods, roots_n);
    if let Ok(mut m) = cache.lock() {
        m.insert((size, target_count), out.clone());
    }
    out
}

fn convolution_crt_ntt_u64<F: PrimeField>(a: &[F], b: &[F], out_len: usize, size: usize) -> Vec<F> {
    // Pick 3 large 64-bit NTT primes (product ~192 bits) to cover large convolution coefficient bounds.
    let (mods, roots_n) = ntt_moduli_roots_for_size_u64(size, /*target_count=*/ 3);
    let k = mods.len();
    assert!(
        k >= 3,
        "CRT+NTT(u64) requires at least 3 NTT primes; got {k} for size={size}"
    );

    // Extract modulus of F as u64 (Frog/BabyBear-in-Frog fits in u64).
    let p_bytes = F::MODULUS.to_bytes_le();
    let mut p_u64: u64 = 0;
    for (i, byte) in p_bytes.iter().enumerate().take(8) {
        p_u64 |= (*byte as u64) << (8 * i);
    }
    assert!(p_u64 > 1, "expected small (<=64-bit) prime modulus");

    let to_u64 = |x: &F| -> u64 {
        let xb = x.into_bigint().to_bytes_le();
        let mut v: u64 = 0;
        for (i, byte) in xb.iter().enumerate().take(8) {
            v |= (*byte as u64) << (8 * i);
        }
        v % p_u64
    };
    let a_u = a.par_iter().map(to_u64).collect::<Vec<_>>();
    let b_u = b.par_iter().map(to_u64).collect::<Vec<_>>();

    // Compute roots and convolve under each modulus.
    // Same memory concern as the u32 CRT+NTT path: do not parallelize across moduli.
    let residues_vec: Vec<Vec<u64>> = (0..k)
        .into_iter()
        .map(|i| {
            let modulus = mods[i];
            let root_n = roots_n[i];
            let aa = a_u.iter().map(|&v| v % modulus).collect::<Vec<_>>();
            let bb = b_u.iter().map(|&v| v % modulus).collect::<Vec<_>>();
            convolution_ntt_mod_u64(&aa, &bb, out_len, size, modulus, root_n)
        })
        .collect();

    // Garner reconstruction mod p using 3 moduli.
    let m0 = mods[0];
    let m1 = mods[1];
    let m2 = mods[2];
    let m0_mod_m1 = (m0 % m1) as u64;
    let inv_m0_m1 = inv_mod_u64(m0_mod_m1, m1);
    let m01_mod_m2 = ((m0 as u128 * m1 as u128) % (m2 as u128)) as u64;
    let inv_m01_m2 = inv_mod_u64(m01_mod_m2, m2);

    let m0_mod_p = (m0 as u128 % (p_u64 as u128)) as u64;
    let m01_mod_p = ((m0 as u128 * m1 as u128) % (p_u64 as u128)) as u64;

    (0..out_len)
        .into_par_iter()
        .map(|idx| {
            let r0 = residues_vec[0][idx];
            let r1 = residues_vec[1][idx];
            let r2 = residues_vec[2][idx];

            // c0 = r0
            let c0 = r0;
            // c1 = (r1 - c0) * inv(m0) mod m1
            let t1 = sub_mod_u64(r1, c0 % m1, m1);
            let c1 = mul_mod_u64(t1, inv_m0_m1, m1);
            // x01 = c0 + m0*c1  (mod m2)
            let x01_mod_m2 = (c0 as u128 + (m0 as u128 * c1 as u128)) % (m2 as u128);
            let t2 = sub_mod_u64(r2, x01_mod_m2 as u64, m2);
            let c2 = mul_mod_u64(t2, inv_m01_m2, m2);

            // x mod p = c0 + m0*c1 + m0*m1*c2  (mod p)
            let mut x_mod_p: u64 = (c0 % p_u64) as u64;
            let add1 = ((m0_mod_p as u128) * (c1 as u128) % (p_u64 as u128)) as u64;
            x_mod_p = (((x_mod_p as u128) + (add1 as u128)) % (p_u64 as u128)) as u64;
            let add2 = ((m01_mod_p as u128) * (c2 as u128) % (p_u64 as u128)) as u64;
            x_mod_p = (((x_mod_p as u128) + (add2 as u128)) % (p_u64 as u128)) as u64;

            F::from_le_bytes_mod_order(&x_mod_p.to_le_bytes())
        })
        .collect()
}

/// Deterministically synthesize `target_count` 32-bit primes p such that `p ≡ 1 (mod size)`,
/// along with a primitive root for each prime.
///
/// This is used for the CRT+NTT convolution fallback when the base field `F` is not FFT-friendly.
fn ntt_moduli_for_size(size: usize, target_count: usize) -> (Vec<u32>, Vec<u32>) {
    assert!(size.is_power_of_two(), "NTT size must be power-of-two");
    assert!(size >= 2, "NTT size too small");

    // Cache per size so RS extrapolation doesn't repeatedly search primes.
    // This function is called from performance-critical code paths.
    static CACHE: OnceLock<Mutex<HashMap<usize, (Vec<u32>, Vec<u32>)>>> = OnceLock::new();
    let cache = CACHE.get_or_init(|| Mutex::new(HashMap::new()));
    if let Some(v) = cache.lock().ok().and_then(|m| m.get(&size).cloned()) {
        return v;
    }

    let mut mods = Vec::<u32>::new();
    let mut roots = Vec::<u32>::new();

    // Search primes of the form p = m*size + 1.
    //
    // NOTE:
    // - `size` is a power-of-two (thus even for all relevant sizes >= 2), so **every** `m`
    //   yields an odd candidate `p` (even * m + 1). Restricting to odd `m` unnecessarily discards
    //   half the valid candidates and can make large sizes (e.g. 2^27) fail to find enough primes.
    // - We intentionally start the search at a multiplier that makes `p >= 2^31` so the primes are
    //   genuinely “32-bit sized”. Using tiny primes can make CRT reconstruction ambiguous for the
    //   coefficient bounds we need in RS extrapolation.
    let min_p: u64 = 1u64 << 31;
    let mut m: u64 = ((min_p.saturating_sub(1)).saturating_div(size as u64)).saturating_add(1);
    if m == 0 {
        m = 1;
    }
    while mods.len() < target_count {
        let p_u64 = m.saturating_mul(size as u64).saturating_add(1);
        if p_u64 >= (u32::MAX as u64) {
            break;
        }
        let p = p_u64 as u32;
        if is_prime_u32(p) {
            let g = primitive_root_u32(p, size as u32)
                .unwrap_or_else(|| panic!("failed to find primitive root for prime p={p}"));
            mods.push(p);
            roots.push(g);
        }
        m += 1;
    }
    assert!(
        !mods.is_empty(),
        "failed to find any NTT primes for size={size}"
    );
    let out = (mods, roots);
    if let Ok(mut m) = cache.lock() {
        m.insert(size, out.clone());
    }
    out
}

fn is_prime_u32(n: u32) -> bool {
    if n < 2 {
        return false;
    }
    if n % 2 == 0 {
        return n == 2;
    }
    if n % 3 == 0 {
        return n == 3;
    }
    // Deterministic Miller-Rabin for 32-bit integers.
    // Bases {2,3,5,7,11} are sufficient for n < 2^32.
    const BASES: [u32; 5] = [2, 3, 5, 7, 11];
    let d = n - 1;
    let s = d.trailing_zeros();
    let d_odd = d >> s;
    'outer: for &a in &BASES {
        if a % n == 0 {
            continue;
        }
        let mut x = pow_mod_u64(a as u64, d_odd as u64, n as u64) as u64;
        if x == 1 || x == (n as u64 - 1) {
            continue;
        }
        for _ in 1..s {
            x = (x * x) % (n as u64);
            if x == (n as u64 - 1) {
                continue 'outer;
            }
        }
        return false;
    }
    true
}

fn pow_mod_u64(mut a: u64, mut e: u64, m: u64) -> u64 {
    let mut r: u64 = 1;
    while e > 0 {
        if (e & 1) == 1 {
            r = (r * a) % m;
        }
        a = (a * a) % m;
        e >>= 1;
    }
    r
}

/// Find a primitive root mod prime `p` (generator of F_p^*), returning it if found.
///
/// `size_pow2` is the target NTT size; we use it to make factoring `p-1` cheap (we know it has a
/// large 2-power factor), but the result is a *full* primitive root.
fn primitive_root_u32(p: u32, size_pow2: u32) -> Option<u32> {
    debug_assert!(p >= 3);
    debug_assert!(((p - 1) as usize) % (size_pow2 as usize) == 0);

    let mut factors: Vec<u32> = Vec::new();
    factors.push(2);

    // Factor the odd part of (p-1). Since p = m*size + 1 with size a large power-of-two,
    // the odd part m is typically small (<= 2^32 / size).
    let mut m = (p - 1) / size_pow2;
    let mut d = 3u32;
    while (d as u64) * (d as u64) <= (m as u64) {
        if m % d == 0 {
            factors.push(d);
            while m % d == 0 {
                m /= d;
            }
        }
        d += 2;
    }
    if m > 1 {
        factors.push(m);
    }
    factors.sort_unstable();
    factors.dedup();

    let p_u64 = p as u64;
    let order = (p - 1) as u64;
    for g in 2u32..(p - 1) {
        let g_u64 = g as u64;
        let mut ok = true;
        for &q in &factors {
            let e = order / (q as u64);
            if pow_mod_u64(g_u64, e, p_u64) == 1 {
                ok = false;
                break;
            }
        }
        if ok {
            return Some(g);
        }
    }
    None
}

/// Extrapolate a degree < n polynomial given its values at 0..n-1 to its values at n..2n-1.
///
/// This is the “next block” extrapolation for consecutive integer points, computed in
/// **O(n log n)** via one convolution over the field (requires `F: FftField`).
///
/// Input: `y[i] = f(i)` for i=0..n-1.
/// Output: `out[s] = f(n+s)` for s=0..n-1.
pub fn extrapolate_consecutive_next_block<F: FftField + PrimeField>(y: &[F]) -> Vec<F> {
    let n = y.len();
    if n == 0 {
        return Vec::new();
    }

    // factorials and inv factorials up to 2n-1
    let mut fact = vec![F::ONE; 2 * n];
    for i in 1..2 * n {
        fact[i] = fact[i - 1] * F::from(i as u64);
    }
    let mut inv_fact = vec![F::ONE; 2 * n];
    inv_fact[2 * n - 1] = fact[2 * n - 1].inverse().expect("fact nonzero");
    for i in (1..2 * n).rev() {
        inv_fact[i - 1] = inv_fact[i] * F::from(i as u64);
    }

    // c_i = y_i * (-1)^{n-1-i} / (i! (n-1-i)!)
    let mut c = vec![F::ZERO; n];
    for i in 0..n {
        let mut ci = y[i] * inv_fact[i] * inv_fact[n - 1 - i];
        if ((n - 1 - i) & 1) == 1 {
            ci = -ci;
        }
        c[i] = ci;
    }

    // b_t = 1/(t+1) for t=0..(2n-2), i.e. inv(1..2n-1)
    let mut b = (1..=(2 * n - 1)).map(|t| F::from(t as u64)).collect::<Vec<_>>();
    batch_inverse_in_place(&mut b);

    // conv = c * b; then S_s = conv[n-1+s]
    let conv = convolution::<F>(&c, &b);

    let mut out = vec![F::ZERO; n];
    for s in 0..n {
        // P(n+s) = (n+s)! / s!
        out[s] = fact[n + s] * inv_fact[s] * conv[n - 1 + s];
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::{Fp64, MontBackend, MontConfig};
    use rand::{RngCore, SeedableRng};
    use rand_chacha::ChaCha20Rng;

    // Frog prime: 15912092521325583641
    #[derive(MontConfig)]
    #[modulus = "15912092521325583641"]
    #[generator = "7"]
    pub struct FrogPrimeConfig;
    type FFrog = Fp64<MontBackend<FrogPrimeConfig, 1>>;

    #[test]
    fn test_convolution_crt_ntt_matches_naive_mod_p() {
        // Moderate sizes: big enough to exercise NTT path, small enough for naive O(n^2).
        let n = 2048usize;
        let m = 1536usize;
        let out_len = n + m - 1;
        let size = out_len.next_power_of_two();

        let mut rng = ChaCha20Rng::seed_from_u64(12345);
        let a = (0..n).map(|_| FFrog::from(rng.next_u64())).collect::<Vec<_>>();
        let b = (0..m).map(|_| FFrog::from(rng.next_u64())).collect::<Vec<_>>();

        let got = convolution_crt_ntt::<FFrog>(&a, &b, out_len, size);

        // Naive integer convolution mod p.
        let p = 15_912_092_521_325_583_641u64;
        let to_u64 = |x: &FFrog| -> u64 {
            let bytes = x.into_bigint().to_bytes_le();
            let mut v: u64 = 0;
            for (i, byte) in bytes.iter().enumerate().take(8) {
                v |= (*byte as u64) << (8 * i);
            }
            v % p
        };
        let a_u = a.iter().map(to_u64).collect::<Vec<_>>();
        let b_u = b.iter().map(to_u64).collect::<Vec<_>>();

        let mut exp_u = vec![0u64; out_len];
        for i in 0..n {
            let ai = a_u[i] as u128;
            for j in 0..m {
                exp_u[i + j] = ((exp_u[i + j] as u128 + ai * (b_u[j] as u128)) % (p as u128)) as u64;
            }
        }
        let exp = exp_u
            .iter()
            .map(|&v| FFrog::from_le_bytes_mod_order(&v.to_le_bytes()))
            .collect::<Vec<_>>();

        assert_eq!(got, exp);
    }
}


