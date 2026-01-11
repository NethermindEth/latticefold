//! WE/DPP gate arithmetization for LatticeFold+ (sparse dR1CS over a prime field).
//!
//! This module is a research/bench frontend: it arithmetizes the *verifier* computation,
//! keeping the relation log-scale in `n` and linear in the verifier-visible message sizes.

use ark_crypto_primitives::sponge::poseidon::PoseidonConfig;
use ark_ff::{BigInteger, Field, PrimeField};
use stark_rings::{psi, unit_monomial, CoeffRing, OverField, PolyRing, Zq};

use crate::recording_transcript::{PoseidonTraceOp as LfPoseidonTraceOp, PoseidonTranscriptTrace};
use crate::we_statement::WeParams;

// Reuse symphony’s sparse dR1CS primitives and Poseidon arithmetizer.
use symphony::dpp_poseidon::{
    merge_sparse_dr1cs_share_one_with_glue, poseidon_sponge_dr1cs_from_trace_with_wiring_and_bytes,
    Constraint, PoseidonByteWiring, PoseidonDr1csWiring, SparseDr1csInstance,
};
use symphony::dpp_sumcheck::Dr1csBuilder;
use symphony::dpp_sumcheck::{sumcheck_verify_degree3, RingVars};

/// Output of WE-gate arithmetization (single merged sparse dR1CS instance).
#[derive(Clone, Debug)]
pub struct WeDr1csOutput<F: PrimeField> {
    pub inst: SparseDr1csInstance<F>,
    pub assignment: Vec<F>,
    /// Number of public variables `l` (prefix of the assignment vector) intended as `x`.
    pub public_len: usize,
}

fn lf_ops_to_symphony_ops<F: PrimeField>(ops: &[LfPoseidonTraceOp<F>]) -> Vec<symphony::transcript::PoseidonTraceOp<F>> {
    ops.iter()
        .map(|op| match op {
            LfPoseidonTraceOp::Absorb(v) => symphony::transcript::PoseidonTraceOp::Absorb(v.clone()),
            LfPoseidonTraceOp::SqueezeField(v) => symphony::transcript::PoseidonTraceOp::SqueezeField(v.clone()),
            LfPoseidonTraceOp::SqueezeBytes { n, out } => symphony::transcript::PoseidonTraceOp::SqueezeBytes { n: *n, out: out.clone() },
        })
        .collect()
}

/// Enforce that each `get_challenge` re-absorb equals the corresponding `SqueezeField` output.
///
/// Our trace transcript records `get_challenge` as:
/// - `SqueezeField(out)`
/// - `Absorb(out)`  (Fiat–Shamir re-absorb)
fn enforce_reabsorb_equals_squeeze<F: PrimeField>(
    inst: &mut SparseDr1csInstance<F>,
    wiring: &PoseidonDr1csWiring,
    ops: &[symphony::transcript::PoseidonTraceOp<F>],
) -> Result<(), String> {
    let mut absorb_idx = 0usize;
    let mut squeeze_idx = 0usize;
    for op in ops {
        match op {
            symphony::transcript::PoseidonTraceOp::Absorb(_) => {
                absorb_idx += 1;
            }
            symphony::transcript::PoseidonTraceOp::SqueezeField(out) => {
                // Next op must be Absorb(out)
                // We enforce equality elementwise: absorb_var == squeeze_var.
                let (sq_start, sq_len) = wiring
                    .squeeze_field_ranges
                    .get(squeeze_idx)
                    .copied()
                    .ok_or("poseidon wiring squeeze_field_ranges oob")?;
                squeeze_idx += 1;
                if sq_len != out.len() {
                    return Err("poseidon squeeze length mismatch".to_string());
                }
                // IMPORTANT: `absorb_idx` tracks how many Absorb ops we've *already processed*.
                // The re-absorb corresponding to this squeeze is the *next* Absorb op in the trace,
                // i.e. it has index `absorb_idx` in `absorb_ranges`. Do NOT increment `absorb_idx`
                // here; it will be incremented when the loop reaches that Absorb op.
                let (ab_start, ab_len) = wiring
                    .absorb_ranges
                    .get(absorb_idx)
                    .copied()
                    .ok_or("poseidon wiring absorb_ranges oob after squeeze")?;
                if ab_len != sq_len {
                    return Err("poseidon reabsorb length mismatch".to_string());
                }
                for j in 0..sq_len {
                    let v_sq = wiring.squeeze_field_vars[sq_start + j];
                    let v_ab = wiring.absorb_vars[ab_start + j];
                    // (v_ab - v_sq) * 1 = 0
                    inst.constraints.push(Constraint {
                        a: vec![(F::ONE, v_ab), (-F::ONE, v_sq)],
                        b: vec![(F::ONE, 0)],
                        c: vec![(F::ZERO, 0)],
                    });
                }
            }
            symphony::transcript::PoseidonTraceOp::SqueezeBytes { .. } => {}
        }
    }
    Ok(())
}

type BF<R> = <<R as PolyRing>::BaseRing as Field>::BasePrimeField;

fn ring_to_ringvars<R>(
    b: &mut Dr1csBuilder<BF<R>>,
    x: &R,
) -> RingVars
where
    R: PolyRing,
    R::BaseRing: Field,
{
    let mut coeffs = Vec::with_capacity(R::dimension());
    for c in x.coeffs() {
        let fp = c.to_base_prime_field_elements().into_iter().next().unwrap();
        let v = b.new_var(fp);
        coeffs.push(v);
    }
    RingVars::new(coeffs)
}

fn scalar_to_ringvars<R>(
    b: &mut Dr1csBuilder<BF<R>>,
    x: BF<R>,
) -> RingVars
where
    R: PolyRing,
    R::BaseRing: Field,
{
    let d = R::dimension();
    let mut coeffs = Vec::with_capacity(d);
    let v0 = b.new_var(x);
    coeffs.push(v0);
    for _ in 1..d {
        let vz = b.new_var(BF::<R>::ZERO);
        b.enforce_var_eq_const(vz, BF::<R>::ZERO);
        coeffs.push(vz);
    }
    RingVars::new(coeffs)
}

fn scalar_var_to_ringvars<R>(
    b: &mut Dr1csBuilder<BF<R>>,
    x0: usize,
) -> RingVars
where
    R: PolyRing,
    R::BaseRing: Field,
{
    let d = R::dimension();
    let mut coeffs = Vec::with_capacity(d);
    coeffs.push(x0);
    for _ in 1..d {
        let vz = b.new_var(BF::<R>::ZERO);
        b.enforce_var_eq_const(vz, BF::<R>::ZERO);
        coeffs.push(vz);
    }
    RingVars::new(coeffs)
}

fn ring_add<F: PrimeField>(b: &mut Dr1csBuilder<F>, x: &RingVars, y: &RingVars) -> RingVars {
    assert_eq!(x.d(), y.d());
    let mut out = Vec::with_capacity(x.d());
    for i in 0..x.d() {
        let val = b.assignment[x.coeffs[i]] + b.assignment[y.coeffs[i]];
        let v = b.new_var(val);
        b.add_constraint(
            vec![(F::ONE, x.coeffs[i]), (F::ONE, y.coeffs[i])],
            vec![(F::ONE, b.one())],
            vec![(F::ONE, v)],
        );
        out.push(v);
    }
    RingVars::new(out)
}

fn ring_sub<F: PrimeField>(b: &mut Dr1csBuilder<F>, x: &RingVars, y: &RingVars) -> RingVars {
    assert_eq!(x.d(), y.d());
    let mut out = Vec::with_capacity(x.d());
    for i in 0..x.d() {
        let val = b.assignment[x.coeffs[i]] - b.assignment[y.coeffs[i]];
        let v = b.new_var(val);
        b.add_constraint(
            vec![(F::ONE, x.coeffs[i]), (-F::ONE, y.coeffs[i])],
            vec![(F::ONE, b.one())],
            vec![(F::ONE, v)],
        );
        out.push(v);
    }
    RingVars::new(out)
}

fn ring_scale<F: PrimeField>(b: &mut Dr1csBuilder<F>, x: &RingVars, s: usize) -> RingVars {
    let mut out = Vec::with_capacity(x.d());
    for i in 0..x.d() {
        let val = b.assignment[x.coeffs[i]] * b.assignment[s];
        let v = b.new_var(val);
        b.enforce_mul(x.coeffs[i], s, v);
        out.push(v);
    }
    RingVars::new(out)
}

fn ring_mul_negacyclic<F: PrimeField>(b: &mut Dr1csBuilder<F>, x: &RingVars, y: &RingVars) -> RingVars {
    // Negacyclic convolution mod (X^d + 1).
    let d = x.d();
    assert_eq!(d, y.d());
    let one = b.one();
    let mut out = Vec::with_capacity(d);
    for k in 0..d {
        // Build coefficient via sum of products; enforce via linear constraints over fresh vars.
        let mut acc_var = b.new_var(F::ZERO);
        b.enforce_var_eq_const(acc_var, F::ZERO);

        for i in 0..d {
            // j = k - i mod d
            let j = if i <= k { k - i } else { d + k - i };
            let sign = if i <= k { F::ONE } else { -F::ONE };
            let prod_val = b.assignment[x.coeffs[i]] * b.assignment[y.coeffs[j]];
            let prod = b.new_var(prod_val);
            b.enforce_mul(x.coeffs[i], y.coeffs[j], prod);
            // acc = acc + sign * prod
            let new_acc = b.new_var(b.assignment[acc_var] + sign * b.assignment[prod]);
            b.add_constraint(
                vec![(F::ONE, acc_var), (sign, prod)],
                vec![(F::ONE, one)],
                vec![(F::ONE, new_acc)],
            );
            acc_var = new_acc;
        }
        out.push(acc_var);
    }
    RingVars::new(out)
}

fn ring_eq<F: PrimeField>(b: &mut Dr1csBuilder<F>, x: &RingVars, y: &RingVars) {
    assert_eq!(x.d(), y.d());
    for i in 0..x.d() {
        b.enforce_lc_times_one_eq_const(vec![(F::ONE, x.coeffs[i]), (-F::ONE, y.coeffs[i])]);
    }
}

fn lc_to_var<F: PrimeField>(b: &mut Dr1csBuilder<F>, lc: Vec<(F, usize)>) -> usize {
    let val = lc
        .iter()
        .fold(F::ZERO, |acc, (c, idx)| acc + (*c * b.assignment[*idx]));
    let v = b.new_var(val);
    // lc * 1 = v
    b.add_constraint(lc, vec![(F::ONE, b.one())], vec![(F::ONE, v)]);
    v
}

fn enforce_lc_eq_var<F: PrimeField>(b: &mut Dr1csBuilder<F>, lc: Vec<(F, usize)>, v: usize) {
    b.add_constraint(lc, vec![(F::ONE, b.one())], vec![(F::ONE, v)]);
}

fn enforce_bool<F: PrimeField>(b: &mut Dr1csBuilder<F>, bit: usize) {
    // bit*(bit-1)=0
    b.add_constraint(
        vec![(F::ONE, bit)],
        vec![(F::ONE, bit), (-F::ONE, b.one())],
        vec![(F::ZERO, b.one())],
    );
}

fn const_var<F: PrimeField>(b: &mut Dr1csBuilder<F>, c: F) -> usize {
    let v = b.new_var(c);
    b.enforce_var_eq_const(v, c);
    v
}

/// Arithmetize one LF+ `short_challenge` coefficient:
///   coeff = (byte % u) - (u/2), where u is a power of two.
///
/// Returns a BF var holding `coeff`.
fn short_challenge_coeff_from_byte<F: PrimeField>(
    b: &mut Dr1csBuilder<F>,
    byte: usize,
    u: u64,
) -> usize {
    debug_assert!(u.is_power_of_two());
    debug_assert!(u <= 256);
    let half = (u / 2) as i64;
    let kbits = (u as f64).log2() as usize;
    debug_assert_eq!(1u64 << kbits, u);

    // Represent byte = r + u*q with:
    // - r = Σ_{i<k} r_i 2^i  (k bits)
    // - q = Σ_{j<8-k} q_j 2^j  ((8-k) bits)
    // Witness assignment for bits is derived from the current byte value.
    let byte_val_u64 = b.assignment[byte]
        .into_bigint()
        .to_bytes_le()
        .get(0)
        .copied()
        .unwrap_or(0) as u64;
    debug_assert!(byte_val_u64 < 256);
    let r_val = (byte_val_u64 % u) as u64;
    let q_val = (byte_val_u64 / u) as u64;

    let mut r_bits = Vec::with_capacity(kbits);
    let mut q_bits = Vec::with_capacity(8 - kbits);
    for i in 0..kbits {
        let bit = (r_val >> i) & 1;
        let v = b.new_var(F::from(bit));
        enforce_bool(b, v);
        r_bits.push(v);
    }
    for j in 0..(8 - kbits) {
        let bit = (q_val >> j) & 1;
        let v = b.new_var(F::from(bit));
        enforce_bool(b, v);
        q_bits.push(v);
    }

    let mut lc = vec![(-F::ONE, byte)];
    // r part
    for (i, &bi) in r_bits.iter().enumerate() {
        lc.push((F::from(1u64 << i), bi));
    }
    // q part
    for (j, &bj) in q_bits.iter().enumerate() {
        let coeff = u * (1u64 << j);
        lc.push((F::from(coeff), bj));
    }
    // Enforce: -byte + r + u*q == 0  => (-byte + ...)*1 = 0
    b.enforce_lc_times_one_eq_const(lc);

    // r = Σ r_i 2^i
    let mut r_lc = Vec::with_capacity(1 + r_bits.len());
    for (i, &bi) in r_bits.iter().enumerate() {
        r_lc.push((F::from(1u64 << i), bi));
    }
    let r = lc_to_var(b, r_lc);

    // coeff = r - half
    let _half = half; // keep for clarity: half = u/2, always positive here
    let coeff_val = b.assignment[r] - F::from((u / 2) as u64);
    let out = b.new_var(coeff_val);
    enforce_lc_eq_var(b, vec![(F::ONE, r), (-F::from((u / 2) as u64), b.one())], out);
    out
}

fn short_challenge_from_bytes<F: PrimeField>(
    b: &mut Dr1csBuilder<F>,
    bytes: &[usize],
    lambda: usize,
    ring_dim: usize,
) -> RingVars {
    debug_assert_eq!(bytes.len(), ring_dim);
    // Matches `utils::short_challenge`: u = 2^(lambda / d).
    let exp = (lambda / ring_dim) as u32;
    let u = 1u64 << exp;
    let mut coeffs = Vec::with_capacity(ring_dim);
    for &by in bytes {
        let c = short_challenge_coeff_from_byte::<F>(b, by, u);
        coeffs.push(c);
    }
    RingVars::new(coeffs)
}

#[derive(Clone, Debug)]
pub struct CmShortChallengeWiring {
    /// Byte variables (one per squeezed byte), in order.
    pub byte_vars: Vec<usize>,
    /// `s[0..3]` short challenges (ring elements as coefficient vars).
    pub s: Vec<RingVars>,
    /// Flattened `s_prime` of length `k*d` (ring elements as coefficient vars).
    pub s_prime_flat: Vec<RingVars>,
}

#[derive(Clone, Debug)]
pub struct CmFieldChallengeWiring {
    pub c0: Vec<usize>,
    pub c1: Vec<usize>,
    pub rc0: usize,
    pub rc1: usize,
    pub sumcheck_r0: Vec<usize>,
    pub sumcheck_r1: Vec<usize>,
}

#[derive(Clone, Debug)]
struct CmChallengeOpWiring {
    /// Poseidon `SqueezeBytes` op indices (in trace order) used for `s` and `s_prime`.
    squeeze_bytes_ops: Vec<usize>,
    /// Poseidon `SqueezeField` op indices (in trace order) used for `c0,c1,rc*,sumcheck r*`.
    squeeze_field_ops: Vec<usize>,
}

fn cm_challenge_op_wiring<R>(
    trace: &PoseidonTranscriptTrace<BF<R>>,
    k: usize,
    log_kappa: usize,
    nvars: usize,
) -> Result<CmChallengeOpWiring, String>
where
    R: PolyRing,
    R::BaseRing: Field,
{
    let d = R::dimension();
    let need_short = 3 + k * d;
    let need_field = 2 * log_kappa + 2 + 2 * nvars;

    let mut squeeze_bytes_ops = Vec::with_capacity(need_short);
    let mut squeeze_field_ops = Vec::with_capacity(need_field);

    let mut seen_first_bytes = false;
    let mut bytes_op_idx = 0usize;
    let mut field_op_idx = 0usize;

    for op in &trace.ops {
        match op {
            LfPoseidonTraceOp::SqueezeBytes { .. } => {
                if !seen_first_bytes {
                    seen_first_bytes = true;
                }
                if seen_first_bytes && squeeze_bytes_ops.len() < need_short {
                    squeeze_bytes_ops.push(bytes_op_idx);
                }
                bytes_op_idx += 1;
            }
            LfPoseidonTraceOp::SqueezeField(_) => {
                // Only start collecting field ops after we've collected all short challenges,
                // since `c0/c1/rc/...` come after `short_challenge` calls in CmProof::verify.
                if seen_first_bytes && squeeze_bytes_ops.len() == need_short && squeeze_field_ops.len() < need_field {
                    squeeze_field_ops.push(field_op_idx);
                }
                field_op_idx += 1;
            }
            _ => {}
        }
    }

    if squeeze_bytes_ops.len() != need_short {
        return Err(format!(
            "cm_challenge_op_wiring: need {} SqueezeBytes ops, saw {}",
            need_short,
            squeeze_bytes_ops.len()
        ));
    }
    if squeeze_field_ops.len() != need_field {
        return Err(format!(
            "cm_challenge_op_wiring: need {} SqueezeField ops, saw {}",
            need_field,
            squeeze_field_ops.len()
        ));
    }
    Ok(CmChallengeOpWiring {
        squeeze_bytes_ops,
        squeeze_field_ops,
    })
}

fn cm_poseidon_challenge_vars<R>(
    pose_wiring: &PoseidonDr1csWiring,
    byte_wiring: &PoseidonByteWiring,
    op_wiring: &CmChallengeOpWiring,
) -> Result<(Vec<usize>, Vec<usize>), String>
where
    R: PolyRing,
    R::BaseRing: Field,
{
    // Flatten byte vars in the order of short_challenges.
    let mut bytes = Vec::new();
    for &op_idx in &op_wiring.squeeze_bytes_ops {
        let (start, len) = *byte_wiring
            .squeeze_byte_ranges
            .get(op_idx)
            .ok_or("poseidon byte wiring squeeze_byte_ranges oob")?;
        bytes.extend_from_slice(&byte_wiring.squeeze_byte_vars[start..start + len]);
    }

    // Flatten field vars in the order we expect.
    let mut fields = Vec::new();
    for &op_idx in &op_wiring.squeeze_field_ops {
        let (start, len) = *pose_wiring
            .squeeze_field_ranges
            .get(op_idx)
            .ok_or("poseidon wiring squeeze_field_ranges oob")?;
        if len != 1 {
            return Err("expected base-field squeeze len=1".to_string());
        }
        fields.push(pose_wiring.squeeze_field_vars[start]);
    }
    Ok((bytes, fields))
}

fn bf_from_base_ring<R>(x: <R as PolyRing>::BaseRing) -> BF<R>
where
    R: PolyRing,
    R::BaseRing: Field,
{
    x.to_base_prime_field_elements()
        .into_iter()
        .next()
        .expect("base ring element has no base prime field elems")
}

fn scalar_one_minus<F: PrimeField>(b: &mut Dr1csBuilder<F>, x: usize) -> usize {
    let one = b.one();
    let v = b.new_var(F::ONE - b.assignment[x]);
    b.add_constraint(
        vec![(F::ONE, one), (-F::ONE, x)],
        vec![(F::ONE, one)],
        vec![(F::ONE, v)],
    );
    v
}

fn scalar_add<F: PrimeField>(b: &mut Dr1csBuilder<F>, x: usize, y: usize) -> usize {
    let v = b.new_var(b.assignment[x] + b.assignment[y]);
    b.add_constraint(
        vec![(F::ONE, x), (F::ONE, y)],
        vec![(F::ONE, b.one())],
        vec![(F::ONE, v)],
    );
    v
}

fn scalar_mul_const<F: PrimeField>(b: &mut Dr1csBuilder<F>, x: usize, c: F) -> usize {
    let v = b.new_var(b.assignment[x] * c);
    b.add_constraint(vec![(c, x)], vec![(F::ONE, b.one())], vec![(F::ONE, v)]);
    v
}

fn scalar_pow_table<F: PrimeField>(b: &mut Dr1csBuilder<F>, base: usize, max_exp: usize) -> Vec<usize> {
    let mut pows = Vec::with_capacity(max_exp + 1);
    let one = b.one();
    let v0 = b.new_var(F::ONE);
    b.enforce_var_eq_const(v0, F::ONE);
    pows.push(v0);
    for i in 0..max_exp {
        let next = b.new_var(b.assignment[pows[i]] * b.assignment[base]);
        b.enforce_mul(pows[i], base, next);
        pows.push(next);
    }
    debug_assert_eq!(pows[0], v0);
    debug_assert_eq!(b.assignment[one], F::ONE);
    pows
}

fn tensor_scalar_vars<F: PrimeField>(b: &mut Dr1csBuilder<F>, c: &[usize]) -> Vec<usize> {
    // Matches utils::tensor ordering: fold tensor_product with [1-c_i, c_i].
    let mut acc: Vec<usize> = vec![const_var(b, F::ONE)];
    for &ci in c {
        let a0 = scalar_one_minus(b, ci);
        let a1 = ci;
        let mut next = Vec::with_capacity(acc.len() * 2);
        for &t in &acc {
            // t*(1-ci)
            let v0 = b.new_var(b.assignment[t] * b.assignment[a0]);
            b.enforce_mul(t, a0, v0);
            next.push(v0);
            // t*ci
            let v1 = b.new_var(b.assignment[t] * b.assignment[a1]);
            b.enforce_mul(t, a1, v1);
            next.push(v1);
        }
        acc = next;
    }
    acc
}

fn eval_small_mle_ring<F: PrimeField>(
    b: &mut Dr1csBuilder<F>,
    evals: &[RingVars],
    r: &[usize],
) -> RingVars {
    // Direct translation of tensor_eval::eval_small_mle (skips zeros not needed here).
    assert!(!evals.is_empty(), "eval_small_mle_ring: empty evals");
    let d = evals[0].d();
    let mut sum = Vec::with_capacity(d);
    for _ in 0..d {
        let vz = b.new_var(F::ZERO);
        b.enforce_var_eq_const(vz, F::ZERO);
        sum.push(vz);
    }
    let mut sum = RingVars::new(sum);

    for (i, ev) in evals.iter().enumerate() {
        debug_assert_eq!(ev.d(), d);
        // eq weight
        let mut w = b.new_var(F::ONE);
        b.enforce_var_eq_const(w, F::ONE);
        for (j, &rj) in r.iter().enumerate() {
            let bit = (i >> j) & 1;
            let term = if bit == 1 {
                rj
            } else {
                scalar_one_minus(b, rj)
            };
            let new_w = b.new_var(b.assignment[w] * b.assignment[term]);
            b.enforce_mul(w, term, new_w);
            w = new_w;
        }
        let scaled = ring_scale(b, ev, w);
        sum = ring_add(b, &sum, &scaled);
    }
    sum
}

fn eval_t_z_optimized_ring<R>(
    b: &mut Dr1csBuilder<BF<R>>,
    c_z_scalars: &[usize], // BF vars, length log_kappa
    s_prime: &[RingVars],  // ring elems, length k*d
    d_prime_powers: &[RingVars], // ring elems, length ell
    x_powers: &[RingVars], // ring elems, length d
    r: &[usize],           // BF vars, length nvars
) -> RingVars
where
    R: PolyRing,
    R::BaseRing: Field,
{
    // tensor(c_z) in BF, then lift to ring scalars.
    let tensor_c = tensor_scalar_vars::<BF<R>>(b, c_z_scalars);
    let tensor_c_ring = tensor_c
        .iter()
        .copied()
        .map(|v| scalar_var_to_ringvars::<R>(b, v))
        .collect::<Vec<_>>();

    let kappa = tensor_c_ring.len();
    let sizes = [x_powers.len(), d_prime_powers.len(), s_prime.len(), kappa];
    // All pow2 is assumed in the protocol/bench; we still follow the optimized path.
    let vars4 = sizes.map(|s| ark_std::log2(s.next_power_of_two()) as usize);
    let tensor_vars = vars4.iter().sum::<usize>();

    // Split r into chunks (innermost to outermost) as in tensor_eval::eval_t_z_optimized.
    let r4 = &r[0..vars4[0]]; // x_powers (lowest bits)
    let r3 = &r[vars4[0]..vars4[0] + vars4[1]];
    let r2 = &r[vars4[0] + vars4[1]..vars4[0] + vars4[1] + vars4[2]];
    let r1 = &r[vars4[0] + vars4[1] + vars4[2]..tensor_vars];

    let v1 = eval_small_mle_ring::<BF<R>>(b, &tensor_c_ring, r1);
    let v2 = eval_small_mle_ring::<BF<R>>(b, s_prime, r2);
    let v3 = eval_small_mle_ring::<BF<R>>(b, d_prime_powers, r3);
    let v4 = eval_small_mle_ring::<BF<R>>(b, x_powers, r4);

    let mut res = ring_mul_negacyclic::<BF<R>>(b, &v1, &v2);
    res = ring_mul_negacyclic::<BF<R>>(b, &res, &v3);
    res = ring_mul_negacyclic::<BF<R>>(b, &res, &v4);

    // Padding factor: Π_{j=tensor_vars..nvars} (1 - r[j]) as BF scalar.
    let mut pad = b.new_var(BF::<R>::ONE);
    b.enforce_var_eq_const(pad, BF::<R>::ONE);
    for &rj in &r[tensor_vars..] {
        let om = scalar_one_minus::<BF<R>>(b, rj);
        let new_pad = b.new_var(b.assignment[pad] * b.assignment[om]);
        b.enforce_mul(pad, om, new_pad);
        pad = new_pad;
    }
    ring_scale::<BF<R>>(b, &res, pad)
}

#[cfg(feature = "we_gate")]
#[derive(Clone, Debug)]
struct CmMathWiring {
    short: CmShortChallengeWiring,
    field: CmFieldChallengeWiring,
    /// Flattened BF variables that must equal Poseidon absorb inputs (non-reabsorb absorbs),
    /// for the CmProof segment starting at `absorb_comh`.
    absorb_flat: Vec<usize>,
    /// Debug-only: cumulative constraint counts inside `cm_inst` after major phases.
    phase_marks: Vec<usize>,
    phase_names: Vec<String>,
}

#[cfg(feature = "we_gate")]
fn cm_verifier_math_dr1cs<R>(
    trace: &PoseidonTranscriptTrace<BF<R>>,
    proof: &crate::cm::CmProof<R>,
    k: usize,
    log_kappa: usize,
    nvars: usize,
    mlen_mats: usize,
) -> Result<
    (
        SparseDr1csInstance<BF<R>>,
        Vec<BF<R>>,
        CmMathWiring,
    ),
    String,
>
where
    R: OverField + CoeffRing + PolyRing,
    R::BaseRing: Zq + Field,
{
    use latticefold::utils::sumcheck::Proof as ScProof;

    let d = R::dimension();
    let l_instances = proof.evals.0.len();
    let ell = proof.dcom.dparams.l;

    if proof.sumcheck_proofs.0.msgs().len() != nvars || proof.sumcheck_proofs.1.msgs().len() != nvars {
        return Err("CmProof: sumcheck proof length mismatch".to_string());
    }

    let mut b = Dr1csBuilder::<BF<R>>::new();
    b.enforce_var_eq_const(b.one(), BF::<R>::ONE);

    fn push_phase<F: PrimeField>(names: &mut Vec<String>, marks: &mut Vec<usize>, name: String, b: &Dr1csBuilder<F>) {
        names.push(name);
        marks.push(b.rows.len());
    }

    let mut phase_marks: Vec<usize> = Vec::new();
    let mut phase_names: Vec<String> = Vec::new();

    // Extract the exact CmProof coin bytes (short_challenge) and field challenges from the trace,
    // so this part's witness assignment matches the Poseidon part (glue constraints).
    let need_short = 3 + k * d;
    let need_bytes = need_short * d;
    let need_field = 2 * log_kappa + 2 + 2 * nvars;

    let mut short_bytes_vals: Vec<u8> = Vec::with_capacity(need_bytes);
    let mut field_vals: Vec<BF<R>> = Vec::with_capacity(need_field);
    let mut seen_bytes_ops = 0usize;
    let mut seen_first_bytes = false;
    for op in &trace.ops {
        match op {
            LfPoseidonTraceOp::SqueezeBytes { out, .. } => {
                if !seen_first_bytes {
                    seen_first_bytes = true;
                }
                if seen_bytes_ops < need_short {
                    short_bytes_vals.extend_from_slice(out);
                    seen_bytes_ops += 1;
                }
            }
            LfPoseidonTraceOp::SqueezeField(v) => {
                if seen_first_bytes && seen_bytes_ops >= need_short && field_vals.len() < need_field {
                    if v.len() != 1 {
                        return Err("cm_verifier_math_dr1cs: expected base-field squeeze len=1".to_string());
                    }
                    field_vals.push(v[0]);
                }
            }
            _ => {}
        }
    }
    if short_bytes_vals.len() < need_bytes {
        return Err("cm_verifier_math_dr1cs: not enough squeeze-bytes for short challenges".to_string());
    }
    short_bytes_vals.truncate(need_bytes);
    if field_vals.len() != need_field {
        return Err("cm_verifier_math_dr1cs: not enough squeeze-field elements for cm challenges".to_string());
    }

    // --- Challenges (allocated locally; caller glues to coin/field wiring) ---
    // short challenges: s (3), s_prime_flat (k*d)
    let mut byte_vars = Vec::new();
    for &by in short_bytes_vals.iter() {
        byte_vars.push(b.new_var(BF::<R>::from(by as u64)));
    }
    let mut rings = Vec::with_capacity(need_short);
    for i in 0..need_short {
        let start = i * d;
        let end = start + d;
        let rv = short_challenge_from_bytes::<BF<R>>(&mut b, &byte_vars[start..end], 128, d);
        rings.push(rv);
    }
    let s = rings[0..3].to_vec();
    let s_prime_flat = rings[3..].to_vec();

    // field challenges: c0,c1,rc0,rc1,sumcheck r0,r1
    let mut cur = 0usize;
    let c0 = (0..log_kappa)
        .map(|_| {
            let v = b.new_var(field_vals[cur]);
            cur += 1;
            v
        })
        .collect::<Vec<_>>();
    let c1 = (0..log_kappa)
        .map(|_| {
            let v = b.new_var(field_vals[cur]);
            cur += 1;
            v
        })
        .collect::<Vec<_>>();
    let rc0 = {
        let v = b.new_var(field_vals[cur]);
        cur += 1;
        v
    };
    let sumcheck_r0 = (0..nvars)
        .map(|_| {
            let v = b.new_var(field_vals[cur]);
            cur += 1;
            v
        })
        .collect::<Vec<_>>();
    let rc1 = {
        let v = b.new_var(field_vals[cur]);
        cur += 1;
        v
    };
    let sumcheck_r1 = (0..nvars)
        .map(|_| {
            let v = b.new_var(field_vals[cur]);
            cur += 1;
            v
        })
        .collect::<Vec<_>>();
    debug_assert_eq!(cur, field_vals.len());

    let short_wiring = CmShortChallengeWiring {
        byte_vars,
        s,
        s_prime_flat,
    };
    let field_wiring = CmFieldChallengeWiring {
        c0,
        c1,
        rc0,
        rc1,
        sumcheck_r0,
        sumcheck_r1,
    };

    // Build the expected absorb surface for the CmProof segment.
    // This excludes all Poseidon-internal reabsorbs performed by `get_challenge`, which we already
    // constrain via `enforce_reabsorb_equals_squeeze` in the Poseidon part.
    let mut absorb_flat: Vec<usize> = Vec::new();

    // --- Witness: commitment surface `comh` (L × κ) ---
    let kappa = proof.comh[0].len();
    if kappa != (1usize << log_kappa) {
        return Err("CmProof: kappa/log_kappa mismatch".to_string());
    }
    if proof.comh.len() != l_instances {
        return Err("CmProof: comh length mismatch".to_string());
    }
    let mut comh_vars: Vec<Vec<RingVars>> = Vec::with_capacity(l_instances);
    for l in 0..l_instances {
        if proof.comh[l].len() != kappa {
            return Err("CmProof: comh inner len mismatch".to_string());
        }
        let mut row = Vec::with_capacity(kappa);
        for j in 0..kappa {
            let rv = ring_to_ringvars::<R>(&mut b, &proof.comh[l][j]);
            // `absorb_comh` absorbs each ring element in coefficient order.
            absorb_flat.extend_from_slice(&rv.coeffs);
            row.push(rv);
        }
        comh_vars.push(row);
    }
    push_phase(&mut phase_names, &mut phase_marks, "after_comh_absorb".to_string(), &b);

    // --- Compute u[l][*] from dcom.out.e and s_prime_flat ---
    // u[l] has length = dcom.out.e.len() (expected 1+Mlen).
    let e_sets = &proof.dcom.out.e;
    let mut u_vars: Vec<Vec<RingVars>> = Vec::with_capacity(l_instances);
    for l in 0..l_instances {
        let mut u_l = Vec::with_capacity(e_sets.len());
        for e_i in e_sets.iter() {
            if e_i.len() < (l + 1) * k {
                return Err("CmProof: dcom.out.e too short for L,k".to_string());
            }
            // Flatten k blocks (each Vec<R> of len d) -> length k*d.
            let mut flat: Vec<RingVars> = Vec::with_capacity(k * d);
            for block in e_i.iter().skip(l * k).take(k) {
                if block.len() != d {
                    return Err("CmProof: dcom.out.e block len != d".to_string());
                }
                for x in block {
                    flat.push(ring_to_ringvars::<R>(&mut b, x));
                }
            }
            if flat.len() != short_wiring.s_prime_flat.len() {
                return Err("CmProof: e_i flatten len mismatch with s_prime_flat".to_string());
            }
            // Σ_j flat[j] * s_prime_flat[j]
            let mut acc = scalar_to_ringvars::<R>(&mut b, BF::<R>::ZERO);
            for (uij, sij) in flat.iter().zip(short_wiring.s_prime_flat.iter()) {
                let prod = ring_mul_negacyclic::<BF<R>>(&mut b, uij, sij);
                acc = ring_add::<BF<R>>(&mut b, &acc, &prod);
            }
            u_l.push(acc);
        }
        u_vars.push(u_l);
    }
    push_phase(&mut phase_names, &mut phase_marks, "after_u_vars".to_string(), &b);

    // --- tensor(c0/c1) and tcch0/tcch1 ---
    let tensor_c0 = tensor_scalar_vars::<BF<R>>(&mut b, &field_wiring.c0);
    let tensor_c1 = tensor_scalar_vars::<BF<R>>(&mut b, &field_wiring.c1);
    if tensor_c0.len() != kappa || tensor_c1.len() != kappa {
        return Err("CmProof: tensor(c) len mismatch".to_string());
    }

    let mut tcch0: Vec<RingVars> = Vec::with_capacity(l_instances);
    let mut tcch1: Vec<RingVars> = Vec::with_capacity(l_instances);
    for l in 0..l_instances {
        let mut acc0 = scalar_to_ringvars::<R>(&mut b, BF::<R>::ZERO);
        let mut acc1 = scalar_to_ringvars::<R>(&mut b, BF::<R>::ZERO);
        for j in 0..kappa {
            let t0 = scalar_var_to_ringvars::<R>(&mut b, tensor_c0[j]);
            let t1 = scalar_var_to_ringvars::<R>(&mut b, tensor_c1[j]);
            let s0 = ring_mul_negacyclic::<BF<R>>(&mut b, &t0, &comh_vars[l][j]);
            let s1 = ring_mul_negacyclic::<BF<R>>(&mut b, &t1, &comh_vars[l][j]);
            acc0 = ring_add::<BF<R>>(&mut b, &acc0, &s0);
            acc1 = ring_add::<BF<R>>(&mut b, &acc1, &s1);
        }
        tcch0.push(acc0);
        tcch1.push(acc1);
    }
    push_phase(&mut phase_names, &mut phase_marks, "after_tcch".to_string(), &b);

    // --- Precompute constants for eval_t_z_optimized ---
    // dpp = [dp^i] as scalar ring elements (length ℓ = dparams.l)
    let dp = (R::dimension() / 2) as u64;
    let mut dpp = Vec::with_capacity(ell);
    let mut pow = BF::<R>::ONE;
    let dp_bf = BF::<R>::from(dp);
    for _ in 0..ell {
        dpp.push(scalar_to_ringvars::<R>(&mut b, pow));
        pow *= dp_bf;
    }
    // xp = unit monomials (length d)
    let mut xp = Vec::with_capacity(d);
    for i in 0..d {
        let mi = stark_rings::unit_monomial::<R>(i);
        xp.push(ring_to_ringvars::<R>(&mut b, &mi));
    }

    // --- Verify the two degree-2 sumchecks + recombination equality ---
    // Helper: parse one sumcheck proof msgs into [[RingVars;3]].
    let parse_sc_msgs = |b: &mut Dr1csBuilder<BF<R>>, p: &ScProof<R>| -> Result<Vec<[RingVars; 3]>, String> {
        let mut out = Vec::with_capacity(nvars);
        for m in p.msgs() {
            if m.evaluations.len() != 3 {
                return Err("CmProof: expected degree-2 evals (len=3)".to_string());
            }
            let e0 = ring_to_ringvars::<R>(b, &m.evaluations[0]);
            let e1 = ring_to_ringvars::<R>(b, &m.evaluations[1]);
            let e2 = ring_to_ringvars::<R>(b, &m.evaluations[2]);
            out.push([e0, e1, e2]);
        }
        Ok(out)
    };

    // Extract eval tables as RingVars: evals[which][l][j][t]
    let extract_evals = |b: &mut Dr1csBuilder<BF<R>>, evals: &[crate::cm::InstanceEvals<R>]| -> Result<Vec<Vec<[RingVars; 4]>>, String> {
        if evals.len() != l_instances {
            return Err("CmProof: evals length mismatch with L".to_string());
        }
        let mut out = Vec::with_capacity(l_instances);
        for l in 0..l_instances {
            let rows = evals[l].rows();
            let mut row = Vec::with_capacity(rows.len());
            for vals in rows {
                let v0 = ring_to_ringvars::<R>(b, &vals[0]);
                let v1 = ring_to_ringvars::<R>(b, &vals[1]);
                let v2 = ring_to_ringvars::<R>(b, &vals[2]);
                let v3 = ring_to_ringvars::<R>(b, &vals[3]);
                row.push([v0, v1, v2, v3]);
            }
            out.push(row);
        }
        Ok(out)
    };

    let sc0 = parse_sc_msgs(&mut b, &proof.sumcheck_proofs.0)?;
    let sc1 = parse_sc_msgs(&mut b, &proof.sumcheck_proofs.1)?;
    let evals0 = extract_evals(&mut b, &proof.evals.0)?;
    let evals1 = extract_evals(&mut b, &proof.evals.1)?;
    push_phase(&mut phase_names, &mut phase_marks, "after_parse_msgs_evals".to_string(), &b);

    // dcom evals for claimed_sum: per l, vectors of len 1+Mlen in (a,b,c)
    let mlen_chunks_usize = mlen_mats;
    let z_idx = l_instances * (4 + 4 * mlen_chunks_usize);
    let max_pow = z_idx + 1;

    // For each of the two sumchecks, compute:
    // - claimed_sum
    // - subclaim_eval via sumcheck_verify_degree2
    // - eval via recombination
    // and enforce equality.
    let do_one = |which: usize,
                  b: &mut Dr1csBuilder<BF<R>>,
                  absorb_flat: &mut Vec<usize>,
                  phase_names: &mut Vec<String>,
                  phase_marks: &mut Vec<usize>,
                  rc: usize,
                  r_sc: &[usize],
                  msgs: &[[RingVars; 3]],
                  evals: &[Vec<[RingVars; 4]>],
                  tcch0: &[RingVars],
                  tcch1: &[RingVars]| -> Result<(), String> {
        let tag = if which == 0 { "cm_sc0" } else { "cm_sc1" };
        // Sumcheck parameter block absorbed by the transcript.
        // NOTE: we assume base field (extension_degree=1), matching our Poseidon wiring usage.
        let v_nvars = const_var(b, BF::<R>::from(nvars as u64));
        let v_deg = const_var(b, BF::<R>::from(2u64));
        absorb_field_elem_as_ring::<R>(b, absorb_flat, v_nvars);
        absorb_field_elem_as_ring::<R>(b, absorb_flat, v_deg);

        // Per-round transcript absorbs:
        // - prover message evaluations (3 ring elems)
        // - then absorbs the sampled randomness scalar r_i
        for (round, m) in msgs.iter().enumerate() {
            absorb_flat.extend_from_slice(&m[0].coeffs);
            absorb_flat.extend_from_slice(&m[1].coeffs);
            absorb_flat.extend_from_slice(&m[2].coeffs);
            absorb_field_elem_as_ring::<R>(b, absorb_flat, r_sc[round]);
        }
        push_phase(phase_names, phase_marks, format!("{tag}:after_absorb_schedule"), b);

        let rc_pows = scalar_pow_table::<BF<R>>(b, rc, max_pow);
        let mut claimed_sum = scalar_to_ringvars::<R>(b, BF::<R>::ZERO);

        for (l, eval) in proof.dcom.evals.iter().enumerate() {
            let l_idx = l * (4 + 4 * mlen_chunks_usize);
            // a terms are scalars in base ring
            let a0 = b.new_var(bf_from_base_ring::<R>(eval.a[0]));
            let a0pow = scalar_mul::<BF<R>>(b, a0, rc_pows[l_idx]);
            let a0t = scalar_var_to_ringvars::<R>(b, a0pow);
            claimed_sum = ring_add::<BF<R>>(b, &claimed_sum, &a0t);

            // b/c are ring
            let b0 = ring_to_ringvars::<R>(b, &eval.b[0]);
            let c0 = ring_to_ringvars::<R>(b, &eval.c[0]);
            let t_b0 = ring_scale::<BF<R>>(b, &b0, rc_pows[l_idx + 1]);
            let t_c0 = ring_scale::<BF<R>>(b, &c0, rc_pows[l_idx + 2]);
            let t_u0 = ring_scale::<BF<R>>(b, &u_vars[l][0], rc_pows[l_idx + 3]);
            claimed_sum = ring_add::<BF<R>>(b, &claimed_sum, &t_b0);
            claimed_sum = ring_add::<BF<R>>(b, &claimed_sum, &t_c0);
            claimed_sum = ring_add::<BF<R>>(b, &claimed_sum, &t_u0);

            for i in 0..mlen_chunks_usize {
                let idx = l_idx + 4 + i * 4;
                let ai = b.new_var(bf_from_base_ring::<R>(eval.a[1 + i]));
                let aipow = scalar_mul::<BF<R>>(b, ai, rc_pows[idx]);
                let ai_t = scalar_var_to_ringvars::<R>(b, aipow);
                claimed_sum = ring_add::<BF<R>>(b, &claimed_sum, &ai_t);

                let bi = ring_to_ringvars::<R>(b, &eval.b[1 + i]);
                let ci = ring_to_ringvars::<R>(b, &eval.c[1 + i]);
                let t_bi = ring_scale::<BF<R>>(b, &bi, rc_pows[idx + 1]);
                claimed_sum = ring_add::<BF<R>>(b, &claimed_sum, &t_bi);
                let t_ci = ring_scale::<BF<R>>(b, &ci, rc_pows[idx + 2]);
                claimed_sum = ring_add::<BF<R>>(b, &claimed_sum, &t_ci);
                let t_ui = ring_scale::<BF<R>>(b, &u_vars[l][1 + i], rc_pows[idx + 3]);
                claimed_sum = ring_add::<BF<R>>(b, &claimed_sum, &t_ui);
            }

            let t_tcch0 = ring_scale::<BF<R>>(b, &tcch0[l], rc_pows[z_idx]);
            claimed_sum = ring_add::<BF<R>>(b, &claimed_sum, &t_tcch0);
            let t_tcch1 = ring_scale::<BF<R>>(b, &tcch1[l], rc_pows[z_idx + 1]);
            claimed_sum = ring_add::<BF<R>>(b, &claimed_sum, &t_tcch1);
        }
        push_phase(phase_names, phase_marks, format!("{tag}:after_claimed_sum"), b);

        let subclaim_eval = sumcheck_verify_degree2::<BF<R>>(b, claimed_sum, msgs, r_sc)?;
        push_phase(phase_names, phase_marks, format!("{tag}:after_sumcheck_verify"), b);

        // t(z) eval at ro (independent of l)
        let t0 = eval_t_z_optimized_ring::<R>(
            b,
            &field_wiring.c0,
            &short_wiring.s_prime_flat,
            &dpp,
            &xp,
            r_sc,
        );
        let t1 = eval_t_z_optimized_ring::<R>(
            b,
            &field_wiring.c1,
            &short_wiring.s_prime_flat,
            &dpp,
            &xp,
            r_sc,
        );

        // eq(r, ro) where r is dcom.out.r (base ring)
        let r_pre = proof
            .dcom
            .out
            .r
            .iter()
            .copied()
            .map(|x| b.new_var(bf_from_base_ring::<R>(x)))
            .collect::<Vec<_>>();
        let eq = eq_eval_vars::<BF<R>>(b, &r_pre, r_sc);
        let mut eval_acc = scalar_to_ringvars::<R>(b, BF::<R>::ZERO);
        push_phase(phase_names, phase_marks, format!("{tag}:after_eq_eval"), b);

        for l in 0..l_instances {
            let l_idx = l * (4 + 4 * mlen_chunks_usize);
            let mut inner = scalar_to_ringvars::<R>(b, BF::<R>::ZERO);
            // First group (tau,m_tau,f,h) is evals[l][0]
            let e00 = &evals[l][0][0];
            let e01 = &evals[l][0][1];
            let e02 = &evals[l][0][2];
            let e03 = &evals[l][0][3];
            let t_e00 = ring_scale::<BF<R>>(b, e00, rc_pows[l_idx]);
            inner = ring_add::<BF<R>>(b, &inner, &t_e00);
            let t_e01 = ring_scale::<BF<R>>(b, e01, rc_pows[l_idx + 1]);
            inner = ring_add::<BF<R>>(b, &inner, &t_e01);
            let t_e02 = ring_scale::<BF<R>>(b, e02, rc_pows[l_idx + 2]);
            inner = ring_add::<BF<R>>(b, &inner, &t_e02);
            let t_e03 = ring_scale::<BF<R>>(b, e03, rc_pows[l_idx + 3]);
            inner = ring_add::<BF<R>>(b, &inner, &t_e03);
            // M chunks
            for i in 0..mlen_chunks_usize {
                let idx = l_idx + 4 + i * 4;
                let Mi = &evals[l][1 + i];
                let t_m0 = ring_scale::<BF<R>>(b, &Mi[0], rc_pows[idx]);
                inner = ring_add::<BF<R>>(b, &inner, &t_m0);
                let t_m1 = ring_scale::<BF<R>>(b, &Mi[1], rc_pows[idx + 1]);
                inner = ring_add::<BF<R>>(b, &inner, &t_m1);
                let t_m2 = ring_scale::<BF<R>>(b, &Mi[2], rc_pows[idx + 2]);
                inner = ring_add::<BF<R>>(b, &inner, &t_m2);
                let t_m3 = ring_scale::<BF<R>>(b, &Mi[3], rc_pows[idx + 3]);
                inner = ring_add::<BF<R>>(b, &inner, &t_m3);
            }
            // eq * inner
            let eq_ring = scalar_var_to_ringvars::<R>(b, eq);
            let eq_inner = ring_mul_negacyclic::<BF<R>>(b, &eq_ring, &inner);
            eval_acc = ring_add::<BF<R>>(b, &eval_acc, &eq_inner);

            // Add t(z) terms (uses el[0][0])
            let t0e = ring_mul_negacyclic::<BF<R>>(b, &t0, e00);
            let t1e = ring_mul_negacyclic::<BF<R>>(b, &t1, e00);
            let t0e_s = ring_scale::<BF<R>>(b, &t0e, rc_pows[z_idx]);
            eval_acc = ring_add::<BF<R>>(b, &eval_acc, &t0e_s);
            let t1e_s = ring_scale::<BF<R>>(b, &t1e, rc_pows[z_idx + 1]);
            eval_acc = ring_add::<BF<R>>(b, &eval_acc, &t1e_s);
        }
        push_phase(phase_names, phase_marks, format!("{tag}:after_recompute_eval"), b);

        ring_eq::<BF<R>>(b, &subclaim_eval, &eval_acc);
        push_phase(phase_names, phase_marks, format!("{tag}:after_final_eq"), b);

        // After sumcheck verification, Cm verifier absorbs the per-instance eval tables.
        // (`absorb_evaluations(evals, transcript)`).
        for l in 0..l_instances {
            for row in &evals[l] {
                // Each row is [R; 4], absorbed in order.
                absorb_flat.extend_from_slice(&row[0].coeffs);
                absorb_flat.extend_from_slice(&row[1].coeffs);
                absorb_flat.extend_from_slice(&row[2].coeffs);
                absorb_flat.extend_from_slice(&row[3].coeffs);
            }
        }
        push_phase(phase_names, phase_marks, format!("{tag}:after_absorb_evals"), b);
        Ok(())
    };

    do_one(
        0,
        &mut b,
        &mut absorb_flat,
        &mut phase_names,
        &mut phase_marks,
        field_wiring.rc0,
        &field_wiring.sumcheck_r0,
        &sc0,
        &evals0,
        &tcch0,
        &tcch1,
    )?;
    push_phase(&mut phase_names, &mut phase_marks, "after_do_one_0".to_string(), &b);
    do_one(
        1,
        &mut b,
        &mut absorb_flat,
        &mut phase_names,
        &mut phase_marks,
        field_wiring.rc1,
        &field_wiring.sumcheck_r1,
        &sc1,
        &evals1,
        &tcch0,
        &tcch1,
    )?;
    push_phase(&mut phase_names, &mut phase_marks, "after_do_one_1".to_string(), &b);

    let (inst, asg) = b.into_instance();
    Ok((
        inst,
        asg,
        CmMathWiring {
            short: short_wiring,
            field: field_wiring,
            absorb_flat,
            phase_marks,
            phase_names,
        },
    ))
}

/// Build a dR1CS part that reconstructs all LF+ `short_challenge(128)` ring elements
/// from Poseidon `SqueezeBytes` outputs (bytes -> coefficients).
///
/// Assumption (holds for current LF+ verifier paths): the only `squeeze_bytes` calls in the
/// verifier transcript are from `short_challenge(128)` within `CmProof::verify`, and each call
/// squeezes exactly `R::dimension()` bytes.
fn cm_short_challenges_dr1cs<R>(
    trace: &PoseidonTranscriptTrace<BF<R>>,
    k: usize,
) -> Result<(SparseDr1csInstance<BF<R>>, Vec<BF<R>>, CmShortChallengeWiring), String>
where
    R: PolyRing,
    R::BaseRing: Field,
{
    let d = R::dimension();
    let lambda = 128usize;

    // Expected total number of short challenges consumed by CmProof::verify:
    // - s: 3
    // - s_prime: k*d
    let need = 3 + k * d;
    let need_bytes = need * d;
    if trace.squeezed_bytes.len() < need_bytes {
        return Err(format!(
            "cm_short_challenges_dr1cs: not enough squeezed bytes: need {}, got {}",
            need_bytes,
            trace.squeezed_bytes.len()
        ));
    }

    let mut b = Dr1csBuilder::<BF<R>>::new();
    b.enforce_var_eq_const(b.one(), BF::<R>::ONE);

    // Allocate byte vars (as field elements) for the needed bytes.
    let mut byte_vars = Vec::with_capacity(need_bytes);
    for &by in trace.squeezed_bytes.iter().take(need_bytes) {
        let v = const_var(&mut b, BF::<R>::from(by as u64));
        byte_vars.push(v);
    }

    // Reconstruct ring elements, chunking by d bytes.
    let mut rings: Vec<RingVars> = Vec::with_capacity(need);
    for i in 0..need {
        let start = i * d;
        let end = start + d;
        let rv = short_challenge_from_bytes::<BF<R>>(&mut b, &byte_vars[start..end], lambda, d);
        rings.push(rv);
    }

    let s = rings[0..3].to_vec();
    let s_prime_flat = rings[3..].to_vec();
    debug_assert_eq!(s_prime_flat.len(), k * d);

    let (inst, asg) = b.into_instance();
    Ok((
        inst,
        asg,
        CmShortChallengeWiring {
            byte_vars,
            s,
            s_prime_flat,
        },
    ))
}

/// Evaluate eq(c, r) where both are vectors of scalar (BF) variables.
fn eq_eval_vars<F: PrimeField>(b: &mut Dr1csBuilder<F>, c: &[usize], r: &[usize]) -> usize {
    assert_eq!(c.len(), r.len());
    let mut acc = b.new_var(F::ONE);
    b.enforce_var_eq_const(acc, F::ONE);
    for (&ci, &ri) in c.iter().zip(r.iter()) {
        let one = b.one();
        let one_minus_ci = b.new_var(F::ONE - b.assignment[ci]);
        b.add_constraint(vec![(F::ONE, one), (-F::ONE, ci)], vec![(F::ONE, one)], vec![(F::ONE, one_minus_ci)]);
        let one_minus_ri = b.new_var(F::ONE - b.assignment[ri]);
        b.add_constraint(vec![(F::ONE, one), (-F::ONE, ri)], vec![(F::ONE, one)], vec![(F::ONE, one_minus_ri)]);
        let ci_ri = b.new_var(b.assignment[ci] * b.assignment[ri]);
        b.enforce_mul(ci, ri, ci_ri);
        let om = b.new_var(b.assignment[one_minus_ci] * b.assignment[one_minus_ri]);
        b.enforce_mul(one_minus_ci, one_minus_ri, om);
        let t = b.new_var(b.assignment[ci_ri] + b.assignment[om]);
        b.add_constraint(vec![(F::ONE, ci_ri), (F::ONE, om)], vec![(F::ONE, one)], vec![(F::ONE, t)]);
        let new_acc = b.new_var(b.assignment[acc] * b.assignment[t]);
        b.enforce_mul(acc, t, new_acc);
        acc = new_acc;
    }
    acc
}

struct ChallengeCursor<'a, F: PrimeField> {
    vals: &'a [F],
    idx: usize,
    vars: Vec<usize>,
}

impl<'a, F: PrimeField> ChallengeCursor<'a, F> {
    fn new(vals: &'a [F]) -> Self {
        Self {
            vals,
            idx: 0,
            vars: Vec::new(),
        }
    }

    fn next(&mut self, b: &mut Dr1csBuilder<F>) -> usize {
        let v = self
            .vals
            .get(self.idx)
            .copied()
            .unwrap_or_else(|| panic!("challenge cursor oob at {}", self.idx));
        self.idx += 1;
        let var = b.new_var(v);
        self.vars.push(var);
        var
    }

    #[allow(dead_code)]
    fn consumed(&self) -> usize {
        self.idx
    }

    fn all_vars(&self) -> &[usize] {
        &self.vars
    }
}

fn comr1cs_verifier_math_dr1cs<R>(
    proof: &crate::r1cs::ComR1CSProof<R>,
    ch: &mut ChallengeCursor<BF<R>>,
) -> Result<(SparseDr1csInstance<BF<R>>, Vec<BF<R>>, Vec<usize>), String>
where
    R: OverField + PolyRing,
    R::BaseRing: Field,
{
    use latticefold::utils::sumcheck::Proof as ScProof;

    let mut b = Dr1csBuilder::<BF<R>>::new();
    b.enforce_var_eq_const(b.one(), BF::<R>::ONE);

    let nvars = proof.nvars;

    // r = transcript.get_challenges(nvars)
    let mut r_pre = Vec::with_capacity(nvars);
    for _ in 0..nvars {
        r_pre.push(ch.next(&mut b));
    }

    // Sumcheck verifier challenges (one per round)
    let mut r_sc = Vec::with_capacity(nvars);
    for _ in 0..nvars {
        r_sc.push(ch.next(&mut b));
    }

    // Sumcheck prover messages: per-round 4 ring elements.
    let msgs: &ScProof<R> = &proof.sumcheck_proof;
    if msgs.msgs().len() != nvars {
        return Err("ComR1CSProof: sumcheck proof length mismatch".to_string());
    }
    let mut msg_vars: Vec<[RingVars; 4]> = Vec::with_capacity(nvars);
    for m in msgs.msgs() {
        if m.evaluations.len() != 4 {
            return Err("ComR1CSProof: expected degree-3 evals (len=4)".to_string());
        }
        let e0 = ring_to_ringvars::<R>(&mut b, &m.evaluations[0]);
        let e1 = ring_to_ringvars::<R>(&mut b, &m.evaluations[1]);
        let e2 = ring_to_ringvars::<R>(&mut b, &m.evaluations[2]);
        let e3 = ring_to_ringvars::<R>(&mut b, &m.evaluations[3]);
        msg_vars.push([e0, e1, e2, e3]);
    }

    // Verify sumcheck with claimed sum = 0.
    let claimed_sum = scalar_to_ringvars::<R>(&mut b, BF::<R>::ZERO);
    let subclaim_eval = sumcheck_verify_degree3::<BF<R>>(&mut b, claimed_sum, &msg_vars, &r_sc)?;

    // Allocate evals absorbed by transcript (we need them for arithmetic check).
    let va = ring_to_ringvars::<R>(&mut b, &proof.va);
    let vb = ring_to_ringvars::<R>(&mut b, &proof.vb);
    let vc = ring_to_ringvars::<R>(&mut b, &proof.vc);

    // e = eq_eval(r_pre, r_sc) (scalar), lifted to a constant-coeff ring element.
    let e = eq_eval_vars::<BF<R>>(&mut b, &r_pre, &r_sc);
    let e_ring = scalar_var_to_ringvars::<R>(&mut b, e);

    // Enforce: e * (va*vb - vc) == subclaim_eval
    let vab = ring_mul_negacyclic::<BF<R>>(&mut b, &va, &vb);
    let diff = ring_sub::<BF<R>>(&mut b, &vab, &vc);
    let lhs = ring_mul_negacyclic::<BF<R>>(&mut b, &e_ring, &diff);
    ring_eq::<BF<R>>(&mut b, &lhs, &subclaim_eval);

    let (inst, asg) = b.into_instance();
    Ok((inst, asg, ch.all_vars().to_vec()))
}

fn scalar_sub_const<F: PrimeField>(b: &mut Dr1csBuilder<F>, r: usize, c: F) -> usize {
    let val = b.assignment[r] - c;
    let v = b.new_var(val);
    // v = r - c
    b.add_constraint(
        vec![(F::ONE, r), (-c, b.one())],
        vec![(F::ONE, b.one())],
        vec![(F::ONE, v)],
    );
    v
}

fn scalar_mul<F: PrimeField>(b: &mut Dr1csBuilder<F>, x: usize, y: usize) -> usize {
    let val = b.assignment[x] * b.assignment[y];
    let v = b.new_var(val);
    b.enforce_mul(x, y, v);
    v
}

fn scalar_sub<F: PrimeField>(b: &mut Dr1csBuilder<F>, x: usize, y: usize) -> usize {
    let val = b.assignment[x] - b.assignment[y];
    let v = b.new_var(val);
    // v = x - y
    b.add_constraint(
        vec![(F::ONE, x), (-F::ONE, y)],
        vec![(F::ONE, b.one())],
        vec![(F::ONE, v)],
    );
    v
}

fn count_squeezed_field_elems_before_first_squeeze_bytes<F: PrimeField>(
    ops: &[LfPoseidonTraceOp<F>],
) -> usize {
    let mut cnt = 0usize;
    for op in ops {
        match op {
            LfPoseidonTraceOp::SqueezeBytes { .. } => break,
            LfPoseidonTraceOp::SqueezeField(v) => cnt += v.len(),
            _ => {}
        }
    }
    cnt
}

fn ring_eval_at_scalar<R>(
    b: &mut Dr1csBuilder<BF<R>>,
    x: &RingVars,
    beta: usize,
) -> usize
where
    R: PolyRing,
    R::BaseRing: Field,
{
    // ev(x, beta) = Σ_{j=0..d-1} x_j * beta^j
    let d = x.d();
    let beta_pows = scalar_pow_table::<BF<R>>(b, beta, d.saturating_sub(1));
    let mut lc: Vec<(BF<R>, usize)> = Vec::with_capacity(d);
    for j in 0..d {
        lc.push((b.assignment[beta_pows[j]], x.coeffs[j]));
    }
    lc_to_var::<BF<R>>(b, lc)
}

fn absorb_field_elem_as_ring<R>(
    b: &mut Dr1csBuilder<BF<R>>,
    absorb_flat: &mut Vec<usize>,
    x0: usize,
) where
    R: PolyRing,
    R::BaseRing: Field,
{
    // Matches latticefold `Transcript::absorb_field_element` default:
    //   absorb(&From::from(*v)) i.e. absorb a constant-coeff ring element.
    let rv = scalar_var_to_ringvars::<R>(b, x0);
    absorb_flat.extend_from_slice(&rv.coeffs);
}

#[cfg(feature = "we_gate")]
#[derive(Clone, Debug)]
struct DcomPrefixMathWiring {
    /// Local vars for all Poseidon `SqueezeField` outputs used in the Dcom prefix, in order.
    squeeze_field_vars: Vec<usize>,
    /// Local vars for the 9 statement-bound params (same order as `WeParams::to_field_vec`).
    params_vars: Vec<usize>,
    /// Local vars for the extra statement-defined public inputs (e.g. SP1 public input digest).
    public_input_vars: Vec<usize>,
    /// Flattened absorb surface for the Dcom prefix excluding reabsorbs:
    /// - `Out::verify` absorbs (sumcheck params, prover msgs, verifier randomness absorbs, and `absorb_evaluations(e,b)`)
    /// - `rgchk::absorb_evaluations(dcom.evals)` absorbs (R::from(a_i) then c_i)
    absorb_flat: Vec<usize>,
}

#[cfg(feature = "we_gate")]
fn ct_psi_mul_ring<R>(b: &mut Dr1csBuilder<BF<R>>, x: &RingVars) -> usize
where
    R: OverField + CoeffRing + PolyRing,
    R::BaseRing: Zq + Field,
{
    // Compute ct(psi * x) as a BF-linear form in the coefficients of x.
    // This avoids implementing full ring multiplication in-circuit.
    let d = R::dimension();
    if x.d() != d {
        panic!("ct_psi_mul_ring: ring dimension mismatch");
    }

    let psi_r = psi::<R>();
    let mut acc = const_var(b, BF::<R>::ZERO);
    for j in 0..d {
        let basis = unit_monomial::<R>(j);
        let w_br = (psi_r * basis).ct();
        let w = bf_from_base_ring::<R>(w_br);
        if w == BF::<R>::ZERO {
            continue;
        }
        let term = scalar_mul_const::<BF<R>>(b, x.coeffs[j], w);
        acc = scalar_add::<BF<R>>(b, acc, term);
    }
    acc
}

#[cfg(feature = "we_gate")]
fn dcom_verifier_math_dr1cs<R>(
    dcom: &crate::rgchk::Dcom<R>,
    trace: &PoseidonTranscriptTrace<BF<R>>,
    params: &WeParams,
    public_inputs: &[BF<R>],
) -> Result<(SparseDr1csInstance<BF<R>>, Vec<BF<R>>, DcomPrefixMathWiring), String>
where
    R: OverField + CoeffRing + PolyRing,
    R::BaseRing: Zq + Field,
{
    use latticefold::utils::sumcheck::Proof as ScProof;

    let out = &dcom.out;
    let nvars = out.nvars;
    let nclaims = out.e[0].len() + out.b.len();
    let has_rc = out.e[0].len() > 1;

    // --- Setchk verifier math (same squeeze schedule as `Out::verify`) ---
    let expected_squeezes = nclaims * (nvars + 2) + if has_rc { 1 } else { 0 } + nvars;
    let prefix_squeezes = count_squeezed_field_elems_before_first_squeeze_bytes::<BF<R>>(&trace.ops);
    if prefix_squeezes != expected_squeezes {
        return Err(format!(
            "dcom/setchk: squeeze_field count mismatch before bytes: expected {}, trace has {}",
            expected_squeezes, prefix_squeezes
        ));
    }
    if trace.squeezed_field.len() < expected_squeezes {
        return Err("dcom/setchk: trace.squeezed_field too short".to_string());
    }

    let mut b = Dr1csBuilder::<BF<R>>::new();
    b.enforce_var_eq_const(b.one(), BF::<R>::ONE);

    // Allocate statement-bound params as local vars and enforce they match the circuit's expected constants.
    let params_vals = params.to_field_vec::<BF<R>>();
    if params_vals.len() != 9 {
        return Err("we params: expected 9 field elements".to_string());
    }
    let mut params_vars = Vec::with_capacity(9);
    for &v in &params_vals {
        params_vars.push(b.new_var(v));
    }
    // Enforce equality to the verifier/gadget constants (so changing public_x alone => UNSAT).
    b.enforce_var_eq_const(params_vars[0], BF::<R>::from(out.nvars as u64)); // nvars_setchk
    b.enforce_var_eq_const(params_vars[1], BF::<R>::from(3u64)); // degree_setchk
    b.enforce_var_eq_const(params_vars[2], BF::<R>::from(out.nvars as u64)); // nvars_cm (same nvars in this bench)
    b.enforce_var_eq_const(params_vars[3], BF::<R>::from(2u64)); // degree_cm
    b.enforce_var_eq_const(params_vars[4], BF::<R>::from(params.kappa)); // kappa (not otherwise used here)
    b.enforce_var_eq_const(params_vars[5], BF::<R>::from(R::dimension() as u64)); // ring_dim_d
    b.enforce_var_eq_const(params_vars[6], BF::<R>::from(dcom.dparams.k as u64)); // k
    b.enforce_var_eq_const(params_vars[7], BF::<R>::from(dcom.dparams.l as u64)); // l
    // mlen is enforced against the builder-provided Mlen in `build_we_dr1cs_for_cm_proof_debug`
    // via the `mlen_mats` argument; here we can at least ensure it is <= u64::MAX (already).
    b.enforce_var_eq_const(params_vars[8], BF::<R>::from(params.mlen)); // mlen

    // Allocate extra statement-defined public inputs as vars (not fixed), and absorb them
    // as field-elements-as-ring at the start of the transcript (proof-agnostic binding).
    let mut public_input_vars = Vec::with_capacity(public_inputs.len());
    for &x in public_inputs {
        public_input_vars.push(b.new_var(x));
    }

    // Allocate local squeeze vars with trace values.
    let mut squeeze_field_vars: Vec<usize> = Vec::with_capacity(expected_squeezes);
    for &v in trace.squeezed_field.iter().take(expected_squeezes) {
        squeeze_field_vars.push(b.new_var(v));
    }
    let mut cur = 0usize;
    let take = |cur: &mut usize, n: usize, xs: &[usize]| -> Vec<usize> {
        let out = xs[*cur..*cur + n].to_vec();
        *cur += n;
        out
    };

    // Parse cba challenges.
    let mut c_vars: Vec<Vec<usize>> = Vec::with_capacity(nclaims);
    let mut beta_vars: Vec<usize> = Vec::with_capacity(nclaims);
    let mut alpha_vars: Vec<usize> = Vec::with_capacity(nclaims);
    for _ in 0..nclaims {
        c_vars.push(take(&mut cur, nvars, &squeeze_field_vars));
        beta_vars.push(take(&mut cur, 1, &squeeze_field_vars)[0]);
        alpha_vars.push(take(&mut cur, 1, &squeeze_field_vars)[0]);
    }
    let rc_var = if has_rc {
        Some(take(&mut cur, 1, &squeeze_field_vars)[0])
    } else {
        None
    };
    let r_point = take(&mut cur, nvars, &squeeze_field_vars);
    debug_assert_eq!(cur, expected_squeezes);

    // Sumcheck prover messages: per-round 4 ring elements.
    let msgs: &ScProof<R> = &out.sumcheck_proof;
    if msgs.msgs().len() != nvars {
        return Err("dcom/setchk: sumcheck proof length mismatch".to_string());
    }

    let mut absorb_flat: Vec<usize> = Vec::new();
    // Public inputs absorbed before verification begins.
    for &v in &public_input_vars {
        absorb_field_elem_as_ring::<R>(&mut b, &mut absorb_flat, v);
    }
    // Sumcheck parameter block absorbs.
    // Use statement-bound params instead of hardcoded constants.
    absorb_field_elem_as_ring::<R>(&mut b, &mut absorb_flat, params_vars[0]); // nvars_setchk
    absorb_field_elem_as_ring::<R>(&mut b, &mut absorb_flat, params_vars[1]); // degree_setchk

    let mut msg_vars: Vec<[RingVars; 4]> = Vec::with_capacity(nvars);
    for (round, m) in msgs.msgs().iter().enumerate() {
        if m.evaluations.len() != 4 {
            return Err("dcom/setchk: expected degree-3 evals (len=4)".to_string());
        }
        let e0 = ring_to_ringvars::<R>(&mut b, &m.evaluations[0]);
        let e1 = ring_to_ringvars::<R>(&mut b, &m.evaluations[1]);
        let e2 = ring_to_ringvars::<R>(&mut b, &m.evaluations[2]);
        let e3 = ring_to_ringvars::<R>(&mut b, &m.evaluations[3]);
        absorb_flat.extend_from_slice(&e0.coeffs);
        absorb_flat.extend_from_slice(&e1.coeffs);
        absorb_flat.extend_from_slice(&e2.coeffs);
        absorb_flat.extend_from_slice(&e3.coeffs);
        absorb_field_elem_as_ring::<R>(&mut b, &mut absorb_flat, r_point[round]);
        msg_vars.push([e0, e1, e2, e3]);
    }

    // Verify sumcheck with claimed sum = 0.
    let z0 = const_var(&mut b, BF::<R>::ZERO);
    let ring_zero = scalar_var_to_ringvars::<R>(&mut b, z0);
    let v = sumcheck_verify_degree3::<BF<R>>(&mut b, ring_zero.clone(), &msg_vars, &r_point)?;

    // Absorb e/b evaluations, and keep their ring vars for rgchk checks.
    let mut out_e_vars: Vec<Vec<Vec<RingVars>>> = Vec::with_capacity(out.e.len());
    for ek in &out.e {
        let mut ek_vars: Vec<Vec<RingVars>> = Vec::with_capacity(ek.len());
        for ej in ek {
            let mut ej_vars: Vec<RingVars> = Vec::with_capacity(ej.len());
            for r in ej {
                let rv = ring_to_ringvars::<R>(&mut b, r);
                absorb_flat.extend_from_slice(&rv.coeffs);
                ej_vars.push(rv);
            }
            ek_vars.push(ej_vars);
        }
        out_e_vars.push(ek_vars);
    }
    let mut out_b_vars: Vec<RingVars> = Vec::with_capacity(out.b.len());
    for bb in &out.b {
        let rv = ring_to_ringvars::<R>(&mut b, bb);
        absorb_flat.extend_from_slice(&rv.coeffs);
        out_b_vars.push(rv);
    }

    // Compute verifier recombination `ver` (scalar-in-ring), and enforce ver == v.
    let rc_pow_base = rc_var.unwrap_or_else(|| const_var(&mut b, BF::<R>::ONE));
    let rc_pows = scalar_pow_table::<BF<R>>(&mut b, rc_pow_base, nclaims.saturating_sub(1));
    let mut ver_scalar = const_var(&mut b, BF::<R>::ZERO);
    b.enforce_var_eq_const(ver_scalar, BF::<R>::ZERO);

    // e[0] claims
    for i in 0..out.e[0].len() {
        let eq = eq_eval_vars::<BF<R>>(&mut b, &c_vars[i], &r_point);
        let beta = beta_vars[i];
        let alpha = alpha_vars[i];
        let beta2 = scalar_mul::<BF<R>>(&mut b, beta, beta);
        let alpha_pows =
            scalar_pow_table::<BF<R>>(&mut b, alpha, out.e[0][i].len().saturating_sub(1));

        let mut e_sum = const_var(&mut b, BF::<R>::ZERO);
        b.enforce_var_eq_const(e_sum, BF::<R>::ZERO);
        for j in 0..out.e[0][i].len() {
            let ejv = &out_e_vars[0][i][j];
            let ev1 = ring_eval_at_scalar::<R>(&mut b, ejv, beta);
            let ev2 = ring_eval_at_scalar::<R>(&mut b, ejv, beta2);
            let ev1_sq = scalar_mul::<BF<R>>(&mut b, ev1, ev1);
            let diff = scalar_sub::<BF<R>>(&mut b, ev1_sq, ev2);
            let term = scalar_mul::<BF<R>>(&mut b, diff, alpha_pows[j]);
            e_sum = scalar_add::<BF<R>>(&mut b, e_sum, term);
        }

        let t = scalar_mul::<BF<R>>(&mut b, eq, e_sum);
        let t = scalar_mul::<BF<R>>(&mut b, t, rc_pows[i]);
        ver_scalar = scalar_add::<BF<R>>(&mut b, ver_scalar, t);
    }

    // b claims
    for i in 0..out.b.len() {
        let offset = out.e[0].len();
        let idx = i + offset;
        let eq = eq_eval_vars::<BF<R>>(&mut b, &c_vars[idx], &r_point);
        let beta = beta_vars[idx];
        let alpha = alpha_vars[idx];
        let beta2 = scalar_mul::<BF<R>>(&mut b, beta, beta);

        let b_ring = &out_b_vars[i];
        let ev1 = ring_eval_at_scalar::<R>(&mut b, b_ring, beta);
        let ev2 = ring_eval_at_scalar::<R>(&mut b, b_ring, beta2);
        let ev1_sq = scalar_mul::<BF<R>>(&mut b, ev1, ev1);
        let b_claim = scalar_sub::<BF<R>>(&mut b, ev1_sq, ev2);

        let t = scalar_mul::<BF<R>>(&mut b, eq, alpha);
        let t = scalar_mul::<BF<R>>(&mut b, t, b_claim);
        let t = scalar_mul::<BF<R>>(&mut b, t, rc_pows[idx]);
        ver_scalar = scalar_add::<BF<R>>(&mut b, ver_scalar, t);
    }

    let ver_ring = scalar_var_to_ringvars::<R>(&mut b, ver_scalar);
    ring_eq::<BF<R>>(&mut b, &ver_ring, &v);

    // --- rgchk::Dcom::verify checks + absorb(dcom.evals) ---
    // Absorb dcom evals in transcript order (R::from(a_i) then c_i),
    // and enforce the psi/ct checks.
    let L = dcom.evals.len();
    let k = dcom.dparams.k;
    let d = R::dimension();
    let d_prime = d / 2;

    // (d')^i constants (BF vars) for u-comb.
    let mut dppow_vars: Vec<usize> = Vec::with_capacity(k);
    for i in 0..k {
        let base = R::BaseRing::from(d_prime as u128);
        let p_br = ark_ff::Field::pow(&base, [i as u64]);
        dppow_vars.push(const_var(&mut b, bf_from_base_ring::<R>(p_br)));
    }

    // Allocate eval witness vars and absorb.
    let mut eval_a_vars: Vec<Vec<usize>> = Vec::with_capacity(L);
    let mut eval_b_vars: Vec<Vec<RingVars>> = Vec::with_capacity(L);
    let mut eval_c_vars: Vec<Vec<RingVars>> = Vec::with_capacity(L);
    let mut eval_v_vars: Vec<Vec<usize>> = Vec::with_capacity(L);

    for eval in &dcom.evals {
        let mut a_l = Vec::with_capacity(eval.a.len());
        for &ai in &eval.a {
            let a_var = b.new_var(bf_from_base_ring::<R>(ai));
            a_l.push(a_var);
            let a_ring = scalar_var_to_ringvars::<R>(&mut b, a_var);
            absorb_flat.extend_from_slice(&a_ring.coeffs);
        }
        eval_a_vars.push(a_l);

        let mut b_l = Vec::with_capacity(eval.b.len());
        for bi in &eval.b {
            b_l.push(ring_to_ringvars::<R>(&mut b, bi));
        }
        eval_b_vars.push(b_l);

        let mut c_l = Vec::with_capacity(eval.c.len());
        for ci in &eval.c {
            let rv = ring_to_ringvars::<R>(&mut b, ci);
            absorb_flat.extend_from_slice(&rv.coeffs);
            c_l.push(rv);
        }
        eval_c_vars.push(c_l);

        let mut v_l = Vec::with_capacity(eval.v.len());
        for &vj in &eval.v {
            v_l.push(b.new_var(bf_from_base_ring::<R>(vj)));
        }
        eval_v_vars.push(v_l);
    }

    // ct(psi * b_i) == a_i
    for l in 0..L {
        if eval_a_vars[l].len() != eval_b_vars[l].len() {
            return Err("dcom/rgchk: eval.a/eval.b length mismatch".to_string());
        }
        for i in 0..eval_a_vars[l].len() {
            let ct = ct_psi_mul_ring::<R>(&mut b, &eval_b_vars[l][i]);
            let diff = scalar_sub::<BF<R>>(&mut b, ct, eval_a_vars[l][i]);
            b.enforce_var_eq_const(diff, BF::<R>::ZERO);
        }
    }

    // ct(psi * (Σ_i (d')^i * u_i[col])) == expected[col]
    for l in 0..L {
        let base = l * k;
        for ni in 0..out.e.len() {
            for col in 0..d {
                let mut acc = ring_zero.clone();
                for i in 0..k {
                    let ui_col = &out_e_vars[ni][base + i][col];
                    let scaled = ring_scale::<BF<R>>(&mut b, ui_col, dppow_vars[i]);
                    acc = ring_add::<BF<R>>(&mut b, &acc, &scaled);
                }
                let ct = ct_psi_mul_ring::<R>(&mut b, &acc);
                let expected = if ni == 0 {
                    *eval_v_vars[l]
                        .get(col)
                        .ok_or("dcom/rgchk: eval.v length mismatch")?
                } else {
                    *eval_c_vars[l]
                        .get(ni)
                        .ok_or("dcom/rgchk: eval.c length mismatch")?
                        .coeffs
                        .get(col)
                        .ok_or("dcom/rgchk: ring coeff index oob")?
                };
                let diff = scalar_sub::<BF<R>>(&mut b, ct, expected);
                b.enforce_var_eq_const(diff, BF::<R>::ZERO);
            }
        }
    }

    let (inst, asg) = b.into_instance();
    Ok((
        inst,
        asg,
        DcomPrefixMathWiring {
            squeeze_field_vars,
            params_vars,
            public_input_vars,
            absorb_flat,
        },
    ))
}

fn lagrange_degree2<F: PrimeField>(b: &mut Dr1csBuilder<F>, r: usize) -> (usize, usize, usize) {
    let inv2 = F::from(2u64).inverse().unwrap();

    // t1=r-1, t2=r-2
    let t1 = scalar_sub_const(b, r, F::ONE);
    let t2 = scalar_sub_const(b, r, F::from(2u64));

    // L0 = (r-1)(r-2)/2
    let p = scalar_mul(b, t1, t2);
    let l0 = b.new_var(b.assignment[p] * inv2);
    b.add_constraint(vec![(inv2, p)], vec![(F::ONE, b.one())], vec![(F::ONE, l0)]);

    // L1 = -r(r-2)
    let p = scalar_mul(b, r, t2);
    let l1 = b.new_var(-b.assignment[p]);
    b.add_constraint(vec![(-F::ONE, p)], vec![(F::ONE, b.one())], vec![(F::ONE, l1)]);

    // L2 = r(r-1)/2
    let p = scalar_mul(b, r, t1);
    let l2 = b.new_var(b.assignment[p] * inv2);
    b.add_constraint(vec![(inv2, p)], vec![(F::ONE, b.one())], vec![(F::ONE, l2)]);

    (l0, l1, l2)
}

pub fn sumcheck_verify_degree2<F: PrimeField>(
    b: &mut Dr1csBuilder<F>,
    mut claimed_sum: RingVars,
    msgs: &[[RingVars; 3]],
    rs: &[usize],
) -> Result<RingVars, String> {
    if msgs.len() != rs.len() {
        return Err("sumcheck_verify_degree2: msgs/rs length mismatch".to_string());
    }
    for (m, &r) in msgs.iter().zip(rs.iter()) {
        // Check g(0)+g(1) == claimed_sum (coefficient-wise).
        let g01 = ring_add(b, &m[0], &m[1]);
        ring_eq(b, &g01, &claimed_sum);

        // Update claim = g(r) by Lagrange interpolation on points 0,1,2.
        let (l0, l1, l2) = lagrange_degree2::<F>(b, r);
        let t0 = ring_scale(b, &m[0], l0);
        let t1 = ring_scale(b, &m[1], l1);
        let t2 = ring_scale(b, &m[2], l2);
        let s01 = ring_add(b, &t0, &t1);
        claimed_sum = ring_add(b, &s01, &t2);
    }
    Ok(claimed_sum)
}

/// WE-gate arithmetization for verifying one `ComR1CSProof` (the Π_lin proof).
///
/// This is a first “end-to-end inside WE” building block: it includes
/// - Poseidon transcript arithmetization
/// - FS re-absorb consistency constraints
/// - statement-bound params prefix (public inputs)
/// - Π_lin verifier arithmetic constraints
/// - glue constraints equating Π_lin challenge variables with Poseidon squeeze outputs
///
/// NOTE: This currently only covers the `ComR1CSProof` verifier path (see `r1cs.rs`).
#[cfg(feature = "we_gate")]
pub fn build_we_dr1cs_for_comr1cs_proof<R>(
    poseidon_cfg: &PoseidonConfig<BF<R>>,
    trace: &PoseidonTranscriptTrace<BF<R>>,
    params: &WeParams,
    proof: &crate::r1cs::ComR1CSProof<R>,
) -> Result<WeDr1csOutput<BF<R>>, String>
where
    R: OverField + CoeffRing + PolyRing,
    R::BaseRing: Zq + Field,
{
    // Poseidon trace -> dR1CS
    let ops = lf_ops_to_symphony_ops::<BF<R>>(&trace.ops);
    let (mut pose_inst, pose_asg, _replay, _byte_wit, wiring, _byte_wiring) =
        poseidon_sponge_dr1cs_from_trace_with_wiring_and_bytes::<BF<R>>(poseidon_cfg, &ops)
            .map_err(|e| format!("poseidon arith failed: {e}"))?;
    enforce_reabsorb_equals_squeeze::<BF<R>>(&mut pose_inst, &wiring, &ops)?;

    // Public statement params prefix (no constraints fixing their value).
    let mut b_params = Dr1csBuilder::<BF<R>>::new();
    b_params.enforce_var_eq_const(b_params.one(), BF::<R>::ONE);
    for &x in &params.to_field_vec::<BF<R>>() {
        let _ = b_params.new_var(x);
    }
    let (params_inst, params_asg) = b_params.into_instance();

    // Π_lin verifier arithmetic.
    let mut ch = ChallengeCursor::<BF<R>>::new(&trace.squeezed_field);
    let (lin_inst, lin_asg, lin_ch_vars) = comr1cs_verifier_math_dr1cs::<R>(proof, &mut ch)?;

    // Glue: each challenge var equals corresponding Poseidon squeeze-field var.
    if wiring.squeeze_field_vars.len() < lin_ch_vars.len() {
        return Err("poseidon wiring: not enough squeeze_field_vars for lin challenges".to_string());
    }
    let mut glue: Vec<(usize, usize, usize, usize)> = Vec::with_capacity(lin_ch_vars.len());
    for (i, &v_lin) in lin_ch_vars.iter().enumerate() {
        let v_pose = wiring.squeeze_field_vars[i];
        glue.push((0, v_pose, 2, v_lin));
    }

    let parts = vec![(pose_inst, pose_asg), (params_inst, params_asg), (lin_inst, lin_asg)];
    let (inst, assignment) =
        merge_sparse_dr1cs_share_one_with_glue(&parts, &glue).map_err(|e| e.to_string())?;

    // Public prefix: [1] + params (fixed 9 scalars)
    let public_len = 1 + 9;
    Ok(WeDr1csOutput {
        inst,
        assignment,
        public_len,
    })
}

/// WE-gate arithmetization for the **short-challenge derivation** portion of `CmProof::verify`.
///
/// Builds:
/// - Poseidon transcript arithmetization (+ byte wiring)
/// - constraints that reconstruct `short_challenge(128)` ring elements from those bytes
/// - glue constraints equating byte variables across the two parts
///
/// Returns both the merged dR1CS output and the wiring for `s` and `s_prime`.
#[cfg(feature = "we_gate")]
pub fn build_we_dr1cs_for_cm_short_challenges<R>(
    poseidon_cfg: &PoseidonConfig<BF<R>>,
    trace: &PoseidonTranscriptTrace<BF<R>>,
    params: &WeParams,
    k: usize,
) -> Result<(WeDr1csOutput<BF<R>>, CmShortChallengeWiring), String>
where
    R: OverField + CoeffRing + PolyRing,
    R::BaseRing: Zq + Field,
{
    // Poseidon trace -> dR1CS (+ wiring with squeeze-field + squeeze-byte var indices).
    let ops = lf_ops_to_symphony_ops::<BF<R>>(&trace.ops);
    let (mut pose_inst, pose_asg, _replay, _byte_wit, wiring, byte_wiring) =
        poseidon_sponge_dr1cs_from_trace_with_wiring_and_bytes::<BF<R>>(poseidon_cfg, &ops)
            .map_err(|e| format!("poseidon arith failed: {e}"))?;
    enforce_reabsorb_equals_squeeze::<BF<R>>(&mut pose_inst, &wiring, &ops)?;

    // Public statement params prefix (no constraints fixing their value).
    let mut b_params = Dr1csBuilder::<BF<R>>::new();
    b_params.enforce_var_eq_const(b_params.one(), BF::<R>::ONE);
    for &x in &params.to_field_vec::<BF<R>>() {
        let _ = b_params.new_var(x);
    }
    let (params_inst, params_asg) = b_params.into_instance();

    // Short-challenge reconstruction part (allocates its own byte vars; we glue them).
    let (coin_inst, coin_asg, coin_wiring) = cm_short_challenges_dr1cs::<R>(trace, k)?;

    // Glue all squeezed bytes in order.
    if byte_wiring.squeeze_byte_vars.len() < coin_wiring.byte_vars.len() {
        return Err("poseidon byte wiring: not enough squeeze_byte_vars".to_string());
    }
    let mut glue: Vec<(usize, usize, usize, usize)> = Vec::with_capacity(coin_wiring.byte_vars.len());
    for i in 0..coin_wiring.byte_vars.len() {
        glue.push((0, byte_wiring.squeeze_byte_vars[i], 2, coin_wiring.byte_vars[i]));
    }

    let parts = vec![(pose_inst, pose_asg), (params_inst, params_asg), (coin_inst, coin_asg)];
    let (inst, assignment) =
        merge_sparse_dr1cs_share_one_with_glue(&parts, &glue).map_err(|e| e.to_string())?;

    let public_len = 1 + 9;
    Ok((WeDr1csOutput { inst, assignment, public_len }, coin_wiring))
}

/// Build a WE dR1CS instance that binds **all transcript coins** used by `CmProof::verify`:
/// - `short_challenge(128)` ring elements (`s`, `s_prime`)
/// - `get_challenge` scalars (`c0,c1`, `rc0,rc1`, and per-round sumcheck `r`s)
///
/// This does **not** yet add the Cm verifier arithmetic constraints; it just provides
/// a properly-wired coin surface to use in subsequent steps.
#[cfg(feature = "we_gate")]
pub fn build_we_dr1cs_for_cm_challenges<R>(
    poseidon_cfg: &PoseidonConfig<BF<R>>,
    trace: &PoseidonTranscriptTrace<BF<R>>,
    params: &WeParams,
    k: usize,
    log_kappa: usize,
    nvars: usize,
) -> Result<(WeDr1csOutput<BF<R>>, CmShortChallengeWiring, CmFieldChallengeWiring), String>
where
    R: OverField + CoeffRing + PolyRing,
    R::BaseRing: Zq + Field,
{
    let ops = lf_ops_to_symphony_ops::<BF<R>>(&trace.ops);
    let (mut pose_inst, pose_asg, _replay, _byte_wit, pose_wiring, byte_wiring) =
        poseidon_sponge_dr1cs_from_trace_with_wiring_and_bytes::<BF<R>>(poseidon_cfg, &ops)
            .map_err(|e| format!("poseidon arith failed: {e}"))?;
    enforce_reabsorb_equals_squeeze::<BF<R>>(&mut pose_inst, &pose_wiring, &ops)?;

    // Public statement params prefix.
    let mut b_params = Dr1csBuilder::<BF<R>>::new();
    b_params.enforce_var_eq_const(b_params.one(), BF::<R>::ONE);
    for &x in &params.to_field_vec::<BF<R>>() {
        let _ = b_params.new_var(x);
    }
    let (params_inst, params_asg) = b_params.into_instance();

    // Short challenges part (bytes -> ring coeffs).
    let (coin_inst, coin_asg, coin_wiring) = cm_short_challenges_dr1cs::<R>(trace, k)?;
    let op_wiring = cm_challenge_op_wiring::<R>(trace, k, log_kappa, nvars)?;
    let (pose_byte_vars, pose_field_vars) =
        cm_poseidon_challenge_vars::<R>(&pose_wiring, &byte_wiring, &op_wiring)?;

    // Glue bytes in the exact order used by short_challenge calls.
    if pose_byte_vars.len() != coin_wiring.byte_vars.len() {
        return Err("poseidon/coin byte length mismatch".to_string());
    }
    let mut glue: Vec<(usize, usize, usize, usize)> = Vec::with_capacity(coin_wiring.byte_vars.len());
    for (pv, lv) in pose_byte_vars.iter().zip(coin_wiring.byte_vars.iter()) {
        glue.push((0, *pv, 2, *lv));
    }

    // Field challenges part: allocate local vars with the expected values from the trace,
    // then glue them to the Poseidon squeeze-field vars selected by op wiring.
    let need_field = 2 * log_kappa + 2 + 2 * nvars;
    if pose_field_vars.len() != need_field {
        return Err("poseidon field var length mismatch".to_string());
    }
    // Extract the matching field values from the trace by scanning SqueezeField ops after short challenges.
    let mut squeezed_field_vals = Vec::with_capacity(need_field);
    let mut seen_first_bytes = false;
    let mut bytes_seen = 0usize;
    for op in &trace.ops {
        match op {
            LfPoseidonTraceOp::SqueezeBytes { .. } => {
                if !seen_first_bytes {
                    seen_first_bytes = true;
                }
                if seen_first_bytes && bytes_seen < (3 + k * R::dimension()) {
                    bytes_seen += 1;
                }
            }
            LfPoseidonTraceOp::SqueezeField(v) => {
                if seen_first_bytes && bytes_seen == (3 + k * R::dimension()) && squeezed_field_vals.len() < need_field {
                    if v.len() != 1 {
                        return Err("expected base-field squeeze len=1".to_string());
                    }
                    squeezed_field_vals.push(v[0]);
                }
            }
            _ => {}
        }
    }
    if squeezed_field_vals.len() != need_field {
        return Err("could not extract enough squeeze_field values for cm".to_string());
    }

    let mut b_fields = Dr1csBuilder::<BF<R>>::new();
    b_fields.enforce_var_eq_const(b_fields.one(), BF::<R>::ONE);
    let mut cur = 0usize;
    let take = |cur: &mut usize, n: usize, vs: &[BF<R>], b: &mut Dr1csBuilder<BF<R>>| -> Vec<usize> {
        let out = vs[*cur..*cur + n].iter().copied().map(|x| b.new_var(x)).collect::<Vec<_>>();
        *cur += n;
        out
    };
    let c0 = take(&mut cur, log_kappa, &squeezed_field_vals, &mut b_fields);
    let c1 = take(&mut cur, log_kappa, &squeezed_field_vals, &mut b_fields);
    let rc0 = take(&mut cur, 1, &squeezed_field_vals, &mut b_fields)[0];
    let sumcheck_r0 = take(&mut cur, nvars, &squeezed_field_vals, &mut b_fields);
    let rc1 = take(&mut cur, 1, &squeezed_field_vals, &mut b_fields)[0];
    let sumcheck_r1 = take(&mut cur, nvars, &squeezed_field_vals, &mut b_fields);
    let (field_inst, field_asg) = b_fields.into_instance();
    let field_wiring_local = CmFieldChallengeWiring { c0, c1, rc0, rc1, sumcheck_r0, sumcheck_r1 };

    // Glue local field vars to selected Poseidon squeeze-field vars in order.
    let mut local_field_vars = Vec::with_capacity(need_field);
    local_field_vars.extend_from_slice(&field_wiring_local.c0);
    local_field_vars.extend_from_slice(&field_wiring_local.c1);
    local_field_vars.push(field_wiring_local.rc0);
    local_field_vars.extend_from_slice(&field_wiring_local.sumcheck_r0);
    local_field_vars.push(field_wiring_local.rc1);
    local_field_vars.extend_from_slice(&field_wiring_local.sumcheck_r1);
    debug_assert_eq!(local_field_vars.len(), pose_field_vars.len());
    for (pv, lv) in pose_field_vars.iter().zip(local_field_vars.iter()) {
        glue.push((0, *pv, 3, *lv));
    }

    let parts = vec![
        (pose_inst, pose_asg),
        (params_inst, params_asg),
        (coin_inst, coin_asg),
        (field_inst, field_asg),
    ];
    let (inst, assignment) =
        merge_sparse_dr1cs_share_one_with_glue(&parts, &glue).map_err(|e| e.to_string())?;

    let public_len = 1 + 9;
    Ok((
        WeDr1csOutput {
            inst,
            assignment,
            public_len,
        },
        coin_wiring,
        field_wiring_local,
    ))
}

#[cfg(feature = "we_gate")]
pub fn build_we_dr1cs_for_cm_proof<R>(
    poseidon_cfg: &PoseidonConfig<BF<R>>,
    trace: &PoseidonTranscriptTrace<BF<R>>,
    params: &WeParams,
    public_inputs: &[BF<R>],
    proof: &crate::cm::CmProof<R>,
    mlen_mats: usize,
) -> Result<WeDr1csOutput<BF<R>>, String>
where
    R: OverField + CoeffRing + PolyRing,
    R::BaseRing: Zq + Field,
{
    let (out, _dbg) =
        build_we_dr1cs_for_cm_proof_debug::<R>(poseidon_cfg, trace, params, public_inputs, proof, mlen_mats)?;
    Ok(out)
}

#[cfg(feature = "we_gate")]
#[derive(Clone, Debug)]
pub struct WeCmBuildDebug {
    pub part_constraints: Vec<usize>,
    pub part_nvars: Vec<usize>,
    pub base_constraints: usize,
    pub glue: Vec<(usize, usize, usize, usize)>,
    pub cm_phase_marks: Vec<usize>,
    pub cm_phase_names: Vec<String>,
}

#[cfg(feature = "we_gate")]
pub fn build_we_dr1cs_for_cm_proof_debug<R>(
    poseidon_cfg: &PoseidonConfig<BF<R>>,
    trace: &PoseidonTranscriptTrace<BF<R>>,
    params: &WeParams,
    public_inputs: &[BF<R>],
    proof: &crate::cm::CmProof<R>,
    mlen_mats: usize,
) -> Result<(WeDr1csOutput<BF<R>>, WeCmBuildDebug), String>
where
    R: OverField + CoeffRing + PolyRing,
    R::BaseRing: Zq + Field,
{
    // Poseidon trace -> dR1CS (+ wiring).
    let ops = lf_ops_to_symphony_ops::<BF<R>>(&trace.ops);
    let (mut pose_inst, pose_asg, _replay, _byte_wit, pose_wiring, byte_wiring) =
        poseidon_sponge_dr1cs_from_trace_with_wiring_and_bytes::<BF<R>>(poseidon_cfg, &ops)
            .map_err(|e| format!("poseidon arith failed: {e}"))?;
    enforce_reabsorb_equals_squeeze::<BF<R>>(&mut pose_inst, &pose_wiring, &ops)?;

    // Public statement params prefix.
    let mut b_params = Dr1csBuilder::<BF<R>>::new();
    b_params.enforce_var_eq_const(b_params.one(), BF::<R>::ONE);
    let mut params_vars = Vec::with_capacity(9);
    for &x in &params.to_field_vec::<BF<R>>() {
        params_vars.push(b_params.new_var(x));
    }
    // Extra statement-defined public inputs.
    let mut pub_input_vars = Vec::with_capacity(public_inputs.len());
    for &x in public_inputs {
        pub_input_vars.push(b_params.new_var(x));
    }
    let (params_inst, params_asg) = b_params.into_instance();

    // Dcom prefix: includes `Out::verify` (setchk) + `rgchk::Dcom::verify` arithmetic checks,
    // and exposes the exact prefix absorb surface (excluding reabsorbs) for transcript gluing.
    let (dcom_inst, dcom_asg, dcom_wiring) =
        dcom_verifier_math_dr1cs::<R>(&proof.dcom, trace, params, public_inputs)?;

    // Cm coin surface (reconstruction gadgets) and glue to Poseidon squeeze outputs.
    let k = params.k as usize;
    let log_kappa = ark_std::log2((params.kappa as usize).next_power_of_two()) as usize;
    let nvars = params.nvars_cm as usize;

    let (coin_inst, coin_asg, coin_wiring) = cm_short_challenges_dr1cs::<R>(trace, k)?;
    let op_wiring = cm_challenge_op_wiring::<R>(trace, k, log_kappa, nvars)?;
    let (pose_byte_vars, pose_field_vars) =
        cm_poseidon_challenge_vars::<R>(&pose_wiring, &byte_wiring, &op_wiring)?;

    if pose_byte_vars.len() != coin_wiring.byte_vars.len() {
        return Err("poseidon/coin byte length mismatch".to_string());
    }
    let mut glue: Vec<(usize, usize, usize, usize)> = Vec::new();
    for (pv, lv) in pose_byte_vars.iter().zip(coin_wiring.byte_vars.iter()) {
        glue.push((0, *pv, 3, *lv));
    }

    // Field challenge local vars (same as in build_we_dr1cs_for_cm_challenges).
    let need_field = 2 * log_kappa + 2 + 2 * nvars;
    if pose_field_vars.len() != need_field {
        return Err("poseidon field var length mismatch".to_string());
    }
    let mut squeezed_field_vals = Vec::with_capacity(need_field);
    let mut seen_first_bytes = false;
    let mut bytes_seen = 0usize;
    for op in &trace.ops {
        match op {
            crate::recording_transcript::PoseidonTraceOp::SqueezeBytes { .. } => {
                if !seen_first_bytes {
                    seen_first_bytes = true;
                }
                if seen_first_bytes && bytes_seen < (3 + k * R::dimension()) {
                    bytes_seen += 1;
                }
            }
            crate::recording_transcript::PoseidonTraceOp::SqueezeField(v) => {
                if seen_first_bytes && bytes_seen == (3 + k * R::dimension()) && squeezed_field_vals.len() < need_field {
                    if v.len() != 1 {
                        return Err("expected base-field squeeze len=1".to_string());
                    }
                    squeezed_field_vals.push(v[0]);
                }
            }
            _ => {}
        }
    }
    if squeezed_field_vals.len() != need_field {
        return Err("could not extract enough squeeze_field values for cm".to_string());
    }

    let mut b_fields = Dr1csBuilder::<BF<R>>::new();
    b_fields.enforce_var_eq_const(b_fields.one(), BF::<R>::ONE);
    let mut cur = 0usize;
    let take = |cur: &mut usize, n: usize, vs: &[BF<R>], b: &mut Dr1csBuilder<BF<R>>| -> Vec<usize> {
        let out = vs[*cur..*cur + n].iter().copied().map(|x| b.new_var(x)).collect::<Vec<_>>();
        *cur += n;
        out
    };
    let c0 = take(&mut cur, log_kappa, &squeezed_field_vals, &mut b_fields);
    let c1 = take(&mut cur, log_kappa, &squeezed_field_vals, &mut b_fields);
    let rc0 = take(&mut cur, 1, &squeezed_field_vals, &mut b_fields)[0];
    let sumcheck_r0 = take(&mut cur, nvars, &squeezed_field_vals, &mut b_fields);
    let rc1 = take(&mut cur, 1, &squeezed_field_vals, &mut b_fields)[0];
    let sumcheck_r1 = take(&mut cur, nvars, &squeezed_field_vals, &mut b_fields);
    let (field_inst, field_asg) = b_fields.into_instance();
    let field_wiring_local = CmFieldChallengeWiring { c0, c1, rc0, rc1, sumcheck_r0, sumcheck_r1 };

    // Glue local field vars to Poseidon squeeze-field vars.
    let mut local_field_vars = Vec::with_capacity(need_field);
    local_field_vars.extend_from_slice(&field_wiring_local.c0);
    local_field_vars.extend_from_slice(&field_wiring_local.c1);
    local_field_vars.push(field_wiring_local.rc0);
    local_field_vars.extend_from_slice(&field_wiring_local.sumcheck_r0);
    local_field_vars.push(field_wiring_local.rc1);
    local_field_vars.extend_from_slice(&field_wiring_local.sumcheck_r1);
    for (pv, lv) in pose_field_vars.iter().zip(local_field_vars.iter()) {
        glue.push((0, *pv, 4, *lv));
    }

    // Cm verifier arithmetic + absorb surface builder.
    let (cm_inst, cm_asg, cm_wiring) =
        cm_verifier_math_dr1cs::<R>(trace, proof, k, log_kappa, nvars, mlen_mats)?;

    // Glue cm_wiring challenges to the coin/field wiring parts (so the math uses the same coins).
    // Bytes:
    for (cv, lv) in cm_wiring.short.byte_vars.iter().zip(coin_wiring.byte_vars.iter()) {
        glue.push((5, *cv, 3, *lv));
    }
    // Field scalars:
    for (cv, lv) in cm_wiring.field.c0.iter().zip(field_wiring_local.c0.iter()) {
        glue.push((5, *cv, 4, *lv));
    }
    for (cv, lv) in cm_wiring.field.c1.iter().zip(field_wiring_local.c1.iter()) {
        glue.push((5, *cv, 4, *lv));
    }
    glue.push((5, cm_wiring.field.rc0, 4, field_wiring_local.rc0));
    glue.push((5, cm_wiring.field.rc1, 4, field_wiring_local.rc1));
    for (cv, lv) in cm_wiring.field.sumcheck_r0.iter().zip(field_wiring_local.sumcheck_r0.iter()) {
        glue.push((5, *cv, 4, *lv));
    }
    for (cv, lv) in cm_wiring.field.sumcheck_r1.iter().zip(field_wiring_local.sumcheck_r1.iter()) {
        glue.push((5, *cv, 4, *lv));
    }

    // Glue Cm absorb surface (non-reabsorb absorbs after first SqueezeBytes) to Poseidon absorb vars.
    // Compute the absorb-op index at which the Cm segment starts (first SqueezeBytes).
    let mut absorb_ops_before_cm = 0usize;
    let mut seen_bytes = false;
    for op in &ops {
        match op {
            symphony::transcript::PoseidonTraceOp::SqueezeBytes { .. } => {
                seen_bytes = true;
                break;
            }
            symphony::transcript::PoseidonTraceOp::Absorb(_) => absorb_ops_before_cm += 1,
            _ => {}
        }
    }
    if !seen_bytes {
        return Err("cm proof: trace has no SqueezeBytes (short_challenge) marker".to_string());
    }

    // Determine which Absorb ops are reabsorbs (immediately following a SqueezeField).
    let mut is_reabsorb = vec![false; pose_wiring.absorb_ranges.len()];
    let mut expect_reabsorb = false;
    let mut absorb_idx = 0usize;
    for op in &ops {
        match op {
            symphony::transcript::PoseidonTraceOp::SqueezeField(_) => expect_reabsorb = true,
            symphony::transcript::PoseidonTraceOp::SqueezeBytes { .. } => {}
            symphony::transcript::PoseidonTraceOp::Absorb(_) => {
                if expect_reabsorb {
                    if absorb_idx < is_reabsorb.len() {
                        is_reabsorb[absorb_idx] = true;
                    }
                    expect_reabsorb = false;
                }
                absorb_idx += 1;
            }
        }
    }

    // Glue Dcom-prefix squeeze-field vars (prefix before first SqueezeBytes) to Poseidon squeeze-field vars.
    if pose_wiring.squeeze_field_vars.len() < dcom_wiring.squeeze_field_vars.len() {
        return Err("poseidon wiring: not enough squeeze_field_vars for dcom prefix".to_string());
    }
    for (i, &sv) in dcom_wiring.squeeze_field_vars.iter().enumerate() {
        glue.push((0, pose_wiring.squeeze_field_vars[i], 2, sv));
    }
    // Glue statement params (public inputs) into the Dcom prefix gadget.
    if params_vars.len() != dcom_wiring.params_vars.len() {
        return Err("params glue length mismatch".to_string());
    }
    for (pv, dv) in params_vars.iter().zip(dcom_wiring.params_vars.iter()) {
        glue.push((1, *pv, 2, *dv));
    }
    // Glue extra public inputs into the Dcom prefix gadget.
    if pub_input_vars.len() != dcom_wiring.public_input_vars.len() {
        return Err("public input glue length mismatch".to_string());
    }
    for (pv, dv) in pub_input_vars.iter().zip(dcom_wiring.public_input_vars.iter()) {
        glue.push((1, *pv, 2, *dv));
    }

    // Flatten Poseidon absorb vars for non-reabsorb absorbs *before* the Cm segment.
    let mut pose_abs_prefix: Vec<usize> = Vec::new();
    for (i, (start, len)) in pose_wiring.absorb_ranges.iter().enumerate() {
        if i >= absorb_ops_before_cm {
            break;
        }
        if is_reabsorb[i] {
            continue;
        }
        pose_abs_prefix.extend_from_slice(&pose_wiring.absorb_vars[*start..*start + *len]);
    }
    if pose_abs_prefix.len() != dcom_wiring.absorb_flat.len() {
        return Err(format!(
            "prefix absorb glue length mismatch: pose={} local={}",
            pose_abs_prefix.len(),
            dcom_wiring.absorb_flat.len()
        ));
    }
    // Glue the entire Dcom-prefix absorb surface.
    for (pv, lv) in pose_abs_prefix.iter().zip(dcom_wiring.absorb_flat.iter()) {
        glue.push((0, *pv, 2, *lv));
    }

    // Flatten Poseidon absorb vars for non-reabsorb absorbs starting at Cm segment.
    let mut pose_abs_flat: Vec<usize> = Vec::new();
    for (i, (start, len)) in pose_wiring.absorb_ranges.iter().enumerate() {
        if i < absorb_ops_before_cm {
            continue;
        }
        if is_reabsorb[i] {
            continue;
        }
        pose_abs_flat.extend_from_slice(&pose_wiring.absorb_vars[*start..*start + *len]);
    }
    if pose_abs_flat.len() != cm_wiring.absorb_flat.len() {
        return Err(format!(
            "cm absorb glue length mismatch: pose={} cm={}",
            pose_abs_flat.len(),
            cm_wiring.absorb_flat.len()
        ));
    }
    for (pv, cv) in pose_abs_flat.iter().zip(cm_wiring.absorb_flat.iter()) {
        glue.push((0, *pv, 5, *cv));
    }

    let parts = vec![
        (pose_inst, pose_asg),   // 0
        (params_inst, params_asg), // 1
        (dcom_inst, dcom_asg),   // 2
        (coin_inst, coin_asg),   // 3
        (field_inst, field_asg), // 4
        (cm_inst, cm_asg),       // 5
    ];
    let part_constraints = parts.iter().map(|(i, _)| i.constraints.len()).collect::<Vec<_>>();
    let part_nvars = parts.iter().map(|(i, _)| i.nvars).collect::<Vec<_>>();
    let base_constraints = part_constraints.iter().sum::<usize>();
    let (inst, assignment) =
        merge_sparse_dr1cs_share_one_with_glue(&parts, &glue).map_err(|e| e.to_string())?;

    let public_len = 1 + 9 + public_inputs.len();
    Ok((
        WeDr1csOutput { inst, assignment, public_len },
        WeCmBuildDebug {
            part_constraints,
            part_nvars,
            base_constraints,
            glue,
            cm_phase_marks: cm_wiring.phase_marks.clone(),
            cm_phase_names: cm_wiring.phase_names.clone(),
        },
    ))
}

#[cfg(all(test, feature = "we_gate"))]
mod tests {
    use super::*;
    use cyclotomic_rings::rings::GoldilocksPoseidonConfig as PC;
    use latticefold::arith::r1cs::R1CS;
    use latticefold::transcript::Transcript;
    use ark_ff::{Fp384, MontBackend, MontConfig};
    use stark_rings::balanced_decomposition::GadgetDecompose;
    use stark_rings::cyclotomic_ring::models::goldilocks::RqPoly as R;
    use stark_rings_linalg::{Matrix, SparseMatrix};

    use crate::lin::Linearize;
    use crate::lin::LinearizedVerify;
    use crate::r1cs::ComR1CS;
    use crate::recording_transcript::TracePoseidonTranscript;
    use crate::rgchk::{DecompParameters, Rg, RgInstance};
    use crate::cm::Cm;

    #[derive(MontConfig)]
    #[modulus = "39402006196394479212279040100143613805079739270465446667948293404245721771496870329047266088258938001861606973112319"]
    #[generator = "2"]
    pub struct Secp384r1Config;
    type FBig = Fp384<MontBackend<Secp384r1Config, 6>>;

    #[derive(Clone)]
    struct ReplayPoseidonTranscript<RR: OverField> {
        idx: usize,
        trace: crate::recording_transcript::PoseidonTranscriptTrace<<RR::BaseRing as Field>::BasePrimeField>,
        scratch: Vec<<RR::BaseRing as Field>::BasePrimeField>,
    }

    impl<RR: OverField> ReplayPoseidonTranscript<RR> {
        fn new(trace: &crate::recording_transcript::PoseidonTranscriptTrace<<RR::BaseRing as Field>::BasePrimeField>) -> Self {
            Self { idx: 0, trace: trace.clone(), scratch: Vec::with_capacity(64) }
        }
        fn advance(&mut self) {
            self.idx += 1;
        }
    }

    impl<RR: OverField> Transcript<RR> for ReplayPoseidonTranscript<RR> {
        type TranscriptConfig = ark_crypto_primitives::sponge::poseidon::PoseidonConfig<<RR::BaseRing as Field>::BasePrimeField>;
        fn new(_config: &Self::TranscriptConfig) -> Self {
            unreachable!("ReplayPoseidonTranscript::new(trace) should be used in tests")
        }

        fn absorb(&mut self, v: &RR) {
            self.scratch.clear();
            for c in v.coeffs() {
                self.scratch.extend(c.to_base_prime_field_elements());
            }
            let op = self.trace.ops.get(self.idx).expect("replay: op index oob").clone();
            match op {
                crate::recording_transcript::PoseidonTraceOp::Absorb(elems) => {
                    assert_eq!(
                        elems.as_slice(),
                        self.scratch.as_slice(),
                        "replay absorb mismatch at op {}",
                        self.idx
                    );
                }
                other => panic!("replay expected Absorb op, got {other:?} at idx {}", self.idx),
            };
            self.advance();
        }

        fn get_challenge(&mut self) -> RR::BaseRing {
            let ext = RR::BaseRing::extension_degree() as usize;
            let op0 = self.trace.ops.get(self.idx).expect("replay: op index oob").clone();
            let c = match op0 {
                crate::recording_transcript::PoseidonTraceOp::SqueezeField(v) => v,
                other => panic!("replay expected SqueezeField op, got {other:?} at idx {}", self.idx),
            };
            assert_eq!(c.len(), ext, "replay squeeze_field length mismatch");
            self.advance();

            // get_challenge reabsorbs the squeezed field elements; the trace records that as Absorb(c).
            let op1 = self.trace.ops.get(self.idx).expect("replay: op index oob").clone();
            match op1 {
                crate::recording_transcript::PoseidonTraceOp::Absorb(v) => {
                    assert_eq!(v.as_slice(), c.as_slice(), "replay reabsorb mismatch");
                }
                other => panic!("replay expected reabsorb Absorb op after SqueezeField, got {other:?}"),
            };
            self.advance();

            <RR::BaseRing as Field>::from_base_prime_field_elems(&c)
                .expect("replay: wrong extension_degree")
        }

        fn squeeze_bytes(&mut self, n: usize) -> Vec<u8> {
            let op = self.trace.ops.get(self.idx).expect("replay: op index oob").clone();
            let out = match op {
                crate::recording_transcript::PoseidonTraceOp::SqueezeBytes { n: nn, out } => {
                    assert_eq!(nn, n, "replay squeeze_bytes n mismatch");
                    out
                }
                other => panic!("replay expected SqueezeBytes op, got {other:?} at idx {}", self.idx),
            };
            self.advance();
            out
        }
    }

    fn identity_cs(n: usize) -> (R1CS<R>, Vec<R>) {
        let r1cs = R1CS::<R> {
            l: 1,
            A: SparseMatrix::identity(n),
            B: SparseMatrix::identity(n),
            C: SparseMatrix::identity(n),
        };
        let z = vec![<R as stark_rings::Ring>::ONE; n];
        (r1cs, z)
    }

    #[test]
    fn test_we_arith_linearized_verify_constraints_satisfy() {
        let n = 1 << 7;
        let k = 4;
        let m = n / k;
        let b = 2;
        let kappa = 2;
        let (mut r1cs, z) = identity_cs(m);
        r1cs.A = r1cs.A.gadget_decompose(b, k);
        r1cs.B = r1cs.B.gadget_decompose(b, k);
        r1cs.C = r1cs.C.gadget_decompose(b, k);

        let A = Matrix::<R>::rand(&mut ark_std::test_rng(), kappa, n);
        let cr1cs = ComR1CS::new(r1cs, z, 1, b, k, &A);

        // Build proof.
        let mut ts = crate::transcript::PoseidonTranscript::empty::<PC>();
        let (_linb, lproof) = cr1cs.linearize(&mut ts);

        // Record the verifier transcript coin stream.
        let mut rec = TracePoseidonTranscript::<R>::empty::<PC>();
        assert!(lproof.verify(&mut rec));
        let trace = rec.trace().clone();

        // Build verifier-math dR1CS and check satisfaction.
        let mut ch = ChallengeCursor::<BF<R>>::new(&trace.squeezed_field);
        let (inst, asg, _ch_vars) = comr1cs_verifier_math_dr1cs::<R>(&lproof, &mut ch).unwrap();
        inst.check(&asg).unwrap();

        // Sanity: we consumed exactly 2*nvars squeeze-field scalars (r_pre and sumcheck rs).
        assert_eq!(ch.consumed(), 2 * lproof.nvars);
    }

    #[test]
    fn test_we_arith_sumcheck_degree2_constraints_satisfy() {
        use latticefold::utils::sumcheck::MLSumcheck;
        use stark_rings_poly::mle::DenseMultilinearExtension;
        use ark_std::UniformRand;

        // Small instance: sumcheck of product of two random MLEs (degree=2).
        let nvars = 6usize;
        let n = 1usize << nvars;
        let mut rng = ark_std::test_rng();
        let evals0 = (0..n).map(|_| R::rand(&mut rng)).collect::<Vec<_>>();
        let evals1 = (0..n).map(|_| R::rand(&mut rng)).collect::<Vec<_>>();
        let mle0 = DenseMultilinearExtension::from_evaluations_vec(nvars, evals0);
        let mle1 = DenseMultilinearExtension::from_evaluations_vec(nvars, evals1);

        let mut ts_p = crate::transcript::PoseidonTranscript::empty::<PC>();
        let (proof, _state) = MLSumcheck::<R, _>::prove_as_subprotocol(
            &mut ts_p,
            vec![mle0, mle1],
            nvars,
            2,
            |vals: &[R]| vals[0] * vals[1],
        );
        // Equivalent to MLSumcheck::extract_sum (avoid transcript type parameter inference).
        let claimed_sum = proof.msgs()[0].evaluations[0] + proof.msgs()[0].evaluations[1];

        // Run verifier to get the real transcript coin stream (r_i).
        let mut rec = crate::recording_transcript::TracePoseidonTranscript::<R>::empty::<PC>();
        let _sub = MLSumcheck::<R, _>::verify_as_subprotocol(
            &mut rec,
            nvars,
            2,
            claimed_sum,
            &proof,
        )
        .unwrap();
        let trace = rec.trace().clone();

        // Build dR1CS for sumcheck verify (standalone, with challenges from trace.squeezed_field).
        type F = BF<R>;
        let mut ch = ChallengeCursor::<F>::new(&trace.squeezed_field);
        let mut b = Dr1csBuilder::<F>::new();
        b.enforce_var_eq_const(b.one(), F::from(1u64));

        // Allocate sumcheck prover msgs (3 evals per round for degree 2).
        let msgs = proof.msgs().to_vec();
        assert_eq!(msgs.len(), nvars);
        let mut msg_vars: Vec<[RingVars; 3]> = Vec::with_capacity(nvars);
        for m in msgs {
            assert_eq!(m.evaluations.len(), 3);
            let e0 = ring_to_ringvars::<R>(&mut b, &m.evaluations[0]);
            let e1 = ring_to_ringvars::<R>(&mut b, &m.evaluations[1]);
            let e2 = ring_to_ringvars::<R>(&mut b, &m.evaluations[2]);
            msg_vars.push([e0, e1, e2]);
        }

        // Sample r_i from trace (matches verify_as_subprotocol schedule).
        let mut r_sc = Vec::with_capacity(nvars);
        for _ in 0..nvars {
            r_sc.push(ch.next(&mut b));
        }

        let claim0 = ring_to_ringvars::<R>(&mut b, &claimed_sum);
        let _final_claim = sumcheck_verify_degree2::<F>(&mut b, claim0, &msg_vars, &r_sc).unwrap();

        let (inst, asg) = b.into_instance();
        inst.check(&asg).unwrap();
    }

    #[test]
    fn test_short_challenge_coeff_from_byte_matches_rust() {
        type F = BF<R>;
        let u = 32u64; // typical for lambda=128,d=24 => floor(128/24)=5 => u=32

        for byte_u8 in [0u8, 1, 2, 15, 16, 31, 32, 33, 63, 64, 127, 128, 200, 255] {
            let mut b = Dr1csBuilder::<F>::new();
            b.enforce_var_eq_const(b.one(), F::from(1u64));
            let byte = const_var(&mut b, F::from(byte_u8 as u64));
            let coeff = short_challenge_coeff_from_byte::<F>(&mut b, byte, u);

            let expected_i64 = ((byte_u8 as u64) % u) as i64 - (u as i64 / 2);
            let expected = if expected_i64 >= 0 {
                F::from(expected_i64 as u64)
            } else {
                -F::from((-expected_i64) as u64)
            };

            assert_eq!(b.assignment[coeff], expected);
            let (inst, asg) = b.into_instance();
            inst.check(&asg).unwrap();
        }
    }

    #[test]
    fn test_short_challenge_from_bytes_matches_rust() {
        use crate::utils::short_challenge;

        type F = BF<R>;
        let d = R::dimension();

        // Get bytes from the real transcript by calling `short_challenge`.
        let mut rec = crate::recording_transcript::TracePoseidonTranscript::<R>::empty::<PC>();
        let r_sc = short_challenge::<R>(128, &mut rec);
        let bytes = match rec.trace().ops.last().unwrap() {
            crate::recording_transcript::PoseidonTraceOp::SqueezeBytes { out, .. } => out.clone(),
            _ => panic!("expected last op to be SqueezeBytes"),
        };
        assert_eq!(bytes.len(), d);

        // Arithmetize from byte vars.
        let mut b = Dr1csBuilder::<F>::new();
        b.enforce_var_eq_const(b.one(), F::from(1u64));
        let byte_vars = bytes
            .iter()
            .map(|&by| const_var(&mut b, F::from(by as u64)))
            .collect::<Vec<_>>();
        let ring = short_challenge_from_bytes::<F>(&mut b, &byte_vars, 128, d);

        // Compare coefficients.
        let expected_coeffs = r_sc.coeffs().to_vec();
        assert_eq!(expected_coeffs.len(), d);
        for i in 0..d {
            let exp_bf = expected_coeffs[i]
                .to_base_prime_field_elements()
                .into_iter()
                .next()
                .unwrap();
            assert_eq!(b.assignment[ring.coeffs[i]], exp_bf);
        }

        let (inst, asg) = b.into_instance();
        inst.check(&asg).unwrap();
    }

    fn part_offset(part: usize, part_nvars: &[usize]) -> usize {
        // Variable 0 is shared "one" across parts. Merge offsets grow by (nvars_i - 1).
        let mut off = 0usize;
        for i in 0..part {
            off += part_nvars[i].saturating_sub(1);
        }
        off
    }

    #[test]
    fn test_we_cm_proof_param_binding_mlen_unsat() {
        // Small-ish instance.
        type PCF = cyclotomic_rings::rings::FrogPoseidonConfig;
        use cyclotomic_rings::rings::GetPoseidonParams;
        use ark_ff::Zero;
        use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as RR;
        use stark_rings::PolyRing;

        let k = 1usize;
        let kappa = 1usize;
        let ell = 32usize;
        let b = 2u128;
        let d = RR::dimension();
        let tau_unpadded_len = kappa * (k * d) * ell * d;
        let n = tau_unpadded_len.next_power_of_two();
        let nvars = ark_std::log2(n) as usize;

        let dparams = DecompParameters { b, k, l: ell };
        let mut rng = ark_std::test_rng();
        let f = vec![RR::from(<RR as PolyRing>::BaseRing::zero()); n];
        let A = Matrix::<RR>::rand(&mut rng, kappa, n);
        let inst = RgInstance::from_f(f, &A, &dparams);
        let rg = Rg { nvars, instances: vec![inst], dparams: dparams.clone() };
        let cm = Cm { rg };
        let M: Vec<SparseMatrix<RR>> = vec![SparseMatrix::identity(n)];

        // Build proof + trace.
        let mut ts = crate::transcript::PoseidonTranscript::empty::<PCF>();
        let (_com, proof) = cm.prove(&M, &mut ts);
        let mut rec = TracePoseidonTranscript::<RR>::empty::<PCF>();
        proof.verify(&M, &mut rec).expect("cm verify");
        let trace = rec.trace().clone();

        let params = WeParams {
            nvars_setchk: nvars as u64,
            degree_setchk: 3,
            nvars_cm: nvars as u64,
            degree_cm: 2,
            kappa: kappa as u64,
            ring_dim_d: RR::dimension() as u64,
            k: k as u64,
            l: ell as u64,
            mlen: M.len() as u64,
        };
        let poseidon_cfg = PCF::get_poseidon_config();

        let (out, dbg) =
            build_we_dr1cs_for_cm_proof_debug::<RR>(&poseidon_cfg, &trace, &params, &[], &proof, M.len())
                .expect("build we dr1cs");
        out.inst.check(&out.assignment).expect("should satisfy");

        // Mutate the public `mlen` parameter only: params part is part index 1, and `mlen` is the 9th param (local var = 9).
        let params_part = 1usize;
        let off = part_offset(params_part, &dbg.part_nvars);
        let mlen_local_var = 9usize;
        let mlen_global = off + mlen_local_var;

        let mut bad = out.assignment.clone();
        bad[mlen_global] += <BF<RR> as ark_ff::Field>::ONE;
        assert!(out.inst.check(&bad).is_err(), "mlen mutation should break satisfaction");
    }

    #[test]
    fn test_we_cm_proof_unsat_on_constraint_var_flip() {
        // Reuse the same construction as the param-binding test, but flip a variable that is
        // guaranteed to appear in some constraint.
        type PCF = cyclotomic_rings::rings::FrogPoseidonConfig;
        use cyclotomic_rings::rings::GetPoseidonParams;
        use ark_ff::Zero;
        use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as RR;
        use stark_rings::PolyRing;

        let k = 1usize;
        let kappa = 1usize;
        let ell = 32usize;
        let b = 2u128;
        let d = RR::dimension();
        let tau_unpadded_len = kappa * (k * d) * ell * d;
        let n = tau_unpadded_len.next_power_of_two();
        let nvars = ark_std::log2(n) as usize;

        let dparams = DecompParameters { b, k, l: ell };
        let mut rng = ark_std::test_rng();
        let f = vec![RR::from(<RR as PolyRing>::BaseRing::zero()); n];
        let A = Matrix::<RR>::rand(&mut rng, kappa, n);
        let inst = RgInstance::from_f(f, &A, &dparams);
        let rg = Rg { nvars, instances: vec![inst], dparams: dparams.clone() };
        let cm = Cm { rg };
        let M: Vec<SparseMatrix<RR>> = vec![SparseMatrix::identity(n)];

        let mut ts = crate::transcript::PoseidonTranscript::empty::<PCF>();
        let (_com, proof) = cm.prove(&M, &mut ts);
        let mut rec = TracePoseidonTranscript::<RR>::empty::<PCF>();
        proof.verify(&M, &mut rec).expect("cm verify");
        let trace = rec.trace().clone();

        let params = WeParams {
            nvars_setchk: nvars as u64,
            degree_setchk: 3,
            nvars_cm: nvars as u64,
            degree_cm: 2,
            kappa: kappa as u64,
            ring_dim_d: RR::dimension() as u64,
            k: k as u64,
            l: ell as u64,
            mlen: M.len() as u64,
        };
        let poseidon_cfg = PCF::get_poseidon_config();

        let (out, _dbg) =
            build_we_dr1cs_for_cm_proof_debug::<RR>(&poseidon_cfg, &trace, &params, &[], &proof, M.len())
                .expect("build we dr1cs");
        out.inst.check(&out.assignment).expect("should satisfy");

        // Pick a var index that appears in the first constraint (and is not the shared ONE=0).
        let c0 = &out.inst.constraints[0];
        let pick = |terms: &Vec<(BF<RR>, usize)>| -> Option<usize> {
            terms.iter().map(|(_, v)| *v).find(|&v| v != 0)
        };
        let v = pick(&c0.a)
            .or_else(|| pick(&c0.b))
            .or_else(|| pick(&c0.c))
            .expect("constraint should mention some non-one variable");

        let mut bad = out.assignment.clone();
        bad[v] += <BF<RR> as ark_ff::Field>::ONE;
        assert!(out.inst.check(&bad).is_err(), "flipping constrained var should break satisfaction");
    }

    #[test]
    fn test_transcript_roundtrip_cm_verify() {
        // Record a trace from a real verifier run, then replay that exact trace and ensure:
        // - all absorbs/squeezes line up
        // - verifier returns Ok
        type PCF = cyclotomic_rings::rings::FrogPoseidonConfig;
        use ark_ff::Zero;
        use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as RR;
        use stark_rings::PolyRing;

        let k = 1usize;
        let kappa = 1usize;
        let ell = 32usize;
        let b = 2u128;
        let d = RR::dimension();
        let tau_unpadded_len = kappa * (k * d) * ell * d;
        let n = tau_unpadded_len.next_power_of_two();
        let nvars = ark_std::log2(n) as usize;

        let dparams = DecompParameters { b, k, l: ell };
        let mut rng = ark_std::test_rng();
        let f = vec![RR::from(<RR as PolyRing>::BaseRing::zero()); n];
        let A = Matrix::<RR>::rand(&mut rng, kappa, n);
        let inst = RgInstance::from_f(f, &A, &dparams);
        let rg = Rg { nvars, instances: vec![inst], dparams: dparams.clone() };
        let cm = Cm { rg };
        let M: Vec<SparseMatrix<RR>> = vec![SparseMatrix::identity(n)];

        let mut ts = crate::transcript::PoseidonTranscript::empty::<PCF>();
        let (_com, proof) = cm.prove(&M, &mut ts);

        // Record.
        let mut rec = TracePoseidonTranscript::<RR>::empty::<PCF>();
        proof.verify(&M, &mut rec).expect("cm verify (record)");
        let trace = rec.trace().clone();

        // Replay.
        let mut replay = ReplayPoseidonTranscript::<RR>::new(&trace);
        proof.verify(&M, &mut replay).expect("cm verify (replay)");
        assert_eq!(replay.idx, replay.trace.ops.len(), "replay should consume full trace");
    }

    #[test]
    fn test_we_cm_proof_public_input_digest_unsat_on_flip() {
        // SP1-style: include exactly one public input digest, absorb it before proving/verifying,
        // and ensure flipping it in `public_x` breaks dR1CS satisfaction.
        type PCF = cyclotomic_rings::rings::FrogPoseidonConfig;
        use ark_ff::Zero;
        use cyclotomic_rings::rings::GetPoseidonParams;
        use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as RR;
        use stark_rings::PolyRing;

        let k = 1usize;
        let kappa = 1usize;
        let ell = 32usize;
        let b = 2u128;
        let d = RR::dimension();
        let tau_unpadded_len = kappa * (k * d) * ell * d;
        let n = tau_unpadded_len.next_power_of_two();
        let nvars = ark_std::log2(n) as usize;

        let dparams = DecompParameters { b, k, l: ell };
        let mut rng = ark_std::test_rng();
        let f = vec![RR::from(<RR as PolyRing>::BaseRing::zero()); n];
        let A = Matrix::<RR>::rand(&mut rng, kappa, n);
        let inst = RgInstance::from_f(f, &A, &dparams);
        let rg = Rg { nvars, instances: vec![inst], dparams: dparams.clone() };
        let cm = Cm { rg };
        let M: Vec<SparseMatrix<RR>> = vec![SparseMatrix::identity(n)];

        // Statement-defined SP1 public input digest (base field element).
        type BF0 = <<RR as PolyRing>::BaseRing as ark_ff::Field>::BasePrimeField;
        let sp1_digest: BF0 = BF0::from(123u64);

        // Prove with digest absorbed before proving.
        let mut ts = crate::transcript::PoseidonTranscript::empty::<PCF>();
        ts.absorb_field_element(&sp1_digest);
        let (_com, proof) = cm.prove(&M, &mut ts);

        // Record verifier trace with the same digest absorbed before verify.
        let mut rec = TracePoseidonTranscript::<RR>::empty::<PCF>();
        rec.absorb_field_element(&sp1_digest);
        proof.verify(&M, &mut rec).expect("cm verify");
        let trace = rec.trace().clone();

        let params = WeParams {
            nvars_setchk: nvars as u64,
            degree_setchk: 3,
            nvars_cm: nvars as u64,
            degree_cm: 2,
            kappa: kappa as u64,
            ring_dim_d: RR::dimension() as u64,
            k: k as u64,
            l: ell as u64,
            mlen: M.len() as u64,
        };
        let poseidon_cfg = PCF::get_poseidon_config();

        let (out, dbg) = build_we_dr1cs_for_cm_proof_debug::<RR>(
            &poseidon_cfg,
            &trace,
            &params,
            &[sp1_digest],
            &proof,
            M.len(),
        )
        .expect("build we dr1cs");
        out.inst.check(&out.assignment).expect("should satisfy");

        // Flip the SP1 digest in the params/public-input part only.
        // Part index 1 is the "params+public_inputs" part; local var index:
        // 0 = ONE, 1..=9 = WeParams, 10 = first extra public input.
        let params_part = 1usize;
        let off = part_offset(params_part, &dbg.part_nvars);
        let digest_local_var = 10usize;
        let digest_global = off + digest_local_var;

        let mut bad = out.assignment.clone();
        bad[digest_global] += <BF<RR> as ark_ff::Field>::ONE;
        assert!(out.inst.check(&bad).is_err(), "digest flip should break satisfaction");
    }

    // Large-scale (n=2^20) end-to-end WE->DPP sanity + digest "before/after" metric.
    // Run with: cargo test --release -p latticefold-plus test_large_trace --features we_gate -- --nocapture --ignored
    #[test]
    #[ignore]
    fn test_large_trace() {
        use cyclotomic_rings::rings::FrogPoseidonConfig as PCF;
        use cyclotomic_rings::rings::GetPoseidonParams;
        use ark_ff::Zero;
        use sha2::{Digest, Sha256};
        use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as RR;
        use stark_rings::PolyRing;

        use crate::we_statement::{digest32_to_field, LFP_WE_GATE_DIGEST_V1, we_statement_hash_lf_plus};
        use dpp::dr1cs_flpcp::Dr1csInstanceSparse as DppInst;
        use dpp::sparse::SparseVec;
        use dpp::pipeline::build_rev2_dpp_sparse_boolean_auto;

        use ark_ff::PrimeField;
        use rand::{rngs::StdRng, SeedableRng};
        #[cfg(feature = "parallel")]
        use rayon::current_num_threads;
        #[cfg(feature = "parallel")]
        use rayon::prelude::*;

        fn lift_to_big<Fs: PrimeField>(x: Fs) -> FBig {
            FBig::from_le_bytes_mod_order(&x.into_bigint().to_bytes_le())
        }

        fn lift_bit_to_big<Fs: PrimeField>(b: Fs) -> FBig {
            if b.is_zero() {
                FBig::from(0u64)
            } else {
                debug_assert_eq!(b, Fs::ONE, "non-boolean bit in proof encoding");
                FBig::from(1u64)
            }
        }

        // Default to the historical large-scale setting (n=2^20), but allow smaller for scaling studies.
        let n_pow: usize = std::env::var("LFP_TRACE_NPOW")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(20);
        let n = 1usize << n_pow;
        let k = 4usize;
        let kappa = 2usize;
        let ell = 32usize;
        let d = RR::dimension();
        let b = (d / 2) as u128;
        let nvars = ark_std::log2(n) as usize;
        let dparams = DecompParameters { b, k, l: ell };

        // Default Mlen=3 like earlier LF+ timings; allow override to inspect Mlen scaling.
        let mlen: usize = std::env::var("LFP_TRACE_MLEN")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(3);
        let M: Vec<SparseMatrix<RR>> = vec![SparseMatrix::identity(n); mlen];

        type FSmall = <<RR as PolyRing>::BaseRing as ark_ff::Field>::BasePrimeField;
        let digest_toy: FSmall = FSmall::from(42u64);
        let digest_hash: FSmall = {
            let d: [u8; 32] = Sha256::digest(b"LFP_SP1_PUBLIC_INPUT_DIGEST_V1").into();
            digest32_to_field::<FSmall>(d)
        };

        let run_one = |label: &str, sp1_digest: FSmall| {
            eprintln!("\n[test_large_trace] case={label} n=2^{n_pow} mlen={mlen}");
            #[cfg(feature = "parallel")]
            eprintln!("[test_large_trace] rayon_threads={}", current_num_threads());
            #[cfg(not(feature = "parallel"))]
            eprintln!("[test_large_trace] rayon_threads=DISABLED(feature=parallel)");
            let mut rng = ark_std::test_rng();
            let f = vec![RR::from(<RR as PolyRing>::BaseRing::zero()); n];
            let A = Matrix::<RR>::rand(&mut rng, kappa, n);
            let inst = RgInstance::from_f(f, &A, &dparams);
            let rg = Rg { nvars, instances: vec![inst], dparams: dparams.clone() };
            let cm = Cm { rg };

            let t0 = std::time::Instant::now();
            let mut ts = crate::transcript::PoseidonTranscript::empty::<PCF>();
            ts.absorb_field_element(&sp1_digest);
            let (_com, proof) = cm.prove(&M, &mut ts);
            eprintln!("[test_large_trace] cm.prove: {:?}", t0.elapsed());

            let t1 = std::time::Instant::now();
            let mut rec = TracePoseidonTranscript::<RR>::empty::<PCF>();
            rec.absorb_field_element(&sp1_digest);
            proof.verify(&M, &mut rec).expect("cm verify");
            let trace = rec.trace().clone();
            eprintln!("[test_large_trace] cm.verify(record): {:?}", t1.elapsed());

            let params = WeParams {
                nvars_setchk: nvars as u64,
                degree_setchk: 3,
                nvars_cm: nvars as u64,
                degree_cm: 2,
                kappa: kappa as u64,
                ring_dim_d: RR::dimension() as u64,
                k: k as u64,
                l: ell as u64,
                mlen: M.len() as u64,
            };
            let poseidon_cfg = PCF::get_poseidon_config();

            let t2 = std::time::Instant::now();
            let (out, _dbg) = build_we_dr1cs_for_cm_proof_debug::<RR>(
                &poseidon_cfg,
                &trace,
                &params,
                &[sp1_digest],
                &proof,
                M.len(),
            )
            .expect("build we dr1cs");
            out.inst.check(&out.assignment).expect("dr1cs sat");
            eprintln!(
                "[test_large_trace] build_we_dr1cs: {:?} (nvars={}, constraints={})",
                t2.elapsed(),
                out.inst.nvars,
                out.inst.constraints.len()
            );

            // DPP verification (single query).
            let t3 = std::time::Instant::now();
            let inst_sparse = DppInst::<FSmall> {
                n: out.inst.nvars,
                a: out.inst.constraints.iter().map(|r| SparseVec::new(r.a.clone())).collect(),
                b: out.inst.constraints.iter().map(|r| SparseVec::new(r.b.clone())).collect(),
                c: out.inst.constraints.iter().map(|r| SparseVec::new(r.c.clone())).collect(),
            };
            eprintln!("[test_large_trace] dr1cs->sparse: {:?}", t3.elapsed());
            let k_rows = inst_sparse.k();
            let ell_rs = 2 * k_rows;
            let l_public = out.public_len;
            let flpcp = dpp::dr1cs_flpcp::RsDr1csNpFlpcpSparse::<FSmall>::new(inst_sparse, l_public, ell_rs);

            let t4 = std::time::Instant::now();
            let x_small = out.assignment[..l_public].to_vec();
            let z_w_small = out.assignment[l_public..].to_vec();
            eprintln!(
                "[test_large_trace] sizes: l_public={} witness_len={} assignment_len={}",
                l_public,
                z_w_small.len(),
                out.assignment.len()
            );
            let pi_field = flpcp.prove(&x_small, &z_w_small);
            eprintln!("[test_large_trace] flpcp.prove: {:?}", t4.elapsed());
            let t5 = std::time::Instant::now();
            let boolized = dpp::BooleanProofFlpcpSparse::<FSmall, _>::new(flpcp.clone());
            let pi_bits = boolized.encode_proof_bits(&pi_field);
            eprintln!(
                "[test_large_trace] booleanize(pi): {:?} (pi_bits_len={})",
                t5.elapsed(),
                pi_bits.len()
            );

            let t6 = std::time::Instant::now();
            let dppv = build_rev2_dpp_sparse_boolean_auto::<FSmall, FBig, _>(
                flpcp,
                dpp::EmbeddingParams { gamma: 2, assume_boolean_proof: true, k_prime: 0 },
            )
            .expect("build dpp");
            eprintln!("[test_large_trace] build_rev2_dpp: {:?}", t6.elapsed());

            let t_xbig = std::time::Instant::now();
            let x_big = x_small.iter().copied().map(lift_to_big::<FSmall>).collect::<Vec<_>>();
            eprintln!("[test_large_trace] lift x_small->x_big: {:?}", t_xbig.elapsed());

            let vk_hash = [1u8; 32];
            let r1cs_digest = [2u8; 32];
            let stmt_digest = we_statement_hash_lf_plus::<RR>(vk_hash, r1cs_digest, LFP_WE_GATE_DIGEST_V1, &[sp1_digest]);
            const ARMER_SEED: [u8; 32] = *b"LFP_ARMER_SEED_V1_00000000000000";
            let lock_j: u64 = 0;
            let coin_seed: [u8; 32] = {
                let mut h = Sha256::new();
                h.update(b"LFP_LOCK_COIN_V1");
                h.update(&ARMER_SEED);
                h.update(&stmt_digest);
                h.update(&lock_j.to_le_bytes());
                h.finalize().into()
            };
            let mut rng = StdRng::from_seed(coin_seed);
            let t7 = std::time::Instant::now();
            let q = dppv.sample_query(&mut rng, &x_big).expect("sample_query");
            eprintln!("[test_large_trace] dpp.sample_query: {:?}", t7.elapsed());

            // `pi_bits` is boolean; lifting via bytes is extremely slow. Use a fast 0/1 lift.
            let t_pibig = std::time::Instant::now();
            #[cfg(feature = "parallel")]
            let pi_big = pi_bits.par_iter().copied().map(lift_bit_to_big::<FSmall>).collect::<Vec<_>>();
            #[cfg(not(feature = "parallel"))]
            let pi_big = pi_bits.iter().copied().map(lift_bit_to_big::<FSmall>).collect::<Vec<_>>();
            eprintln!("[test_large_trace] lift pi_bits->pi_big: {:?}", t_pibig.elapsed());

            let t8 = std::time::Instant::now();
            let ok = dppv.verify_with_query(&x_big, &pi_big, &q).expect("verify");
            eprintln!("[test_large_trace] dpp verify_with_query: {:?} ok={ok}", t8.elapsed());
        };

        run_one("toy_42", digest_toy);
        run_one("sha256->field", digest_hash);
    }
}

