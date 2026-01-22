//! Poseidon-permutation arithmetization (prime-field dR1CS skeleton).
//!
//! This is an incremental step toward arithmetizing the full WE-facing verifier relation:
//! we first arithmetize the Poseidon permutation(s) used by the transcript.

use ark_ff::{BigInteger, PrimeField};
use rayon::prelude::*;

use crate::poseidon_trace::{permute_with_round_trace, PoseidonReplayError};
use crate::poseidon_trace::{replay_ops, PoseidonReplayError as ReplayErr, PoseidonSpongeReplayResult};
use crate::transcript::PoseidonTraceOp;

#[derive(Clone, Debug)]
pub struct Constraint<F: PrimeField> {
    pub a: Vec<(F, usize)>,
    pub b: Vec<(F, usize)>,
    pub c: Vec<(F, usize)>,
}

#[derive(Clone, Debug)]
pub struct SparseDr1csInstance<F: PrimeField> {
    pub nvars: usize,
    pub constraints: Vec<Constraint<F>>,
}

impl<F: PrimeField> SparseDr1csInstance<F> {
    pub fn eval_lc(terms: &[(F, usize)], assignment: &[F]) -> F {
        if terms.len() > 64 {
            terms
                .par_iter()
                .map(|(c, idx)| *c * assignment[*idx])
                .reduce(|| F::ZERO, |a, b| a + b)
        } else {
        terms
            .iter()
            .fold(F::ZERO, |acc, (c, idx)| acc + (*c * assignment[*idx]))
        }
    }

    pub fn check(&self, assignment: &[F]) -> Result<(), String> {
        if assignment.len() != self.nvars {
            return Err(format!(
                "assignment length mismatch: expected {}, got {}",
                self.nvars,
                assignment.len()
            ));
        }

        let failed = self
            .constraints
            .par_iter()
            .enumerate()
            .find_any(|(_, row)| {
                let a = Self::eval_lc(&row.a, assignment);
                let b = Self::eval_lc(&row.b, assignment);
                let c = Self::eval_lc(&row.c, assignment);
                a * b != c
            });

        if let Some((i, _)) = failed {
            return Err(format!("constraint {i} failed"));
        }
        Ok(())
    }

    /// Convert to the dense dR1CS format used by the prototype RS FLPCP.
    pub fn to_dense(&self) -> dpp::dr1cs_flpcp::Dr1csInstance<F> {
        let n = self.nvars;
        
        let rows: Vec<(Vec<F>, Vec<F>, Vec<F>)> = self.constraints.par_iter()
            .map(|row| {
            let mut ra = vec![F::ZERO; n];
            let mut rb = vec![F::ZERO; n];
            let mut rc = vec![F::ZERO; n];
            for (coeff, idx) in &row.a {
                ra[*idx] += *coeff;
            }
            for (coeff, idx) in &row.b {
                rb[*idx] += *coeff;
            }
            for (coeff, idx) in &row.c {
                rc[*idx] += *coeff;
            }
                (ra, rb, rc)
            })
            .collect();
        
        let (a, b, c): (Vec<_>, Vec<_>, Vec<_>) = rows.into_iter()
            .map(|(ra, rb, rc)| (ra, rb, rc))
            .fold((Vec::new(), Vec::new(), Vec::new()), |(mut a, mut b, mut c), (ra, rb, rc)| {
            a.push(ra);
            b.push(rb);
            c.push(rc);
                (a, b, c)
            });
        
        dpp::dr1cs_flpcp::Dr1csInstance { a, b, c }
    }
}

/// Merge multiple sparse dR1CS instances into one, sharing variable 0 as the constant-1 slot.
///
/// Each part is assumed to have `assignment[0] = 1`. In the merged instance:
/// - var 0 is shared across all parts
/// - all other variables are appended, and constraints are re-indexed accordingly
pub fn merge_sparse_dr1cs_share_one<F: PrimeField>(
    parts: &[(SparseDr1csInstance<F>, Vec<F>)],
) -> Result<(SparseDr1csInstance<F>, Vec<F>), String> {
    if parts.is_empty() {
        return Err("merge_sparse_dr1cs_share_one: empty parts".to_string());
    }

    let mut merged_assignment: Vec<F> = vec![F::ONE];
    let mut merged_constraints: Vec<Constraint<F>> = Vec::new();

    for (inst, asg) in parts {
        if asg.is_empty() || asg[0] != F::ONE {
            return Err("merge_sparse_dr1cs_share_one: each part must have assignment[0]=1".to_string());
        }
        if inst.nvars != asg.len() {
            return Err("merge_sparse_dr1cs_share_one: inst/assignment length mismatch".to_string());
        }

        // Map part var0 -> merged var0, and shift the rest by current offset.
        let offset = merged_assignment.len() - 1;
        let remap_idx = |idx: usize| -> usize { if idx == 0 { 0 } else { idx + offset } };

        for row in &inst.constraints {
            let remap_lc = |lc: &[(F, usize)]| -> Vec<(F, usize)> {
                lc.iter().map(|(c, i)| (*c, remap_idx(*i))).collect()
            };
            merged_constraints.push(Constraint {
                a: remap_lc(&row.a),
                b: remap_lc(&row.b),
                c: remap_lc(&row.c),
            });
        }

        // Append assignment sans constant slot.
        merged_assignment.extend_from_slice(&asg[1..]);
    }

    Ok((
        SparseDr1csInstance {
            nvars: merged_assignment.len(),
            constraints: merged_constraints,
        },
        merged_assignment,
    ))
}

/// Merge multiple sparse dR1CS instances (sharing var 0) and add *glue* equality constraints
/// between variables belonging to different parts.
///
/// `glue` entries are `(part_a, var_a, part_b, var_b)` in **local** indices.
/// The merged instance enforces `var_a == var_b` for each glue entry.
pub fn merge_sparse_dr1cs_share_one_with_glue<F: PrimeField>(
    parts: &[(SparseDr1csInstance<F>, Vec<F>)],
    glue: &[(usize, usize, usize, usize)],
) -> Result<(SparseDr1csInstance<F>, Vec<F>), String> {
    if parts.is_empty() {
        return Err("merge_sparse_dr1cs_share_one_with_glue: empty parts".to_string());
    }

    // Compute offsets for each part (how much its non-const vars are shifted by in merged space).
    let mut offsets: Vec<usize> = Vec::with_capacity(parts.len());
    let mut merged_assignment: Vec<F> = vec![F::ONE];
    for (inst, asg) in parts {
        if asg.is_empty() || asg[0] != F::ONE {
            return Err("merge_sparse_dr1cs_share_one_with_glue: each part must have assignment[0]=1".to_string());
        }
        if inst.nvars != asg.len() {
            return Err("merge_sparse_dr1cs_share_one_with_glue: inst/assignment length mismatch".to_string());
        }
        offsets.push(merged_assignment.len() - 1);
        merged_assignment.extend_from_slice(&asg[1..]);
    }

    let remap_global = |part_idx: usize, local: usize, offsets: &[usize]| -> usize {
        if local == 0 { 0 } else { local + offsets[part_idx] }
    };

    let total_constraints: usize =
        parts.iter().map(|(inst, _)| inst.constraints.len()).sum::<usize>() + glue.len();
    let mut merged_constraints: Vec<Constraint<F>> = Vec::with_capacity(total_constraints);

    // Merge constraints with remapped indices.
    //
    // Default is a single-pass, allocation-friendly merge.
    //
    // For very large gates (multi-million constraints), this merge can become a wall-time
    // bottleneck and appear “single-core bound” in system monitors. In that case, we allow a
    // parallel remap that preserves constraint order and does not allocate *extra* constraints
    // beyond the final merged list (it still must allocate the final `merged_constraints`).
    let use_parallel_merge =
        total_constraints >= 2_000_000;

    if use_parallel_merge {
        let remapped_parts: Vec<Vec<Constraint<F>>> = parts
            .par_iter()
            .enumerate()
            .map(|(part_idx, (inst, _asg))| {
                let offset = offsets[part_idx];
                let remap_idx = |idx: usize| -> usize { if idx == 0 { 0 } else { idx + offset } };
                let remap_lc = |lc: &[(F, usize)]| -> Vec<(F, usize)> {
                    let mut out = Vec::with_capacity(lc.len());
                    for (c, i) in lc {
                        out.push((*c, remap_idx(*i)));
                    }
                    out
                };
                inst.constraints
                    .par_iter()
                    .map(|row| Constraint {
                        a: remap_lc(&row.a),
                        b: remap_lc(&row.b),
                        c: remap_lc(&row.c),
                    })
                    .collect::<Vec<_>>()
            })
            .collect();

        for v in remapped_parts {
            merged_constraints.extend(v);
        }
    } else {
        for (part_idx, (inst, _asg)) in parts.iter().enumerate() {
            let offset = offsets[part_idx];
            let remap_idx = |idx: usize| -> usize { if idx == 0 { 0 } else { idx + offset } };
            for row in &inst.constraints {
                let remap_lc = |lc: &[(F, usize)]| -> Vec<(F, usize)> {
                    let mut out = Vec::with_capacity(lc.len());
                    for (c, i) in lc {
                        out.push((*c, remap_idx(*i)));
                    }
                    out
                };
                merged_constraints.push(Constraint {
                    a: remap_lc(&row.a),
                    b: remap_lc(&row.b),
                    c: remap_lc(&row.c),
                });
            }
        }
    }

    // Add glue constraints: (x - y) * 1 = 0.
    for &(pa, xa, pb, xb) in glue {
        if pa >= parts.len() || pb >= parts.len() {
            return Err("merge_sparse_dr1cs_share_one_with_glue: glue part idx out of range".to_string());
        }
        let ga = remap_global(pa, xa, &offsets);
        let gb = remap_global(pb, xb, &offsets);
        merged_constraints.push(Constraint {
            a: vec![(F::ONE, ga), (-F::ONE, gb)],
            b: vec![(F::ONE, 0)],
            c: vec![(F::ZERO, 0)],
        });
    }

    Ok((
        SparseDr1csInstance {
            nvars: merged_assignment.len(),
            constraints: merged_constraints,
        },
        merged_assignment,
    ))
}

#[derive(Clone, Debug)]
struct Dr1csBuilder<F: PrimeField> {
    assignment: Vec<F>,
    rows: Vec<Constraint<F>>,
}

impl<F: PrimeField> Dr1csBuilder<F> {
    fn new() -> Self {
        // var 0 is the constant-1 slot
        Self {
            assignment: vec![F::ONE],
            rows: Vec::new(),
        }
    }

    fn one(&self) -> usize {
        0
    }

    fn new_var(&mut self, value: F) -> usize {
        let idx = self.assignment.len();
        self.assignment.push(value);
        idx
    }

    fn add_constraint(&mut self, a: Vec<(F, usize)>, b: Vec<(F, usize)>, c: Vec<(F, usize)>) {
        self.rows.push(Constraint { a, b, c });
    }

    fn enforce_lc_times_one_eq_var(&mut self, lc: Vec<(F, usize)>, out: usize) {
        self.add_constraint(lc, vec![(F::ONE, self.one())], vec![(F::ONE, out)]);
    }

    fn enforce_var_eq_const(&mut self, x: usize, c: F) {
        // x * 1 = c
        self.add_constraint(vec![(F::ONE, x)], vec![(F::ONE, self.one())], vec![(c, self.one())]);
    }

    fn enforce_var_eq_var(&mut self, x: usize, y: usize) {
        self.enforce_lc_times_one_eq_var(vec![(F::ONE, x)], y);
    }

    fn enforce_mul(&mut self, x: usize, y: usize, out: usize) {
        self.add_constraint(vec![(F::ONE, x)], vec![(F::ONE, y)], vec![(F::ONE, out)]);
    }

    fn enforce_pow_u64(&mut self, base: usize, alpha: u64) -> usize {
        // Build constraints for base^alpha over variables.
        //
        // IMPORTANT:
        // This gadget is used heavily by Poseidon (S-box exponentiation). For common small exponents,
        // we use short addition chains to minimize multiplication constraints.
        if alpha == 0 {
            return self.new_var(F::ONE);
        }
        if alpha == 1 {
            return base;
        }
        // Common Poseidon alphas.
        // - alpha=3: x^3 = (x^2)*x (2 muls)
        // - alpha=5: x^5 = (x^2)^2 * x (3 muls)
        // - alpha=7: x^7 = (x^2)^2 * x^2 * x (4 muls)
        if alpha == 3 {
            let x2_val = self.assignment[base] * self.assignment[base];
            let x2 = self.new_var(x2_val);
            self.enforce_mul(base, base, x2);
            let x3_val = x2_val * self.assignment[base];
            let x3 = self.new_var(x3_val);
            self.enforce_mul(x2, base, x3);
            return x3;
        }
        if alpha == 5 {
            let x = base;
            let x_val = self.assignment[x];
            let x2_val = x_val * x_val;
            let x2 = self.new_var(x2_val);
            self.enforce_mul(x, x, x2);
            let x4_val = x2_val * x2_val;
            let x4 = self.new_var(x4_val);
            self.enforce_mul(x2, x2, x4);
            let x5_val = x4_val * x_val;
            let x5 = self.new_var(x5_val);
            self.enforce_mul(x4, x, x5);
            return x5;
        }
        if alpha == 7 {
            let x = base;
            let x_val = self.assignment[x];
            let x2_val = x_val * x_val;
            let x2 = self.new_var(x2_val);
            self.enforce_mul(x, x, x2);
            let x4_val = x2_val * x2_val;
            let x4 = self.new_var(x4_val);
            self.enforce_mul(x2, x2, x4);
            let x6_val = x4_val * x2_val;
            let x6 = self.new_var(x6_val);
            self.enforce_mul(x4, x2, x6);
            let x7_val = x6_val * x_val;
            let x7 = self.new_var(x7_val);
            self.enforce_mul(x6, x, x7);
            return x7;
        }

        // current = base
        let mut cur_var = base;
        let mut cur_val = self.assignment[base];

        // acc = 1
        let mut acc_var = self.new_var(F::ONE);
        let mut acc_val = F::ONE;

        let mut e = alpha;
        while e > 0 {
            if (e & 1) == 1 {
                let out_val = acc_val * cur_val;
                let out_var = self.new_var(out_val);
                self.enforce_mul(acc_var, cur_var, out_var);
                acc_var = out_var;
                acc_val = out_val;
            }
            e >>= 1;
            if e == 0 {
                break;
            }
            // square cur
            let sq_val = cur_val * cur_val;
            let sq_var = self.new_var(sq_val);
            self.enforce_mul(cur_var, cur_var, sq_var);
            cur_var = sq_var;
            cur_val = sq_val;
        }

        acc_var
    }

    fn into_sparse_instance(self) -> (SparseDr1csInstance<F>, Vec<F>) {
        let nvars = self.assignment.len();
        let inst = SparseDr1csInstance {
            nvars,
            constraints: self.rows,
        };
        (inst, self.assignment)
    }
}

/// Build a sparse dR1CS instance for a single Poseidon permutation, given an input state.
///
/// Returns `(instance, assignment, out_state_var_indices)`.
pub fn poseidon_permutation_dr1cs<F: PrimeField>(
    cfg: &ark_crypto_primitives::sponge::poseidon::PoseidonConfig<F>,
    before_state: &[F],
) -> Result<(SparseDr1csInstance<F>, Vec<F>, Vec<usize>), PoseidonReplayError> {
    let t = cfg.rate + cfg.capacity;
    if before_state.len() != t {
        return Err(PoseidonReplayError::Invalid(format!(
            "bad state length: expected {}, got {}",
            t,
            before_state.len()
        )));
    }

    // Compute round states (ground truth for witness materialization).
    let (after_state, round_states) = permute_with_round_trace(cfg, before_state)?;

    let mut b = Dr1csBuilder::<F>::new();

    // Materialize the input state variables.
    let mut state_vars: Vec<usize> = Vec::with_capacity(t);
    for &s in before_state {
        state_vars.push(b.new_var(s));
    }

    let full_rounds_over_2 = cfg.full_rounds / 2;
    let total_rounds = cfg.full_rounds + cfg.partial_rounds;
    assert_eq!(round_states.len(), total_rounds);

    // For each round: ARK -> SBOX -> MDS.
    for r in 0..total_rounds {
        let is_full = r < full_rounds_over_2 || r >= (full_rounds_over_2 + cfg.partial_rounds);

        // ARK (affine shift).
        let mut ark_vars: Vec<usize> = Vec::with_capacity(t);
        for i in 0..t {
            let val = b.assignment[state_vars[i]] + cfg.ark[r][i];
            let v = b.new_var(val);
            // (state_i + ark*r,i*one) * 1 = v
            b.enforce_lc_times_one_eq_var(
                vec![(F::ONE, state_vars[i]), (cfg.ark[r][i], b.one())],
                v,
            );
            ark_vars.push(v);
        }

        // S-box.
        let mut sbox_vars: Vec<usize> = Vec::with_capacity(t);
        for i in 0..t {
            if is_full || i == 0 {
                let out_var = b.enforce_pow_u64(ark_vars[i], cfg.alpha);
                sbox_vars.push(out_var);
            } else {
                // identity: materialize as a fresh var with equality constraint for uniformity
                let v = b.new_var(b.assignment[ark_vars[i]]);
                b.enforce_var_eq_var(ark_vars[i], v);
                sbox_vars.push(v);
            }
        }

        // MDS.
        let mut next_vars: Vec<usize> = Vec::with_capacity(t);
        for i in 0..t {
            let mut lc: Vec<(F, usize)> = Vec::with_capacity(t);
            let mut val = F::ZERO;
            for j in 0..t {
                let coeff = cfg.mds[i][j];
                lc.push((coeff, sbox_vars[j]));
                val += coeff * b.assignment[sbox_vars[j]];
            }
            let v = b.new_var(val);
            b.enforce_lc_times_one_eq_var(lc, v);
            next_vars.push(v);
        }

        // Sanity: witness should match ground-truth round state.
        // (Recorded round state is after ARK+SBOX+MDS).
        let expected = &round_states[r];
        for i in 0..t {
            if b.assignment[next_vars[i]] != expected[i] {
                return Err(PoseidonReplayError::Mismatch(format!(
                    "round {r} state mismatch at i={i}"
                )));
            }
        }

        state_vars = next_vars;
    }

    // Final state must match.
    for i in 0..t {
        if b.assignment[state_vars[i]] != after_state[i] {
            return Err(PoseidonReplayError::Mismatch(format!(
                "final state mismatch at i={i}"
            )));
        }
    }

    let (inst, assignment) = b.into_sparse_instance();
    Ok((inst, assignment, state_vars))
}

/// Byte squeeze witness info (for later byte-decomposition constraints).
#[derive(Clone, Debug)]
pub struct ByteSqueezeWitness {
    pub n: usize,
    pub usable_bytes: usize,
    pub src_elems: Vec<usize>, // variable indices of squeezed field elements
    pub out: Vec<u8>,          // recorded bytes
}

/// Wiring information for a Poseidon transcript dR1CS instance.
///
/// This is intended for *higher-level verifier arithmetizations* that need to reference
/// specific squeezed field elements (e.g., Fiat–Shamir challenges) as variables.
#[derive(Clone, Debug, Default)]
pub struct PoseidonDr1csWiring {
    /// Flattened variable indices for all absorbed field elements (BF variables) in trace order.
    pub absorb_vars: Vec<usize>,
    /// For each `Absorb` op, `(start, len)` into `absorb_vars`.
    pub absorb_ranges: Vec<(usize, usize)>,
    /// Flattened variable indices for all `SqueezeField` outputs in trace order.
    pub squeeze_field_vars: Vec<usize>,
    /// For each `SqueezeField` op, `(start, len)` into `squeeze_field_vars`.
    pub squeeze_field_ranges: Vec<(usize, usize)>,
}

/// Wiring for `SqueezeBytes` outputs.
#[derive(Clone, Debug, Default)]
pub struct PoseidonByteWiring {
    /// Flattened variable indices for all squeezed bytes (BF variables) in trace order.
    pub squeeze_byte_vars: Vec<usize>,
    /// For each `SqueezeBytes` op, `(start, len)` into `squeeze_byte_vars`.
    pub squeeze_byte_ranges: Vec<(usize, usize)>,
}

/// Build a dR1CS instance for the *entire* Poseidon sponge transcript trace, including:
/// - permutation constraints,
/// - absorb updates (linear constraints),
/// - squeeze-field outputs (linear constraints).
///
/// For `SqueezeBytes`, we constrain the squeezed **field elements** and return `ByteSqueezeWitness`
/// entries so byte-decomposition constraints can be added in a later step.
pub fn poseidon_sponge_dr1cs_from_trace<F: PrimeField>(
    cfg: &ark_crypto_primitives::sponge::poseidon::PoseidonConfig<F>,
    ops: &[PoseidonTraceOp<F>],
) -> Result<
    (
        SparseDr1csInstance<F>,
        Vec<F>,
        PoseidonSpongeReplayResult<F>,
        Vec<ByteSqueezeWitness>,
    ),
    ReplayErr,
> {
    poseidon_sponge_dr1cs_from_trace_with_wiring(cfg, ops)
        .map(|(inst, asg, replay, bytes, _wiring)| (inst, asg, replay, bytes))
}

/// Same as `poseidon_sponge_dr1cs_from_trace`, but also returns `PoseidonDr1csWiring` describing
/// where each `SqueezeField` output element lives in the dR1CS assignment vector.
pub fn poseidon_sponge_dr1cs_from_trace_with_wiring<F: PrimeField>(
    cfg: &ark_crypto_primitives::sponge::poseidon::PoseidonConfig<F>,
    ops: &[PoseidonTraceOp<F>],
) -> Result<
    (
        SparseDr1csInstance<F>,
        Vec<F>,
        PoseidonSpongeReplayResult<F>,
        Vec<ByteSqueezeWitness>,
        PoseidonDr1csWiring,
    ),
    ReplayErr,
> {
    poseidon_sponge_dr1cs_from_trace_impl(cfg, ops, false).map(|(inst, asg, replay, bytes, wiring, _bw)| {
        (inst, asg, replay, bytes, wiring)
    })
}

/// Like `poseidon_sponge_dr1cs_from_trace_with_wiring`, but also **arithmetizes `SqueezeBytes`**:
/// - allocates byte variables,
/// - constrains each byte is 8-bit,
/// - links bytes to the underlying squeezed field elements via radix-256 decomposition.
pub fn poseidon_sponge_dr1cs_from_trace_with_wiring_and_bytes<F: PrimeField>(
    cfg: &ark_crypto_primitives::sponge::poseidon::PoseidonConfig<F>,
    ops: &[PoseidonTraceOp<F>],
) -> Result<
    (
        SparseDr1csInstance<F>,
        Vec<F>,
        PoseidonSpongeReplayResult<F>,
        Vec<ByteSqueezeWitness>,
        PoseidonDr1csWiring,
        PoseidonByteWiring,
    ),
    ReplayErr,
> {
    poseidon_sponge_dr1cs_from_trace_impl(cfg, ops, true)
}

fn poseidon_sponge_dr1cs_from_trace_impl<F: PrimeField>(
    cfg: &ark_crypto_primitives::sponge::poseidon::PoseidonConfig<F>,
    ops: &[PoseidonTraceOp<F>],
    with_bytes: bool,
) -> Result<
    (
        SparseDr1csInstance<F>,
        Vec<F>,
        PoseidonSpongeReplayResult<F>,
        Vec<ByteSqueezeWitness>,
        PoseidonDr1csWiring,
        PoseidonByteWiring,
    ),
    ReplayErr,
> {
    // First replay to ensure the ops are consistent and to get permute boundaries.
    let replay = replay_ops(cfg, ops)?;

    let t = cfg.rate + cfg.capacity;
    let mut b = Dr1csBuilder::<F>::new();
    let one = b.one();

    // Initial state is all zeros.
    let mut state_vars = Vec::with_capacity(t);
    for _ in 0..t {
        let v = b.new_var(F::ZERO);
        b.enforce_var_eq_const(v, F::ZERO);
        state_vars.push(v);
    }

    let mut mode = ark_crypto_primitives::sponge::DuplexSpongeMode::Absorbing {
        next_absorb_index: 0,
    };

    let mut byte_witnesses: Vec<ByteSqueezeWitness> = Vec::new();
    let mut wiring = PoseidonDr1csWiring::default();
    let mut byte_wiring = PoseidonByteWiring::default();

    // Helper: apply a Poseidon permutation to the current `state_vars`.
    let mut permute_ptr: usize = 0;
    let mut apply_perm = |b: &mut Dr1csBuilder<F>, state_vars: &mut Vec<usize>| -> Result<(), ReplayErr> {
        if permute_ptr >= replay.permutes.len() {
            return Err(ReplayErr::Invalid("permute ptr out of range".to_string()));
        }
        let before = replay.permutes[permute_ptr].before.clone();
        let (after_state, round_states) = permute_with_round_trace(cfg, &before)?;

        // Ensure the current state witness matches `before` (sanity; not a constraint).
        for i in 0..t {
            if b.assignment[state_vars[i]] != before[i] {
                return Err(ReplayErr::Mismatch(format!(
                    "state mismatch before permute #{permute_ptr} at i={i}"
                )));
            }
        }

        let full_rounds_over_2 = cfg.full_rounds / 2;
        let total_rounds = cfg.full_rounds + cfg.partial_rounds;
        assert_eq!(round_states.len(), total_rounds);

        for r in 0..total_rounds {
            let is_full = r < full_rounds_over_2 || r >= (full_rounds_over_2 + cfg.partial_rounds);

            let mut next_vars: Vec<usize> = Vec::with_capacity(t);
            if is_full {
                // Full rounds: ARK -> SBOX on all lanes -> MDS.
                //
                // Optimization (safe): avoid materializing `ark_vars[i] = state[i] + ark[r][i]`
                // as separate linear constraints/variables. Instead, feed the affine form
                // `(state[i] + const)` directly into the S-box exponentiation constraints as
                // a linear combination.
                //
                // This reduces ~t linear constraints per full round (and their vars), which is
                // a meaningful fraction of Poseidon cost for large permute counts.

                #[inline]
                fn eval_lc<F: PrimeField>(b: &Dr1csBuilder<F>, lc: &[(F, usize)]) -> F {
                    lc.iter()
                        .fold(F::ZERO, |acc, (c, idx)| acc + (*c * b.assignment[*idx]))
                }

                #[inline]
                fn enforce_mul_lc_lc<F: PrimeField>(
                    b: &mut Dr1csBuilder<F>,
                    a: Vec<(F, usize)>,
                    bb: Vec<(F, usize)>,
                ) -> usize {
                    let out_val = eval_lc::<F>(b, &a) * eval_lc::<F>(b, &bb);
                    let out = b.new_var(out_val);
                    b.add_constraint(a, bb, vec![(F::ONE, out)]);
                    out
                }

                #[inline]
                fn enforce_mul_var_lc<F: PrimeField>(
                    b: &mut Dr1csBuilder<F>,
                    x: usize,
                    lc: Vec<(F, usize)>,
                ) -> usize {
                    let out_val = b.assignment[x] * eval_lc::<F>(b, &lc);
                    let out = b.new_var(out_val);
                    b.add_constraint(vec![(F::ONE, x)], lc, vec![(F::ONE, out)]);
                    out
                }

                // S-box outputs per lane.
                let mut sbox_vars: Vec<usize> = Vec::with_capacity(t);
                for i in 0..t {
                    let ark = cfg.ark[r][i];
                    // lc_in = state[i] + ark
                    let lc_in = if ark.is_zero() {
                        vec![(F::ONE, state_vars[i])]
                    } else {
                        vec![(F::ONE, state_vars[i]), (ark, one)]
                    };

                    let out = match cfg.alpha {
                        3 => {
                            // x^3 = (x^2)*x
                            let x2 = enforce_mul_lc_lc::<F>(b, lc_in.clone(), lc_in.clone());
                            enforce_mul_var_lc::<F>(b, x2, lc_in)
                        }
                        5 => {
                            // x^5 = (x^2)^2 * x
                            let x2 = enforce_mul_lc_lc::<F>(b, lc_in.clone(), lc_in.clone());
                            let x4_val = b.assignment[x2] * b.assignment[x2];
                            let x4 = b.new_var(x4_val);
                            b.enforce_mul(x2, x2, x4);
                            enforce_mul_var_lc::<F>(b, x4, lc_in)
                        }
                        7 => {
                            // x^7 = (x^2)^2 * x^2 * x
                            let x2 = enforce_mul_lc_lc::<F>(b, lc_in.clone(), lc_in.clone());
                            let x4_val = b.assignment[x2] * b.assignment[x2];
                            let x4 = b.new_var(x4_val);
                            b.enforce_mul(x2, x2, x4);
                            let x6_val = x4_val * b.assignment[x2];
                            let x6 = b.new_var(x6_val);
                            b.enforce_mul(x4, x2, x6);
                            enforce_mul_var_lc::<F>(b, x6, lc_in)
                        }
                        _ => {
                            // Fallback: materialize ARK as a var and use the generic exponentiation gadget.
                            let val = b.assignment[state_vars[i]] + ark;
                            let v = b.new_var(val);
                            b.enforce_lc_times_one_eq_var(vec![(F::ONE, state_vars[i]), (ark, one)], v);
                            b.enforce_pow_u64(v, cfg.alpha)
                        }
                    };
                    sbox_vars.push(out);
                }

                for i in 0..t {
                    let mut lc: Vec<(F, usize)> = Vec::with_capacity(t);
                    let mut val = F::ZERO;
                    for j in 0..t {
                        let coeff = cfg.mds[i][j];
                        lc.push((coeff, sbox_vars[j]));
                        val += coeff * b.assignment[sbox_vars[j]];
                    }
                    let v = b.new_var(val);
                    b.enforce_lc_times_one_eq_var(lc, v);
                    next_vars.push(v);
                }
            } else {
                // Partial rounds: only lane 0 is S-boxed.
                //
                // Optimization: avoid materializing ARK vars (and equality vars) for the identity lanes.
                // We fold their ARK constants directly into the MDS linear constraints:
                //   input[j] = state[j] + ark[r][j] for j>0, and input[0] = SBOX(state[0] + ark[r][0]).
                let ark0_val = b.assignment[state_vars[0]] + cfg.ark[r][0];
                let ark0 = b.new_var(ark0_val);
                b.enforce_lc_times_one_eq_var(vec![(F::ONE, state_vars[0]), (cfg.ark[r][0], one)], ark0);
                let sbox0 = b.enforce_pow_u64(ark0, cfg.alpha);

                for i in 0..t {
                    // y_i = Σ_j mds[i][j] * input[j]
                    let mut lc: Vec<(F, usize)> = Vec::with_capacity(t + 1);
                    let mut val = F::ZERO;
                    // j=0 nonlinear lane
                    let c0 = cfg.mds[i][0];
                    lc.push((c0, sbox0));
                    val += c0 * b.assignment[sbox0];

                    // j>0 linear lanes: use state var directly, and fold constants into the LC.
                    let mut const_term = F::ZERO;
                    for j in 1..t {
                        let coeff = cfg.mds[i][j];
                        lc.push((coeff, state_vars[j]));
                        val += coeff * b.assignment[state_vars[j]];
                        const_term += coeff * cfg.ark[r][j];
                    }
                    if !const_term.is_zero() {
                        lc.push((const_term, one));
                        val += const_term;
                    }
                    let v = b.new_var(val);
                    b.enforce_lc_times_one_eq_var(lc, v);
                    next_vars.push(v);
                }
            }

            // Sanity against traced round state.
            let expected = &round_states[r];
            for i in 0..t {
                if b.assignment[next_vars[i]] != expected[i] {
                    return Err(ReplayErr::Mismatch(format!(
                        "round {r} state mismatch at i={i} for permute #{permute_ptr}"
                    )));
                }
            }
            *state_vars = next_vars;
        }

        // Final sanity against permute-after.
        for i in 0..t {
            if b.assignment[state_vars[i]] != after_state[i] {
                return Err(ReplayErr::Mismatch(format!(
                    "after state mismatch at i={i} for permute #{permute_ptr}"
                )));
            }
        }
        permute_ptr += 1;
        Ok(())
    };

    for op in ops {
        match op {
            PoseidonTraceOp::Absorb(elems) => {
                if elems.is_empty() {
                    continue;
                }
                let range_start = wiring.absorb_vars.len();
                let mut range_len = 0usize;
                for &e in elems {
                    // If we were squeezing, permute first.
                    if matches!(mode, ark_crypto_primitives::sponge::DuplexSpongeMode::Squeezing { .. }) {
                        apply_perm(&mut b, &mut state_vars)?;
                        mode = ark_crypto_primitives::sponge::DuplexSpongeMode::Absorbing { next_absorb_index: 0 };
                    }
                    let mut absorb_index = match mode {
                        ark_crypto_primitives::sponge::DuplexSpongeMode::Absorbing { next_absorb_index } => next_absorb_index,
                        _ => unreachable!(),
                    };
                    if absorb_index == cfg.rate {
                        apply_perm(&mut b, &mut state_vars)?;
                        absorb_index = 0;
                    }

                    // Materialize element as a fixed var (so it can be part of the witness vector).
                    let e_var = b.new_var(e);
                    b.enforce_var_eq_const(e_var, e);
                    wiring.absorb_vars.push(e_var);
                    range_len += 1;

                    // Update one state slot: state[cap + absorb_index] += e
                    let pos = cfg.capacity + absorb_index;
                    let new_val = b.assignment[state_vars[pos]] + b.assignment[e_var];
                    let new_var = b.new_var(new_val);
                    b.enforce_lc_times_one_eq_var(vec![(F::ONE, state_vars[pos]), (F::ONE, e_var)], new_var);
                    state_vars[pos] = new_var;

                    mode = ark_crypto_primitives::sponge::DuplexSpongeMode::Absorbing {
                        next_absorb_index: absorb_index + 1,
                    };
                }
                wiring.absorb_ranges.push((range_start, range_len));
            }
            PoseidonTraceOp::SqueezeField(out) => {
                if out.is_empty() {
                    continue;
                }
                // If we were absorbing, permute first.
                let mut squeeze_index = match mode {
                    ark_crypto_primitives::sponge::DuplexSpongeMode::Absorbing { .. } => {
                        apply_perm(&mut b, &mut state_vars)?;
                        0
                    }
                    ark_crypto_primitives::sponge::DuplexSpongeMode::Squeezing { next_squeeze_index } => next_squeeze_index,
                };
                if squeeze_index == cfg.rate {
                    apply_perm(&mut b, &mut state_vars)?;
                    squeeze_index = 0;
                }

                let mut produced = 0usize;
                let range_start = wiring.squeeze_field_vars.len();
                while produced < out.len() {
                    let take = core::cmp::min(cfg.rate - squeeze_index, out.len() - produced);
                    for j in 0..take {
                        let pos = cfg.capacity + squeeze_index + j;
                        let expected = out[produced + j];
                        let v = b.new_var(expected);
                        b.enforce_var_eq_const(v, expected);
                        // v == state[pos]
                        b.enforce_var_eq_var(state_vars[pos], v);
                        wiring.squeeze_field_vars.push(v);
                    }
                    produced += take;
                    squeeze_index += take;
                    if produced < out.len() && squeeze_index == cfg.rate {
                        apply_perm(&mut b, &mut state_vars)?;
                        squeeze_index = 0;
                    }
                }
                wiring.squeeze_field_ranges.push((range_start, out.len()));

                mode = ark_crypto_primitives::sponge::DuplexSpongeMode::Squeezing {
                    next_squeeze_index: squeeze_index,
                };
            }
            PoseidonTraceOp::SqueezeBytes { n, out } => {
                let usable_bytes = ((F::MODULUS_BIT_SIZE - 1) / 8) as usize;
                let num_elements = (*n + usable_bytes - 1) / usable_bytes;

                // Squeeze native field elements and constrain them (like SqueezeField),
                // then check the bytes in Rust and return them for later constraints.
                let mut squeeze_index = match mode {
                    ark_crypto_primitives::sponge::DuplexSpongeMode::Absorbing { .. } => {
                        apply_perm(&mut b, &mut state_vars)?;
                        0
                    }
                    ark_crypto_primitives::sponge::DuplexSpongeMode::Squeezing { next_squeeze_index } => next_squeeze_index,
                };
                if squeeze_index == cfg.rate {
                    apply_perm(&mut b, &mut state_vars)?;
                    squeeze_index = 0;
                }

                let mut src_vars: Vec<usize> = Vec::with_capacity(num_elements);
                let mut produced = 0usize;
                while produced < num_elements {
                    let take = core::cmp::min(cfg.rate - squeeze_index, num_elements - produced);
                    for j in 0..take {
                        let pos = cfg.capacity + squeeze_index + j;
                        // Use the state var directly as the squeezed element source.
                        src_vars.push(state_vars[pos]);
                    }
                    produced += take;
                    squeeze_index += take;
                    if produced < num_elements && squeeze_index == cfg.rate {
                        apply_perm(&mut b, &mut state_vars)?;
                        squeeze_index = 0;
                    }
                }

                // Compute bytes from witness values and check they match recorded.
                let mut bytes: Vec<u8> = Vec::with_capacity(usable_bytes * num_elements);
                for &v in &src_vars {
                    let elem_bytes = b.assignment[v].into_bigint().to_bytes_le();
                    bytes.extend_from_slice(&elem_bytes[..usable_bytes]);
                }
                bytes.truncate(*n);
                if &bytes != out {
                    return Err(ReplayErr::Mismatch("SqueezeBytes bytes mismatch".to_string()));
                }

                byte_witnesses.push(ByteSqueezeWitness {
                    n: *n,
                    usable_bytes,
                    src_elems: src_vars.clone(),
                    out: out.clone(),
                });

                if with_bytes {
                    // Allocate and constrain byte vars + link to src elements.
                    //
                    // We allocate the **full** usable_bytes*num_elements bytes (even if n truncates
                    // mid-element) so that every src element is linked to its low bytes.
                    let range_start = byte_wiring.squeeze_byte_vars.len();
                    let full_len = usable_bytes * num_elements;

                    // Build full bytes from current witness values (untruncated).
                    let mut full_bytes: Vec<u8> = Vec::with_capacity(full_len);
                    for &v in &src_vars {
                        let elem_bytes = b.assignment[v].into_bigint().to_bytes_le();
                        full_bytes.extend_from_slice(&elem_bytes[..usable_bytes]);
                    }
                    debug_assert_eq!(full_bytes.len(), full_len);

                    // Allocate byte variables in the same order as `full_bytes`.
                    // We will expose the first `n` bytes as transcript output vars (linked to state).
                    let mut op_byte_vars: Vec<usize> = Vec::with_capacity(full_len);

                    // Constants for radix-256 recomposition.
                    let mut pow256: Vec<F> = Vec::with_capacity(usable_bytes);
                    let mut acc = F::ONE;
                    let base = F::from(256u64);
                    for _ in 0..usable_bytes {
                        pow256.push(acc);
                        acc *= base;
                    }
                    let pow256k = acc; // 256^{usable_bytes}

                    for e in 0..num_elements {
                        // Allocate byte vars for this element.
                        //
                        // IMPORTANT: these bytes are already fixed by the recorded trace
                        // (`enforce_var_eq_const`), so additional bit-decomposition constraints
                        // are unnecessary for soundness of this arithmetization and are very costly.
                        let mut byte_vars: Vec<usize> = Vec::with_capacity(usable_bytes);
                        for i in 0..usable_bytes {
                            let bval = F::from(full_bytes[e * usable_bytes + i] as u64);
                            let bv = b.new_var(bval);
                            b.enforce_var_eq_const(bv, bval);
                            byte_vars.push(bv);
                            op_byte_vars.push(bv);
                        }

                        // Link src element: src = Σ 256^i * byte_i + 256^k * high
                        let src = src_vars[e];
                        let mut low = F::ZERO;
                        for i in 0..usable_bytes {
                            low += pow256[i] * b.assignment[byte_vars[i]];
                        }
                        let high_val = (b.assignment[src] - low) * pow256k.inverse().unwrap();
                        let high = b.new_var(high_val);
                        // src - Σ 256^i*byte_i - 256^k*high = 0
                        let mut lc: Vec<(F, usize)> = Vec::with_capacity(2 + usable_bytes);
                        lc.push((F::ONE, src));
                        for i in 0..usable_bytes {
                            lc.push((-pow256[i], byte_vars[i]));
                        }
                        lc.push((-pow256k, high));
                        let z = b.new_var(F::ZERO);
                        b.enforce_var_eq_const(z, F::ZERO);
                        b.enforce_lc_times_one_eq_var(lc, z);
                    }

                    // Expose only the first n bytes in trace order (truncated), as transcript output vars.
                    // These are already linked to the squeezed field elements through the constraints above.
                    for i in 0..*n {
                        byte_wiring.squeeze_byte_vars.push(op_byte_vars[i]);
                    }
                    byte_wiring.squeeze_byte_ranges.push((range_start, *n));
                }

                mode = ark_crypto_primitives::sponge::DuplexSpongeMode::Squeezing {
                    next_squeeze_index: squeeze_index,
                };
            }
        }
    }

    if permute_ptr != replay.permutes.len() {
        return Err(ReplayErr::Invalid(format!(
            "permute count mismatch: used {permute_ptr}, replay has {}",
            replay.permutes.len()
        )));
    }

    let (inst, assignment) = b.into_sparse_instance();
    Ok((inst, assignment, replay, byte_witnesses, wiring, byte_wiring))
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::UniformRand;
    use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
    use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as R;
    use stark_rings::PolyRing;

    #[test]
    fn test_poseidon_perm_dr1cs_satisfies() {
        type BF = <<R as PolyRing>::BaseRing as ark_ff::Field>::BasePrimeField;
        let cfg = <PC as cyclotomic_rings::rings::GetPoseidonParams<BF>>::get_poseidon_config();

        let t = cfg.rate + cfg.capacity;
        let mut rng = ark_std::test_rng();
        let before = (0..t).map(|_| BF::rand(&mut rng)).collect::<Vec<_>>();

        let (inst, assignment, _out_vars) =
            poseidon_permutation_dr1cs::<BF>(&cfg, &before).unwrap();
        inst.check(&assignment).unwrap();
    }

    #[test]
    fn test_poseidon_sponge_dr1cs_from_real_trace_satisfies_constraints() {
        use ark_crypto_primitives::sponge::{
            poseidon::PoseidonSponge, CryptographicSponge, FieldBasedCryptographicSponge,
        };
        use stark_rings::PolyRing;

        type BF = <<R as PolyRing>::BaseRing as ark_ff::Field>::BasePrimeField;
        let cfg = <PC as cyclotomic_rings::rings::GetPoseidonParams<BF>>::get_poseidon_config();

        // Build a small synthetic ops trace directly from PoseidonSponge.
        let mut rng = ark_std::test_rng();
        let mut sponge = PoseidonSponge::<BF>::new(&cfg);
        let mut ops: Vec<PoseidonTraceOp<BF>> = Vec::new();

        let absorb1 = (0..(cfg.rate + 3)).map(|_| BF::rand(&mut rng)).collect::<Vec<_>>();
        sponge.absorb(&absorb1);
        ops.push(PoseidonTraceOp::Absorb(absorb1.clone()));

        let out1 = sponge.squeeze_field_elements::<BF>(5);
        ops.push(PoseidonTraceOp::SqueezeField(out1.clone()));
        sponge.absorb(&out1);
        ops.push(PoseidonTraceOp::Absorb(out1));

        let bytes = sponge.squeeze_bytes(17);
        ops.push(PoseidonTraceOp::SqueezeBytes { n: 17, out: bytes });

        let absorb2 = (0..7).map(|_| BF::rand(&mut rng)).collect::<Vec<_>>();
        sponge.absorb(&absorb2);
        ops.push(PoseidonTraceOp::Absorb(absorb2));

        let out2 = sponge.squeeze_native_field_elements(3);
        ops.push(PoseidonTraceOp::SqueezeField(out2.clone()));

        let (inst, assignment, _replay, _bytes) =
            poseidon_sponge_dr1cs_from_trace::<BF>(&cfg, &ops).expect("build dr1cs failed");
        inst.check(&assignment).unwrap();
    }
}

