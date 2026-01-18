//! Sumcheck verifier arithmetization (sparse dR1CS over a prime field).
//!
//! This module is an intermediate building block toward arithmetizing the full Π_fold verifier.
//! It encodes **verification** of the multilinear sumcheck proof messages as constraints.
//!
//! Current scope:
//! - degree-3 prover messages (4 evaluations at points 0,1,2,3), matching the usage in Symphony.
//! - shared-randomness schedule for two sumchecks (same `r_i` per round).
//!
//! Note: this does **not** (yet) encode the post-sumcheck algebraic checks (Eq(26), monomial
//! recomputation, Step-5). Those are the next layer on top of the verified subclaims.

use ark_ff::PrimeField;

use crate::dpp_poseidon::{Constraint, SparseDr1csInstance};

#[derive(Clone, Debug)]
pub struct Dr1csBuilder<F: PrimeField> {
    pub assignment: Vec<F>,
    pub rows: Vec<Constraint<F>>,
}

impl<F: PrimeField> Dr1csBuilder<F> {
    pub fn new() -> Self {
        Self { assignment: vec![F::ONE], rows: Vec::new() }
    }
    pub fn one(&self) -> usize { 0 }
    pub fn new_var(&mut self, value: F) -> usize {
        let idx = self.assignment.len();
        self.assignment.push(value);
        idx
    }
    pub fn add_constraint(&mut self, a: Vec<(F, usize)>, b: Vec<(F, usize)>, c: Vec<(F, usize)>) {
        self.rows.push(Constraint { a, b, c });
    }
    pub fn enforce_lc_times_one_eq_const(&mut self, lc: Vec<(F, usize)>) {
        self.add_constraint(lc, vec![(F::ONE, self.one())], vec![(F::ZERO, self.one())]);
    }
    pub fn enforce_var_eq_const(&mut self, x: usize, c: F) {
        self.add_constraint(vec![(F::ONE, x)], vec![(F::ONE, self.one())], vec![(c, self.one())]);
    }
    pub fn enforce_mul(&mut self, x: usize, y: usize, out: usize) {
        self.add_constraint(vec![(F::ONE, x)], vec![(F::ONE, y)], vec![(F::ONE, out)]);
    }
    pub fn into_instance(self) -> (SparseDr1csInstance<F>, Vec<F>) {
        let inst = SparseDr1csInstance { nvars: self.assignment.len(), constraints: self.rows };
        (inst, self.assignment)
    }
}

/// A "ring element" represented as `d` prime-field variables (coefficients).
#[derive(Clone, Debug)]
pub struct RingVars {
    pub coeffs: Vec<usize>,
}

impl RingVars {
    pub fn new(coeffs: Vec<usize>) -> Self { Self { coeffs } }
    pub fn d(&self) -> usize { self.coeffs.len() }
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

fn scalar_sub_const<F: PrimeField>(b: &mut Dr1csBuilder<F>, r: usize, c: F) -> usize {
    let val = b.assignment[r] - c;
    let v = b.new_var(val);
    b.add_constraint(vec![(F::ONE, r), (-c, b.one())], vec![(F::ONE, b.one())], vec![(F::ONE, v)]);
    v
}

fn scalar_mul<F: PrimeField>(b: &mut Dr1csBuilder<F>, x: usize, y: usize) -> usize {
    let val = b.assignment[x] * b.assignment[y];
    let v = b.new_var(val);
    b.enforce_mul(x, y, v);
    v
}

/// Degree-3 Lagrange interpolation at `r` for points 0,1,2,3.
///
/// Returns scalar variables `(L0(r), L1(r), L2(r), L3(r))`.
fn lagrange_degree3<F: PrimeField>(b: &mut Dr1csBuilder<F>, r: usize) -> (usize, usize, usize, usize) {
    let inv2 = F::from(2u64).inverse().unwrap();
    let inv6 = F::from(6u64).inverse().unwrap();

    // t1 = r-1, t2 = r-2, t3 = r-3
    let t1 = scalar_sub_const(b, r, F::ONE);
    let t2 = scalar_sub_const(b, r, F::from(2u64));
    let t3 = scalar_sub_const(b, r, F::from(3u64));

    // L0 = -(t1*t2*t3)/6
    let p12 = scalar_mul(b, t1, t2);
    let p123 = scalar_mul(b, p12, t3);
    let l0_val = -(b.assignment[p123] * inv6);
    let l0 = b.new_var(l0_val);
    // l0 = (-inv6)*p123
    b.add_constraint(vec![(-inv6, p123)], vec![(F::ONE, b.one())], vec![(F::ONE, l0)]);

    // L1 = r*(r-2)*(r-3)/2
    let p = scalar_mul(b, r, t2);
    let p = scalar_mul(b, p, t3);
    let l1_val = b.assignment[p] * inv2;
    let l1 = b.new_var(l1_val);
    b.add_constraint(vec![(inv2, p)], vec![(F::ONE, b.one())], vec![(F::ONE, l1)]);

    // L2 = -r*(r-1)*(r-3)/2
    let p = scalar_mul(b, r, t1);
    let p = scalar_mul(b, p, t3);
    let l2_val = -(b.assignment[p] * inv2);
    let l2 = b.new_var(l2_val);
    b.add_constraint(vec![(-inv2, p)], vec![(F::ONE, b.one())], vec![(F::ONE, l2)]);

    // L3 = r*(r-1)*(r-2)/6
    let p = scalar_mul(b, r, t1);
    let p = scalar_mul(b, p, t2);
    let l3_val = b.assignment[p] * inv6;
    let l3 = b.new_var(l3_val);
    b.add_constraint(vec![(inv6, p)], vec![(F::ONE, b.one())], vec![(F::ONE, l3)]);

    (l0, l1, l2, l3)
}

/// Verify one degree-3 sumcheck over "ring elements" represented coefficient-wise.
///
/// Inputs:
/// - `claimed_sum`: current claim (ring vars)
/// - `msgs[i][t]`: per-round, per-evaluation point ring vars (t in 0..4)
/// - `rs[i]`: per-round verifier challenge scalar vars
///
/// Returns the final subclaim value (ring vars) after all rounds.
pub fn sumcheck_verify_degree3<F: PrimeField>(
    b: &mut Dr1csBuilder<F>,
    mut claimed_sum: RingVars,
    msgs: &[ [RingVars; 4] ],
    rs: &[usize],
) -> Result<RingVars, String> {
    if msgs.len() != rs.len() {
        return Err("sumcheck_verify_degree3: msgs/rs length mismatch".to_string());
    }
    for (round, (m, &r)) in msgs.iter().zip(rs.iter()).enumerate() {
        // Check g(0)+g(1) == claimed_sum (coefficient-wise).
        let g01 = ring_add(b, &m[0], &m[1]);
        // enforce g01 == claimed_sum
        for i in 0..claimed_sum.d() {
            b.enforce_lc_times_one_eq_const(vec![
                (F::ONE, g01.coeffs[i]),
                (-F::ONE, claimed_sum.coeffs[i]),
            ]);
        }

        // Update claim = g(r) by Lagrange interpolation.
        let (l0, l1, l2, l3) = lagrange_degree3::<F>(b, r);
        let t0 = ring_scale(b, &m[0], l0);
        let t1 = ring_scale(b, &m[1], l1);
        let t2 = ring_scale(b, &m[2], l2);
        let t3 = ring_scale(b, &m[3], l3);
        let s01 = ring_add(b, &t0, &t1);
        let s23 = ring_add(b, &t2, &t3);
        let new_claim = ring_add(b, &s01, &s23);
        claimed_sum = new_claim;

        let _ = round;
    }
    Ok(claimed_sum)
}

