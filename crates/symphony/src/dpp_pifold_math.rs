//! Π_fold verifier arithmetic arithmetization (sparse dR1CS over the transcript prime field).
//!
//! This is the next layer above:
//! - Poseidon transcript arithmetization (`dpp_poseidon`)
//! - Sumcheck verifier arithmetization (`dpp_sumcheck`)
//!
//! Goal: encode the arithmetic checks performed by
//! `verify_pi_fold_batched_and_fold_outputs_with_openings_and_aux_hetero_m`:
//! - shared-randomness sumcheck verification (had + mon)
//! - Eq.(26) check (Hadamard linkage)
//! - monomial recomputation check
//! - Step-5 folded check
//!
//! ## Important limitation (current state)
//! We treat all prover messages / aux values as witness variables, but we do **not** yet
//! enforce the *Poseidon<->challenge* binding here; instead callers should glue the
//! challenge variables to Poseidon `SqueezeField` outputs via
//! `merge_sparse_dr1cs_share_one_with_glue` (see `dpp_poseidon`).
//!
//! Also: `derive_beta_chi` uses `SqueezeBytes`; the byte->beta constraints are not
//! arithmetized yet. For now, callers provide `beta_cts` as witness scalars.

use ark_ff::{Field, PrimeField};
use stark_rings::{PolyRing, Zq};

use crate::dpp_poseidon::SparseDr1csInstance;
use crate::dpp_sumcheck::{sumcheck_verify_degree3, Dr1csBuilder, RingVars};
use crate::symphony_pifold_batched::{PiFoldAuxWitness, PiFoldBatchedProof};

#[derive(Clone, Debug)]
pub struct PiFoldMathWiring {
    pub beta_cts: Vec<usize>,
    pub s_base: Vec<usize>,
    pub alpha_base: usize,
    /// Flattened per-instance-per-digit vectors.
    pub c_all: Vec<Vec<usize>>,
    pub beta_i_all: Vec<usize>,
    pub alpha_i_all: Vec<usize>,
    pub rc_all: Vec<Option<usize>>,
    pub rhos: Vec<usize>,
    pub rs_shared: Vec<usize>,
    /// Step-5 folded scalar LHS per digit: lhs_ct[dig] = ct(psi * u_folded[dig]).
    ///
    /// This is the exact scalar that the Π_fold verifier checks against `dot(v_digits_folded[dig], ts_s)`.
    /// Exposing it enables batching / binding checks in higher layers (e.g. batchlin PCS).
    pub step5_lhs_ct: Vec<usize>,
}

type BF<R> = <<R as PolyRing>::BaseRing as Field>::BasePrimeField;

/// Evaluate the multilinear eq polynomial at `(a,b)`:
/// eq(a,b) = Π_i (a_i*b_i + (1-a_i)*(1-b_i)).
fn eq_eval_vars<F: PrimeField>(b: &mut Dr1csBuilder<F>, a: &[usize], r: &[usize]) -> usize {
    assert_eq!(a.len(), r.len());
    let mut acc = b.new_var(F::ONE);
    b.enforce_var_eq_const(acc, F::ONE);
    for (&ai, &ri) in a.iter().zip(r.iter()) {
        // t = ai*ri + (1-ai)*(1-ri)
        let one = b.one();
        let one_minus_ai = b.new_var(F::ONE - b.assignment[ai]);
        b.add_constraint(vec![(F::ONE, one), (-F::ONE, ai)], vec![(F::ONE, one)], vec![(F::ONE, one_minus_ai)]);
        let one_minus_ri = b.new_var(F::ONE - b.assignment[ri]);
        b.add_constraint(vec![(F::ONE, one), (-F::ONE, ri)], vec![(F::ONE, one)], vec![(F::ONE, one_minus_ri)]);
        let ai_ri = b.new_var(b.assignment[ai] * b.assignment[ri]);
        b.enforce_mul(ai, ri, ai_ri);
        let om = b.new_var(b.assignment[one_minus_ai] * b.assignment[one_minus_ri]);
        b.enforce_mul(one_minus_ai, one_minus_ri, om);
        let t = b.new_var(b.assignment[ai_ri] + b.assignment[om]);
        b.add_constraint(vec![(F::ONE, ai_ri), (F::ONE, om)], vec![(F::ONE, one)], vec![(F::ONE, t)]);

        // acc *= t
        let new_acc = b.new_var(b.assignment[acc] * b.assignment[t]);
        b.enforce_mul(acc, t, new_acc);
        acc = new_acc;
    }
    acc
}

fn base_to_bf<R: PolyRing>(x: R::BaseRing) -> <R::BaseRing as Field>::BasePrimeField
where
    R::BaseRing: Field,
{
    // We require extension degree 1 throughout Symphony's current base field.
    x.to_base_prime_field_elements()
        .into_iter()
        .next()
        .expect("base_to_bf expects extension degree 1")
}

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
        // treat as witness var: no eq_const
        coeffs.push(v);
    }
    RingVars::new(coeffs)
}

/// Build a dR1CS instance that enforces the Π_fold verifier arithmetic checks, given:
/// - the proof object (`PiFoldBatchedProof`)
/// - auxiliary messages (`aux`) supplying `had_u` and `mon_b`
/// - explicit challenge scalars (in transcript order, as BF values)
///
/// Returns `(inst, assignment, challenge_vars)` where `challenge_vars[i]` is the variable index
/// corresponding to `challenges[i]` in the assignment.
pub fn pifold_verifier_math_dr1cs<R>(
    proof: &PiFoldBatchedProof<R>,
    aux: &PiFoldAuxWitness<R>,
    beta_cts: &[BF<R>],
    s_base: &[BF<R>],
    alpha_base: BF<R>,
    // Per-instance per-digit: (c vector, beta, alpha)
    cba_all: &[Vec<(Vec<BF<R>>, BF<R>, BF<R>)>],
    rc_all: &[Option<BF<R>>],
    rhos: &[BF<R>],
    rs_shared: &[BF<R>],
) -> Result<
    (
        SparseDr1csInstance<BF<R>>,
        Vec<BF<R>>,
        PiFoldMathWiring,
    ),
    String,
>
where
    R: PolyRing,
    R: stark_rings::OverField,
    R::BaseRing: Zq + Field + PrimeField,
{
    let d = R::dimension();
    let ell = aux.had_u.len();
    if ell == 0 {
        return Err("pifold math dr1cs: empty batch".to_string());
    }
    if aux.mon_b.len() != ell {
        return Err("pifold math dr1cs: aux length mismatch".to_string());
    }
    if beta_cts.len() != ell {
        return Err("pifold math dr1cs: beta_cts length mismatch".to_string());
    }

    let mut b = Dr1csBuilder::<BF<R>>::new();

    let log_m = (proof.m as f64).log2() as usize;
    let g_len = proof.m * d;
    let g_nvars = (g_len as f64).log2() as usize;
    let k_g = proof.rg_params.k_g;

    if s_base.len() != log_m {
        return Err("pifold math dr1cs: s_base length mismatch".to_string());
    }
    if cba_all.len() != ell || rc_all.len() != ell || rhos.len() != ell {
        return Err("pifold math dr1cs: per-instance challenge shapes mismatch".to_string());
    }
    if rs_shared.len() != g_nvars {
        return Err("pifold math dr1cs: rs_shared length mismatch".to_string());
    }

    // Allocate witness vars for coin values. Sequential allocation required.
    let beta_vars: Vec<_> = beta_cts.iter().copied().map(|x| b.new_var(x)).collect();
    let s_vars: Vec<_> = s_base.iter().copied().map(|x| b.new_var(x)).collect();
    let alpha_base_var = b.new_var(alpha_base);

    let mut c_all: Vec<Vec<usize>> = Vec::with_capacity(ell * k_g);
    let mut beta_i_all: Vec<usize> = Vec::with_capacity(ell * k_g);
    let mut alpha_i_all: Vec<usize> = Vec::with_capacity(ell * k_g);
    let mut rc_vars: Vec<Option<usize>> = Vec::with_capacity(ell);
    for inst_idx in 0..ell {
        if cba_all[inst_idx].len() != k_g {
            return Err("pifold math dr1cs: cba_all[inst] k_g mismatch".to_string());
        }
        for dig in 0..k_g {
            let (c_vec, beta_i, alpha_i) = &cba_all[inst_idx][dig];
            if c_vec.len() != g_nvars {
                return Err("pifold math dr1cs: c vector length mismatch".to_string());
            }
            let c_vars = c_vec.iter().copied().map(|x| b.new_var(x)).collect::<Vec<_>>();
            c_all.push(c_vars);
            beta_i_all.push(b.new_var(*beta_i));
            alpha_i_all.push(b.new_var(*alpha_i));
        }
        let rc = rc_all[inst_idx].map(|x| b.new_var(x));
        rc_vars.push(rc);
    }
    let rho_vars = rhos.iter().copied().map(|x| b.new_var(x)).collect::<Vec<_>>();
    let rs_vars = rs_shared.iter().copied().map(|x| b.new_var(x)).collect::<Vec<_>>();

    // Step5 uses s_chals = mon_point[log_m..] where mon_point = rs_shared (length g_nvars)
    let s_chals = &rs_vars[log_m..];

    let mut step5_lhs_ct_vars: Vec<usize> = Vec::with_capacity(k_g);

    // -------------------------------------------------------------------------
    // Sumcheck verification (degree 3): enforce transcript consistency g(0)+g(1)=claim and update.
    // -------------------------------------------------------------------------
    // Prepare prover message vars from proofs.
    // had proof has log_m rounds, mon proof has g_nvars rounds.
    if proof.had_sumcheck.msgs().len() != log_m || proof.mon_sumcheck.msgs().len() != g_nvars {
        return Err("pifold math dr1cs: proof round length mismatch".to_string());
    }

    let mut had_msgs: Vec<[RingVars; 4]> = Vec::with_capacity(log_m);
    for pm in proof.had_sumcheck.msgs() {
        if pm.evaluations.len() != 4 {
            return Err("pifold math dr1cs: had msg eval len != 4".to_string());
        }
        let m0 = ring_to_ringvars::<R>(&mut b, &pm.evaluations[0]);
        let m1 = ring_to_ringvars::<R>(&mut b, &pm.evaluations[1]);
        let m2 = ring_to_ringvars::<R>(&mut b, &pm.evaluations[2]);
        let m3 = ring_to_ringvars::<R>(&mut b, &pm.evaluations[3]);
        had_msgs.push([m0, m1, m2, m3]);
    }

    let mut mon_msgs: Vec<[RingVars; 4]> = Vec::with_capacity(g_nvars);
    for pm in proof.mon_sumcheck.msgs() {
        if pm.evaluations.len() != 4 {
            return Err("pifold math dr1cs: mon msg eval len != 4".to_string());
        }
        let m0 = ring_to_ringvars::<R>(&mut b, &pm.evaluations[0]);
        let m1 = ring_to_ringvars::<R>(&mut b, &pm.evaluations[1]);
        let m2 = ring_to_ringvars::<R>(&mut b, &pm.evaluations[2]);
        let m3 = ring_to_ringvars::<R>(&mut b, &pm.evaluations[3]);
        mon_msgs.push([m0, m1, m2, m3]);
    }

    // claimed sums are ZERO ring elements.
    let mut zero_coeffs = Vec::with_capacity(d);
    for _ in 0..d {
        let v = b.new_var(BF::<R>::ZERO);
        b.enforce_var_eq_const(v, BF::<R>::ZERO);
        zero_coeffs.push(v);
    }
    let claim0 = RingVars::new(zero_coeffs.clone());
    let had_eval = sumcheck_verify_degree3(&mut b, claim0.clone(), &had_msgs, &rs_vars[..log_m])?;
    let mon_eval = sumcheck_verify_degree3(&mut b, claim0.clone(), &mon_msgs, &rs_vars)?;

    // -------------------------------------------------------------------------
    // Eq(26) check (Hadamard linkage): enforce only the constant-term coefficient (SP1 regime).
    // -------------------------------------------------------------------------
    let eq_sr = eq_eval_vars(&mut b, &s_vars, &rs_vars[..log_m]);

    // alpha_pows[j] are scalars in base ring for our constant-term check.
    let mut alpha_pows: Vec<usize> = Vec::with_capacity(d);
    // alpha^0 = 1
    let a0 = b.new_var(BF::<R>::ONE);
    b.enforce_var_eq_const(a0, BF::<R>::ONE);
    alpha_pows.push(a0);
    for j in 1..d {
        let prev = alpha_pows[j - 1];
        let v = b.new_var(b.assignment[prev] * b.assignment[alpha_base_var]);
        b.enforce_mul(prev, alpha_base_var, v);
        alpha_pows.push(v);
    }

    // had_u are base scalars; interpret as scalars.
    let mut lhs = b.new_var(BF::<R>::ZERO);
    b.enforce_var_eq_const(lhs, BF::<R>::ZERO);
    for inst_idx in 0..ell {
        // acc_i = Σ_j alpha^j * (u1*u2 - u3) at lane j.
        let mut acc_i = b.new_var(BF::<R>::ZERO);
        b.enforce_var_eq_const(acc_i, BF::<R>::ZERO);
        for j in 0..d {
            let u1 = base_to_bf::<R>(aux.had_u[inst_idx][0][j].into());
            let u2 = base_to_bf::<R>(aux.had_u[inst_idx][1][j].into());
            let u3 = base_to_bf::<R>(aux.had_u[inst_idx][2][j].into());
            let vu1 = b.new_var(u1);
            let vu2 = b.new_var(u2);
            let vu3 = b.new_var(u3);
            // witness vars (no eq const)
            let prod = b.new_var(b.assignment[vu1] * b.assignment[vu2]);
            b.enforce_mul(vu1, vu2, prod);
            let diff = b.new_var(b.assignment[prod] - b.assignment[vu3]);
            b.add_constraint(vec![(BF::<R>::ONE, prod), (-BF::<R>::ONE, vu3)], vec![(BF::<R>::ONE, b.one())], vec![(BF::<R>::ONE, diff)]);
            let term = b.new_var(b.assignment[diff] * b.assignment[alpha_pows[j]]);
            b.enforce_mul(diff, alpha_pows[j], term);
            let new_acc = b.new_var(b.assignment[acc_i] + b.assignment[term]);
            b.add_constraint(vec![(BF::<R>::ONE, acc_i), (BF::<R>::ONE, term)], vec![(BF::<R>::ONE, b.one())], vec![(BF::<R>::ONE, new_acc)]);
            acc_i = new_acc;
        }
        // inst_term = rho_i * eq_sr * acc_i
        let t = b.new_var(b.assignment[eq_sr] * b.assignment[acc_i]);
        b.enforce_mul(eq_sr, acc_i, t);
        let t2 = b.new_var(b.assignment[t] * b.assignment[rho_vars[inst_idx]]);
        b.enforce_mul(t, rho_vars[inst_idx], t2);
        let new_lhs = b.new_var(b.assignment[lhs] + b.assignment[t2]);
        b.add_constraint(vec![(BF::<R>::ONE, lhs), (BF::<R>::ONE, t2)], vec![(BF::<R>::ONE, b.one())], vec![(BF::<R>::ONE, new_lhs)]);
        lhs = new_lhs;
    }
    // Enforce had_eval.ct == lhs (ct is coeff 0), and other coeffs == 0.
    b.enforce_lc_times_one_eq_const(vec![(BF::<R>::ONE, had_eval.coeffs[0]), (-BF::<R>::ONE, lhs)]);
    for j in 1..d {
        b.enforce_lc_times_one_eq_const(vec![(BF::<R>::ONE, had_eval.coeffs[j])]);
    }

    // -------------------------------------------------------------------------
    // Monomial recomputation check (constant-term only) and Step-5 (psi-ct linear form).
    // -------------------------------------------------------------------------
    // Build ver scalar as in Rust, but using BF vars and `ev` over witness ring elems (aux.mon_b).
    let mut ver = b.new_var(BF::<R>::ZERO);
    b.enforce_var_eq_const(ver, BF::<R>::ZERO);

    for inst_idx in 0..ell {
        let mut inst_acc = b.new_var(BF::<R>::ZERO);
        b.enforce_var_eq_const(inst_acc, BF::<R>::ZERO);
        let mut rc_pow = b.new_var(BF::<R>::ONE);
        b.enforce_var_eq_const(rc_pow, BF::<R>::ONE);
        for dig in 0..k_g {
            let c = &c_all[inst_idx * k_g + dig];
            let eq = eq_eval_vars(&mut b, c, &rs_vars);
            // b_i: represent ring element coeffs as BF vars.
            let b_i = ring_to_ringvars::<R>(&mut b, &aux.mon_b[inst_idx][dig]);

            // ev1 = Σ coeff[k]*beta^k
            let beta_i = beta_i_all[inst_idx * k_g + dig];
            // compute powers
            let mut pow = b.new_var(BF::<R>::ONE);
            b.enforce_var_eq_const(pow, BF::<R>::ONE);
            let mut ev1 = b.new_var(BF::<R>::ZERO);
            b.enforce_var_eq_const(ev1, BF::<R>::ZERO);
            for k in 0..d {
                let term = b.new_var(b.assignment[b_i.coeffs[k]] * b.assignment[pow]);
                b.enforce_mul(b_i.coeffs[k], pow, term);
                let new_ev1 = b.new_var(b.assignment[ev1] + b.assignment[term]);
                b.add_constraint(vec![(BF::<R>::ONE, ev1), (BF::<R>::ONE, term)], vec![(BF::<R>::ONE, b.one())], vec![(BF::<R>::ONE, new_ev1)]);
                ev1 = new_ev1;
                if k + 1 < d {
                    let new_pow = b.new_var(b.assignment[pow] * b.assignment[beta_i]);
                    b.enforce_mul(pow, beta_i, new_pow);
                    pow = new_pow;
                }
            }
            // ev2 at beta^2
            let beta2 = b.new_var(b.assignment[beta_i] * b.assignment[beta_i]);
            b.enforce_mul(beta_i, beta_i, beta2);
            let mut pow2 = b.new_var(BF::<R>::ONE);
            b.enforce_var_eq_const(pow2, BF::<R>::ONE);
            let mut ev2 = b.new_var(BF::<R>::ZERO);
            b.enforce_var_eq_const(ev2, BF::<R>::ZERO);
            for k in 0..d {
                let term = b.new_var(b.assignment[b_i.coeffs[k]] * b.assignment[pow2]);
                b.enforce_mul(b_i.coeffs[k], pow2, term);
                let new_ev2 = b.new_var(b.assignment[ev2] + b.assignment[term]);
                b.add_constraint(vec![(BF::<R>::ONE, ev2), (BF::<R>::ONE, term)], vec![(BF::<R>::ONE, b.one())], vec![(BF::<R>::ONE, new_ev2)]);
                ev2 = new_ev2;
                if k + 1 < d {
                    let new_pow = b.new_var(b.assignment[pow2] * b.assignment[beta2]);
                    b.enforce_mul(pow2, beta2, new_pow);
                    pow2 = new_pow;
                }
            }
            // b_claim = ev1^2 - ev2
            let ev1sq = b.new_var(b.assignment[ev1] * b.assignment[ev1]);
            b.enforce_mul(ev1, ev1, ev1sq);
            let b_claim = b.new_var(b.assignment[ev1sq] - b.assignment[ev2]);
            b.add_constraint(vec![(BF::<R>::ONE, ev1sq), (-BF::<R>::ONE, ev2)], vec![(BF::<R>::ONE, b.one())], vec![(BF::<R>::ONE, b_claim)]);

            let alpha_i = alpha_i_all[inst_idx * k_g + dig];
            let mut term = b.new_var(b.assignment[eq] * b.assignment[alpha_i]);
            b.enforce_mul(eq, alpha_i, term);
            let term0 = term;
            term = b.new_var(b.assignment[term0] * b.assignment[b_claim]);
            b.enforce_mul(term0, b_claim, term);
            if let Some(rc) = rc_vars[inst_idx] {
                let term1 = term;
                term = b.new_var(b.assignment[term1] * b.assignment[rc_pow]);
                b.enforce_mul(term1, rc_pow, term);
                let new_rc_pow = b.new_var(b.assignment[rc_pow] * b.assignment[rc]);
                b.enforce_mul(rc_pow, rc, new_rc_pow);
                rc_pow = new_rc_pow;
            }
            let new_inst = b.new_var(b.assignment[inst_acc] + b.assignment[term]);
            b.add_constraint(vec![(BF::<R>::ONE, inst_acc), (BF::<R>::ONE, term)], vec![(BF::<R>::ONE, b.one())], vec![(BF::<R>::ONE, new_inst)]);
            inst_acc = new_inst;
        }
        let inst_term = b.new_var(b.assignment[rho_vars[inst_idx]] * b.assignment[inst_acc]);
        b.enforce_mul(rho_vars[inst_idx], inst_acc, inst_term);
        let new_ver = b.new_var(b.assignment[ver] + b.assignment[inst_term]);
        b.add_constraint(vec![(BF::<R>::ONE, ver), (BF::<R>::ONE, inst_term)], vec![(BF::<R>::ONE, b.one())], vec![(BF::<R>::ONE, new_ver)]);
        ver = new_ver;
    }

    // Enforce mon_eval.ct == ver and other coeffs 0.
    b.enforce_lc_times_one_eq_const(vec![(BF::<R>::ONE, mon_eval.coeffs[0]), (-BF::<R>::ONE, ver)]);
    for j in 1..d {
        b.enforce_lc_times_one_eq_const(vec![(BF::<R>::ONE, mon_eval.coeffs[j])]);
    }

    // Step-5: enforce (psi*u_folded[dig]).ct == dot(v_digits_folded[dig], ts_s)
    //
    // Use the manual negacyclic constant-term formula pinned in `rp_rgchk` tests.
    let psi_poly = stark_rings::psi::<R>();
    let psi_coeffs = psi_poly.coeffs();

    // Compute ts_s weights for log_d = log2(d).
    let log_d = (d as f64).log2() as usize;
    if s_chals.len() != log_d {
        return Err("pifold math dr1cs: unexpected s_chals length".to_string());
    }
    // Compute multilinear weights of length d from s_chals.
    let mut weights: Vec<usize> = vec![b.new_var(BF::<R>::ONE)];
    b.enforce_var_eq_const(weights[0], BF::<R>::ONE);
    for &r in s_chals {
        let mut next = Vec::with_capacity(weights.len() * 2);
        for &w in &weights {
            // w*(1-r)
            let one_minus_r = b.new_var(BF::<R>::ONE - b.assignment[r]);
            b.add_constraint(vec![(BF::<R>::ONE, b.one()), (-BF::<R>::ONE, r)], vec![(BF::<R>::ONE, b.one())], vec![(BF::<R>::ONE, one_minus_r)]);
            let a = b.new_var(b.assignment[w] * b.assignment[one_minus_r]);
            b.enforce_mul(w, one_minus_r, a);
            // w*r
            let c = b.new_var(b.assignment[w] * b.assignment[r]);
            b.enforce_mul(w, r, c);
            next.push(a);
            next.push(c);
        }
        weights = next;
    }

    // Fold u*: u_folded[dig] = Σ_i beta_cts[i] * mon_b[i][dig]
    for dig in 0..k_g {
        // u_folded ring vars
        let mut u_coeffs = vec![b.new_var(BF::<R>::ZERO); d];
        for j in 0..d { b.enforce_var_eq_const(u_coeffs[j], BF::<R>::ZERO); }
        for inst_idx in 0..ell {
            let b_i = ring_to_ringvars::<R>(&mut b, &aux.mon_b[inst_idx][dig]);
            for j in 0..d {
                let scaled = b.new_var(b.assignment[b_i.coeffs[j]] * b.assignment[beta_vars[inst_idx]]);
                b.enforce_mul(b_i.coeffs[j], beta_vars[inst_idx], scaled);
                let new_u = b.new_var(b.assignment[u_coeffs[j]] + b.assignment[scaled]);
                b.add_constraint(vec![(BF::<R>::ONE, u_coeffs[j]), (BF::<R>::ONE, scaled)], vec![(BF::<R>::ONE, b.one())], vec![(BF::<R>::ONE, new_u)]);
                u_coeffs[j] = new_u;
            }
        }
        // lhs_ct = psi0*u0 - Σ psi_i * u_{d-i}
        let mut lhs_ct = b.new_var(BF::<R>::ZERO);
        b.enforce_var_eq_const(lhs_ct, BF::<R>::ZERO);
        // psi0*u0
    let psi0 = base_to_bf::<R>(psi_coeffs[0]);
        let psi0v = b.new_var(psi0);
        let t0 = b.new_var(b.assignment[psi0v] * b.assignment[u_coeffs[0]]);
        b.enforce_mul(psi0v, u_coeffs[0], t0);
        lhs_ct = t0;
        for i_c in 1..d {
            let psi_i = base_to_bf::<R>(psi_coeffs[i_c]);
            let psiiv = b.new_var(psi_i);
            let ui = u_coeffs[d - i_c];
            let prod = b.new_var(b.assignment[psiiv] * b.assignment[ui]);
            b.enforce_mul(psiiv, ui, prod);
            let new_lhs = b.new_var(b.assignment[lhs_ct] - b.assignment[prod]);
            b.add_constraint(vec![(BF::<R>::ONE, lhs_ct), (-BF::<R>::ONE, prod)], vec![(BF::<R>::ONE, b.one())], vec![(BF::<R>::ONE, new_lhs)]);
            lhs_ct = new_lhs;
        }

        // rhs = dot(v_digits_folded[dig], weights)
        if proof.v_digits_folded.len() != k_g || proof.v_digits_folded[dig].len() != d {
            return Err("pifold math dr1cs: v_digits_folded shape mismatch".to_string());
        }
        let mut rhs = b.new_var(BF::<R>::ZERO);
        b.enforce_var_eq_const(rhs, BF::<R>::ZERO);
        for j in 0..d {
            let vij = base_to_bf::<R>(proof.v_digits_folded[dig][j]);
            let vv = b.new_var(vij);
            let prod = b.new_var(b.assignment[vv] * b.assignment[weights[j]]);
            b.enforce_mul(vv, weights[j], prod);
            let new_rhs = b.new_var(b.assignment[rhs] + b.assignment[prod]);
            b.add_constraint(vec![(BF::<R>::ONE, rhs), (BF::<R>::ONE, prod)], vec![(BF::<R>::ONE, b.one())], vec![(BF::<R>::ONE, new_rhs)]);
            rhs = new_rhs;
        }

        // lhs_ct == rhs
        b.enforce_lc_times_one_eq_const(vec![(BF::<R>::ONE, lhs_ct), (-BF::<R>::ONE, rhs)]);
        step5_lhs_ct_vars.push(lhs_ct);
    }

    let wiring = PiFoldMathWiring {
        beta_cts: beta_vars.clone(),
        s_base: s_vars.clone(),
        alpha_base: alpha_base_var,
        c_all: c_all.clone(),
        beta_i_all: beta_i_all.clone(),
        alpha_i_all: alpha_i_all.clone(),
        rc_all: rc_vars.clone(),
        rhos: rho_vars.clone(),
        rs_shared: rs_vars.clone(),
        step5_lhs_ct: step5_lhs_ct_vars,
    };

    let (inst, assignment) = b.into_instance();
    Ok((inst, assignment, wiring))
}

