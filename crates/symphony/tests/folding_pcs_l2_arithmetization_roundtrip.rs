use ark_ff::{Field, PrimeField, UniformRand};
use ark_std::test_rng;
use stark_rings::{cyclotomic_ring::models::frog_ring::RqPoly as R, PolyRing};

use symphony::dpp_sumcheck::Dr1csBuilder;
use symphony::pcs::dpp_folding_pcs_l2::{folding_pcs_l2_params, folding_pcs_l2_verify_dr1cs_with_c_bits};
use symphony::pcs::folding_pcs_l2::{
    gadget_apply_digits, kron_ikn_xt_mul, verify_folding_pcs_l2_with_c_matrices, BinMatrix,
    DenseMatrix, FoldingPcsL2ProofCore,
};

// Match the rest of Symphony's DPP tests: arithmetize over Frog's base prime field.
type F = <<R as PolyRing>::BaseRing as ark_ff::Field>::BasePrimeField;

fn rand_small_signed_vec<Ff: PrimeField>(len: usize, beta: u64, rng: &mut impl rand::RngCore) -> Vec<Ff> {
    (0..len)
        .map(|_| {
            let mag = u64::rand(rng) % (beta + 1);
            let sign: bool = (u32::rand(rng) & 1) == 1;
            let v = Ff::from(mag);
            if sign { -v } else { v }
        })
        .collect()
}

fn rand_bin_matrix<Ff: PrimeField>(rows: usize, cols: usize, rng: &mut impl rand::RngCore) -> BinMatrix<Ff> {
    let mut data = Vec::with_capacity(rows * cols);
    for _ in 0..rows * cols {
        let b: bool = (u32::rand(rng) & 1) == 1;
        data.push(if b { Ff::ONE } else { Ff::ZERO });
    }
    BinMatrix { rows, cols, data }
}

fn alloc_c_bits<Ff: PrimeField>(b: &mut Dr1csBuilder<Ff>, c: &BinMatrix<Ff>) -> Vec<usize> {
    let mut out = Vec::with_capacity(c.rows * c.cols);
    for &x in &c.data {
        let v = b.new_var(x);
        b.enforce_var_eq_const(v, x);
        out.push(v);
    }
    out
}

fn enforce_eq_const<Ff: PrimeField>(b: &mut Dr1csBuilder<Ff>, var: usize, c: Ff) {
    let cv = b.new_var(c);
    b.enforce_var_eq_const(cv, c);
    b.enforce_lc_times_one_eq_const(vec![(Ff::ONE, var), (-Ff::ONE, cv)]);
}

#[test]
fn folding_pcs_l2_dr1cs_roundtrip_small() {
    let mut rng = test_rng();

    // Tiny parameters for fast test.
    // Pick a tiny shape that avoids solving linear systems:
    // r=1, alpha=1, and A=I_n, so (I_k ⊗ A) is the identity map on each block.
    let r = 1usize;
    let kappa = 2usize;
    let n = 4usize;
    let delta = 4u64;
    let alpha = 1usize;

    // Choose a small per-coordinate bound for y0, and set y1/y2 bounds large enough for 0/1 mixing.
    let beta0 = 1u64 << 10;
    let beta1 = 2 * beta0; // sum of <= (r*kappa)=2 terms
    let beta2 = 2 * beta1; // another sum of <=2 terms

    // A = I_n.
    let mut a_data = vec![F::ZERO; n * (r * n * alpha)];
    for i in 0..n {
        a_data[i * (r * n * alpha) + i] = F::ONE;
    }
    let a = DenseMatrix::new(n, r * n * alpha, a_data);
    let p = folding_pcs_l2_params(r, kappa, n, delta, alpha, beta0, beta1, beta2, a);

    // Public inputs.
    let x0 = vec![F::ONE; r];
    let x1 = vec![F::ONE; r];
    let x2 = vec![F::ONE; r];

    // Short y0 sampled from [-beta0, beta0].
    let y0 = rand_small_signed_vec::<F>(p.y0_len(), beta0, &mut rng);

    // Random binary challenges (as witness bits for now).
    let c1 = rand_bin_matrix::<F>(r * kappa, kappa, &mut rng);
    let c2 = rand_bin_matrix::<F>(r * kappa, kappa, &mut rng);

    // With A=I and alpha=1 (G=identity), the y-chain is:
    // y1 = (C1^T ⊗ I_n) * y0
    // y2 = (C2^T ⊗ I_n) * y1
    let y1 = symphony::pcs::folding_pcs_l2::kron_ct_in_mul(&c1, n, &y0);
    let y2 = symphony::pcs::folding_pcs_l2::kron_ct_in_mul(&c2, n, &y1);

    // Compute t from (I_k ⊗ A) y0.
    let t = symphony::pcs::folding_pcs_l2::kron_i_a_mul(&p.a, p.kappa, p.r * p.n * p.alpha, &y0);

    let mut delta_pows = Vec::with_capacity(alpha);
    let mut acc = F::ONE;
    let delta_f = F::from(delta);
    for _ in 0..alpha {
        delta_pows.push(acc);
        acc *= delta_f;
    }
    let v0 = gadget_apply_digits(&delta_pows, r * kappa * n, &y0);
    let v1 = gadget_apply_digits(&delta_pows, r * kappa * n, &y1);
    let v2 = gadget_apply_digits(&delta_pows, r * kappa * n, &y2);
    let u = kron_ikn_xt_mul(&x2, kappa, n, &v0);

    let core = FoldingPcsL2ProofCore { y0, v0, y1, v1, y2, v2 };

    // Plain verifier (sanity).
    verify_folding_pcs_l2_with_c_matrices(&p, &t, &x0, &x1, &x2, &u, &core, &c1, &c2).unwrap();

    // Arithmetize verifier into dR1CS and check satisfiable.
    let mut b = Dr1csBuilder::<F>::new();
    let c1_bits = alloc_c_bits(&mut b, &c1);
    let c2_bits = alloc_c_bits(&mut b, &c2);
    let wiring =
        folding_pcs_l2_verify_dr1cs_with_c_bits(&mut b, &p, &t, &x0, &x1, &x2, &core, &c1_bits, &c2_bits)
            .unwrap();
    // Bind u_re to expected u constants.
    assert_eq!(wiring.u_re_vars.len(), u.len());
    for (v, c) in wiring.u_re_vars.iter().copied().zip(u.iter().copied()) {
        enforce_eq_const(&mut b, v, c);
    }
    let (inst, asg) = b.into_instance();
    inst.check(&asg).unwrap();

    // Tamper: flip one bit in C1 and ensure constraints break.
    let mut b2 = Dr1csBuilder::<F>::new();
    let mut c1_bad = c1.clone();
    c1_bad.data[0] = if c1_bad.data[0] == F::ZERO { F::ONE } else { F::ZERO };
    let c1_bits_bad = alloc_c_bits(&mut b2, &c1_bad);
    let c2_bits_bad = alloc_c_bits(&mut b2, &c2);
    let wiring_bad =
        folding_pcs_l2_verify_dr1cs_with_c_bits(&mut b2, &p, &t, &x0, &x1, &x2, &core, &c1_bits_bad, &c2_bits_bad)
            .unwrap();
    for (v, c) in wiring_bad.u_re_vars.iter().copied().zip(u.iter().copied()) {
        enforce_eq_const(&mut b2, v, c);
    }
    let (inst2, asg2) = b2.into_instance();
    assert!(inst2.check(&asg2).is_err());
}

#[test]
fn folding_pcs_l2_dr1cs_roundtrip_r2_nontrivial() {
    let mut rng = test_rng();

    // Nontrivial r=2 test:
    // - alpha=1 so G is identity (no gadget digits), but r>1 exercises (I_{kappa*n} ⊗ x^T).
    // - A is chosen as [I_n | 0] so we can construct y1,y2 preimages without solving.
    let r = 2usize;
    let kappa = 2usize;
    let n = 4usize;
    let delta = 4u64;
    let alpha = 1usize;

    // Bound growth: rhs1 sums <= (r*kappa)=4 terms; rhs2 sums <=4 terms of y1.
    let beta0 = 1u64 << 10;
    let beta1 = 4 * beta0;
    let beta2 = 4 * beta1;

    // A = [I_n | 0] ∈ F^{n × (r*n)}.
    let cols = r * n * alpha; // = 2n
    let mut a_data = vec![F::ZERO; n * cols];
    for i in 0..n {
        a_data[i * cols + i] = F::ONE;
    }
    let a = DenseMatrix::new(n, cols, a_data);
    let p = folding_pcs_l2_params(r, kappa, n, delta, alpha, beta0, beta1, beta2, a);

    // Choose x0=x1=x2=[1,0], so (I_{kappa*n} ⊗ x^T) selects the "first lane" of each r-block.
    let x0 = vec![F::ONE, F::ZERO];
    let x1 = vec![F::ONE, F::ZERO];
    let x2 = vec![F::ONE, F::ZERO];

    // Short y0 sampled from [-beta0, beta0].
    let y0 = rand_small_signed_vec::<F>(p.y0_len(), beta0, &mut rng);

    // Random binary challenges (as witness bits for now).
    let c1 = rand_bin_matrix::<F>(r * kappa, kappa, &mut rng);
    let c2 = rand_bin_matrix::<F>(r * kappa, kappa, &mut rng);

    // Compute RHS1 = (C1^T ⊗ I_n) y0  (since G is identity when alpha=1).
    let rhs1 = symphony::pcs::folding_pcs_l2::kron_ct_in_mul(&c1, n, &y0); // length kappa*n

    // Build y1 (length kappa*r*n) such that (I_kappa ⊗ A) y1 = rhs1.
    // With A=[I|0], we can set first n coords to rhs1 block and second n coords to 0.
    let mut y1 = vec![F::ZERO; p.y1_len()];
    for blk in 0..kappa {
        let rhs_blk = &rhs1[blk * n..(blk + 1) * n];
        let dst = &mut y1[blk * (r * n)..blk * (r * n) + n];
        dst.copy_from_slice(rhs_blk);
    }

    // Compute RHS2 = (C2^T ⊗ I_n) y1.
    let rhs2 = symphony::pcs::folding_pcs_l2::kron_ct_in_mul(&c2, n, &y1); // length kappa*n

    // Build y2 similarly so (I_kappa ⊗ A) y2 = rhs2.
    let mut y2 = vec![F::ZERO; p.y2_len()];
    for blk in 0..kappa {
        let rhs_blk = &rhs2[blk * n..(blk + 1) * n];
        let dst = &mut y2[blk * (r * n)..blk * (r * n) + n];
        dst.copy_from_slice(rhs_blk);
    }

    // Compute t from (I_k ⊗ A) y0.
    let t = symphony::pcs::folding_pcs_l2::kron_i_a_mul(&p.a, p.kappa, p.r * p.n * p.alpha, &y0);

    let mut delta_pows = Vec::with_capacity(alpha);
    let mut acc = F::ONE;
    let delta_f = F::from(delta);
    for _ in 0..alpha {
        delta_pows.push(acc);
        acc *= delta_f;
    }
    let v0 = gadget_apply_digits(&delta_pows, r * kappa * n, &y0);
    let v1 = gadget_apply_digits(&delta_pows, r * kappa * n, &y1);
    let v2 = gadget_apply_digits(&delta_pows, r * kappa * n, &y2);

    // u = (I_{kappa*n} ⊗ x2^T) v0.
    let u = kron_ikn_xt_mul(&x2, kappa, n, &v0);

    let core = FoldingPcsL2ProofCore { y0, v0, y1, v1, y2, v2 };

    // Plain verifier (sanity).
    verify_folding_pcs_l2_with_c_matrices(&p, &t, &x0, &x1, &x2, &u, &core, &c1, &c2).unwrap();

    // Arithmetize verifier into dR1CS and check satisfiable.
    let mut b = Dr1csBuilder::<F>::new();
    let c1_bits = alloc_c_bits(&mut b, &c1);
    let c2_bits = alloc_c_bits(&mut b, &c2);
    let wiring =
        folding_pcs_l2_verify_dr1cs_with_c_bits(&mut b, &p, &t, &x0, &x1, &x2, &core, &c1_bits, &c2_bits)
            .unwrap();
    assert_eq!(wiring.u_re_vars.len(), u.len());
    for (v, c) in wiring.u_re_vars.iter().copied().zip(u.iter().copied()) {
        enforce_eq_const(&mut b, v, c);
    }
    let (inst, asg) = b.into_instance();
    inst.check(&asg).unwrap();

    // Tamper: flip one bit in C2 and ensure constraints break.
    let mut b2 = Dr1csBuilder::<F>::new();
    let mut c2_bad = c2.clone();
    c2_bad.data[0] = if c2_bad.data[0] == F::ZERO { F::ONE } else { F::ZERO };
    let c1_bits_bad = alloc_c_bits(&mut b2, &c1);
    let c2_bits_bad = alloc_c_bits(&mut b2, &c2_bad);
    let wiring_bad =
        folding_pcs_l2_verify_dr1cs_with_c_bits(&mut b2, &p, &t, &x0, &x1, &x2, &core, &c1_bits_bad, &c2_bits_bad)
            .unwrap();
    for (v, c) in wiring_bad.u_re_vars.iter().copied().zip(u.iter().copied()) {
        enforce_eq_const(&mut b2, v, c);
    }
    let (inst2, asg2) = b2.into_instance();
    assert!(inst2.check(&asg2).is_err());
}

#[test]
fn folding_pcs_l2_dr1cs_norm_violation_fails() {
    // Tiny shape (r=1) but with a very small beta so we can force a violation.
    let r = 1usize;
    let kappa = 1usize;
    let n = 2usize;
    let delta = 4u64;
    let alpha = 1usize;
    let beta0 = 1u64; // force failure with value 2
    let beta1 = 1u64 << 20;
    let beta2 = 1u64 << 20;

    // A = I_n.
    let mut a_data = vec![F::ZERO; n * (r * n * alpha)];
    for i in 0..n {
        a_data[i * (r * n * alpha) + i] = F::ONE;
    }
    let a = DenseMatrix::new(n, r * n * alpha, a_data);
    let p = folding_pcs_l2_params(r, kappa, n, delta, alpha, beta0, beta1, beta2, a);

    let x0 = vec![F::ONE; r];
    let x1 = vec![F::ONE; r];
    let x2 = vec![F::ONE; r];

    // Deterministic violation: y0 = [2,0]
    let mut y0 = vec![F::ZERO; p.y0_len()];
    y0[0] = F::from(2u64);

    let mut rng = test_rng();
    let c1 = rand_bin_matrix::<F>(r * kappa, kappa, &mut rng);
    let c2 = rand_bin_matrix::<F>(r * kappa, kappa, &mut rng);
    let y1 = symphony::pcs::folding_pcs_l2::kron_ct_in_mul(&c1, n, &y0);
    let y2 = symphony::pcs::folding_pcs_l2::kron_ct_in_mul(&c2, n, &y1);

    let t = symphony::pcs::folding_pcs_l2::kron_i_a_mul(&p.a, p.kappa, p.r * p.n * p.alpha, &y0);
    let mut delta_pows = Vec::with_capacity(alpha);
    let mut acc = F::ONE;
    let delta_f = F::from(delta);
    for _ in 0..alpha {
        delta_pows.push(acc);
        acc *= delta_f;
    }
    let v0 = gadget_apply_digits(&delta_pows, r * kappa * n, &y0);
    let v1 = gadget_apply_digits(&delta_pows, r * kappa * n, &y1);
    let v2 = gadget_apply_digits(&delta_pows, r * kappa * n, &y2);
    let u = kron_ikn_xt_mul(&x2, kappa, n, &v0);

    let core = FoldingPcsL2ProofCore { y0, v0, y1, v1, y2, v2 };

    // Native verifier should fail on norm bound.
    assert!(verify_folding_pcs_l2_with_c_matrices(&p, &t, &x0, &x1, &x2, &u, &core, &c1, &c2).is_err());

    // dR1CS should be unsatisfiable.
    let mut b = Dr1csBuilder::<F>::new();
    let c1_bits = alloc_c_bits(&mut b, &c1);
    let c2_bits = alloc_c_bits(&mut b, &c2);
    let wiring =
        folding_pcs_l2_verify_dr1cs_with_c_bits(&mut b, &p, &t, &x0, &x1, &x2, &core, &c1_bits, &c2_bits)
            .unwrap();
    for (v, c) in wiring.u_re_vars.iter().copied().zip(u.iter().copied()) {
        enforce_eq_const(&mut b, v, c);
    }
    let (inst, asg) = b.into_instance();
    assert!(inst.check(&asg).is_err());
}

