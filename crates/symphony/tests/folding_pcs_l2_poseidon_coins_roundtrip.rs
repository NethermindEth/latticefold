use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;
use ark_crypto_primitives::sponge::CryptographicSponge;
use ark_ff::{Field, UniformRand};
use ark_std::test_rng;
use cyclotomic_rings::rings::{FrogPoseidonConfig as PC, GetPoseidonParams};
use stark_rings::{cyclotomic_ring::models::frog_ring::RqPoly as R, PolyRing};

use symphony::dpp_poseidon::{
    merge_sparse_dr1cs_share_one_with_glue, poseidon_sponge_dr1cs_from_trace_with_wiring_and_bytes,
};
use symphony::dpp_sumcheck::Dr1csBuilder;
use symphony::pcs::dpp_folding_pcs_l2::folding_pcs_l2_verify_dr1cs_with_c_bytes;
use symphony::pcs::dpp_folding_pcs_l2::folding_pcs_l2_params;
use symphony::pcs::folding_pcs_l2::{
    gadget_apply_digits, kron_ct_in_mul, kron_i_a_mul, kron_ikn_xt_mul, BinMatrix, DenseMatrix,
    FoldingPcsL2ProofCore,
    verify_folding_pcs_l2_with_c_matrices,
};
use symphony::transcript::PoseidonTraceOp;

type BF = <<R as PolyRing>::BaseRing as Field>::BasePrimeField;

fn bits_le_from_bytes(bytes: &[u8]) -> Vec<bool> {
    let mut out = Vec::with_capacity(bytes.len() * 8);
    for &b in bytes {
        for i in 0..8 {
            out.push(((b >> i) & 1) == 1);
        }
    }
    out
}

fn bin_matrix_from_bits(rows: usize, cols: usize, bits: &[bool]) -> BinMatrix<BF> {
    let mut data = Vec::with_capacity(rows * cols);
    for i in 0..rows * cols {
        data.push(if bits[i] { BF::ONE } else { BF::ZERO });
    }
    BinMatrix { rows, cols, data }
}

fn rand_small_signed_vec(len: usize, beta: u64, rng: &mut impl rand::RngCore) -> Vec<BF> {
    (0..len)
        .map(|_| {
            let mag = u64::rand(rng) % (beta + 1);
            let sign: bool = (u32::rand(rng) & 1) == 1;
            let v = BF::from(mag);
            if sign { -v } else { v }
        })
        .collect()
}

fn enforce_eq_const(b: &mut Dr1csBuilder<BF>, var: usize, c: BF) {
    let cv = b.new_var(c);
    b.enforce_var_eq_const(cv, c);
    b.enforce_lc_times_one_eq_const(vec![(BF::ONE, var), (-BF::ONE, cv)]);
}

#[test]
fn folding_pcs_l2_poseidon_derived_c_matrices_bind() {
    let mut rng = test_rng();

    // 1) Build a tiny Poseidon trace that squeezes enough bytes for C1/C2.
    let cfg = <PC as GetPoseidonParams<BF>>::get_poseidon_config();
    let mut sponge = PoseidonSponge::<BF>::new(&cfg);
    let absorb = (0..8).map(|_| BF::rand(&mut rng)).collect::<Vec<_>>();
    sponge.absorb(&absorb);
    let c_bytes = sponge.squeeze_bytes(2); // 16 bits

    let ops: Vec<PoseidonTraceOp<BF>> = vec![
        PoseidonTraceOp::Absorb(absorb.clone()),
        PoseidonTraceOp::SqueezeBytes { n: 2, out: c_bytes.clone() },
    ];

    let (pose_inst, pose_asg, _replay, _byte_wit, _wiring, byte_wiring) =
        poseidon_sponge_dr1cs_from_trace_with_wiring_and_bytes::<BF>(&cfg, &ops)
            .expect("poseidon dr1cs build failed");
    assert_eq!(byte_wiring.squeeze_byte_ranges.len(), 1);
    let (byte_start, byte_len) = byte_wiring.squeeze_byte_ranges[0];
    assert_eq!(byte_len, 2);
    let pose_byte_vars = &byte_wiring.squeeze_byte_vars[byte_start..byte_start + byte_len];

    // 2) Derive C1/C2 bits (same convention as dR1CS: little-endian bits per byte).
    let bits = bits_le_from_bytes(&c_bytes);
    // r=1,kappa=2 => rows=2, cols=2 => 4 bits per C
    let c1 = bin_matrix_from_bits(2, 2, &bits[0..4]);
    let c2 = bin_matrix_from_bits(2, 2, &bits[4..8]);

    // 3) Build a tiny FoldingPCS instance whose correctness depends on C1/C2.
    let r = 1usize;
    let kappa = 2usize;
    let n = 4usize;
    let delta = 4u64;
    let alpha = 1usize;
    let beta0 = 1u64 << 10;
    let beta1 = 2 * beta0;
    let beta2 = 2 * beta1;

    // A = I_n.
    let mut a_data = vec![BF::ZERO; n * (r * n * alpha)];
    for i in 0..n {
        a_data[i * (r * n * alpha) + i] = BF::ONE;
    }
    let a = DenseMatrix::new(n, r * n * alpha, a_data);
    let p = folding_pcs_l2_params(r, kappa, n, delta, alpha, beta0, beta1, beta2, a);

    let x0 = vec![BF::ONE; r];
    let x1 = vec![BF::ONE; r];
    let x2 = vec![BF::ONE; r];

    let y0 = rand_small_signed_vec(p.y0_len(), beta0, &mut rng);
    let y1 = kron_ct_in_mul(&c1, n, &y0);
    let y2 = kron_ct_in_mul(&c2, n, &y1);

    let t = kron_i_a_mul(&p.a, p.kappa, p.r * p.n * p.alpha, &y0);
    let mut delta_pows = Vec::with_capacity(alpha);
    let mut acc = BF::ONE;
    let delta_f = BF::from(delta);
    for _ in 0..alpha {
        delta_pows.push(acc);
        acc *= delta_f;
    }
    let v0 = gadget_apply_digits(&delta_pows, r * kappa * n, &y0);
    let v1 = gadget_apply_digits(&delta_pows, r * kappa * n, &y1);
    let v2 = gadget_apply_digits(&delta_pows, r * kappa * n, &y2);
    let u = kron_ikn_xt_mul(&x2, kappa, n, &v0);

    let core = FoldingPcsL2ProofCore { y0, v0, y1, v1, y2, v2 };
    // Native sanity check with explicit C matrices.
    verify_folding_pcs_l2_with_c_matrices(&p, &t, &x0, &x1, &x2, &u, &core, &c1, &c2)
        .expect("native verify failed");

    // 4) Build PCS verifier dR1CS that derives C1/C2 from byte-vars.
    let mut b_pcs = Dr1csBuilder::<BF>::new();
    let pcs_byte_vars: Vec<usize> = c_bytes
        .iter()
        .map(|&by| b_pcs.new_var(BF::from(by as u64)))
        .collect();

    let pcs_wiring = folding_pcs_l2_verify_dr1cs_with_c_bytes(
        &mut b_pcs,
        &p,
        &t,
        &x0,
        &x1,
        &x2,
        &core,
        &pcs_byte_vars,
    )
    .expect("pcs dr1cs build failed");
    // Bind u_re to the expected u constants.
    assert_eq!(pcs_wiring.u_re_vars.len(), u.len());
    for (v, c) in pcs_wiring.u_re_vars.iter().copied().zip(u.iter().copied()) {
        enforce_eq_const(&mut b_pcs, v, c);
    }
    let (pcs_inst, pcs_asg) = b_pcs.into_instance();

    // 5) Merge Poseidon + PCS instances and glue the byte vars together.
    let glue: Vec<(usize, usize, usize, usize)> = pose_byte_vars
        .iter()
        .zip(pcs_byte_vars.iter())
        .map(|(&pv, &cv)| (0usize, pv, 1usize, cv))
        .collect();

    let (merged_inst, merged_asg) = merge_sparse_dr1cs_share_one_with_glue::<BF>(
        &[(pose_inst, pose_asg), (pcs_inst, pcs_asg)],
        &glue,
    )
    .expect("merge failed");
    merged_inst.check(&merged_asg).expect("merged instance should satisfy");

    // Tamper: flip one Poseidon-derived byte var; merged should fail (both Poseidon and glue binding).
    let mut bad = merged_asg.clone();
    bad[pose_byte_vars[0]] += BF::ONE;
    assert!(merged_inst.check(&bad).is_err());
}

