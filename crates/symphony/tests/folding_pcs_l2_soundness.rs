use ark_ff::{Field, UniformRand};
use ark_std::test_rng;
use stark_rings::{cyclotomic_ring::models::frog_ring::RqPoly as R, PolyRing};

use symphony::pcs::dpp_folding_pcs_l2::folding_pcs_l2_params;
use symphony::pcs::folding_pcs_l2::{
    commit, open, verify_folding_pcs_l2_with_c_matrices, BinMatrix, DenseMatrix,
};

// Match the rest of Symphony's DPP tests: run PCS over Frog's base prime field.
type F = <<R as PolyRing>::BaseRing as Field>::BasePrimeField;

fn rand_bin_matrix(rows: usize, cols: usize, rng: &mut impl rand::RngCore) -> BinMatrix<F> {
    let mut data = Vec::with_capacity(rows * cols);
    for _ in 0..rows * cols {
        let b: bool = (u32::rand(rng) & 1) == 1;
        data.push(if b { F::ONE } else { F::ZERO });
    }
    BinMatrix { rows, cols, data }
}

#[test]
fn folding_pcs_commit_rejects_truncating_gadget_params() {
    // delta^alpha must cover the modulus. Use a deliberately too-small gadget.
    let r = 1usize;
    let kappa = 1usize;
    let n = 1usize;
    let delta = 4u64;
    let alpha = 1usize; // delta^alpha = 4 << modulus
    let beta0 = 1u64 << 10;
    let beta1 = beta0;
    let beta2 = beta0;

    let a = DenseMatrix::new(n, r * n * alpha, vec![F::ONE]);
    let p = folding_pcs_l2_params(r, kappa, n, delta, alpha, beta0, beta1, beta2, a);

    let f = vec![F::ZERO; p.f_len()];
    assert!(commit::<F>(&p, &f).is_err());
}

#[test]
fn folding_pcs_native_verify_fails_on_tampered_messages() {
    let mut rng = test_rng();

    // Small but nontrivial params with exact gadget: (2^32)^2 = 2^64 >= Frog modulus.
    let r = 2usize;
    let kappa = 1usize;
    let n = 1usize;
    let delta = 1u64 << 32;
    let alpha = 2usize;
    let beta0 = 1u64 << 63;
    let beta1 = beta0;
    let beta2 = beta0;

    let cols = r * n * alpha;
    let a_data = (0..(n * cols)).map(|_| F::rand(&mut rng)).collect::<Vec<_>>();
    let a = DenseMatrix::new(n, cols, a_data);
    let p = folding_pcs_l2_params(r, kappa, n, delta, alpha, beta0, beta1, beta2, a);

    let f = (0..p.f_len()).map(|_| F::rand(&mut rng)).collect::<Vec<_>>();
    let (t, s) = commit::<F>(&p, &f).expect("commit");

    let x0 = (0..r).map(|_| F::rand(&mut rng)).collect::<Vec<_>>();
    let x1 = (0..r).map(|_| F::rand(&mut rng)).collect::<Vec<_>>();
    let x2 = (0..r).map(|_| F::rand(&mut rng)).collect::<Vec<_>>();
    let c1 = rand_bin_matrix(r * kappa, kappa, &mut rng);
    let c2 = rand_bin_matrix(r * kappa, kappa, &mut rng);

    let (u, mut core) = open::<F>(&p, &f, &s, &x0, &x1, &x2, &c1, &c2).expect("open");
    verify_folding_pcs_l2_with_c_matrices(&p, &t, &x0, &x1, &x2, &u, &core, &c1, &c2)
        .expect("native verify should pass");

    // Mutate v2: should fail via G(y2)=v2.
    core.v2[0] += F::ONE;
    assert!(verify_folding_pcs_l2_with_c_matrices(&p, &t, &x0, &x1, &x2, &u, &core, &c1, &c2).is_err());
}

