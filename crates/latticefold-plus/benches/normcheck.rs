//! Microbenchmark for the norm-check path: ell-infinity (monomial range check)
//! vs ell-2 (JL projection + exact shortening) of ePrint 2026/721.
//!
//! Scoped to the norm-check prover only, at matched parameters (same witness
//! length `m`, same `L`, same `kappa`). The expected direction is that the ell-2
//! dominant term `L * m * 256` beats the ell-infinity `L * m * N * log_N(beta)`
//! (here `n = m`, single block). This produces the comparison number cited in
//! the PR; it is not the full end-to-end benchmark.

#![allow(non_snake_case)]

use ark_ff::PrimeField;
use ark_std::Zero;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
use latticefold_plus::{
    conjugation::l2_norm_sq,
    normcheck::{L2NormCheck, L2Params, LinfNormCheck, LinfParams, NormCheck},
    rgchk::DecompParameters,
    transcript::PoseidonTranscript,
};
use rand::Rng;
use stark_rings::{cyclotomic_ring::models::frog_ring::RqPoly as R, PolyRing};
use stark_rings_linalg::Matrix;

fn small_witness(rng: &mut impl Rng, l: usize, m: usize) -> Vec<Vec<R>> {
    (0..l)
        .map(|_| {
            (0..m)
                .map(|_| {
                    let mut c = R::zero();
                    for k in 0..R::dimension() {
                        let v: i64 = rng.gen_range(-1..=1);
                        c.coeffs_mut()[k] = if v < 0 {
                            -<R as PolyRing>::BaseRing::from(1u128)
                        } else {
                            <R as PolyRing>::BaseRing::from(v as u128)
                        };
                    }
                    c
                })
                .collect()
        })
        .collect()
}

fn linf_params() -> LinfParams {
    let b = (R::dimension() / 2) as u128;
    let l = ((<<R as PolyRing>::BaseRing>::MODULUS.0[0] as f64).ln()
        / ((R::dimension() / 2) as f64).ln())
    .ceil() as usize;
    LinfParams {
        kappa: 2,
        decomp: DecompParameters { b, k: 2, l },
    }
}

fn bench_normcheck(c: &mut Criterion) {
    let mut group = c.benchmark_group("normcheck_prover");
    group.sample_size(10);
    let L = 2usize;

    // The ell-infinity path requires a large witness (commits the gadget-
    // decomposed witness), so the matched-parameter comparison starts at 2^15.
    // Both are multiples of b = 256 * L for the ell-2 block structure.
    for &m in &[1usize << 15, 1usize << 16] {
        let mut rng = ark_std::test_rng();
        let w = small_witness(&mut rng, L, m);
        let A = Matrix::<R>::rand(&mut ark_std::test_rng(), 2, m);
        let beta_sq = w
            .iter()
            .map(|w_j| l2_norm_sq(w_j).into_bigint().0[0] as u128)
            .max()
            .unwrap();

        let l2_params = L2Params { m, beta_sq };
        let lp = linf_params();

        group.bench_with_input(BenchmarkId::new("l2", m), &m, |bencher, _| {
            bencher.iter(|| {
                let mut ts = PoseidonTranscript::empty::<PC>();
                L2NormCheck.prove(&w, &A, &l2_params, &mut ts)
            });
        });

        group.bench_with_input(BenchmarkId::new("linf", m), &m, |bencher, _| {
            bencher.iter(|| {
                let mut ts = PoseidonTranscript::empty::<PC>();
                LinfNormCheck.prove(&w, &A, &lp, &mut ts)
            });
        });
    }

    group.finish();
}

criterion_group!(benches, bench_normcheck);
criterion_main!(benches);
