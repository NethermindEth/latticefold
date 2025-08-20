//! LF+ double commitments

#![allow(missing_docs)]
#![allow(non_snake_case)]

use ark_ff::PrimeField;
use ark_std::time::Duration;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use latticefold::arith::r1cs::R1CS;
use latticefold_plus::{
    lin::LinParameters,
    r1cs::{r1cs_decomposed_square, ComR1CS},
    rgchk::{DecompParameters, RgInstance},
    utils::estimate_bound,
};
use rand::prelude::*;
use stark_rings::{cyclotomic_ring::models::frog_ring::RqPoly as R, PolyRing, Ring};
use stark_rings_linalg::{Matrix, SparseMatrix};

fn criterion_benchmark(c: &mut Criterion) {
    let n = 1 << 16;
    let sop = R::dimension() * 128; // S inf-norm = 128
    let L = 3;
    let k = 4;
    let d = R::dimension();
    let b = (R::dimension() / 2) as u128;
    let B = estimate_bound(sop, L, d, k) / 2; // + 1;
    let m = n / k;
    let kappa = 2;
    // log_d' (q)
    let l = ((<<R as PolyRing>::BaseRing>::MODULUS.0[0] as f64).ln()
        / ((R::dimension() / 2) as f64).ln())
    .ceil() as usize;
    let dparams = DecompParameters { b, k, l };
    let params = LinParameters {
        kappa,
        decomp: dparams.clone(),
    };

    let mut rng = ark_std::test_rng();
    let pop = [R::ZERO, R::ONE];
    let z: Vec<R> = (0..m).map(|_| *pop.choose(&mut rng).unwrap()).collect();

    let r1cs = r1cs_decomposed_square(
        R1CS::<R> {
            l: 1,
            A: SparseMatrix::identity(m),
            B: SparseMatrix::identity(m),
            C: SparseMatrix::identity(m),
        },
        n,
        B,
        k,
    );

    let A = Matrix::<R>::rand(&mut ark_std::test_rng(), params.kappa, n);

    let cr1cs = ComR1CS::new(r1cs, z, 1, B, k, &A);

    c.bench_function("RgInstance::from_f", |b| {
        b.iter_batched(
            || cr1cs.f.clone(),
            |f| {
                black_box(RgInstance::from_f(f, &A, &dparams));
            },
            criterion::BatchSize::SmallInput,
        )
    });
}

criterion_group!(
    name=benches;
    config = Criterion::default().sample_size(10).measurement_time(Duration::from_secs(50)).warm_up_time(Duration::from_secs(1));
    targets = criterion_benchmark);
criterion_main!(benches);
