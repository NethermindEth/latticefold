//! LF+ E2E (All protocols) prove and verify

#![allow(missing_docs)]
#![allow(non_snake_case)]

use ark_ff::PrimeField;
use ark_std::time::Duration;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
use latticefold::arith::r1cs::R1CS;
use latticefold_plus::{
    lin::LinParameters,
    plus::{PlusParameters, PlusProver, PlusVerifier},
    r1cs::{r1cs_decomposed_square, ComR1CS},
    rgchk::DecompParameters,
    transcript::PoseidonTranscript,
    utils::estimate_bound,
};
use rand::prelude::*;
use stark_rings::{cyclotomic_ring::models::frog_ring::RqPoly as R, PolyRing, Ring};
use stark_rings_linalg::{Matrix, SparseMatrix};

fn criterion_benchmark(c: &mut Criterion) {
    let n = 1 << 17;
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
    let params = LinParameters {
        kappa,
        decomp: DecompParameters { b, k, l },
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

    let M = cr1cs.x.matrices();
    let pparams = PlusParameters { lin: params, B };

    // Prover / Fold
    c.bench_function("prove", |b| {
        b.iter_batched(
            || {
                let ts = PoseidonTranscript::empty::<PC>();
                // L=3 (equal) instances are folded here
                // TODO Do accumulated instance (2) + one online (1)
                (
                    PlusProver::init(A.clone(), M.clone(), 1, pparams.clone(), ts),
                    [cr1cs.clone(), cr1cs.clone(), cr1cs.clone()],
                )
            },
            |(mut prover, c1rcs)| {
                black_box(prover.prove(&c1rcs));
            },
            criterion::BatchSize::SmallInput,
        )
    });

    // Verifier
    let ts = PoseidonTranscript::empty::<PC>();
    let mut prover = PlusProver::init(A.clone(), M.clone(), 1, pparams.clone(), ts);
    let proof = prover.prove(&[cr1cs.clone(), cr1cs.clone(), cr1cs.clone()]);

    c.bench_function("verify", |b| {
        b.iter_batched(
            || {
                PlusVerifier::init(
                    A.clone(),
                    M.clone(),
                    pparams.clone(),
                    PoseidonTranscript::empty::<PC>(),
                )
            },
            |mut verifier| {
                black_box(verifier.verify(&proof));
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
