//! LF+ E2E (All protocols) prove and verify

#![allow(missing_docs)]
#![allow(non_snake_case)]

use ark_ff::PrimeField;
use ark_std::time::Duration;
use criterion::{criterion_group, criterion_main, Criterion};
use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
use latticefold::{arith::r1cs::R1CS, transcript::poseidon::PoseidonTS};
use latticefold_plus::{
    lin::{LinParameters, Linearize, Verify},
    mlin::Mlin,
    plus::{PlusParameters, PlusProver, PlusVerifier},
    r1cs::{r1cs_decomposed_square, ComR1CS},
    rgchk::DecompParameters,
    utils::estimate_bound,
};
use rand::prelude::*;
use stark_rings::{
    balanced_decomposition::GadgetDecompose, cyclotomic_ring::models::frog_ring::RqPoly as R,
    PolyRing, Ring,
};
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

    let mut rng = rand::thread_rng();
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

    let A = Matrix::<R>::rand(&mut rand::thread_rng(), params.kappa, n);

    let cr1cs = ComR1CS::new(r1cs, z, 1, B, k, &A);

    let M = cr1cs.x.matrices();

    // Prover / Fold
    c.bench_function("prove", |b| {
        b.iter_batched(
            || PoseidonTS::default::<PC>(),
            |mut ts| {
                let (linb, _lproof) = cr1cs.linearize(&mut ts);
                // L=3 (equal) instances are folded here
                // TODO Do accumulated instance (2) + one online (1)
                let mlin = Mlin {
                    lins: vec![linb.clone(), linb.clone(), linb],
                    params: params.clone(),
                };
                let prover = PlusProver {
                    instances: mlin,
                    A: A.clone(),
                    params: PlusParameters {
                        lin: params.clone(),
                        B,
                    },
                };
                let (_acc, _proof) = prover.prove(&M, &mut ts);
            },
            criterion::BatchSize::SmallInput,
        )
    });

    // Verifier
    let mut ts = PoseidonTS::default::<PC>();
    let (linb, lproof) = cr1cs.linearize(&mut ts);
    let mlin = Mlin {
        lins: vec![linb.clone(), linb.clone(), linb],
        params: params.clone(),
    };
    let prover = PlusProver {
        instances: mlin,
        A: A.clone(),
        params: PlusParameters {
            lin: params.clone(),
            B,
        },
    };
    let (_acc, proof) = prover.prove(&M, &mut ts);
    c.bench_function("verify", |b| {
        b.iter_batched(
            || PoseidonTS::default::<PC>(),
            |mut ts_v| {
                let verifier = PlusVerifier {
                    A: A.clone(),
                    params: PlusParameters {
                        lin: params.clone(),
                        B,
                    },
                };
                lproof.verify(&mut ts_v);
                verifier.verify(&proof, &M, &mut ts_v);
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
