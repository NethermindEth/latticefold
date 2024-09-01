use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use rand::thread_rng;
use std::time::Duration;

use latticefold::{
    arith::{r1cs::R1CS, Arith, Witness, CCCS, CCS},
    commitment::AjtaiCommitmentScheme,
    nifs::linearization::{
        LFLinearizationProver, LFLinearizationVerifier, LinearizationProver, LinearizationVerifier,
    },
    parameters::{
        DecompositionParamData, DecompositionParams, DilithiumTestParams, Pow2_57TestParams,
        Pow2_59TestParams, DILITHIUM_PRIME, POW2_57_PRIME, POW2_59_PRIME,
    },
    transcript::{self, poseidon::PoseidonTranscript},
};
use lattirust_arithmetic::{
    challenge_set::latticefold_challenge_set::BinarySmallSet,
    linear_algebra::SparseMatrix,
    ring::{Pow2CyclotomicPolyRingNTT, Ring},
};

fn prover_linearization_benchmark<
    const Q: u64,
    const N: usize,
    const C: usize,
    const W: usize,
    P: DecompositionParams,
>(
    c: &mut Criterion,
    p: P,
    prime_name: &str,
) {
    let ccs = get_test_ccs::<Pow2CyclotomicPolyRingNTT<Q, N>, W>();
    let (one, x_ccs, w_ccs) = get_test_z_split::<Pow2CyclotomicPolyRingNTT<Q, N>, W>();
    let mut z = vec![one];
    z.extend(&x_ccs);
    z.extend(&w_ccs);
    match ccs.check_relation(&z) {
        Ok(_) => println!("R1CS valid!"),
        Err(_) => println!("R1CS invalid"),
    }

    let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());

    let wit: Witness<Pow2CyclotomicPolyRingNTT<Q, N>> =
        Witness::from_w_ccs::<Pow2CyclotomicPolyRingNTT<Q, N>, P>(&w_ccs);
    let cm_i: CCCS<C, Pow2CyclotomicPolyRingNTT<Q, N>> = CCCS {
        cm: wit
            .commit::<C, W, Pow2CyclotomicPolyRingNTT<Q, N>, P>(&scheme)
            .unwrap(),
        x_ccs,
    };
    
    c.bench_with_input(
        BenchmarkId::new(format!("Linearization Prover {}", prime_name), DecompositionParamData::from(p)),
        &(cm_i, wit, ccs),
        |b, (cm_i, wit, ccs)| {
            b.iter(|| {
                println!();
                let mut transcript = PoseidonTranscript::<Pow2CyclotomicPolyRingNTT<Q, N>, BinarySmallSet<Q, N>>::default();
                let _ = LFLinearizationProver::<
                    _,
                    PoseidonTranscript<Pow2CyclotomicPolyRingNTT<Q, N>, BinarySmallSet<Q, N>>,
                >::prove(cm_i, wit, &mut transcript, ccs);
            })
        },
    );
}

fn verifier_linearization_benchmark<
    const Q: u64,
    const N: usize,
    const C: usize,
    const W: usize,
    P: DecompositionParams,
>(
    c: &mut Criterion,
    p: P,
    prime_name: &str,
) {
    let ccs = get_test_ccs::<Pow2CyclotomicPolyRingNTT<Q, N>, W>();
    let (one, x_ccs, w_ccs) = get_test_z_split::<Pow2CyclotomicPolyRingNTT<Q, N>, W>();
    let mut z = vec![one];
    z.extend(&x_ccs);
    z.extend(&w_ccs);
    match ccs.check_relation(&z) {
        Ok(_) => println!("R1CS valid!"),
        Err(_) => println!("R1CS invalid"),
    }

    let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());

    let wit: Witness<Pow2CyclotomicPolyRingNTT<Q, N>> =
        Witness::from_w_ccs::<Pow2CyclotomicPolyRingNTT<Q, N>, P>(&w_ccs);
    let cm_i: CCCS<C, Pow2CyclotomicPolyRingNTT<Q, N>> = CCCS {
        cm: wit
            .commit::<C, W, Pow2CyclotomicPolyRingNTT<Q, N>, P>(&scheme)
            .unwrap(),
        x_ccs,
    };
    let mut transcript =
        PoseidonTranscript::<Pow2CyclotomicPolyRingNTT<Q, N>, BinarySmallSet<Q, N>>::default();
    let res = LFLinearizationProver::<
        _,
        PoseidonTranscript<Pow2CyclotomicPolyRingNTT<Q, N>, BinarySmallSet<Q, N>>,
    >::prove(&cm_i, &wit, &mut transcript, &ccs);

    let mut transcript =
        PoseidonTranscript::<Pow2CyclotomicPolyRingNTT<Q, N>, BinarySmallSet<Q, N>>::default();

    c.bench_with_input(
        BenchmarkId::new(format!("Linearization Verifier {}", prime_name), DecompositionParamData::from(p)),
        &(cm_i, res.unwrap().1, ccs),
        |b, (cm_i, proof, ccs)| {
            b.iter(|| {
                let _ = LFLinearizationVerifier::<
                    _,
                    PoseidonTranscript<Pow2CyclotomicPolyRingNTT<Q, N>, BinarySmallSet<Q, N>>,
                >::verify(&cm_i, &proof, &mut transcript, &ccs);
            })
        },
    );
}

fn linearization_benchmarks(c: &mut Criterion) {
    prover_linearization_benchmark::<DILITHIUM_PRIME, 256, 4, 4, _>(
        c,
        DilithiumTestParams,
        "Dilithium prime",
    );
    // prover_linearization_benchmark::<POW2_59_PRIME, 256, 9, { 1 << 15 }, _>(c, Pow2_59TestParams);
    // prover_linearization_benchmark::<POW2_57_PRIME, 256, 9, { 1 << 15 }, _>(c, Pow2_57TestParams);
    verifier_linearization_benchmark::<DILITHIUM_PRIME, 256, 9, { 1 << 15 }, _>(
        c,
        DilithiumTestParams,
        "Dilithium prime"
    );
    // verifier_linearization_benchmark::<POW2_59_PRIME, 256, 9, { 1 << 15 }, _>(c, Pow2_59TestParams);
    // verifier_linearization_benchmark::<POW2_57_PRIME, 256, 9, { 1 << 15 }, _>(c, Pow2_57TestParams);
}

fn benchmarks_main(c: &mut Criterion) {
    linearization_benchmarks(c);
}

pub fn get_test_z_split<R: Ring, const W: usize>() -> (R, Vec<R>, Vec<R>) {
    // z = (1, io, w)
    (
        R::one(),
        to_F_vec(vec![
            1, // io
        ]),
        to_F_vec(vec![1; W/2]) // This should be the witness size but is failing
    )
}
pub fn get_test_ccs<R: Ring, const W: usize>() -> CCS<R> {
    let r1cs = get_test_r1cs::<R, W>();
    CCS::<R>::from_r1cs(r1cs)
}
pub fn get_test_r1cs<R: Ring, const W: usize>() -> R1CS<R> {
    let A = to_F_matrix::<R>(create_identity_matrix(W));
    let B = A.clone();
    let C = A.clone();

    R1CS::<R> { l: 1, A, B, C }
}
pub fn to_F_matrix<R: Ring>(M: Vec<Vec<usize>>) -> SparseMatrix<R> {
    to_F_dense_matrix::<R>(M).as_slice().into()
}
pub fn to_F_dense_matrix<R: Ring>(M: Vec<Vec<usize>>) -> Vec<Vec<R>> {
    M.iter()
        .map(|m| m.iter().map(|r| R::from(*r as u64)).collect())
        .collect()
}
pub fn to_F_vec<R: Ring>(z: Vec<usize>) -> Vec<R> {
    z.iter().map(|c| R::from(*c as u64)).collect()
}

fn create_identity_matrix(size: usize) -> Vec<Vec<usize>> {
    let mut matrix = Vec::with_capacity(size);
    
    for i in 0..size {
        let mut row = vec![0; size];
        row[i] = 1;
        matrix.push(row);
    }
    
    matrix
}

criterion_group!(
    name=benches;
    config = Criterion::default().sample_size(10).measurement_time(Duration::from_secs(50)).warm_up_time(Duration::from_secs(1));
    targets = benchmarks_main);
criterion_main!(benches);
