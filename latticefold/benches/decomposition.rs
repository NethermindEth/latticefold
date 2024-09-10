use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use rand::thread_rng;
use std::time::Duration;

use latticefold::{
    arith::{r1cs::R1CS, Witness, CCCS, CCS},
    commitment::AjtaiCommitmentScheme,
    nifs::{
        decomposition::{
            DecompositionProver, DecompositionVerifier, LFDecompositionProver,
            LFDecompositionVerifier,
        },
        linearization::{
            LFLinearizationProver, LFLinearizationVerifier, LinearizationProver,
            LinearizationVerifier,
        },
    },
    parameters::{
        DecompositionParamData, DecompositionParams, DilithiumTestParams, Pow2_57TestParams,
        Pow2_59TestParams, DILITHIUM_PRIME, POW2_57_PRIME, POW2_59_PRIME,
    },
    transcript::poseidon::PoseidonTranscript,
};
use lattirust_arithmetic::{
    challenge_set::latticefold_challenge_set::BinarySmallSet,
    linear_algebra::SparseMatrix,
    ring::{Pow2CyclotomicPolyRing, Pow2CyclotomicPolyRingNTT, Ring, Zq},
};

fn prover_decomposition_benchmark<
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
    let (_, x_ccs, w_ccs) = get_test_z_split::<Pow2CyclotomicPolyRingNTT<Q, N>, W>();
    let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());
    let wit: Witness<Pow2CyclotomicPolyRingNTT<Q, N>> =
        Witness::from_w_ccs::<Pow2CyclotomicPolyRing<Zq<Q>, N>, P>(&w_ccs);
    let cm_i: CCCS<C, Pow2CyclotomicPolyRingNTT<Q, N>> = CCCS {
        cm: wit
            .commit::<C, W, Pow2CyclotomicPolyRing<Zq<Q>, N>, P>(&scheme)
            .unwrap(),
        x_ccs,
    };

    let mut prover_transcript =
        PoseidonTranscript::<Pow2CyclotomicPolyRingNTT<Q, N>, BinarySmallSet<Q, N>>::default();
    let mut verifier_transcript =
        PoseidonTranscript::<Pow2CyclotomicPolyRingNTT<Q, N>, BinarySmallSet<Q, N>>::default();

    let (_, linearization_proof) = LFLinearizationProver::<
        _,
        PoseidonTranscript<Pow2CyclotomicPolyRingNTT<Q, N>, BinarySmallSet<Q, N>>,
    >::prove(&cm_i, &wit, &mut prover_transcript, &ccs)
    .unwrap();

    let lcccs = LFLinearizationVerifier::<
        _,
        PoseidonTranscript<Pow2CyclotomicPolyRingNTT<Q, N>, BinarySmallSet<Q, N>>,
    >::verify(&cm_i, &linearization_proof, &mut verifier_transcript, &ccs)
    .unwrap();

    c.bench_with_input(
        BenchmarkId::new(
            format!("Decomposition Prover {}", prime_name),
            DecompositionParamData::from(p),
        ),
        &(lcccs, wit, ccs),
        |b, (lcccs, wit, ccs)| {
            b.iter(|| {
                let (_, _, _) = LFDecompositionProver::<
                    _,
                    PoseidonTranscript<Pow2CyclotomicPolyRingNTT<Q, N>, BinarySmallSet<Q, N>>,
                >::prove::<W, C, Pow2CyclotomicPolyRing<Zq<Q>, N>, P>(
                    lcccs,
                    &wit,
                    &mut prover_transcript,
                    &ccs,
                    &scheme,
                )
                .unwrap();
            })
        },
    );
}

fn verifier_decomposition_benchmark<
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
    let (_, x_ccs, w_ccs) = get_test_z_split::<Pow2CyclotomicPolyRingNTT<Q, N>, W>();
    let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());
    let wit: Witness<Pow2CyclotomicPolyRingNTT<Q, N>> =
        Witness::from_w_ccs::<Pow2CyclotomicPolyRing<Zq<Q>, N>, P>(&w_ccs);
    let cm_i: CCCS<C, Pow2CyclotomicPolyRingNTT<Q, N>> = CCCS {
        cm: wit
            .commit::<C, W, Pow2CyclotomicPolyRing<Zq<Q>, N>, P>(&scheme)
            .unwrap(),
        x_ccs,
    };

    let mut prover_transcript =
        PoseidonTranscript::<Pow2CyclotomicPolyRingNTT<Q, N>, BinarySmallSet<Q, N>>::default();
    let mut verifier_transcript =
        PoseidonTranscript::<Pow2CyclotomicPolyRingNTT<Q, N>, BinarySmallSet<Q, N>>::default();

    let (_, linearization_proof) = LFLinearizationProver::<
        _,
        PoseidonTranscript<Pow2CyclotomicPolyRingNTT<Q, N>, BinarySmallSet<Q, N>>,
    >::prove(&cm_i, &wit, &mut prover_transcript, &ccs)
    .unwrap();

    let lcccs = LFLinearizationVerifier::<
        _,
        PoseidonTranscript<Pow2CyclotomicPolyRingNTT<Q, N>, BinarySmallSet<Q, N>>,
    >::verify(&cm_i, &linearization_proof, &mut verifier_transcript, &ccs)
    .unwrap();

    let (_, _, decomposition_proof) = LFDecompositionProver::<
        _,
        PoseidonTranscript<Pow2CyclotomicPolyRingNTT<Q, N>, BinarySmallSet<Q, N>>,
    >::prove::<W, C, Pow2CyclotomicPolyRing<Zq<Q>, N>, P>(
        &lcccs,
        &wit,
        &mut prover_transcript,
        &ccs,
        &scheme,
    )
    .unwrap();
    c.bench_with_input(
        BenchmarkId::new(
            format!("Decomposition Verifier {}", prime_name),
            DecompositionParamData::from(p),
        ),
        &(lcccs, decomposition_proof, ccs),
        |b, (lcccs, proof, ccs)| {
            b.iter(|| {
                let _ = LFDecompositionVerifier::<
                    _,
                    PoseidonTranscript<Pow2CyclotomicPolyRingNTT<Q, N>, BinarySmallSet<Q, N>>,
                >::verify::<C, P>(
                    lcccs, proof, &mut verifier_transcript, &ccs
                );
            })
        },
    );
}
fn decomposition_benchmarks(c: &mut Criterion) {
    prover_decomposition_benchmark::<DILITHIUM_PRIME, 256, 9, { 1 << 15 }, _>(
        c,
        DilithiumTestParams,
        "Dilithium prime",
    );
    prover_decomposition_benchmark::<POW2_59_PRIME, 256, 9, { 1 << 15 }, _>(
        c,
        Pow2_59TestParams,
        "p = 27 * (1 << 59) + 1",
    );
    prover_decomposition_benchmark::<POW2_57_PRIME, 256, 9, { 1 << 15 }, _>(
        c,
        Pow2_57TestParams,
        "p = 71 * (1 << 57) + 1",
    );
    verifier_decomposition_benchmark::<DILITHIUM_PRIME, 256, 9, { 1 << 15 }, _>(
        c,
        DilithiumTestParams,
        "Dilithium prime",
    );
    verifier_decomposition_benchmark::<POW2_59_PRIME, 256, 9, { 1 << 15 }, _>(
        c,
        Pow2_59TestParams,
        "p = 27 * (1 << 59) + 1",
    );
    verifier_decomposition_benchmark::<POW2_57_PRIME, 256, 9, { 1 << 15 }, _>(
        c,
        Pow2_57TestParams,
        "p = 71 * (1 << 57) + 1",
    );
}

fn benchmarks_main(c: &mut Criterion) {
    decomposition_benchmarks(c);
}

pub fn get_test_z_split<R: Ring, const W: usize>() -> (R, Vec<R>, Vec<R>) {
    // z = (1, io, w)
    (
        R::one(),
        to_F_vec(vec![
            1, // io
        ]),
        to_F_vec(vec![1; W / 2]), // This should be the witness size but is failing
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
