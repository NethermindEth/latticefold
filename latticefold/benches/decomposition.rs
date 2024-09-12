use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use cyclotomic_rings::{
    challenge_set::{BinarySmallSet, LatticefoldChallengeSet},
    SuitableRing,
};
use lattirust_ring::Pow2CyclotomicPolyRingNTT;
use rand::thread_rng;
mod utils;
use std::time::Duration;
use utils::{get_test_dummy_ccs, get_test_dummy_z_split};

use latticefold::{
    arith::{Arith, Witness, CCCS, CCS},
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
        DecompositionParamData, DecompositionParams, SomeFermatTestParams, SOME_FERMAT_PRIME,
    },
    transcript::poseidon::PoseidonTranscript,
};

fn wit_and_ccs_gen<
    const IO: usize,
    const C: usize,
    const WIT_LEN: usize,
    const W: usize,
    P: DecompositionParams,
    R: SuitableRing,
>(
    r1cs_rows: usize,
) -> (
    CCCS<C, R>,
    Witness<R>,
    CCS<R>,
    AjtaiCommitmentScheme<C, W, R>,
) {
    let ccs = get_test_dummy_ccs::<R, IO, WIT_LEN>(r1cs_rows);
    let (one, x_ccs, w_ccs) = get_test_dummy_z_split::<R, IO, WIT_LEN>();
    let mut z = vec![one];
    z.extend(&x_ccs);
    z.extend(&w_ccs);
    match ccs.check_relation(&z) {
        Ok(_) => println!("R1CS valid!"),
        Err(_) => println!("R1CS invalid"),
    }

    let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());

    let wit: Witness<R> = Witness::from_w_ccs::<P>(&w_ccs);
    let cm_i: CCCS<C, R> = CCCS {
        cm: wit.commit::<C, W, P>(&scheme).unwrap(),
        x_ccs,
    };

    (cm_i, wit, ccs, scheme)
}

fn prover_decomposition_benchmark<
    const C: usize,
    const W: usize,
    P: DecompositionParams,
    R: SuitableRing,
    CS: LatticefoldChallengeSet<R>,
>(
    c: &mut Criterion,
    p: P,
    prime_name: &str,
    cm_i: &CCCS<C, R>,
    wit: &Witness<R>,
    ccs: &CCS<R>,
    scheme: &AjtaiCommitmentScheme<C, W, R>,
) {
    let mut prover_transcript = PoseidonTranscript::<R, CS>::default();
    let mut verifier_transcript = PoseidonTranscript::<R, CS>::default();

    let (_, linearization_proof) = LFLinearizationProver::<_, PoseidonTranscript<R, CS>>::prove(
        &cm_i,
        &wit,
        &mut prover_transcript,
        &ccs,
    )
    .unwrap();

    let lcccs = LFLinearizationVerifier::<_, PoseidonTranscript<R, CS>>::verify(
        &cm_i,
        &linearization_proof,
        &mut verifier_transcript,
        &ccs,
    )
    .unwrap();

    c.bench_with_input(
        BenchmarkId::new(
            format!("Decomposition Prover {}", prime_name),
            DecompositionParamData::from(p),
        ),
        &(lcccs, wit, ccs),
        |b, (lcccs, wit, ccs)| {
            b.iter(|| {
                let (_, _, _) = LFDecompositionProver::<_, PoseidonTranscript<R, CS>>::prove::<
                    W,
                    C,
                    P,
                >(lcccs, &wit, &mut prover_transcript, &ccs, scheme)
                .unwrap();
            })
        },
    );
}

fn verifier_decomposition_benchmark<
    const C: usize,
    const W: usize,
    P: DecompositionParams,
    R: SuitableRing,
    CS: LatticefoldChallengeSet<R>,
>(
    c: &mut Criterion,
    p: P,
    prime_name: &str,
    cm_i: &CCCS<C, R>,
    wit: &Witness<R>,
    ccs: &CCS<R>,
    scheme: &AjtaiCommitmentScheme<C, W, R>,
) {
    let mut prover_transcript = PoseidonTranscript::<R, CS>::default();
    let mut verifier_transcript = PoseidonTranscript::<R, CS>::default();

    let (_, linearization_proof) = LFLinearizationProver::<_, PoseidonTranscript<R, CS>>::prove(
        &cm_i,
        &wit,
        &mut prover_transcript,
        &ccs,
    )
    .unwrap();

    let lcccs = LFLinearizationVerifier::<_, PoseidonTranscript<R, CS>>::verify(
        &cm_i,
        &linearization_proof,
        &mut verifier_transcript,
        &ccs,
    )
    .unwrap();

    let (_, _, decomposition_proof) =
        LFDecompositionProver::<_, PoseidonTranscript<R, CS>>::prove::<W, C, P>(
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
                let _ = LFDecompositionVerifier::<_, PoseidonTranscript<R, CS>>::verify::<C, P>(
                    lcccs,
                    proof,
                    &mut verifier_transcript,
                    &ccs,
                );
            })
        },
    );
}

const IO: usize = 1;
const C: usize = 10;
const WIT_LEN: usize = 1 << 10;
fn decomposition_benchmarks(c: &mut Criterion) {
    const W: usize = WIT_LEN * SomeFermatTestParams::L;
    let r1cs_rows = 4;
    let (cm_i, wit, ccs, scheme) = wit_and_ccs_gen::<
        IO,
        C,
        WIT_LEN,
        W,
        SomeFermatTestParams,
        Pow2CyclotomicPolyRingNTT<{ SOME_FERMAT_PRIME }, 16>,
    >(r1cs_rows);

    prover_decomposition_benchmark::<
        C,
        W,
        _,
        Pow2CyclotomicPolyRingNTT<SOME_FERMAT_PRIME, 16>,
        BinarySmallSet<SOME_FERMAT_PRIME, 16>,
    >(
        c,
        SomeFermatTestParams,
        "Some fermat(2^16 + 1) prime",
        &cm_i,
        &wit,
        &ccs,
        &scheme,
    );

    verifier_decomposition_benchmark::<
        C,
        W,
        _,
        Pow2CyclotomicPolyRingNTT<SOME_FERMAT_PRIME, 16>,
        BinarySmallSet<SOME_FERMAT_PRIME, 16>,
    >(
        c,
        SomeFermatTestParams,
        "Some fermat(2^16 + 1) prime",
        &cm_i,
        &wit,
        &ccs,
        &scheme,
    );
}

fn benchmarks_main(c: &mut Criterion) {
    decomposition_benchmarks(c);
}

criterion_group!(
    name=benches;
    config = Criterion::default().sample_size(10).measurement_time(Duration::from_secs(50)).warm_up_time(Duration::from_secs(1));
    targets = benchmarks_main);
criterion_main!(benches);
