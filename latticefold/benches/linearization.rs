use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use cyclotomic_rings::{
    challenge_set::{BinarySmallSet, LatticefoldChallengeSet},
    SuitableRing,
};
use lattirust_ring::cyclotomic_ring::models::pow2_debug::Pow2CyclotomicPolyRingNTT;
use rand::thread_rng;
mod utils;
use std::time::Duration;
use utils::{get_test_dummy_ccs, get_test_dummy_z_split};

use latticefold::{
    arith::{Arith, Witness, CCCS, CCS, LCCCS},
    commitment::AjtaiCommitmentScheme,
    nifs::linearization::{
        LFLinearizationProver, LFLinearizationVerifier, LinearizationProof, LinearizationProver,
        LinearizationVerifier,
    },
    parameters::{
        DecompositionParamData, DecompositionParams
    },
    transcript::poseidon::PoseidonTranscript,
};
/*
fn wit_and_ccs_gen<
    const IO: usize,
    const C: usize,
    const WIT_LEN: usize,
    const W: usize,
    P: DecompositionParams,
    R: SuitableRing,
>(
    r1cs_rows: usize,
) -> (CCCS<C, R>, Witness<R>, CCS<R>) {
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
    (cm_i, wit, ccs)
}

fn prover_linearization_benchmark<
    const C: usize,
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
) -> (LCCCS<C, R>, LinearizationProof<R>) {
    let mut transcript = PoseidonTranscript::<R, CS>::default();
    let res = LFLinearizationProver::<_, PoseidonTranscript<R, CS>>::prove(
        &cm_i,
        &wit,
        &mut transcript,
        &ccs,
    );
    c.bench_with_input(
        BenchmarkId::new(
            format!("Linearization Prover {}", prime_name),
            DecompositionParamData::from(p),
        ),
        &(cm_i, wit, ccs),
        |b, (cm_i, wit, ccs)| {
            b.iter(|| {
                let _ = LFLinearizationProver::<_, PoseidonTranscript<R, CS>>::prove(
                    cm_i,
                    wit,
                    &mut transcript,
                    ccs,
                );
            })
        },
    );
    res.unwrap()
}

fn verifier_linearization_benchmark<
    const C: usize,
    P: DecompositionParams,
    R: SuitableRing,
    CS: LatticefoldChallengeSet<R>,
>(
    c: &mut Criterion,
    p: P,
    prime_name: &str,
    cm_i: &CCCS<C, R>,
    ccs: &CCS<R>,
    proof: (LCCCS<C, R>, LinearizationProof<R>),
) {
    c.bench_with_input(
        BenchmarkId::new(
            format!("Linearization Verifier {}", prime_name),
            DecompositionParamData::from(p),
        ),
        &(cm_i, proof.1, ccs),
        |b, (cm_i, proof, ccs)| {
            b.iter(|| {
                let mut transcript = PoseidonTranscript::<R, CS>::default();
                let _ = LFLinearizationVerifier::<_, PoseidonTranscript<R, CS>>::verify(
                    &cm_i,
                    &proof,
                    &mut transcript,
                    &ccs,
                );
            })
        },
    );
}

const IO: usize = 1;
const C: usize = 10;
const WIT_LEN: usize = 1 << 10;
fn linearization_benchmarks(c: &mut Criterion) {
    const W: usize = WIT_LEN * SomeFermatTestParams::L;
    let r1cs_rows = 5;
    let (cm_i, wit, ccs) = wit_and_ccs_gen::<
        IO,
        C,
        WIT_LEN,
        W,
        SomeFermatTestParams,
        Pow2CyclotomicPolyRingNTT<{ SOME_FERMAT_PRIME }, 16>,
    >(r1cs_rows);

    let proof = prover_linearization_benchmark::<
        C,
        _,
        Pow2CyclotomicPolyRingNTT<{ SOME_FERMAT_PRIME }, 16>,
        BinarySmallSet<{ SOME_FERMAT_PRIME }, 16>,
    >(
        c,
        SomeFermatTestParams,
        "Some fermat prime(2^16 + 1)",
        &cm_i,
        &wit,
        &ccs,
    );
    verifier_linearization_benchmark::<
        C,
        _,
        Pow2CyclotomicPolyRingNTT<{ SOME_FERMAT_PRIME }, 16>,
        BinarySmallSet<{ SOME_FERMAT_PRIME }, 16>,
    >(
        c,
        SomeFermatTestParams,
        "Some fermat prime(2^16 + 1)",
        &cm_i,
        &ccs,
        proof,
    );
}

fn benchmarks_main(c: &mut Criterion) {
    linearization_benchmarks(c);
}

criterion_group!(
    name=benches;
    config = Criterion::default().sample_size(10).measurement_time(Duration::from_secs(50)).warm_up_time(Duration::from_secs(1));
    targets = benchmarks_main);
criterion_main!(benches);
*/