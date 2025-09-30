#![allow(dead_code)]

use std::{fmt::Debug, time::Instant};

use ark_serialize::{CanonicalSerialize, Compress};
use ark_std::{vec::Vec, UniformRand};
use cyclotomic_rings::{challenge_set::LatticefoldChallengeSet, rings::SuitableRing};
use latticefold::{
    arith::{
        ccs::get_test_dummy_degree_three_ccs_non_scalar, r1cs::get_test_dummy_z_split_ntt, Arith,
        Witness, CCCS, CCS, LCCCS,
    },
    commitment::AjtaiCommitmentScheme,
    nifs::{
        linearization::{LFLinearizationProver, LinearizationProver},
        NIFSProver, NIFSVerifier,
    },
    transcript::poseidon::PoseidonTranscript,
};

include!(concat!(env!("OUT_DIR"), "/examples_generated.rs"));

#[allow(dead_code)]
pub fn wit_and_ccs_gen_degree_three_non_scalar<
    P: DecompositionParams,
    R: Clone + UniformRand + Debug + SuitableRing,
>(
    x_len: usize,
    n: usize,
    wit_len: usize,
    r1cs_rows: usize,
    kappa: usize,
) -> (CCCS<R>, Witness<R>, CCS<R>, AjtaiCommitmentScheme<R>) {
    let mut rng = ark_std::test_rng();

    let new_r1cs_rows = if P::L == 1 && (wit_len > 0 && (wit_len & (wit_len - 1)) == 0) {
        r1cs_rows - 2
    } else {
        r1cs_rows // This makes a square matrix but is too much memory
    };
    let (one, x_ccs, w_ccs) = get_test_dummy_z_split_ntt::<R>(x_len, wit_len);

    let mut z = vec![one];
    z.extend(&x_ccs);
    z.extend(&w_ccs);
    let ccs: CCS<R> =
        get_test_dummy_degree_three_ccs_non_scalar::<R>(&z, x_len, n, wit_len, P::L, new_r1cs_rows);
    ccs.check_relation(&z).expect("R1CS invalid!");

    let scheme: AjtaiCommitmentScheme<R> = AjtaiCommitmentScheme::rand(kappa, n, &mut rng);
    let wit: Witness<R> = Witness::from_w_ccs::<P>(w_ccs);

    let cm_i: CCCS<R> = CCCS {
        cm: wit.commit::<P>(&scheme).unwrap(),
        x_ccs,
    };

    (cm_i, wit, ccs, scheme)
}

#[allow(clippy::type_complexity)]
fn setup_example_environment<
    RqNTT: SuitableRing,
    DP: DecompositionParams,
    CS: LatticefoldChallengeSet<RqNTT>,
>() -> (
    LCCCS<RqNTT>,
    Witness<RqNTT>,
    CCCS<RqNTT>,
    Witness<RqNTT>,
    CCS<RqNTT>,
    AjtaiCommitmentScheme<RqNTT>,
) {
    let r1cs_rows = X_LEN + WIT_LEN + 1;

    let (cm_i, wit, ccs, scheme) =
        wit_and_ccs_gen_degree_three_non_scalar::<DP, RqNTT>(X_LEN, N, WIT_LEN, r1cs_rows, KAPPA);

    let rand_w_ccs: Vec<RqNTT> = (0..WIT_LEN).map(|i| RqNTT::from(i as u64)).collect();
    let wit_acc = Witness::from_w_ccs::<DP>(rand_w_ccs);

    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();

    let (acc, _) = LFLinearizationProver::<_, PoseidonTranscript<RqNTT, CS>>::prove(
        &cm_i,
        &wit_acc,
        &mut transcript,
        &ccs,
    )
    .expect("Failed to generate linearization proof");

    (acc, wit_acc, cm_i, wit, ccs, scheme)
}

type T = PoseidonTranscript<RqNTT, CS>;

fn main() {
    println!("Setting up example environment...");

    println!("Decomposition parameters:");
    println!("\tB: {}", DP::B);
    println!("\tL: {}", DP::L);
    println!("\tB_SMALL: {}", DP::B_SMALL);
    println!("\tK: {}", DP::K);

    let (acc, wit_acc, cm_i, wit_i, ccs, scheme) = setup_example_environment::<RqNTT, DP, CS>();

    let mut prover_transcript = PoseidonTranscript::<RqNTT, CS>::default();
    let mut verifier_transcript = PoseidonTranscript::<RqNTT, CS>::default();
    println!("Generating proof...");
    let start = Instant::now();

    let (_, _, proof) = NIFSProver::<RqNTT, DP, T>::prove(
        &acc,
        &wit_acc,
        &cm_i,
        &wit_i,
        &mut prover_transcript,
        &ccs,
        &scheme,
    )
    .unwrap();
    let duration = start.elapsed();
    println!("Proof generated in {:?}", duration);

    let mut serialized_proof = Vec::new();

    println!("Serializing proof (with compression)...");
    proof
        .serialize_with_mode(&mut serialized_proof, Compress::Yes)
        .unwrap();
    let compressed_size = serialized_proof.len();
    println!(
        "Proof size (with compression) size: {}",
        humansize::format_size(compressed_size, humansize::BINARY)
    );

    println!("Serializing proof (without compression)...");
    proof
        .serialize_with_mode(&mut serialized_proof, Compress::No)
        .unwrap();
    let uncompressed_size = serialized_proof.len();
    println!(
        "Proof (without compression) size: {}",
        humansize::format_size(uncompressed_size, humansize::BINARY)
    );

    println!("Verifying proof");
    let start = Instant::now();
    NIFSVerifier::<RqNTT, DP, T>::verify(&acc, &cm_i, &proof, &mut verifier_transcript, &ccs)
        .unwrap();
    let duration = start.elapsed();
    println!("Proof verified in {:?}", duration);
}
