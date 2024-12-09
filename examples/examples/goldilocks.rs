use ark_serialize::{CanonicalSerialize, Compress};
use ark_std::rand::Rng;
use ark_std::test_rng;
use cyclotomic_rings::challenge_set::LatticefoldChallengeSet;
use cyclotomic_rings::rings::{GoldilocksChallengeSet, GoldilocksRingNTT, SuitableRing};
use latticefold::arith::r1cs::{get_test_z_split, to_F_matrix, R1CS};
use latticefold::arith::{Witness, CCCS, CCS, LCCCS};
use latticefold::commitment::AjtaiCommitmentScheme;
use latticefold::nifs::linearization::{LFLinearizationProver, LinearizationProver};
use latticefold::nifs::{NIFSProver, NIFSVerifier};
use latticefold::transcript::poseidon::PoseidonTranscript;
use lattirust_ring::Ring;
use std::time::Instant;

include!(concat!(env!("OUT_DIR"), "/generated.rs"));

/// Generates a test R1CS matrix for the equation x^3 + x + 5 = y
/// This will be changed to be easier to customized
pub fn get_test_r1cs<R: Ring>() -> R1CS<R> {
    let a = to_F_matrix::<R>(vec![
        vec![1, 0, 0, 0, 0, 0],
        vec![0, 0, 0, 1, 0, 0],
        vec![1, 0, 0, 0, 1, 0],
        vec![0, 5, 0, 0, 0, 1],
    ]);
    let b = to_F_matrix::<R>(vec![
        vec![1, 0, 0, 0, 0, 0],
        vec![1, 0, 0, 0, 0, 0],
        vec![0, 1, 0, 0, 0, 0],
        vec![0, 1, 0, 0, 0, 0],
    ]);
    let c = to_F_matrix::<R>(vec![
        vec![0, 0, 0, 1, 0, 0],
        vec![0, 0, 0, 0, 1, 0],
        vec![0, 0, 0, 0, 0, 1],
        vec![0, 0, 1, 0, 0, 0],
    ]);

    R1CS {
        l: 1,
        A: a,
        B: b,
        C: c,
    }
}

/// Generates a CCS from a test R1CS matrix with specified parameters
pub fn get_test_ccs<R: Ring>(w: usize, l: usize) -> CCS<R> {
    let r1cs = get_test_r1cs::<R>();
    CCS::<R>::from_r1cs_padded(r1cs, w, l)
}

#[allow(clippy::type_complexity)]
fn setup_example_environment<
    const C: usize,
    RqNTT: SuitableRing,
    DP: DecompositionParams,
    const W: usize,
    const WIT_LEN: usize,
    CS: LatticefoldChallengeSet<RqNTT>,
>() -> (
    LCCCS<C, RqNTT>,
    Witness<RqNTT>,
    CCCS<C, RqNTT>,
    Witness<RqNTT>,
    CCS<RqNTT>,
    AjtaiCommitmentScheme<C, W, RqNTT>,
) {
    let ccs = get_test_ccs::<RqNTT>(W, DP::L);
    let mut rng = test_rng();
    let (_, x_ccs, w_ccs) = get_test_z_split::<RqNTT>(rng.gen_range(0..64));
    let scheme = AjtaiCommitmentScheme::rand(&mut rng);

    let wit_i = Witness::from_w_ccs::<DP>(w_ccs);
    let cm_i = CCCS {
        cm: wit_i.commit::<C, W, DP>(&scheme).unwrap(),
        x_ccs: x_ccs.clone(),
    };

    let rand_w_ccs: Vec<RqNTT> = (0..WIT_LEN).map(|i| RqNTT::from(i as u64)).collect();
    let wit_acc = Witness::from_w_ccs::<DP>(rand_w_ccs);

    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();
    let (acc, _) = LFLinearizationProver::<_, PoseidonTranscript<RqNTT, CS>>::prove(
        &cm_i,
        &wit_acc,
        &mut transcript,
        &ccs,
    )
    .unwrap();

    (acc, wit_acc, cm_i, wit_i, ccs, scheme)
}

type RqNTT = GoldilocksRingNTT;
type CS = GoldilocksChallengeSet;
type T = PoseidonTranscript<RqNTT, CS>;

fn main() {
    println!("Setting up example environment...");
    let (acc, wit_acc, cm_i, wit_i, ccs, scheme) =
        setup_example_environment::<C, RqNTT, ExampleDP, W, WIT_LEN, CS>();

    let mut prover_transcript = PoseidonTranscript::<RqNTT, CS>::default();
    let mut verifier_transcript = PoseidonTranscript::<RqNTT, CS>::default();
    println!("Generating proof...");
    let start = Instant::now();

    let (_, _, proof) = NIFSProver::<C, W, RqNTT, ExampleDP, T>::prove(
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
    NIFSVerifier::<C, RqNTT, ExampleDP, T>::verify(
        &acc,
        &cm_i,
        &proof,
        &mut verifier_transcript,
        &ccs,
    )
    .unwrap();
    let duration = start.elapsed();
    println!("Proof verified in {:?}", duration);
}
