//! Sequential IVC-style folding of 1024 Ajtai-committed steps.
//!
//! Each "step" is a degree-three CCS instance whose witness is committed with
//! the Ajtai commitment scheme. Starting from a linearized accumulator, we fold
//! the incoming committed instance into the accumulator 1024 times in sequence,
//! threading a single Fiat-Shamir transcript through the whole chain (exactly as
//! an IVC prover would). After each fold the corresponding `NIFSVerifier` step is
//! run so the chain is verified end-to-end.
//!
//! This is the native LatticeFold verification path (there is no on-chain
//! decider in this codebase — the final artifact is a lattice proof verified in
//! Rust).
//!
//! Run with:
//!   cargo run --release --example fold_1024

use std::{fmt::Debug, time::Instant};

use ark_serialize::{CanonicalSerialize, Compress};
use ark_std::UniformRand;
use cyclotomic_rings::rings::{GoldilocksChallengeSet, GoldilocksRingNTT, SuitableRing};
use latticefold::{
    arith::{
        ccs::get_test_dummy_degree_three_ccs_non_scalar, r1cs::get_test_dummy_z_split_ntt, Arith,
        Witness, CCCS, CCS,
    },
    commitment::AjtaiCommitmentScheme,
    decomposition_parameters::DecompositionParams,
    nifs::{
        linearization::{LFLinearizationProver, LinearizationProver},
        NIFSProver, NIFSVerifier,
    },
    transcript::poseidon::PoseidonTranscript,
};

// ---- Concrete instantiation ------------------------------------------------

type RqNTT = GoldilocksRingNTT;
type CS = GoldilocksChallengeSet;
type T = PoseidonTranscript<RqNTT, CS>;

/// Decomposition parameters (Goldilocks defaults from examples/README.md).
#[derive(Clone)]
struct DP {}
impl DecompositionParams for DP {
    const B: u128 = 1 << 15;
    const L: usize = 5;
    const B_SMALL: usize = 2;
    const K: usize = 15;
}

const X_LEN: usize = 1;
const WIT_LEN: usize = 4;
const KAPPA: usize = 4;
const N: usize = WIT_LEN * DP::L;

/// Number of sequential folding steps.
const STEPS: usize = 1024;

// ---- Step-circuit / instance generation ------------------------------------

/// Build one Ajtai-committed CCS step instance (the incoming `cm_i` + witness),
/// plus the fixed CCS and commitment scheme shared by every step.
fn gen_step<P: DecompositionParams, R: Clone + UniformRand + Debug + SuitableRing>(
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
        r1cs_rows
    };

    let (one, x_ccs, w_ccs) = get_test_dummy_z_split_ntt::<R>(x_len, wit_len);

    let mut z = vec![one];
    z.extend(&x_ccs);
    z.extend(&w_ccs);

    let ccs: CCS<R> =
        get_test_dummy_degree_three_ccs_non_scalar::<R>(&z, x_len, n, wit_len, P::L, new_r1cs_rows);
    ccs.check_relation(&z).expect("step CCS relation invalid!");

    let scheme: AjtaiCommitmentScheme<R> = AjtaiCommitmentScheme::rand(kappa, n, &mut rng);
    let wit: Witness<R> = Witness::from_w_ccs::<P>(w_ccs);

    let cm_i: CCCS<R> = CCCS {
        cm: wit.commit::<P>(&scheme).unwrap(),
        x_ccs,
    };

    (cm_i, wit, ccs, scheme)
}

fn main() {
    println!("LatticeFold sequential folding — {STEPS} Ajtai-committed steps");
    println!("Ring: Goldilocks | KAPPA={KAPPA} WIT_LEN={WIT_LEN} N={N}");
    println!(
        "Decomposition: B={} L={} B_SMALL={} K={}",
        DP::B,
        DP::L,
        DP::B_SMALL,
        DP::K
    );

    let r1cs_rows = X_LEN + WIT_LEN + 1;

    // Fixed step circuit + commitment scheme, and the (repeated) incoming instance.
    let (cm_i, wit_i, ccs, scheme) =
        gen_step::<DP, RqNTT>(X_LEN, N, WIT_LEN, r1cs_rows, KAPPA);

    // Bootstrap the accumulator by linearizing an initial committed instance.
    let init_w: Vec<RqNTT> = (0..WIT_LEN).map(|i| RqNTT::from(i as u64)).collect();
    let mut w_acc = Witness::from_w_ccs::<DP>(init_w);

    let mut bootstrap_transcript = PoseidonTranscript::<RqNTT, CS>::default();
    let (mut acc, _) = LFLinearizationProver::<_, T>::prove(
        &cm_i,
        &w_acc,
        &mut bootstrap_transcript,
        &ccs,
    )
    .expect("failed to bootstrap accumulator");

    // One transcript per party, threaded through the entire IVC chain.
    let mut prover_transcript = PoseidonTranscript::<RqNTT, CS>::default();
    let mut verifier_transcript = PoseidonTranscript::<RqNTT, CS>::default();

    println!("\nFolding {STEPS} steps...");
    let start = Instant::now();
    let mut last_proof = None;

    for step in 0..STEPS {
        let (new_acc, new_w_acc, proof) = NIFSProver::<RqNTT, DP, T>::prove(
            &acc,
            &w_acc,
            &cm_i,
            &wit_i,
            &mut prover_transcript,
            &ccs,
            &scheme,
        )
        .expect("folding prover failed");

        // Verify this fold step against the same accumulator/instance.
        let verified_acc = NIFSVerifier::<RqNTT, DP, T>::verify(
            &acc,
            &cm_i,
            &proof,
            &mut verifier_transcript,
            &ccs,
        )
        .expect("folding verifier failed");

        // Prover and verifier must agree on the folded accumulator.
        assert_eq!(
            new_acc, verified_acc,
            "prover/verifier accumulator mismatch at step {step}"
        );

        acc = new_acc;
        w_acc = new_w_acc;

        if (step + 1) % 128 == 0 {
            println!("  step {:>4}/{STEPS} folded & verified", step + 1);
        }
        last_proof = Some(proof);
    }

    let elapsed = start.elapsed();
    println!("\nAll {STEPS} steps folded and verified in {elapsed:?}");
    println!("Average per step: {:?}", elapsed / STEPS as u32);

    if let Some(proof) = last_proof {
        let mut buf = Vec::new();
        proof.serialize_with_mode(&mut buf, Compress::Yes).unwrap();
        println!(
            "Per-step fold proof size (compressed): {}",
            humansize::format_size(buf.len(), humansize::BINARY)
        );
    }

    println!("\nFinal accumulator is the single folded LCCCS attesting all {STEPS} steps.");
}
