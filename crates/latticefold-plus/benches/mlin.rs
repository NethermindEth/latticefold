//! Benchmarks for Construction 5.2: Multilinear folding protocol.
//!
//! The multilinear folding protocol is the main LatticeFold+ operation, batching
//! multiple LinB instances (L>1) into a single LinB2 instance with improved
//! amortization. This enables efficient recursive composition for IVC/PCD schemes.
//!
//! ## Protocol Overview
//!
//! The multilinear folding operates in three phases:
//! 1. **Instance batching**: Takes L committed R1CS instances and linearizes each
//!    to produce L LinB instances
//! 2. **Multilinear aggregation**: Uses sumcheck protocols to aggregate all L
//!    instances into a single folded constraint system
//! 3. **Commitment transformation**: Converts L commitments into one LinB2 commitment
//!    via the commitment transformation protocol
//!
//! The key insight is that amortization improves with larger L, as many operations
//! (especially sumcheck rounds) can be shared across instances.
//!
//! ## Benchmarked Operations
//!
//! - **Prover**: Multilinear aggregation and commitment transformation across L instances
//! - **Verifier**: Sumcheck-based verification of the L×LinB → LinB2 transformation
//! - **Folding arity (L) scaling**: Performance across L ∈ [2,3,4,5,6,7,8]
//! - **Decomposition width (k) scaling**: Impact of k ∈ [2,3,4] on folding cost
//! - **Large witness scaling**: Performance on witnesses up to 512K elements
//! - **Security parameter (κ) scaling**: Impact of κ ∈ [2,3,4,5] on commitment operations
//!
//! ## Paper Reference
//! Construction 5.2, Section 5.2 of the LatticeFold+ paper

#![allow(non_snake_case)]

use criterion::Criterion;
use latticefold_plus::{
    lin::{LinParameters, Linearize},
    mlin::{LinB2, Mlin},
    r1cs::ComR1CS,
};
use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as R;
use stark_rings_linalg::{Matrix, SparseMatrix};

#[path = "utils/mod.rs"]
mod utils;
use utils::{
    helpers::{
        bench_prover_protocol, bench_rng, bench_verifier_protocol, create_ajtai_matrix,
        create_test_m_matrix, create_transcript, get_validated_decomp_params, ProverBenchmark,
        R1CSBuilder, VerifierBenchmark, WitnessPattern,
    },
    multilinear_fold,
};

// ============================================================================
// Setup Functions
// ============================================================================

/// Creates multilinear fold input for prover benchmarks.
///
/// Generates L committed R1CS instances, linearizes each to produce LinB instances,
/// and packages them into an `Mlin` structure. Uses varied witness patterns across
/// instances to avoid unrealistic uniformity in benchmarks.
fn setup_input(L: usize, n: usize, k: usize, kappa: usize, B: usize) -> Mlin<R> {
    let mut rng = bench_rng();
    let dparams = get_validated_decomp_params(k, kappa, n);

    let params = LinParameters {
        kappa,
        decomp: dparams,
    };

    let r1cs = R1CSBuilder::new(n, k, B as u128).build_decomposed_square();
    let A = create_ajtai_matrix(kappa, n, &mut rng);

    let mut lins = Vec::with_capacity(L);
    let mut ts = create_transcript();

    for i in 0..L {
        let z = if i == 0 {
            WitnessPattern::AllOnes.generate(n / k, &mut rng)
        } else {
            WitnessPattern::Custom(Box::new(move |idx| {
                if idx == 0 {
                    R::from((i % 10) as u128)
                } else {
                    R::from(1u128)
                }
            }))
            .generate(n / k, &mut rng)
        };

        let cr1cs = ComR1CS::new(r1cs.clone(), z, 1, B as u128, k, &A);
        let (linb, _lproof) = cr1cs.linearize(&mut ts);
        lins.push(linb);
    }

    Mlin { lins, params }
}

/// Generates a valid multilinear fold proof for verifier benchmarks.
///
/// Creates an `Mlin` input with L instances, executes the mlin protocol to
/// generate a `CmProof`, and validates it before returning. This ensures the
/// verifier benchmarks measure only verification time, not error handling.
fn setup_proof(
    L: usize,
    n: usize,
    k: usize,
    kappa: usize,
    B: usize,
) -> (Mlin<R>, latticefold_plus::cm::CmProof<R>) {
    let mlin = setup_input(L, n, k, kappa, B);
    let mut rng = bench_rng();
    let M = create_test_m_matrix(n);

    let A = create_ajtai_matrix(kappa, n, &mut rng);
    let mut ts = create_transcript();
    let (linb2, proof) = mlin.mlin(&A, &M, &mut ts);

    let mut verify_ts = create_transcript();
    proof
        .verify(&M, &mut verify_ts)
        .expect("Generated multilinear folding proof should be valid");

    (mlin, proof)
}

// ============================================================================
// Trait Implementations
// ============================================================================

/// Prover benchmark for standard multilinear folding protocol.
///
/// Measures L×LinB → LinB2 folding time across varying folding arity
/// L ∈ [2,3,4,5,6,7,8] while keeping n=65536, k=2, κ=2, and B=50 fixed.
/// Demonstrates batching efficiency as L increases.
struct MultilinearFoldProver;

impl ProverBenchmark for MultilinearFoldProver {
    type Input = (Mlin<R>, Matrix<R>, Vec<SparseMatrix<R>>);
    type Output = (LinB2<R>, latticefold_plus::cm::CmProof<R>);
    type Params = (usize, usize, usize, usize, usize);

    fn group_name() -> &'static str {
        "MultilinearFolding-Prover"
    }

    fn setup_input((L, n, k, kappa, B): Self::Params) -> Self::Input {
        let mut rng = bench_rng();
        let mlin = setup_input(L, n, k, kappa, B);
        let M = create_test_m_matrix(n);
        let A = create_ajtai_matrix(kappa, n, &mut rng);
        (mlin, A, M)
    }

    fn param_label((L, n, k, kappa, B): Self::Params) -> String {
        format!("L={}_n={}_k={}_κ={}_B={}", L, n, k, kappa, B)
    }

    fn throughput((L, n, _, _, _): Self::Params) -> u64 {
        (L * n) as u64
    }

    fn run_prover((mlin, A, M): Self::Input) -> Self::Output {
        let mut ts = create_transcript();
        mlin.mlin(&A, &M, &mut ts)
    }
}

/// Verifier benchmark for standard multilinear folding protocol.
///
/// Measures sumcheck verification time for multilinear folding proofs.
/// Tests how verification scales with L, demonstrating sublinear cost
/// relative to the prover's linear work in L.
struct MultilinearFoldVerifier;

impl VerifierBenchmark for MultilinearFoldVerifier {
    type Input = Mlin<R>;
    type Proof = latticefold_plus::cm::CmProof<R>;
    type Params = (usize, usize, usize, usize, usize);

    fn group_name() -> &'static str {
        "MultilinearFolding-Verifier"
    }

    fn setup_proof((L, n, k, kappa, B): Self::Params) -> (Self::Input, Self::Proof) {
        setup_proof(L, n, k, kappa, B)
    }

    fn param_label((L, n, k, kappa, B): Self::Params) -> String {
        format!("L={}_n={}_k={}_κ={}_B={}", L, n, k, kappa, B)
    }

    fn throughput((L, n, _, _, _): Self::Params) -> u64 {
        (L * n) as u64
    }

    fn run_verifier(_input: &Self::Input, proof: &Self::Proof) {
        let n = 65536;
        let M = create_test_m_matrix(n);
        let mut ts = create_transcript();
        proof.verify(&M, &mut ts).unwrap();
    }
}

/// Prover benchmark measuring decomposition width (k) scaling.
///
/// Tests how prover performance scales with k ∈ [2,3,4] while keeping
/// L=4, κ=2, and B=50 fixed. Note that witness size n is adjusted
/// proportionally with k to maintain parameter constraints.
struct MultilinearFoldKScaling;

impl ProverBenchmark for MultilinearFoldKScaling {
    type Input = (Mlin<R>, Matrix<R>, Vec<SparseMatrix<R>>);
    type Output = (LinB2<R>, latticefold_plus::cm::CmProof<R>);
    type Params = (usize, usize, usize, usize, usize);

    fn group_name() -> &'static str {
        "MultilinearFolding-KScaling"
    }

    fn setup_input((L, n, k, kappa, B): Self::Params) -> Self::Input {
        let mut rng = bench_rng();
        let mlin = setup_input(L, n, k, kappa, B);
        let M = create_test_m_matrix(n);
        let A = create_ajtai_matrix(kappa, n, &mut rng);
        (mlin, A, M)
    }

    fn param_label((_L, n, k, _kappa, _B): Self::Params) -> String {
        format!("k={}_n={}", k, n)
    }

    fn throughput((L, n, _, _, _): Self::Params) -> u64 {
        (L * n) as u64
    }

    fn run_prover((mlin, A, M): Self::Input) -> Self::Output {
        let mut ts = create_transcript();
        mlin.mlin(&A, &M, &mut ts)
    }
}

/// Prover benchmark measuring large witness scaling.
///
/// Tests performance on very large witnesses from 128K to 512K elements
/// while keeping L=4, k=2, κ=2, and B=50 fixed. Demonstrates scalability
/// to real-world proof sizes.
struct MultilinearFoldLargeWitness;

impl ProverBenchmark for MultilinearFoldLargeWitness {
    type Input = (Mlin<R>, Matrix<R>, Vec<SparseMatrix<R>>);
    type Output = (LinB2<R>, latticefold_plus::cm::CmProof<R>);
    type Params = (usize, usize, usize, usize, usize);

    fn group_name() -> &'static str {
        "MultilinearFolding-LargeWitness"
    }

    fn setup_input((L, n, k, kappa, B): Self::Params) -> Self::Input {
        let mut rng = bench_rng();
        let mlin = setup_input(L, n, k, kappa, B);
        let M = create_test_m_matrix(n);
        let A = create_ajtai_matrix(kappa, n, &mut rng);
        (mlin, A, M)
    }

    fn param_label((_L, n, _k, _kappa, _B): Self::Params) -> String {
        format!("n={}", n)
    }

    fn throughput((L, n, _, _, _): Self::Params) -> u64 {
        (L * n) as u64
    }

    fn run_prover((mlin, A, M): Self::Input) -> Self::Output {
        let mut ts = create_transcript();
        mlin.mlin(&A, &M, &mut ts)
    }
}

/// Prover benchmark measuring security parameter (κ) scaling.
///
/// Tests how prover performance scales with κ ∈ [2,3,4,5] while keeping
/// L=4, n=65536, k=2, and B=50 fixed. Higher κ values increase commitment
/// matrix dimensions and thus affect both commitment operations and proof size.
struct MultilinearFoldKappaScaling;

impl ProverBenchmark for MultilinearFoldKappaScaling {
    type Input = (Mlin<R>, Matrix<R>, Vec<SparseMatrix<R>>);
    type Output = (LinB2<R>, latticefold_plus::cm::CmProof<R>);
    type Params = (usize, usize, usize, usize, usize);

    fn group_name() -> &'static str {
        "MultilinearFolding-KappaScaling"
    }

    fn setup_input((L, n, k, kappa, B): Self::Params) -> Self::Input {
        let mut rng = bench_rng();
        let mlin = setup_input(L, n, k, kappa, B);
        let M = create_test_m_matrix(n);
        let A = create_ajtai_matrix(kappa, n, &mut rng);
        (mlin, A, M)
    }

    fn param_label((_L, _n, _k, kappa, _B): Self::Params) -> String {
        format!("{}", kappa)
    }

    fn throughput((L, n, _, _, _): Self::Params) -> u64 {
        (L * n) as u64
    }

    fn run_prover((mlin, A, M): Self::Input) -> Self::Output {
        let mut ts = create_transcript();
        mlin.mlin(&A, &M, &mut ts)
    }
}

// ============================================================================
// Benchmark Entry Points
// ============================================================================

/// Benchmark entry point for multilinear folding prover with L scaling.
fn bench_mlin_prover(c: &mut Criterion) {
    bench_prover_protocol::<MultilinearFoldProver>(c, multilinear_fold::FOLDING_ARITY);
}

/// Benchmark entry point for multilinear folding verifier with L scaling.
fn bench_mlin_verifier(c: &mut Criterion) {
    bench_verifier_protocol::<MultilinearFoldVerifier>(c, multilinear_fold::FOLDING_ARITY);
}

/// Benchmark entry point for multilinear folding prover with k scaling.
fn bench_mlin_k_scaling(c: &mut Criterion) {
    bench_prover_protocol::<MultilinearFoldKScaling>(c, multilinear_fold::K_SCALING);
}

/// Benchmark entry point for multilinear folding prover with large witness scaling.
fn bench_mlin_large_witness(c: &mut Criterion) {
    bench_prover_protocol::<MultilinearFoldLargeWitness>(c, multilinear_fold::LARGE_WITNESS);
}

/// Benchmark entry point for multilinear folding prover with κ scaling.
fn bench_mlin_kappa_scaling(c: &mut Criterion) {
    bench_prover_protocol::<MultilinearFoldKappaScaling>(c, multilinear_fold::KAPPA_SCALING);
}

criterion::criterion_group!(
    benches,
    bench_mlin_prover,
    bench_mlin_verifier,
    bench_mlin_k_scaling,
    bench_mlin_large_witness,
    bench_mlin_kappa_scaling,
);
criterion::criterion_main!(benches);
