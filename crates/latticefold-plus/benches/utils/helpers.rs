//! Shared utility functions and traits for LatticeFold+ benchmarks.
//!
//! This module provides common infrastructure used across all protocol benchmarks:
//! - Trait-based benchmark abstractions to eliminate code duplication
//! - Helper functions for deterministic setup and configuration
//! - Builders for common cryptographic structures

use std::time::Duration;

use ark_ff::{PrimeField, Zero};
use criterion::{BatchSize, BenchmarkGroup, BenchmarkId, Criterion, Throughput};
use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
use latticefold::arith::r1cs::R1CS;
use latticefold_plus::{rgchk::DecompParameters, transcript::PoseidonTranscript};
use rand::{CryptoRng, Rng};
use stark_rings::{cyclotomic_ring::models::frog_ring::RqPoly as R, PolyRing};
use stark_rings_linalg::{Matrix, SparseMatrix};

// ============================================================================
// Core Types
// ============================================================================

/// Decomposition parameters derived from ring dimension and modulus.
///
/// These parameters determine the base and length for gadget decomposition,
/// which is fundamental to double commitment schemes used throughout the protocols.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct DecompParams {
    /// Decomposition base: d/2 where d is the ring dimension
    pub b: u128,
    /// Decomposition length: ⌈log_b(q)⌉ where q is the base ring modulus
    pub l: usize,
}

impl DecompParams {
    /// Computes decomposition parameters for the current ring configuration.
    ///
    /// Uses FrogRing dimension and modulus to calculate optimal base and length
    /// for gadget decomposition in double commitment schemes.
    pub fn compute() -> Self {
        let d = R::dimension();
        let b = (d / 2) as u128;

        let modulus = <<R as PolyRing>::BaseRing as PrimeField>::MODULUS.0[0] as f64;
        let l = (modulus.ln() / (b as f64).ln()).ceil() as usize;

        Self { b, l }
    }

    /// Returns base multiplied by 2 for balanced decomposition contexts.
    #[inline]
    pub fn base_pow2(&self) -> u128 {
        self.b * 2
    }
}

/// Witness generation strategy for benchmark input construction.
///
/// Defines different patterns for creating witness vectors with specific
/// properties required by various protocols.
pub enum WitnessPattern {
    /// All ring elements set to the multiplicative identity.
    AllOnes,
    /// Small random coefficients in [0, 10) placed in the first coefficient position.
    /// Other coefficient positions are zero.
    SmallRandom,
    /// Binary witness: randomly chooses between 0 and 1 for each element.
    /// Commonly used in boolean circuits and binary constraint systems.
    BinaryChoice,
    /// Custom generation function for specialized witness patterns.
    Custom(Box<dyn Fn(usize) -> R>),
}

impl WitnessPattern {
    /// Generates a witness vector of the specified size following this pattern.
    ///
    /// # Arguments
    /// * `size` - Number of ring elements to generate
    /// * `rng` - Random number generator (used only for random patterns)
    pub fn generate(&self, size: usize, rng: &mut impl Rng) -> Vec<R> {
        use rand::seq::SliceRandom;
        use stark_rings::Ring;

        match self {
            WitnessPattern::AllOnes => vec![R::from(1u128); size],
            WitnessPattern::SmallRandom => (0..size)
                .map(|_| {
                    let mut r = R::zero();
                    r.coeffs_mut()[0] = ((rng.gen::<u32>() % 10) as u64).into();
                    r
                })
                .collect(),
            WitnessPattern::BinaryChoice => {
                let choices = [R::ZERO, R::ONE];
                (0..size).map(|_| *choices.choose(rng).unwrap()).collect()
            }
            WitnessPattern::Custom(f) => (0..size).map(|i| f(i)).collect(),
        }
    }
}

/// Fluent builder for R1CS instances with common benchmark configurations.
///
/// Provides convenient construction of R1CS constraints with identity matrices
/// and optional modifications commonly used in benchmark protocols.
pub struct R1CSBuilder {
    n: usize,
    k: usize,
    b: u128,
    modify_first: bool,
}

impl R1CSBuilder {
    /// Creates a new R1CS builder with the given parameters.
    ///
    /// # Arguments
    /// * `n` - Witness size (total number of variables)
    /// * `k` - Decomposition width
    /// * `b` - Norm bound parameter for constraint satisfaction
    pub fn new(n: usize, k: usize, b: u128) -> Self {
        Self {
            n,
            k,
            b,
            modify_first: true,
        }
    }

    /// Controls whether to modify the first coefficient of A and C matrices.
    ///
    /// When enabled (default), sets `A[0][0] = 2` and `C[0][0] = 2` to create
    /// non-trivial constraints for testing.
    pub fn with_first_element_modification(mut self, modify: bool) -> Self {
        self.modify_first = modify;
        self
    }

    /// Builds an R1CS instance with decomposed square transformation.
    ///
    /// Creates identity matrices for A, B, C, applies `r1cs_decomposed_square`
    /// transformation, and optionally modifies first elements for non-trivial constraints.
    pub fn build_decomposed_square(self) -> R1CS<R> {
        use latticefold_plus::r1cs::r1cs_decomposed_square;

        let mut r1cs = r1cs_decomposed_square(
            R1CS::<R> {
                l: 1,
                A: SparseMatrix::identity(self.n / self.k),
                B: SparseMatrix::identity(self.n / self.k),
                C: SparseMatrix::identity(self.n / self.k),
            },
            self.n,
            self.b,
            self.k,
        );

        if self.modify_first {
            r1cs.A.coeffs[0][0].0 = 2u128.into();
            r1cs.C.coeffs[0][0].0 = 2u128.into();
        }

        r1cs
    }

    /// Builds a basic R1CS instance without decomposed square transformation.
    ///
    /// Creates identity matrices for A, B, C with optional first element modification.
    pub fn build_basic(self) -> R1CS<R> {
        let mut r1cs = R1CS::<R> {
            l: 1,
            A: SparseMatrix::identity(self.n / self.k),
            B: SparseMatrix::identity(self.n / self.k),
            C: SparseMatrix::identity(self.n / self.k),
        };

        if self.modify_first {
            r1cs.A.coeffs[0][0].0 = 2u128.into();
            r1cs.C.coeffs[0][0].0 = 2u128.into();
        }

        r1cs
    }
}

// ============================================================================
// Benchmark Traits
// ============================================================================

/// Trait defining the structure of a prover benchmark for a cryptographic protocol.
///
/// Implementors specify how to set up inputs, format parameter labels, calculate
/// throughput, and execute the prover. The generic `bench_prover_protocol` function
/// handles all criterion orchestration (group creation, iteration, timing).
///
/// # Type Parameters
/// - `Input`: The protocol's input type (may be a tuple for complex setups)
/// - `Output`: The protocol's output/proof type
/// - `Params`: Parameter tuple type (must be `Copy` for iteration)
pub trait ProverBenchmark {
    /// Input type for the prover
    type Input;

    /// Output/proof type from the prover
    type Output;

    /// Parameter tuple type (must be Copy for iteration)
    type Params: Copy;

    /// Returns the benchmark group name displayed in criterion output.
    fn group_name() -> &'static str;

    /// Creates protocol input from the given parameters.
    ///
    /// This is called inside `iter_batched` setup phase and should perform
    /// all necessary initialization for a single benchmark iteration.
    fn setup_input(params: Self::Params) -> Self::Input;

    /// Formats parameters into a human-readable label for criterion reports.
    ///
    /// Should include all relevant parameter values in a compact format,
    /// e.g., "w=65536_k=2_κ=2".
    fn param_label(params: Self::Params) -> String;

    /// Calculates throughput measurement for the given parameters.
    ///
    /// Return value represents the number of "elements" processed, which
    /// criterion uses to compute elements/second metrics.
    fn throughput(params: Self::Params) -> u64;

    /// Executes the prover with the given input.
    ///
    /// This is the timed portion of the benchmark. Should perform the actual
    /// cryptographic operation without any setup or validation overhead.
    fn run_prover(input: Self::Input) -> Self::Output;
}

/// Trait defining the structure of a verifier benchmark for a cryptographic protocol.
///
/// Similar to `ProverBenchmark` but for verification operations. Implementors
/// specify how to generate valid proofs and execute verification.
///
/// # Type Parameters
/// - `Input`: The protocol's input type
/// - `Proof`: The proof type to be verified
/// - `Params`: Parameter tuple type (must be `Copy` for iteration)
pub trait VerifierBenchmark {
    /// Input type (often same as prover input)
    type Input;

    /// Proof type to be verified
    type Proof;

    /// Parameter tuple type (must be Copy for iteration)
    type Params: Copy;

    /// Returns the benchmark group name displayed in criterion output.
    fn group_name() -> &'static str;

    /// Generates a valid proof for benchmarking the verifier.
    ///
    /// This is called once per parameter set before timing begins. The generated
    /// proof should be valid to avoid measuring error handling overhead.
    ///
    /// Returns (input, proof) tuple where input may be needed for verification context.
    fn setup_proof(params: Self::Params) -> (Self::Input, Self::Proof);

    /// Formats parameters into a human-readable label for criterion reports.
    fn param_label(params: Self::Params) -> String;

    /// Calculates throughput measurement for the given parameters.
    fn throughput(params: Self::Params) -> u64;

    /// Executes verification with the given proof.
    ///
    /// This is the timed portion of the benchmark. Should perform only the
    /// verification operation without proof generation overhead.
    fn run_verifier(input: &Self::Input, proof: &Self::Proof);
}

// ============================================================================
// Generic Benchmark Runners
// ============================================================================

/// Generic benchmark runner for prover protocols.
///
/// Eliminates boilerplate by handling benchmark group creation, parameter
/// iteration, throughput measurement, and timing orchestration. The protocol-specific
/// behavior is defined via the `ProverBenchmark` trait implementation.
///
/// # Arguments
/// * `c` - Criterion instance
/// * `param_sets` - Slice of parameter tuples to benchmark
///
/// # Type Parameters
/// * `P` - Protocol implementing `ProverBenchmark`
///
/// # Examples
/// ```rust,ignore
/// struct MyProver;
/// impl ProverBenchmark for MyProver { /* ... */ }
///
/// fn bench_my_prover(c: &mut Criterion) {
///     bench_prover_protocol::<MyProver>(c, MY_PARAMS);
/// }
/// ```
pub fn bench_prover_protocol<P: ProverBenchmark>(c: &mut Criterion, param_sets: &[P::Params]) {
    let mut group = c.benchmark_group(P::group_name());
    configure_benchmark_group(&mut group);

    for &params in param_sets {
        group.throughput(Throughput::Elements(P::throughput(params)));

        group.bench_with_input(
            BenchmarkId::from_parameter(&P::param_label(params)),
            &params,
            |bencher, &params| {
                bencher.iter_batched(
                    || P::setup_input(params),
                    P::run_prover,
                    BatchSize::SmallInput,
                );
            },
        );
    }

    group.finish();
}

/// Generic benchmark runner for verifier protocols.
///
/// Similar to `bench_prover_protocol` but for verification. Handles proof generation
/// in the setup phase so only verification time is measured.
///
/// # Arguments
/// * `c` - Criterion instance
/// * `param_sets` - Slice of parameter tuples to benchmark
///
/// # Type Parameters
/// * `V` - Protocol implementing `VerifierBenchmark`
pub fn bench_verifier_protocol<V: VerifierBenchmark>(c: &mut Criterion, param_sets: &[V::Params]) {
    let mut group = c.benchmark_group(V::group_name());
    configure_benchmark_group(&mut group);

    for &params in param_sets {
        group.throughput(Throughput::Elements(V::throughput(params)));

        group.bench_with_input(
            BenchmarkId::from_parameter(&V::param_label(params)),
            &params,
            |bencher, &params| {
                let (input, proof) = V::setup_proof(params);
                bencher.iter(|| V::run_verifier(&input, &proof));
            },
        );
    }

    group.finish();
}

// ============================================================================
// Helper Functions
// ============================================================================

/// Creates a deterministic RNG for reproducible benchmarks.
///
/// Uses a fixed seed (0x42424242) to ensure benchmark results are comparable
/// across runs and different machines. This RNG implements both `Rng` and
/// `CryptoRng` traits for use in cryptographic setup.
#[inline]
pub fn bench_rng() -> impl Rng + CryptoRng {
    use rand::SeedableRng;
    rand::rngs::StdRng::seed_from_u64(0x42424242)
}

/// Creates a modified identity matrix for protocol testing.
///
/// Generates an n×n identity matrix with the first diagonal entry set to 2
/// instead of 1, creating non-trivial constraints for R1CS benchmarks.
///
/// Returns a single-element vector containing the modified sparse matrix.
pub fn create_test_m_matrix(n: usize) -> Vec<SparseMatrix<R>> {
    let mut m = SparseMatrix::identity(n);
    m.coeffs[0][0].0 = 2u128.into();
    vec![m]
}

/// Creates a random Ajtai commitment matrix.
///
/// Generates a kappa × n matrix with random ring elements for use in
/// lattice-based commitment schemes.
///
/// # Arguments
/// * `kappa` - Number of rows (security parameter)
/// * `n` - Number of columns (witness size)
/// * `rng` - Random number generator
pub fn create_ajtai_matrix(kappa: usize, n: usize, rng: &mut impl Rng) -> Matrix<R> {
    Matrix::<R>::rand(rng, kappa, n)
}

/// Creates an empty Fiat-Shamir transcript for proof generation.
///
/// Returns a Poseidon transcript configured for the FrogRing, used to make
/// protocols non-interactive via the Fiat-Shamir transform.
#[inline]
pub fn create_transcript() -> PoseidonTranscript<R> {
    PoseidonTranscript::empty::<PC>()
}

/// Configures a criterion benchmark group with standard settings.
///
/// Applies consistent parameters across all protocol benchmarks:
/// - Sample size: 10 iterations per parameter set
/// - Measurement time: 10 seconds total
/// - Warm-up time: 3 seconds before measurements begin
pub fn configure_benchmark_group(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(10));
    group.warm_up_time(Duration::from_secs(3));
}

/// Validates that witness size satisfies decomposition parameter constraints.
///
/// Ensures the fundamental constraint: `witness_size >= κ × k × d × l × d`
/// where κ (kappa) is the security parameter, k is decomposition width,
/// d is ring dimension, and l is decomposition length.
///
/// # Panics
/// Panics with a detailed error message if the constraint is violated.
pub fn validate_witness_params(witness_size: usize, k: usize, kappa: usize) {
    let decomp = DecompParams::compute();
    let d = R::dimension();

    let min_witness_size = kappa * k * d * decomp.l * d;
    assert!(
        witness_size >= min_witness_size,
        "Invalid parameters: witness_size ({}) must be >= κ × k × d × l × d = {} × {} × {} × {} × {} = {}",
        witness_size, kappa, k, d, decomp.l, d, min_witness_size
    );
}

/// Computes and validates decomposition parameters for the given configuration.
///
/// Combines parameter validation with decomposition parameter computation,
/// returning `DecompParameters` ready for use in protocol construction.
///
/// # Arguments
/// * `k` - Decomposition width
/// * `kappa` - Security parameter (number of commitment rows)
/// * `witness_size` - Size of the witness vector
///
/// # Panics
/// Panics if witness size constraint is violated (see `validate_witness_params`).
pub fn get_validated_decomp_params(
    k: usize,
    kappa: usize,
    witness_size: usize,
) -> DecompParameters {
    validate_witness_params(witness_size, k, kappa);
    let decomp = DecompParams::compute();
    DecompParameters {
        b: decomp.b,
        k,
        l: decomp.l,
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bench_rng_deterministic() {
        let mut rng1 = bench_rng();
        let mut rng2 = bench_rng();

        let val1: u64 = rng1.gen();
        let val2: u64 = rng2.gen();

        assert_eq!(val1, val2, "RNG should be deterministic");
    }

    #[test]
    fn test_decomp_params_compute() {
        let params = DecompParams::compute();

        let d = R::dimension();
        assert_eq!(params.b, (d / 2) as u128);
        assert!(params.l > 0, "Decomposition width should be positive");
    }

    #[test]
    fn test_decomp_params_base_pow2() {
        let params = DecompParams::compute();
        assert_eq!(params.base_pow2(), params.b * 2);
    }

    #[test]
    fn test_create_test_m_matrix() {
        let n = 16;
        let matrices = create_test_m_matrix(n);

        assert_eq!(matrices.len(), 1);
        assert_eq!(matrices[0].nrows, n);
        assert_eq!(matrices[0].ncols, n);

        let expected_val: <R as PolyRing>::ElemTy = 2u128.into();
        assert_eq!(matrices[0].coeffs[0][0].0, expected_val);
    }

    #[test]
    fn test_create_ajtai_matrix() {
        let mut rng = bench_rng();
        let kappa = 4;
        let n = 256;

        let matrix = create_ajtai_matrix(kappa, n, &mut rng);

        assert_eq!(matrix.nrows, kappa);
        assert_eq!(matrix.ncols, n);
    }

    #[test]
    fn test_create_transcript() {
        let _ts = create_transcript();
    }
}
