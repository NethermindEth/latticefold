# LatticeFold+ Benchmarks

Performance benchmarks for all protocols in the LatticeFold+ folding scheme.

## Architecture

The benchmark suite uses a trait-based architecture to eliminate code duplication while maintaining type safety and clarity:

```
benches/
├── utils/
│   ├── helpers.rs          # Core benchmark infrastructure
│   └── mod.rs              # Parameter sets for all protocols
│
├── rgchk.rs                # Construction 4.3-4.4: Range check
├── setchk.rs               # Construction 4.2: Set check (monomial verification)
├── split.rs                # Construction 4.1: Gadget decomposition
├── double_commitment.rs    # Section 4.1: Double commitment (RgInstance creation)
├── cm.rs                   # Construction 4.5: Commitment transformation
├── decomp.rs               # Construction 5.3: Decomposition (LinB2 → 2×LinB)
├── lin.rs                  # Construction 5.1: Single instance folding (L=1)
├── mlin.rs                 # Construction 5.2: Multilinear folding (L>1)
└── e2e.rs                  # Section 5: End-to-end LatticeFold+ protocol
```

## Trait-Based Design

All benchmarks implement one of two core traits:

### ProverBenchmark

Defines how to set up inputs and execute the prover:

```rust
pub trait ProverBenchmark {
    type Input;
    type Output;
    type Params: Copy;

    fn group_name() -> &'static str;
    fn setup_input(params: Self::Params) -> Self::Input;
    fn param_label(params: Self::Params) -> String;
    fn throughput(params: Self::Params) -> u64;
    fn run_prover(input: Self::Input) -> Self::Output;
}
```

### VerifierBenchmark

Defines how to generate proofs and execute verification:

```rust
pub trait VerifierBenchmark {
    type Input;
    type Proof;
    type Params: Copy;

    fn group_name() -> &'static str;
    fn setup_proof(params: Self::Params) -> (Self::Input, Self::Proof);
    fn param_label(params: Self::Params) -> String;
    fn throughput(params: Self::Params) -> u64;
    fn run_verifier(input: &Self::Input, proof: &Self::Proof);
}
```

### Generic Runners

The generic functions `bench_prover_protocol<P>()` and `bench_verifier_protocol<V>()` handle all criterion orchestration:

- Benchmark group creation and configuration
- Parameter iteration
- Throughput measurement
- Timing setup with `iter_batched`

## Protocol Benchmarks

Each protocol benchmark file follows a consistent structure:

1. **Module documentation** - Protocol overview, paper references
2. **Setup functions** - `setup_input()` and `setup_proof()` helpers
3. **Trait implementations** - One struct per benchmark variant
4. **Entry point functions** - Thin wrappers calling generic runners
5. **Criterion registration** - `criterion_group!` and `criterion_main!`

### Range Check (Construction 4.3-4.4)

**File**: `rgchk.rs`
**Protocol**: Verifies witness coefficients lie within range [-B, B]
**Benchmarks**:
- `RangeCheckProver` - Witness scaling (32K-512K)
- `RangeCheckVerifier` - Verification scaling
- `RangeCheckKScaling` - Decomposition width k ∈ [2,3,4,5]
- `RangeCheckKappaScaling` - Security parameter κ ∈ [2,3,4,5]

### Set Check (Construction 4.2)

**File**: `setchk.rs`
**Protocol**: Verifies matrices contain monomials (one non-zero per row/column)
**Benchmarks**:
- `SetCheckProver` - Set size scaling (256, 512, 1024)
- `SetCheckVerifier` - Verification scaling
- `SetCheckBatching` - Batching efficiency (1, 2, 4, 8, 16 sets)

### Split Function (Construction 4.1)

**File**: `split.rs`
**Protocol**: Gadget decomposition on double commitments
**Benchmarks**:
- `SplitVaryingParams` - Combined parameter scaling
- `SplitScalingKFirst` - First decomposition width k_first ∈ [2,4,6,8]
- `SplitScalingKappa` - Security parameter κ ∈ [1,2,3,4]

### Double Commitment (Section 4.1-4.3)

**File**: `double_commitment.rs`
**Protocol**: Creates RgInstance structures from witness vectors using double commitment
**Benchmarks**:
- `DoubleCommitmentProver` - Witness scaling (32K, 64K, 128K)
- `DoubleCommitmentKScaling` - Decomposition width k ∈ [2,4]

### Commitment Transformation (Construction 4.5)

**File**: `cm.rs`
**Protocol**: Converts double commitments to folded commitments
**Benchmarks**:
- `CommitmentTransformProver` - Folding arity L ∈ [2,3,4,5,6,7,8]
- `CommitmentTransformVerifier` - Verification scaling

### Decomposition (Construction 5.3)

**File**: `decomp.rs`
**Protocol**: Splits LinB2 (norm B²) → 2×LinB (norm B each)
**Benchmarks**:
- `DecompositionProver` - Witness scaling (32K, 64K, 128K)
- `DecompositionVerifier` - Verification scaling
- `bench_decomp_roundtrip` - Complete fold→decompose cycle

### Single Instance Folding (Construction 5.1)

**File**: `lin.rs`
**Protocol**: Baseline folding (L=1): LinB → LinB2
**Benchmarks**:
- `SingleInstanceFoldProver` - Witness scaling (32K, 64K, 128K)
- `SingleInstanceFoldVerifier` - Verification scaling

### Multilinear Folding (Construction 5.2)

**File**: `mlin.rs`
**Protocol**: Main folding protocol: L×LinB → LinB2
**Benchmarks**:
- `MultilinearFoldProver` - Folding arity L ∈ [2,3,4,5,6,7,8]
- `MultilinearFoldVerifier` - Verification scaling
- `MultilinearFoldKScaling` - Decomposition width k ∈ [2,3,4]
- `MultilinearFoldLargeWitness` - Large witnesses (128K, 256K, 512K)
- `MultilinearFoldKappaScaling` - Security parameter κ ∈ [2,3,4,5]

### End-to-End Protocol (Section 5)

**File**: `e2e.rs`
**Protocol**: Complete LatticeFold+ stack (R1CS commitment → linearization → range check → commitment transformation → multilinear folding → verification)
**Benchmarks**:
- `E2EProver` - Protocol scaling (64K-128K witness, L ∈ [2,3])
- `E2EVerifier` - Verification scaling
- `E2EFoldingArity` - Folding arity L ∈ [2,3,4,5]

## Running Benchmarks

### Run all benchmarks:
```bash
cargo bench
```

### Run specific protocol:
```bash
cargo bench --bench rgchk
cargo bench --bench double_commitment
cargo bench --bench mlin
cargo bench --bench e2e
```

### Run with test mode (faster, for validation):
```bash
cargo bench -- --test
```

### Run specific benchmark within a file:
```bash
cargo bench --bench mlin -- "Prover"
cargo bench --bench rgchk -- "KScaling"
```

## Benchmark Configuration

All benchmarks use consistent criterion settings (defined in `helpers.rs::configure_benchmark_group`):

- **Sample size**: 10 iterations
- **Measurement time**: 10 seconds
- **Warm-up time**: 3 seconds
- **Batch size**: `SmallInput` (setup outside timing, execution inside)

## Parameter Sets

Parameter sets are defined in `utils/mod.rs` and follow naming conventions:

- `WITNESS_SCALING` - Varies witness/input size
- `K_SCALING` - Varies decomposition width k
- `KAPPA_SCALING` - Varies security parameter κ
- `FOLDING_ARITY` - Varies number of batched instances L

Each parameter set is documented with:
- What varies and the range
- What remains fixed
- Parameter tuple format


## Adding New Benchmarks

To add a new protocol benchmark:

1. Add parameter set to `utils/mod.rs`
2. Create `protocol_name.rs` with:
   ```rust
   // 1. Module docs
   //! Benchmarks for Construction X.Y: Protocol Name

   // 2. Setup functions
   fn setup_input(...) -> ProtocolInput { /* ... */ }
   fn setup_proof(...) -> (ProtocolInput, ProtocolProof) { /* ... */ }

   // 3. Trait implementation
   struct ProtocolProver;
   impl ProverBenchmark for ProtocolProver { /* ... */ }

   // 4. Entry point
   fn bench_protocol_prover(c: &mut Criterion) {
       bench_prover_protocol::<ProtocolProver>(c, params::PARAM_SET);
   }

   // 5. Criterion registration
   criterion_group!(benches, bench_protocol_prover);
   criterion_main!(benches);
   ```

## Implementation Notes

- **Deterministic RNG**: All benchmarks use `bench_rng()` with fixed seed (0x42424242) for reproducibility
- **Proof validation**: Verifier benchmarks validate generated proofs in setup to ensure correctness
- **Tuple inputs**: Complex protocols use tuple `Input` types (e.g., `(Mlin<R>, Matrix<R>, Vec<SparseMatrix<R>>)`) to handle multi-component setups
- **Zero-cost abstractions**: Trait-based design compiles to identical code as hand-written benchmarks
- **Witness patterns**: `WitnessPattern::BinaryChoice` introduced for binary (0/1) witness generation, used by `double_commitment.rs` and `e2e.rs`
- **Parameter constraints**: Some protocols enforce minimum witness sizes (e.g., e2e requires n ≥ 45056 due to tau_unpadded constraint)
- **Decomposition constraints**: Witness size must be divisible by decomposition width k

## Paper References

All constructions reference the LatticeFold+ paper:

- **Section 4.1**: Split function (gadget decomposition)
- **Section 4.1-4.3**: Double commitment (RgInstance creation)
- **Section 4.2**: Set check (monomial verification)
- **Section 4.3-4.4**: Range check
- **Section 4.5**: Commitment transformation
- **Section 5**: End-to-end protocol (complete LatticeFold+ stack)
- **Section 5.1**: Single instance folding
- **Section 5.2**: Multilinear folding (main protocol)
- **Section 5.3**: Decomposition protocol

## Benchmark Results

Results are saved to `target/criterion/` with:
- HTML reports in `target/criterion/report/index.html`
- CSV data for each benchmark
- Historical comparison data

View results:
```bash
open target/criterion/report/index.html
```
