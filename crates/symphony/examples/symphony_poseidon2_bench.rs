//! Symphony Poseidon2 Benchmark — SP1 Verifier Component GM-1B
//!
//! This implements the SP1 Poseidon2 hash function as a Hadamard relation (M1·f ⊙ M2·f = M3·f)
//! and folds it using Symphony's Π_fold to measure performance.
//!
//! SP1 Poseidon2 parameters:
//! - WIDTH = 16
//! - NUM_EXTERNAL_ROUNDS = 8 (4 at start, 4 at end)
//! - NUM_INTERNAL_ROUNDS = 13
//! - S-box: x^7 (computed as x² → x³ → x⁶ → x⁷, 4 mult constraints per S-box)
//!
//! Run:
//! - `cargo run -p symphony --example symphony_poseidon2_bench --release

#![allow(non_local_definitions)]

use std::sync::Arc;
use std::time::Instant;
use ark_ff::{Field, Fp384, MontBackend, MontConfig, PrimeField};
use rayon::ThreadPoolBuilder;
use dpp::BoundedFlpcpSparse;

fn main() {
    use std::panic::{catch_unwind, AssertUnwindSafe};
    use stark_rings::PolyRing;
    use cyclotomic_rings::rings::FrogPoseidonConfig;
    use stark_rings::cyclotomic_ring::models::frog_ring;

    struct PanicHookGuard(Option<Box<dyn Fn(&std::panic::PanicHookInfo<'_>) + Sync + Send + 'static>>);
    impl PanicHookGuard {
        fn silence_expected_panics() -> Self {
            let prev = std::panic::take_hook();
            std::panic::set_hook(Box::new(|_| {
                // Intentionally silent: the benchmark uses `catch_unwind` for expected failures
                // (e.g. too-small k_g in balanced decomposition).
            }));
            Self(Some(prev))
        }
    }
    impl Drop for PanicHookGuard {
        fn drop(&mut self) {
            if let Some(prev) = self.0.take() {
                std::panic::set_hook(prev);
            }
        }
    }

    fn run_for_ring<
        R: stark_rings::CoeffRing,
        PC: cyclotomic_rings::rings::GetPoseidonParams<<<R>::BaseRing as ark_ff::Field>::BasePrimeField>,
    >(
        ring_name: &str,
    ) where
        R::BaseRing: stark_rings::Zq + stark_rings::balanced_decomposition::Decompose + ark_ff::PrimeField,
    {
        const WIDTH: usize = 16;
        const NUM_EXTERNAL_ROUNDS: usize = 8;
        const NUM_INTERNAL_ROUNDS: usize = 13;

        println!("=========================================================");
        println!("Symphony Poseidon2 Benchmark — {ring_name}");
        println!("=========================================================");
        println!();
        println!("SP1 Poseidon2 parameters:");
        println!("  WIDTH = {WIDTH}");
        println!("  NUM_EXTERNAL_ROUNDS = {NUM_EXTERNAL_ROUNDS} (4 + 4)");
        println!("  NUM_INTERNAL_ROUNDS = {NUM_INTERNAL_ROUNDS}");
        println!("  S-box = x^7 (4 mult constraints per S-box)");
        println!();

        println!("Expected constraints per permutation:");
        println!("  External: 8 rounds × 16 S-boxes × 4 mults = 512");
        println!("  Internal: 13 rounds × 1 S-box × 4 mults = 52");
        println!("  Total S-box: ~564 constraints");
        println!();

        // Digit-sizing hint based on modulus size.
        let d_prime: u32 = (R::dimension() as u32) / 2;
        let bits_per_digit = (d_prime as f64).log2();
        let mod_bits = <<R as PolyRing>::BaseRing as ark_ff::Field>::BasePrimeField::MODULUS_BIT_SIZE;
        let k_suggest = ((mod_bits as f64) / bits_per_digit).ceil() as usize;
        println!(
            "Digit-sizing hint: dim={} => d'={d_prime} (~{bits_per_digit:.2} bits/digit), modulus bits={mod_bits} => suggest k_g≈{k_suggest}",
            R::dimension()
        );

        // Sweep permutation counts.
        //
        // IMPORTANT:
        // This benchmark stresses witness magnitudes; for Frog (~64-bit modulus) and small ring dim,
        // a tiny k_g (like 3) will often overflow balanced decomposition. We therefore default to
        // the computed `k_suggest` unless overridden via `K_G`.
        let k_list = if let Ok(s) = std::env::var("K_G") {
            vec![s.parse::<usize>().expect("K_G must be a usize")]
        } else {
            vec![k_suggest.max(3)]
        };
        let perms_list: Vec<usize> = vec![1, 2, 4, 8];
        for &k_g in &k_list {
            for &num_perms in &perms_list {
                // Silence panic-hook output only for the `catch_unwind` region so expected failures
                // don't print scary backtraces.
                let _hook_guard = PanicHookGuard::silence_expected_panics();
                let r = catch_unwind(AssertUnwindSafe(|| {
                    test_symphony_poseidon2_folding::<R, PC>(num_perms, k_g);
                }));
                if let Err(e) = r {
                    // Print panic payload to distinguish decomposition overflow from other issues.
                    let msg = if let Some(s) = e.downcast_ref::<&str>() {
                        *s
                    } else if let Some(s) = e.downcast_ref::<String>() {
                        s.as_str()
                    } else {
                        "<non-string panic payload>"
                    };
                    println!("  !! PANIC for perms={num_perms}, k_g={k_g}: {msg}");
                }
            }
        }
    }

    // NOTE: Frog is ~64-bit modulus and is very expensive once k_g is raised high enough.
    // For B1 we focus on BabyBear modulus first.
    // Default to Frog (matches the paper's instantiation modulus family and has d=16 in this repo).
    // Note: the built-in BabyBear ring in this repo has d=72, which is incompatible with our Π_fold
    // MLE domain requirement that m*d is a power-of-two.
    run_for_ring::<frog_ring::RqPoly, FrogPoseidonConfig>("Frog ring (64-bit prime, d=16)");

    // Constraint count per permutation:
    // - 8 external rounds × 16 S-boxes × 4 mults = 512 constraints
    // - 13 internal rounds × 1 S-box × 4 mults = 52 constraints
    // Total: ~564 mult constraints per permutation

    // All printing is done inside `run_for_ring`.
}

/// Helper alias: Poseidon's base prime field for a given ring `R`.
type BF<R> = <<R as stark_rings::PolyRing>::BaseRing as Field>::BasePrimeField;

// Large prime field for the large-field embedding/packing path (Rev2).
//
// IMPORTANT:
// With Frog (~64-bit modulus) and multi-million-length witness vectors, the Rev2 packing step
// (Construction 5.21) can require a modulus larger than 256 bits under the conservative bounds.
// A 256-bit field (like secp256k1) can therefore fail with `PackingError::ModulusTooSmall`.
//
// We use the NIST P-384 prime here:
//   p = 2^384 - 2^128 - 2^96 + 2^32 - 1
#[derive(MontConfig)]
#[modulus = "39402006196394479212279040100143613805079739270465446667948293404245721771496870329047266088258938001861606973112319"]
#[generator = "2"]
pub struct Secp384r1Config;
type FLarge = Fp384<MontBackend<Secp384r1Config, 6>>;

/// Build Poseidon2 Hadamard relation matrices (M1, M2, M3) for Symphony.
///
/// For the Hadamard relation M1·f ⊙ M2·f = M3·f:
/// - Each S-box x^7 requires 4 mult constraints
/// - Witness f contains: [state_0..state_{WIDTH-1}, aux_vars...]
fn build_poseidon2_hadamard_matrices<R: stark_rings::PolyRing>(
    num_permutations: usize,
) -> (
    stark_rings_linalg::SparseMatrix<R>,
    stark_rings_linalg::SparseMatrix<R>,
    stark_rings_linalg::SparseMatrix<R>,
    Vec<R>,
    usize,
) {
    use stark_rings_linalg::SparseMatrix;

    const WIDTH: usize = 16;
    const NUM_EXTERNAL_ROUNDS: usize = 8;
    const NUM_INTERNAL_ROUNDS: usize = 13;

    // Count constraints and variables first
    let sboxes_per_perm = NUM_EXTERNAL_ROUNDS * WIDTH + NUM_INTERNAL_ROUNDS;
    let constraints_per_perm = sboxes_per_perm * 4; // 4 mult constraints per S-box
    let total_constraints = constraints_per_perm * num_permutations;

    // Variables: initial state + intermediate values
    // Each S-box creates 5 values: rc_add + x² + x³ + x⁶ + x⁷
    // External round: 16 sboxes × 5 + 16 state update = 96 per round
    // Internal round: 1 sbox × 5 + 16 state update = 21 per round
    let aux_per_round_external = WIDTH * 5 + WIDTH; // 5 per sbox + state update
    let aux_per_round_internal = 5 + WIDTH; // 1 sbox + state update
    let aux_per_perm = NUM_EXTERNAL_ROUNDS * aux_per_round_external / 2 * 2
        + NUM_INTERNAL_ROUNDS * aux_per_round_internal;

    let total_aux = aux_per_perm * num_permutations;
    let n = WIDTH + total_aux; // witness length

    // Round to power of 2 for folding (with extra headroom)
    let n_padded = (n * 2).next_power_of_two();
    
    // m must be >= m_J where m_J = (n / l_h) * lambda_pj
    // With l_h=64, lambda_pj=32: m_J = n_padded / 64 * 32 = n_padded / 2
    // So ensure m >= n_padded / 2 and m is a multiple of m_J
    let m_j_estimate = (n_padded / 64) * 32;
    let m_min = total_constraints.max(m_j_estimate);
    let m = m_min.next_power_of_two();

    // Build sparse matrices - format is (value, column_index) per row
    let mut m1_coeffs: Vec<Vec<(R, usize)>> = vec![Vec::new(); m];
    let mut m2_coeffs: Vec<Vec<(R, usize)>> = vec![Vec::new(); m];
    let mut m3_coeffs: Vec<Vec<(R, usize)>> = vec![Vec::new(); m];

    // Build witness
    let mut witness = vec![R::ZERO; n_padded];

    // Initialize state.
    //
    // IMPORTANT:
    // This benchmark only constrains the Poseidon2 S-box *multiplication* structure (x^7),
    // not the full set of linear/additive constraints of a real Poseidon2 permutation.
    //
    // We still use it to stress-test the *magnitude* of witness values for Π_rg decomposition
    // by actually computing x^7 chains (this is what broke and forced `ONE`-witnesses).
    let mut state_indices: Vec<usize> = (0..WIDTH).collect();
    let mut state_vals: Vec<R> = (0..WIDTH).map(|_| R::from(5u128)).collect(); // 5^7 = 78125 > 4096
    for i in 0..WIDTH {
        witness[i] = state_vals[i];
    }

    let mut next_aux = WIDTH;
    let mut constraint_idx = 0;
    let one = R::ONE;

    for _perm in 0..num_permutations {
        // === External rounds (first half) ===
        for _round in 0..NUM_EXTERNAL_ROUNDS / 2 {
            let mut new_state_indices = Vec::with_capacity(WIDTH);
            let mut new_state_vals = Vec::with_capacity(WIDTH);

            for i in 0..WIDTH {
                let with_rc_idx = next_aux;
                let with_rc_val = state_vals[i] + R::ONE;
                witness[with_rc_idx] = with_rc_val;
                next_aux += 1;

                let t1_idx = next_aux;
                let t1_val = with_rc_val * with_rc_val;
                witness[t1_idx] = t1_val;
                next_aux += 1;
                m1_coeffs[constraint_idx].push((one, with_rc_idx));
                m2_coeffs[constraint_idx].push((one, with_rc_idx));
                m3_coeffs[constraint_idx].push((one, t1_idx));
                constraint_idx += 1;

                let t2_idx = next_aux;
                let t2_val = t1_val * with_rc_val;
                witness[t2_idx] = t2_val;
                next_aux += 1;
                m1_coeffs[constraint_idx].push((one, t1_idx));
                m2_coeffs[constraint_idx].push((one, with_rc_idx));
                m3_coeffs[constraint_idx].push((one, t2_idx));
                constraint_idx += 1;

                let t3_idx = next_aux;
                let t3_val = t2_val * t2_val;
                witness[t3_idx] = t3_val;
                next_aux += 1;
                m1_coeffs[constraint_idx].push((one, t2_idx));
                m2_coeffs[constraint_idx].push((one, t2_idx));
                m3_coeffs[constraint_idx].push((one, t3_idx));
                constraint_idx += 1;

                let t4_idx = next_aux;
                let t4_val = t3_val * with_rc_val;
                witness[t4_idx] = t4_val;
                next_aux += 1;
                m1_coeffs[constraint_idx].push((one, t3_idx));
                m2_coeffs[constraint_idx].push((one, with_rc_idx));
                m3_coeffs[constraint_idx].push((one, t4_idx));
                constraint_idx += 1;

                new_state_indices.push(t4_idx);
                new_state_vals.push(t4_val);
            }

            // Linear layer placeholder: sum and add (like the R1CS toy circuit).
            // Not constrained here, but keeps subsequent S-box inputs growing realistically.
            let sum_val = new_state_vals.iter().copied().fold(R::ZERO, |a, b| a + b);
            state_indices.clear();
            state_vals.clear();
            for i in 0..WIDTH {
                let new_idx = next_aux;
                let new_val = new_state_vals[i] + sum_val;
                witness[new_idx] = new_val;
                next_aux += 1;
                state_indices.push(new_idx);
                state_vals.push(new_val);
            }
        }

        // === Internal rounds ===
        for _round in 0..NUM_INTERNAL_ROUNDS {
            // Round constant (+1)
            let s0_with_rc_idx = next_aux;
            let s0_with_rc_val = state_vals[0] + R::ONE;
            witness[s0_with_rc_idx] = s0_with_rc_val;
            next_aux += 1;

            // S-box for first element (4 constraints)
            let t1_idx = next_aux;
            let t1_val = s0_with_rc_val * s0_with_rc_val;
            witness[t1_idx] = t1_val;
            next_aux += 1;
            m1_coeffs[constraint_idx].push((one, s0_with_rc_idx));
            m2_coeffs[constraint_idx].push((one, s0_with_rc_idx));
            m3_coeffs[constraint_idx].push((one, t1_idx));
            constraint_idx += 1;

            let t2_idx = next_aux;
            let t2_val = t1_val * s0_with_rc_val;
            witness[t2_idx] = t2_val;
            next_aux += 1;
            m1_coeffs[constraint_idx].push((one, t1_idx));
            m2_coeffs[constraint_idx].push((one, s0_with_rc_idx));
            m3_coeffs[constraint_idx].push((one, t2_idx));
            constraint_idx += 1;

            let t3_idx = next_aux;
            let t3_val = t2_val * t2_val;
            witness[t3_idx] = t3_val;
            next_aux += 1;
            m1_coeffs[constraint_idx].push((one, t2_idx));
            m2_coeffs[constraint_idx].push((one, t2_idx));
            m3_coeffs[constraint_idx].push((one, t3_idx));
            constraint_idx += 1;

            let t4_idx = next_aux;
            let t4_val = t3_val * s0_with_rc_val;
            witness[t4_idx] = t4_val;
            next_aux += 1;
            m1_coeffs[constraint_idx].push((one, t3_idx));
            m2_coeffs[constraint_idx].push((one, s0_with_rc_idx));
            m3_coeffs[constraint_idx].push((one, t4_idx));
            constraint_idx += 1;

            // State update (all bounded)
            state_indices[0] = t4_idx;
            state_vals[0] = t4_val;
            let sum_val = state_vals.iter().copied().fold(R::ZERO, |a, b| a + b);
            for i in 0..WIDTH {
                let new_idx = next_aux;
                let new_val = state_vals[i] + sum_val;
                witness[new_idx] = new_val;
                next_aux += 1;
                state_indices[i] = new_idx;
                state_vals[i] = new_val;
            }
        }

        // === External rounds (second half) ===
        for _round in NUM_EXTERNAL_ROUNDS / 2..NUM_EXTERNAL_ROUNDS {
            let mut new_state_indices = Vec::with_capacity(WIDTH);
            let mut new_state_vals = Vec::with_capacity(WIDTH);

            for i in 0..WIDTH {
                let with_rc_idx = next_aux;
                let with_rc_val = state_vals[i] + R::ONE;
                witness[with_rc_idx] = with_rc_val;
                next_aux += 1;

                // S-box (4 constraints)
                let t1_idx = next_aux;
                let t1_val = with_rc_val * with_rc_val;
                witness[t1_idx] = t1_val;
                next_aux += 1;
                m1_coeffs[constraint_idx].push((one, with_rc_idx));
                m2_coeffs[constraint_idx].push((one, with_rc_idx));
                m3_coeffs[constraint_idx].push((one, t1_idx));
                constraint_idx += 1;

                let t2_idx = next_aux;
                let t2_val = t1_val * with_rc_val;
                witness[t2_idx] = t2_val;
                next_aux += 1;
                m1_coeffs[constraint_idx].push((one, t1_idx));
                m2_coeffs[constraint_idx].push((one, with_rc_idx));
                m3_coeffs[constraint_idx].push((one, t2_idx));
                constraint_idx += 1;

                let t3_idx = next_aux;
                let t3_val = t2_val * t2_val;
                witness[t3_idx] = t3_val;
                next_aux += 1;
                m1_coeffs[constraint_idx].push((one, t2_idx));
                m2_coeffs[constraint_idx].push((one, t2_idx));
                m3_coeffs[constraint_idx].push((one, t3_idx));
                constraint_idx += 1;

                let t4_idx = next_aux;
                let t4_val = t3_val * with_rc_val;
                witness[t4_idx] = t4_val;
                next_aux += 1;
                m1_coeffs[constraint_idx].push((one, t3_idx));
                m2_coeffs[constraint_idx].push((one, with_rc_idx));
                m3_coeffs[constraint_idx].push((one, t4_idx));
                constraint_idx += 1;

                new_state_indices.push(t4_idx);
                new_state_vals.push(t4_val);
            }

            // Linear layer
            let sum_val = new_state_vals.iter().copied().fold(R::ZERO, |a, b| a + b);
            state_indices.clear();
            state_vals.clear();
            for i in 0..WIDTH {
                let new_idx = next_aux;
                let new_val = new_state_vals[i] + sum_val;
                witness[new_idx] = new_val;
                next_aux += 1;
                state_indices.push(new_idx);
                state_vals.push(new_val);
            }
        }
    }

    // Construct SparseMatrix directly
    let m1 = SparseMatrix {
        nrows: m,
        ncols: n_padded,
        coeffs: m1_coeffs,
    };
    let m2 = SparseMatrix {
        nrows: m,
        ncols: n_padded,
        coeffs: m2_coeffs,
    };
    let m3 = SparseMatrix {
        nrows: m,
        ncols: n_padded,
        coeffs: m3_coeffs,
    };

    (m1, m2, m3, witness, constraint_idx)
}

fn r1cs_decompose_witness_and_expand_matrices<R: stark_rings::CoeffRing>(
    m1: &stark_rings_linalg::SparseMatrix<R>,
    m2: &stark_rings_linalg::SparseMatrix<R>,
    m3: &stark_rings_linalg::SparseMatrix<R>,
    witness: &[R],
    k_cs: usize,
    base_b: u128,
) -> (
    stark_rings_linalg::SparseMatrix<R>,
    stark_rings_linalg::SparseMatrix<R>,
    stark_rings_linalg::SparseMatrix<R>,
    Vec<R>,
)
where
    R::BaseRing: stark_rings::Zq + ark_ff::PrimeField,
{
    use ark_ff::Field as _;
    use ark_ff::PrimeField;
    use stark_rings_linalg::SparseMatrix;

    let d = R::dimension();
    let n = witness.len();

    // Read modulus into u128 (works for p <= 2^128).
    let modulus_big = <R::BaseRing as PrimeField>::MODULUS;
    let mut q: u128 = 0;
    for (i, limb) in modulus_big.as_ref().iter().enumerate() {
        q |= (*limb as u128) << (64 * i);
    }
    let q_half = q / 2;
    let b = base_b as i128;
    let b_half = (base_b / 2) as i128;

    // Decompose each ring element coefficient-wise into k_cs digits base b.
    // New witness layout: for original index j, digit t is at (j*k_cs + t).
    let mut out = vec![R::ZERO; n * k_cs];
    for (j, r) in witness.iter().enumerate() {
        let coeffs = r.coeffs();
        debug_assert_eq!(coeffs.len(), d);

        // digits[t][k] is the digit at position t of coefficient k.
        let mut digits: Vec<Vec<R::BaseRing>> = vec![vec![R::BaseRing::ZERO; d]; k_cs];
        for (k, c) in coeffs.iter().enumerate() {
            // Convert coefficient to centered i128 in (-q/2, q/2].
            let big = c.into_bigint();
            let mut u: u128 = 0;
            for (i, limb) in big.as_ref().iter().enumerate() {
                u |= (*limb as u128) << (64 * i);
            }
            let mut x: i128 = if u > q_half { u as i128 - q as i128 } else { u as i128 };

            for t in 0..k_cs {
                // Balanced remainder in [-b/2, b/2].
                let mut r0 = x.rem_euclid(b);
                if r0 > b_half {
                    r0 -= b;
                }
                x = (x - r0) / b;
                let fe = if r0 >= 0 {
                    R::BaseRing::from(r0 as u64)
                } else {
                    -R::BaseRing::from((-r0) as u64)
                };
                digits[t][k] = fe;
            }
        }
        for t in 0..k_cs {
            out[j * k_cs + t] = R::from(digits[t].clone());
        }
    }

    // Expand sparse matrices: each original column j becomes k_cs columns with weights b^t.
    fn expand<R: stark_rings::CoeffRing>(
        m: &SparseMatrix<R>,
        k_cs: usize,
        base_b: u128,
    ) -> SparseMatrix<R>
    where
        R::BaseRing: stark_rings::Zq,
    {
        let mut coeffs = vec![Vec::<(R, usize)>::new(); m.nrows];
        for (row, row_terms) in m.coeffs.iter().enumerate() {
            let mut new_row = Vec::with_capacity(row_terms.len() * k_cs);
            for (a, col) in row_terms.iter() {
                let mut bpow: u128 = 1;
                for t in 0..k_cs {
                    let w = R::from(bpow);
                    new_row.push((*a * w, col * k_cs + t));
                    bpow = bpow.saturating_mul(base_b);
                }
            }
            coeffs[row] = new_row;
        }
        SparseMatrix::<R> {
            nrows: m.nrows,
            ncols: m.ncols * k_cs,
            coeffs,
        }
    }

    (expand(m1, k_cs, base_b), expand(m2, k_cs, base_b), expand(m3, k_cs, base_b), out)
}

fn test_symphony_poseidon2_folding<
    R: stark_rings::CoeffRing,
    PC: cyclotomic_rings::rings::GetPoseidonParams<<<R>::BaseRing as ark_ff::Field>::BasePrimeField>,
>(
    num_permutations: usize,
    k_g: usize,
)
where
    R::BaseRing: stark_rings::Zq + stark_rings::balanced_decomposition::Decompose + ark_ff::PrimeField,
{
    use ark_ff::{BigInteger, Field, PrimeField};
    use latticefold::commitment::AjtaiCommitmentScheme;
    use symphony::{
        pcs::folding_pcs_l2::{
            BinMatrix, verify_folding_pcs_l2_with_c_matrices,
        },
        rp_rgchk::RPParams,
        symphony_open::MultiAjtaiOpenVerifier,
        symphony_pifold_batched::{verify_pi_fold_cp_poseidon_fs, PiFoldMatrices},
        symphony_pifold_streaming::{prove_pi_fold_poseidon_fs, PiFoldStreamingConfig},
        transcript::PoseidonTraceOp,
        poseidon_trace::find_squeeze_bytes_idx_after_absorb_marker,
        we_gate_arith::WeGateDr1csBuilder,
    };
    use cyclotomic_rings::rings::GetPoseidonParams;

    println!("\n=== Testing Symphony Poseidon2 with {} permutation(s), k_g={} ===", num_permutations, k_g);

    let build_start = Instant::now();
    let (m1, m2, m3, witness, actual_constraints) =
        build_poseidon2_hadamard_matrices::<R>(num_permutations);
    let build_time = build_start.elapsed();

    // Paper-style R1CS witness decomposition (k_cs digits base 2^4), expanding the witness and the
    // linear maps so the Hadamard relation is unchanged. This is the Table 1 instantiation setup
    // and is enabled unconditionally in this benchmark.
    // R1CS witness decomposition always enabled in this benchmark
    println!("  R1CS witness decomposition: ENABLED (k_cs=16, b=2^4)");
    let (m1, m2, m3, witness) =
        r1cs_decompose_witness_and_expand_matrices::<R>(&m1, &m2, &m3, &witness, 16, 16);

    // After witness decomposition, the witness length grows by k_cs, which increases m_J.
    // Π_fold requires: m_J <= m and m multiple of m_J, and also requires m*d is a power-of-two.
    // This benchmark pads the matrices with empty rows to satisfy these shape constraints.
    let mut m1 = m1;
    let mut m2 = m2;
    let mut m3 = m3;

    let n = witness.len();
    let mut m = m1.nrows;

    println!("  Hadamard matrices built:");
    println!("    Permutations: {num_permutations}");
    println!("    Actual constraints: {actual_constraints}");
    println!("    Matrix rows (m): {m}");
    println!("    Witness length (n): {n}");
    println!("    Build time: {build_time:?}");

    // Build `rg_params` early so we can compute the required m_J for this run.
    let rg_params = if k_g == 3 {
        // Paper-style prototype mode:
        // - use d' := d-2 in Eq.(29) / (33)
        // - keep m_J small by using lambda_pj=1 (toy bench knob; paper uses 2^8)
        RPParams {
            l_h: 64,
            lambda_pj: 1,
            k_g,
            d_prime: (R::dimension() as u128) - 2,
        }
    } else {
        RPParams {
            l_h: 64,
            lambda_pj: 32,
            k_g,
            d_prime: (R::dimension() as u128) / 2,
        }
    };

    let blocks = n / rg_params.l_h;
    let m_j = blocks
        .checked_mul(rg_params.lambda_pj)
        .expect("m_J overflow");
    if m < m_j || m % m_j != 0 {
        let m_new = (m.max(m_j)).next_power_of_two();
        println!(
            "  !! Padding matrices: m={} -> {} to satisfy m_J={} (l_h={}, lambda_pj={})",
            m, m_new, m_j, rg_params.l_h, rg_params.lambda_pj
        );
        // Extend coefficient rows with empties (zero rows).
        m1.coeffs.resize(m_new, Vec::new());
        m2.coeffs.resize(m_new, Vec::new());
        m3.coeffs.resize(m_new, Vec::new());
        m1.nrows = m_new;
        m2.nrows = m_new;
        m3.nrows = m_new;
        m = m_new;
        println!("    Matrix rows (m): {m} (padded)");
    }

    // Witness magnitude proxy: maximum bit-length among base-prime-field limbs.
    //
    // This is not "signed magnitude" in Zq, but it is a useful indicator of whether values
    // are staying tiny vs spanning the full modulus.
    let mut max_bits = 0u32;
    for r in &witness {
        for c in r.coeffs() {
            for limb in c.to_base_prime_field_elements() {
                max_bits = max_bits.max(limb.into_bigint().num_bits());
            }
        }
    }
    println!("    Witness max limb bits: {max_bits}");

    // Verify Hadamard relation holds: M1·f ⊙ M2·f = M3·f
    let y1 = m1.try_mul_vec(&witness).expect("M1*f failed");
    let y2 = m2.try_mul_vec(&witness).expect("M2*f failed");
    let y3 = m3.try_mul_vec(&witness).expect("M3*f failed");
    let mut hadamard_ok = true;
    for row in 0..m {
        if y1[row] * y2[row] != y3[row] {
            hadamard_ok = false;
            println!("    Hadamard FAILED at row {row}");
            break;
        }
    }
    if hadamard_ok {
        println!("  Hadamard relation: PASSED ✓");
    } else {
        println!("  Hadamard relation: FAILED ✗");
        return;
    }

    // Use Arcs so we can run multiple provers without cloning huge structures.
    let m1 = Arc::new(m1);
    let m2 = Arc::new(m2);
    let m3 = Arc::new(m3);
    let witness = Arc::new(witness);

    // Commitment setup (cm_f is PCS now).
    //
    // We commit to the **flattened witness coefficients** (Option A) using FoldingPCS(ℓ=2),
    // and treat the commitment surface `t` as the public `cm_f` object that is absorbed into Π_fold.
    let kappa = 8; // PCS commitment length (t_len)
    let flat_witness: Vec<BF<R>> = witness
        .iter()
        .flat_map(|re| {
            re.coeffs()
                .iter()
                .map(|c| c.to_base_prime_field_elements().into_iter().next().expect("bf limb"))
        })
        .collect();
    let pcs_params_f = symphony::pcs::cmf_pcs::cmf_pcs_params_for_flat_len::<BF<R>>(
        flat_witness.len(),
        kappa,
    )
    .expect("cm_f pcs params");
    let f_pcs_f = symphony::pcs::cmf_pcs::pad_flat_message(&pcs_params_f, &flat_witness);
    let (t_pcs_f, s_pcs_f) = symphony::pcs::folding_pcs_l2::commit(&pcs_params_f, &f_pcs_f)
        .expect("cm_f pcs commit failed");
    // Pack the PCS commitment surface into ring elements for Π_fold plumbing (fills coefficients).
    let cm = symphony::pcs::cmf_pcs::pack_t_as_ring::<R>(&t_pcs_f);

    // Symphony Π_rg parameters (defined above, after decomposition so it matches padded shapes).

    // CP commitment schemes for aux messages
    const MASTER_SEED: [u8; 32] = *b"SYMPHONY_AJTAI_SEED_V1_000000000";
    let scheme_had =
        AjtaiCommitmentScheme::<R>::seeded(b"cfs_had_u", MASTER_SEED, kappa, 3 * R::dimension());
    let scheme_mon = AjtaiCommitmentScheme::<R>::seeded(b"cfs_mon_b", MASTER_SEED, kappa, rg_params.k_g);
    let scheme_g = AjtaiCommitmentScheme::<R>::seeded(b"cm_g", MASTER_SEED, kappa, m * R::dimension());

    // Public inputs (statement binding)
    let public_inputs: Vec<R::BaseRing> = vec![
        R::BaseRing::from(0x5350315fu128), // "SP1_"
        R::BaseRing::from(num_permutations as u128),
        R::BaseRing::from(k_g as u128),
    ];

    // Open verifier for CP commitments (cloneable, so we can reuse across runs).
    let open = MultiAjtaiOpenVerifier::<R>::new()
        .with_scheme("cfs_had_u", scheme_had.clone())
        .with_scheme("cfs_mon_b", scheme_mon.clone());

    // Helper: verify R_cp for a produced output
    let verify = |out: &symphony::symphony_pifold_batched::PiFoldProverOutput<R>| -> Result<(), String> {
        let attempt = verify_pi_fold_cp_poseidon_fs::<R, PC>(
            PiFoldMatrices::Shared([m1.as_ref(), m2.as_ref(), m3.as_ref()]),
            &[cm.clone()],
            &out.proof,
            &open,
            &out.cfs_had_u,
            &out.cfs_mon_b,
            &out.aux,
            &public_inputs,
        );
        let _ = attempt.result?;
        Ok(())
    };

    // Streaming prover (canonical path)
    println!("  Π_fold prove (streaming): START...");

    // Build Ms: single instance
    let ms = vec![[m1.clone(), m2.clone(), m3.clone()]];

    let cfg = PiFoldStreamingConfig::default();
    let prove_start = Instant::now();
    let out = prove_pi_fold_poseidon_fs::<R, PC>(
        ms.as_slice(),
        &[cm.clone()],
        &[witness.clone()],
        &public_inputs,
        Some(&scheme_had),
        Some(&scheme_mon),
        &scheme_g,
        rg_params.clone(),
        &cfg,
    );
    let prove_time = prove_start.elapsed();
    let out = match out {
        Ok(o) => o,
        Err(e) => {
            println!("  Π_fold prove (streaming) FAILED: {e}");
            return;
        }
    };
    println!("  Π_fold prove (streaming): PASSED ✓ ({prove_time:?})");

    println!("  R_cp verify (streaming): START...");
    let verify_start = Instant::now();
    match verify(&out) {
        Ok(_) => {
            let verify_time = verify_start.elapsed();
            println!("  R_cp verify (streaming): PASSED ✓ ({verify_time:?})");
        }
        Err(e) => {
            println!("  R_cp verify (streaming) FAILED: {e}");
            return;
        }
    }

    // Gate dR1CS size report (R_cp vs full with PCS delta), using the *real verifier trace*.
    let poseidon_cfg = <PC as GetPoseidonParams<BF<R>>>::get_poseidon_config();
    let attempt = verify_pi_fold_cp_poseidon_fs::<R, PC>(
        PiFoldMatrices::Shared([m1.as_ref(), m2.as_ref(), m3.as_ref()]),
        &[cm.clone()],
        &out.proof,
        &open,
        &out.cfs_had_u,
        &out.cfs_mon_b,
        &out.aux,
        &public_inputs,
    );
    let _ = attempt
        .result
        .expect("unexpected verify failure when extracting trace for dr1cs sizing");
    let trace = attempt.trace;

    let (rcp, rcp_asg) = WeGateDr1csBuilder::r_cp_poseidon_pifold_math_and_cfs_openings::<R>(
        &poseidon_cfg,
        &trace.ops,
        &[cm.clone()],
        &out.proof,
        &scheme_had,
        &scheme_mon,
        &out.aux,
        &out.cfs_had_u,
        &out.cfs_mon_b,
    )
    .expect("build r_cp dr1cs failed");
    rcp.check(&rcp_asg).expect("r_cp dr1cs unsat");

    let squeeze_bytes: Vec<Vec<u8>> = trace
        .ops
        .iter()
        .filter_map(|op| match op {
            PoseidonTraceOp::SqueezeBytes { out, .. } => Some(out.clone()),
            _ => None,
        })
        .collect();
    if squeeze_bytes.is_empty() {
        println!(
            "    gate dr1cs: r_cp(nvars={}, constraints={}) (no SqueezeBytes => skip full pcs)",
            rcp.nvars,
            rcp.constraints.len()
        );
    } else {
        // Choose which Poseidon `SqueezeBytes` output seeds each PCS verifier.
        //
        // - PCS_f: keep the bench default (first SqueezeBytes) for now.
        // - PCS_g (batchlin): use the *squeeze immediately after* the mandatory
        // PCS#1 coin source: SqueezeBytes right after `CMF_PCS_DOMAIN_SEP`.
        let pcs_coin_squeeze_idx = find_squeeze_bytes_idx_after_absorb_marker(
            &trace.ops,
            <BF<R> as From<u128>>::from(symphony::pcs::cmf_pcs::CMF_PCS_DOMAIN_SEP),
        )
        .expect("missing SqueezeBytes after CMF_PCS_DOMAIN_SEP");
        let c_bytes = &squeeze_bytes[pcs_coin_squeeze_idx];
        let mut bits = Vec::with_capacity(c_bytes.len() * 8);
        for &b in c_bytes {
            for i in 0..8 {
                bits.push(((b >> i) & 1) == 1);
            }
        }

        // PCS#1 (cm_f PCS): MUST match the actual cm_f commitment surface absorbed into Π_fold.
        // So we use the prover's `pcs_params_f` and commitment `t_pcs_f`.
        let pcs_params = pcs_params_f.clone();
        let r = pcs_params.r;
        let kappa_pcs = pcs_params.kappa;
        let c1 = BinMatrix {
            rows: r * kappa_pcs,
            cols: kappa_pcs,
            data: (0..(r * kappa_pcs * kappa_pcs))
                .map(|i| if bits[i] { <BF<R> as Field>::ONE } else { <BF<R> as Field>::ZERO })
                .collect(),
        };
        let c2 = BinMatrix {
            rows: r * kappa_pcs,
            cols: kappa_pcs,
            data: (0..(r * kappa_pcs * kappa_pcs))
                .map(|i| {
                    if bits[(r * kappa_pcs * kappa_pcs) + i] {
                        <BF<R> as Field>::ONE
                    } else {
                        <BF<R> as Field>::ZERO
                    }
                })
                .collect(),
        };
        let x0 = vec![<BF<R> as Field>::ONE; r];
        let x1 = vec![<BF<R> as Field>::ONE; r];
        let x2 = vec![<BF<R> as Field>::ONE; r];
        let (u_pcs, pcs_core) = symphony::pcs::folding_pcs_l2::open(
            &pcs_params,
            &f_pcs_f,
            &s_pcs_f,
            &x0,
            &x1,
            &x2,
            &c1,
            &c2,
        )
        .expect("cm_f pcs open failed");
        verify_folding_pcs_l2_with_c_matrices(&pcs_params, &t_pcs_f, &x0, &x1, &x2, &u_pcs, &pcs_core, &c1, &c2)
            .expect("native folding pcs sanity failed");

        // use prover-produced commitment surface + proof core.
        use symphony::pcs::batchlin_pcs::{batchlin_scalar_pcs_params, BATCHLIN_PCS_DOMAIN_SEP};
        let log_n = ((out.proof.m * R::dimension()) as f64).log2() as usize;
        let pcs_params_g = batchlin_scalar_pcs_params::<BF<R>>(log_n).expect("batchlin pcs params");
        let t_pcs_g = &out.proof.batchlin_pcs_t[0];
        let x0_g = &out.proof.batchlin_pcs_x0;
        let x1_g = &out.proof.batchlin_pcs_x1;
        let x2_g = &out.proof.batchlin_pcs_x2;
        let pcs_core_g = &out.proof.batchlin_pcs_core;
        let pcs_coin_squeeze_idx2 = find_squeeze_bytes_idx_after_absorb_marker(
            &trace.ops,
            <BF<R> as From<u128>>::from(BATCHLIN_PCS_DOMAIN_SEP),
        )
        .expect("missing SqueezeBytes after BATCHLIN_PCS_DOMAIN_SEP");
        let (full, full_asg) = WeGateDr1csBuilder::poseidon_plus_pifold_plus_cfs_plus_pcs::<R>(
            &poseidon_cfg,
            &trace.ops,
            &[cm.clone()],
            &out.proof,
            &scheme_had,
            &scheme_mon,
            &out.aux,
            &out.cfs_had_u,
            &out.cfs_mon_b,
            &pcs_params,
            &t_pcs_f,
            &x0,
            &x1,
            &x2,
            &pcs_core,
            pcs_coin_squeeze_idx,
            &pcs_params_g,
            t_pcs_g,
            x0_g,
            x1_g,
            x2_g,
            pcs_core_g,
            pcs_coin_squeeze_idx2,
        )
        .expect("build full gate dr1cs failed");
        full.check(&full_asg).expect("full gate dr1cs unsat");

        println!(
            "    gate dr1cs: r_cp(nvars={}, constraints={})  full(nvars={}, constraints={})  delta(nvars={}, constraints={})",
            rcp.nvars,
            rcp.constraints.len(),
            full.nvars,
            full.constraints.len(),
            full.nvars.saturating_sub(rcp.nvars),
            full.constraints.len().saturating_sub(rcp.constraints.len()),
        );

        // Optional: run the *actual DPP proving + verify-with-query* over the full gate dR1CS.
        //
        // Enable on a fast machine:
        //   DPP=1 cargo run -p symphony --example symphony_poseidon2_bench --release
        //
        // This is intentionally not CI-friendly (can be very slow for large gates).
        if std::env::var("DPP").ok().as_deref() == Some("1") {
            use dpp::{
                dr1cs_flpcp::Dr1csInstanceSparse as DppDr1csInstanceSparse,
                dr1cs_flpcp::RsDr1csNpFlpcpSparse,
                embedding::EmbeddingParams,
                pipeline::build_rev2_dpp_sparse_boolean_auto,
                sparse::SparseVec,
                BooleanProofFlpcpSparse,
            };
            use sha2::{Digest, Sha256};
            use rand::{rngs::StdRng, RngCore, SeedableRng};
            use rayon::prelude::*;
            use symphony::we_statement::we_statement_hash_hetero_m;

            println!("    DPP: building sparse Dr1csInstance for full gate...");
            let t0 = Instant::now();
            let k = full.constraints.len();
            let nvars = full.nvars;
            let a_rows: Vec<SparseVec<BF<R>>> = full
                .constraints
                .par_iter()
                .map(|row| SparseVec::new(row.a.clone()))
                .collect();
            let b_rows: Vec<SparseVec<BF<R>>> = full
                .constraints
                .par_iter()
                .map(|row| SparseVec::new(row.b.clone()))
                .collect();
            let c_rows: Vec<SparseVec<BF<R>>> = full
                .constraints
                .par_iter()
                .map(|row| SparseVec::new(row.c.clone()))
                .collect();
            let dr1cs = DppDr1csInstanceSparse::<BF<R>> { n: nvars, a: a_rows, b: b_rows, c: c_rows };
            println!("    DPP: dr1cs build: {:?} (nvars={}, k={})", t0.elapsed(), nvars, k);

            // NP-style RS FLPCP over BF: x=[] (public), witness z_w = full assignment.
            let ell = 2 * k;
            let flpcp = RsDr1csNpFlpcpSparse::<BF<R>>::new(dr1cs, 0, ell);
            let x_small: Vec<BF<R>> = vec![];

            println!("    DPP: RS-FLPCP prove: START...");
            let t1 = Instant::now();
            // Run inside an explicit pool to ensure we actually use all cores for the heavy steps.
            let avail = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
            let pool = ThreadPoolBuilder::new()
                .num_threads(avail)
                .build()
                .expect("build rayon threadpool");
            let (pi_field, cw) = pool.install(|| flpcp.prove_with_codewords(&x_small, &full_asg));
            let t1e = t1.elapsed();
            println!(
                "    DPP: RS-FLPCP prove: {:?} (pi_field_len={})",
                t1e,
                pi_field.len()
            );

            println!("    DPP: booleanize + embed/pack: START...");
            let t2 = Instant::now();
            let boolized = BooleanProofFlpcpSparse::<BF<R>, _>::new(flpcp.clone());
            // Use bitpacked Boolean proof to avoid multi-GB allocations.
            let pi_bits_packed = pool.install(|| boolized.encode_proof_bits_packed(&pi_field));
            let dpp = build_rev2_dpp_sparse_boolean_auto::<BF<R>, FLarge, _>(
                flpcp,
                EmbeddingParams { gamma: 2, assume_boolean_proof: true, k_prime: 0 },
            )
            .expect("build_rev2_dpp_sparse_boolean_auto");
            let t2e = t2.elapsed();
            println!(
                "    DPP: booleanize+build: {:?} (pi_bits_len={}, packed_bytes={})",
                t2e,
                boolized.m_bits(),
                pi_bits_packed.len()
            );

            // Lift x into FLarge (empty for NP mode).
            let _x_large: Vec<FLarge> = vec![];
            let _pi_bits_len = boolized.m_bits();

            // Sample query and verify.
            println!("    DPP: verify_with_query: START...");
            println!("    DPP: rayon available_parallelism={avail}");
            let _pool = pool; // keep pool alive for the remainder of the DPP section

            let t4 = Instant::now();
            // WE-friendly / statement-only coins: derive a deterministic RNG seed from the
            // (vk_hash, r1cs_digest, gate_digest, public_inputs) statement hash.
            //
            // NOTE: In this benchmark we don't have a real SP1 vk hash nor a canonical gate digest
            // yet, so we use placeholders plus a shape-dependent gate digest derived from the
            // benchmark parameters. This is sufficient to demonstrate the intended structure.
            // For this bench (no SP1), we still model statement binding by using:
            // - vk_hash: a fixed "program id" digest
            // - r1cs_digest: a digest of the Poseidon2/Hadamard instance parameters
            // - gate_digest: a digest of the *actual built* WE gate constraint system shape
            let vk_hash = {
                let mut h = Sha256::new();
                h.update(b"SYMPHONY_BENCH_VK_V1::POSEIDON2_HADAMARD");
                h.finalize().into()
            };
            let r1cs_digest = {
                let mut h = Sha256::new();
                h.update(b"SYMPHONY_BENCH_R1CS_V1::POSEIDON2_HADAMARD");
                h.update(&(num_permutations as u64).to_le_bytes());
                h.update(&(k_g as u64).to_le_bytes());
                h.update(&(n as u64).to_le_bytes());
                h.update(&(m as u64).to_le_bytes());
                h.finalize().into()
            };
            let gate_digest = {
                let mut h = Sha256::new();
                h.update(b"SYMPHONY_GATE_DIGEST_BENCH_V1");
                h.update(&(full.nvars as u64).to_le_bytes());
                h.update(&(full.constraints.len() as u64).to_le_bytes());
                // Hash the sparse constraints deterministically. (This is O(|constraints|) work;
                // in production we'd treat the gate digest as a precomputed per-shape constant.)
                for row in &full.constraints {
                    h.update(&(row.a.len() as u64).to_le_bytes());
                    for (c, idx) in &row.a {
                        h.update(&(*idx as u64).to_le_bytes());
                        h.update(c.into_bigint().to_bytes_le());
                    }
                    h.update(&(row.b.len() as u64).to_le_bytes());
                    for (c, idx) in &row.b {
                        h.update(&(*idx as u64).to_le_bytes());
                        h.update(c.into_bigint().to_bytes_le());
                    }
                    h.update(&(row.c.len() as u64).to_le_bytes());
                    for (c, idx) in &row.c {
                        h.update(&(*idx as u64).to_le_bytes());
                        h.update(c.into_bigint().to_bytes_le());
                    }
                }
                h.finalize().into()
            };
            let stmt_digest =
                we_statement_hash_hetero_m::<R>(vk_hash, r1cs_digest, gate_digest, &public_inputs);

            // Model WE arming: an armer derives per-lock public coins from (armer_seed, stmt_digest, j).
            // Those public coins define the query randomness (idx, λ) and packing weights.
            const ARMER_SEED: [u8; 32] = *b"SYMPHONY_ARMER_SEED_V1_000000000";
            let lock_j: u64 = 0;
            let coin_j: [u8; 32] = {
                let mut h = Sha256::new();
                h.update(b"SYMPHONY_LOCK_COIN_V1");
                h.update(&ARMER_SEED);
                h.update(&stmt_digest);
                h.update(&lock_j.to_le_bytes());
                h.finalize().into()
            };
            let mut rng = StdRng::from_seed(coin_j);
            // NOTE: Avoid building booleanized query vectors (which explode to 100M+ terms).
            // For lock evaluation, we can sample the *base* 3-query RS-FLPCP query over BF<R>,
            // compute its 3 answers against the BF witness/proof, then pack those answers with `w`.
            //
            // This models the intended WE flow: locks are published as coins, not expanded vectors.
            let b = dpp.flpcp.bounds_b();
            let w = dpp::packing::sample_packing_weights::<FLarge>(&mut rng, dpp.params.ell, &b)
                .expect("sample_packing_weights");
            let pred = dpp::packing::FlpcpPredicate::MulEqModP {
                p_small: num_bigint::BigInt::from_bytes_le(
                    num_bigint::Sign::Plus,
                    &BF::<R>::MODULUS.to_bytes_le(),
                ),
            };

            // Sample RS-FLPCP verifier coins (idx, lambda) directly.
            let k = cw.y_a.len();
            let ell = 2 * k;
            let idx = (rng.next_u64() as usize) % ell;
            let lambda_small = BF::<R>::from(rng.next_u64());

            // Split verification timing into:
            // - compute k answers a_i = <q_i, v>
            // - pack a = <a, w> (integers) and reduce to field
            // - bounded decoding + predicate
            let (a, t_dot) = {
                let t = Instant::now();
                // Evaluate the 3-query RS-FLPCP answers in coin form, by indexing codewords.
                let (a_small, b_small, c_small) = {
                    if idx < k {
                        let a = cw.y_a[idx];
                        let b = cw.y_b[idx];
                        let wv = cw.w[idx];
                        let cx_minus = cw.y_c[idx] - wv;
                        let c = wv + lambda_small * cx_minus;
                        (a, b, c)
                    } else {
                        let j = idx - k;
                        let a = cw.y_a_tail[j];
                        let b = cw.y_b_tail[j];
                        let wv = cw.w[idx];
                        // For α in the "tail" half (k+1..2k), the verifier's q3 uses only the w-part:
                        // the C-part coefficient vector is zero because λ_{2k}(α) is the basis vector at idx>=k
                        // and `sample_queries_and_predicate_sparse` uses only `lam_2k[..k]` for the C-part.
                        // So the correct answer is just w(α) = a*b.
                        let c = wv;
                        (a, b, c)
                    }
                };

                let ans_field: [FLarge; 3] = [
                    FLarge::from_le_bytes_mod_order(&a_small.into_bigint().to_bytes_le()),
                    FLarge::from_le_bytes_mod_order(&b_small.into_bigint().to_bytes_le()),
                    FLarge::from_le_bytes_mod_order(&c_small.into_bigint().to_bytes_le()),
                ];

                // Pack into one integer a_int = Σ w_i * [ans_i]_centered, then reduce to field.
                let mut a_int = num_bigint::BigInt::from(0);
                for (wi, ai) in w.iter().zip(ans_field.iter()) {
                    let ai_int = dpp::packing::field_to_centered_bigint::<FLarge>(ai);
                    a_int += wi * ai_int;
                }
                let a = dpp::packing::centered_bigint_to_field::<FLarge>(&a_int);
                (a, t.elapsed())
            };
            println!("    DPP: unpacked-dot+pack: {:?}", t_dot);

            let (ok, t_dec) = {
                let t = Instant::now();
                // Use the same decoding + predicate as the packed DPP verifier.
                let q_meta = dpp::packing::PackedDppQuerySparse {
                    q: dpp::sparse::SparseVec::default(),
                    w,
                    b,
                    pred,
                };
                let ok = dpp.verify_packed_answer(&a, &q_meta);
                (ok, t.elapsed())
            };
            let ok = ok.expect("verify_packed_answer");
            println!("    DPP: decode+pred: {:?}", t_dec);
            let t4e = t4.elapsed();
            println!("    DPP: verify_with_query: {:?} (ok={})", t4e, ok);
        }
    }

    println!("\n  === RESULTS ===");
    println!("  Permutations: {num_permutations}");
    println!("  Constraints: {actual_constraints} (padded to {m})");
    println!("  Witness size: {n}");
    println!("  Streaming: prove={prove_time:?}, coins_bytes={}", out.proof.coins.bytes.len());
}
