use ark_ff::{One, PrimeField, Zero};
use rand::prelude::*;
use stark_rings::{
    cyclotomic_ring::models::frog_ring::RqPoly as R,
    PolyRing,
};
use stark_rings_linalg::{Matrix, SparseMatrix};
use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
use latticefold_plus::{
    rgchk::{DecompParameters, Rg, RgInstance},
    setchk::{In, MonomialSet, Out},
    transcript::PoseidonTranscript,
};

pub fn setup_split_input(k_first: usize, kappa: usize) -> Matrix<R> {
    let mut rng = bench_rng();
    let d = R::dimension(); // 16 for FrogRing
    let cols = k_first * d;
    let mat = Matrix::<R>::rand(&mut rng, kappa, cols);

    // Validation: ensure matrix has correct dimensions
    assert_eq!(mat.nrows, kappa, "Matrix should have kappa rows");
    assert_eq!(mat.ncols, cols, "Matrix should have k_first * d columns");

    mat
}

pub fn bench_rng() -> impl rand::Rng {
    use rand::SeedableRng;
    rand::rngs::StdRng::seed_from_u64(0x42424242)
}

/// Generate test input for set check benchmarks
///
/// Creates monomial sets (identity matrices) which are guaranteed to satisfy
/// the monomial property.
///
/// # Arguments
/// * `set_size` - Size of each monomial set
/// * `num_batches` - Number of sets to check (tests batching)
///
/// # Returns
/// `In<R>` structure containing valid monomial sets
///
/// # Validation
/// - All matrices are identity matrices (guaranteed monomials)
/// - nvars is correctly computed as log2(set_size)
pub fn setup_setchk_input(set_size: usize, num_batches: usize) -> In<R> {
    let mut sets = Vec::with_capacity(num_batches);

    for _ in 0..num_batches {
        // Identity matrix guarantees monomials (one 1 per row/column)
        let m = SparseMatrix::<R>::identity(set_size);
        sets.push(MonomialSet::Matrix(m));
    }

    // nvars = log2(set_size) rounded up
    let nvars = (set_size as f64).log2().ceil() as usize;

    In { sets, nvars }
}

/// Generate set check proof for verification benchmarks
///
/// Generates a complete valid proof by running the prover.
///
/// # Arguments
/// * `set_size` - Size of monomial set
/// * `num_batches` - Number of batched sets
///
/// # Returns
/// Tuple of (input, output/proof)
///
/// # Validation
/// - Proof is verified to be valid before returning
pub fn setup_setchk_proof(set_size: usize, num_batches: usize) -> (In<R>, Out<R>) {
    let input = setup_setchk_input(set_size, num_batches);
    let mut ts = PoseidonTranscript::empty::<PC>();

    // Generate proof
    let output = input.set_check(&[], &mut ts);

    // Cryptographic validation: verify proof is valid
    let mut verify_ts = PoseidonTranscript::empty::<PC>();
    output
        .verify(&mut verify_ts)
        .expect("Generated set check proof should be valid");

    (input, output)
}

/// Generate test input for range check benchmarks
///
/// Creates a valid range check instance following the exact pattern from test cases.
/// Decomposes witness vector f and creates monomial matrices M_f for range checking.
///
/// # Arguments
/// * `witness_size` - Length of witness vector (must satisfy: witness_size >= kappa * k * d * l * d)
/// * `k` - Decomposition width for first decomposition (determines range B = (d/2)^k)
/// * `kappa` - Number of commitment rows (security parameter)
///
/// # Returns
/// `Rg<R>` structure ready for range checking
///
/// # Panics
/// Panics if witness_size violates the constraint: witness_size >= kappa * k * d * l * d
pub fn setup_rgchk_input(witness_size: usize, k: usize, kappa: usize) -> Rg<R> {
    let mut rng = bench_rng();

    // Compute decomposition parameters from ring parameters (NOT hardcoded!)
    let d = R::dimension();
    let b = (d / 2) as u128;

    // Compute l = ⌈log_b(q)⌉ where q is the base ring modulus
    let l = ((<<R as PolyRing>::BaseRing>::MODULUS.0[0] as f64).ln()
        / (b as f64).ln())
    .ceil() as usize;

    // Validate constraint: witness_size >= kappa * k * d * l * d
    let min_witness_size = kappa * k * d * l * d;
    assert!(
        witness_size >= min_witness_size,
        "Invalid parameters: witness_size ({}) must be >= kappa * k * d * l * d = {} * {} * {} * {} * {} = {}",
        witness_size, kappa, k, d, l, d, min_witness_size
    );

    // Generate witness vector with small coefficients
    let f: Vec<R> = (0..witness_size)
        .map(|_| {
            let mut r = R::zero();
            // Use small coefficients like the test cases (range [0, 10))
            r.coeffs_mut()[0] = ((rng.next_u32() % 10) as u64).into();
            r
        })
        .collect();

    // Generate Ajtai commitment matrix (size: kappa × witness_size)
    let A = Matrix::<R>::rand(&mut rng, kappa, witness_size);

    let dparams = DecompParameters { b, k, l };

    // Create range check instance using from_f
    let instance = RgInstance::from_f(f, &A, &dparams);

    let nvars = (witness_size as f64).log2().ceil() as usize;

    Rg {
        nvars,
        instances: vec![instance],
        dparams,
    }
}

/// Generate range check proof for verification benchmarks
///
/// Creates a complete valid proof by running the prover.
///
/// # Arguments
/// * `witness_size` - Length of witness vector
/// * `k` - Decomposition width
/// * `kappa` - Number of commitment rows
///
/// # Returns
/// Tuple of (input, output/proof)
///
/// # Validation
/// Proof is verified to be valid before returning
pub fn setup_rgchk_proof(
    witness_size: usize,
    k: usize,
    kappa: usize,
) -> (Rg<R>, latticefold_plus::rgchk::Dcom<R>) {
    let rg = setup_rgchk_input(witness_size, k, kappa);
    let mut ts = PoseidonTranscript::empty::<PC>();

    // Generate proof via range check
    let dcom = rg.range_check(&[], &mut ts);

    // Verify proof is valid
    let mut verify_ts = PoseidonTranscript::empty::<PC>();
    dcom.verify(&mut verify_ts)
        .expect("Generated range check proof should be valid");

    (rg, dcom)
}