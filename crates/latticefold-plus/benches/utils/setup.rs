use stark_rings::{
    cyclotomic_ring::models::frog_ring::RqPoly as R,
    PolyRing,
};
use stark_rings_linalg::{Matrix, SparseMatrix};
use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
use latticefold_plus::{
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