use ark_ff::{One, PrimeField, Zero};
use rand::prelude::*;
use stark_rings::{
    cyclotomic_ring::models::frog_ring::RqPoly as R,
    PolyRing,
};
use stark_rings_linalg::{Matrix, SparseMatrix};
use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
use latticefold::arith::r1cs::R1CS;
use latticefold_plus::{
    decomp::{Decomp, DecompProof},
    lin::{LinB, Linearize, LinParameters, LinearizedVerify},
    mlin::{LinB2, Mlin},
    r1cs::{r1cs_decomposed_square, ComR1CS},
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

/// Generate test input for commitment transformation benchmarks
///
/// # Arguments
/// * `L` - Number of instances to transform/fold
/// * `witness_size` - Length of witness vector (must satisfy: witness_size >= kappa * k * d * l * d)
/// * `k` - Decomposition width for range check
/// * `kappa` - Number of commitment rows (security parameter)
///
/// # Returns
/// `Cm<R>` structure ready for commitment transformation
///
/// # Panics
/// Panics if witness_size violates the constraint: witness_size >= kappa * k * d * l * d
pub fn setup_cm_input(L: usize, witness_size: usize, k: usize, kappa: usize) -> latticefold_plus::cm::Cm<R> {
    let mut rng = bench_rng();

    // Compute decomposition parameters from ring parameters
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

    let dparams = DecompParameters { b, k, l };

    // Create L instances (for folding)
    let instances: Vec<RgInstance<R>> = (0..L)
        .map(|_| {
            // Generate witness vector with small coefficients
            let f: Vec<R> = (0..witness_size)
                .map(|_| {
                    let mut r = R::zero();
                    // Use small coefficients (range [0, 10))
                    r.coeffs_mut()[0] = ((rng.next_u32() % 10) as u64).into();
                    r
                })
                .collect();

            // Generate Ajtai commitment matrix (size: kappa × witness_size)
            let A = Matrix::<R>::rand(&mut rng, kappa, witness_size);

            // Create range check instance using from_f
            RgInstance::from_f(f, &A, &dparams)
        })
        .collect();

    let nvars = (witness_size as f64).log2().ceil() as usize;

    let rg = Rg {
        nvars,
        instances,
        dparams,
    };

    latticefold_plus::cm::Cm { rg }
}

/// Generate commitment transformation proof for verification benchmarks
///
/// # Arguments
/// * `L` - Number of instances to transform/fold
/// * `witness_size` - Length of witness vector
/// * `k` - Decomposition width
/// * `kappa` - Number of commitment rows
///
/// # Returns
/// Tuple of (input, proof)
///
/// # Validation
/// Proof is verified to be valid before returning
pub fn setup_cm_proof(
    L: usize,
    witness_size: usize,
    k: usize,
    kappa: usize,
) -> (latticefold_plus::cm::Cm<R>, latticefold_plus::cm::CmProof<R>) {
    let cm = setup_cm_input(L, witness_size, k, kappa);

    // Create modified identity matrix M
    let mut m = SparseMatrix::identity(witness_size);
    m.coeffs[0][0].0 = 2u128.into();
    let M = vec![m];

    let mut ts = PoseidonTranscript::empty::<PC>();

    // Generate proof via commitment transformation
    let (_com, proof) = cm.prove(&M, &mut ts);

    // Verify proof is valid
    let mut verify_ts = PoseidonTranscript::empty::<PC>();
    proof
        .verify(&M, &mut verify_ts)
        .expect("Generated commitment transformation proof should be valid");

    (cm, proof)
}

/// Generate test input for multilinear folding benchmarks
///
/// Creates valid Mlin<R> input.
/// This function creates L linearized R1CS instances ready for multilinear folding.
///
/// # Arguments
/// * `L` - Number of instances to fold (higher L = better amortization)
/// * `n` - Witness size (length of witness vector after decomposition)
/// * `k` - Decomposition width
/// * `kappa` - Number of commitment rows
/// * `B` - Norm bound parameter
///
/// # Returns
/// `Mlin<R>`
///
/// # Panics
/// Panics if n violates the constraint: n >= kappa * k * d * l * d
pub fn setup_mlin_input(L: usize, n: usize, k: usize, kappa: usize, B: usize) -> Mlin<R> {
    let mut rng = bench_rng();

    // Compute decomposition parameters from ring parameters
    let d = R::dimension();
    let b = (d / 2) as u128;

    // Compute l = ⌈log_b(q)⌉ where q is the base ring modulus
    let l = ((<<R as PolyRing>::BaseRing>::MODULUS.0[0] as f64).ln()
        / (b as f64).ln())
    .ceil() as usize;

    // Validate constraint: n >= kappa * k * d * l * d
    let min_witness_size = kappa * k * d * l * d;
    assert!(
        n >= min_witness_size,
        "Invalid parameters: n ({}) must be >= kappa * k * d * l * d = {} * {} * {} * {} * {} = {}",
        n, kappa, k, d, l, d, min_witness_size
    );

    let params = LinParameters {
        kappa,
        decomp: DecompParameters { b, k, l },
    };

    // Create R1CS with modified first entry
    let mut r1cs = r1cs_decomposed_square(
        R1CS::<R> {
            l: 1,
            A: SparseMatrix::identity(n / k),
            B: SparseMatrix::identity(n / k),
            C: SparseMatrix::identity(n / k),
        },
        n,
        B as u128,
        k,
    );

    r1cs.A.coeffs[0][0].0 = 2u128.into();
    r1cs.C.coeffs[0][0].0 = 2u128.into();

    // Generate Ajtai commitment matrix (size: kappa × n)
    let A = Matrix::<R>::rand(&mut rng, kappa, n);

    // Create L different witnesses with small variations
    let mut lins = Vec::with_capacity(L);
    let mut ts = PoseidonTranscript::empty::<PC>();

    for i in 0..L {
        // Create witness: mostly ones, with small variations
        let mut z = vec![R::one(); n / k];
        if i > 0 {
            // Vary first element for non-first instances
            z[0] = R::from((i % 10) as u128);
        }

        // Create ComR1CS instance
        let cr1cs = ComR1CS::new(r1cs.clone(), z, 1, B as u128, k, &A);

        // Linearize to get LinB instance
        let (linb, _lproof) = cr1cs.linearize(&mut ts);
        lins.push(linb);
    }

    Mlin { lins, params }
}

/// Generate multilinear folding proof for verification benchmarks
///
/// Creates a complete valid proof by running the prover.
///
/// # Arguments
/// * `L` - Number of instances to fold
/// * `n` - Witness size (length of witness vector after decomposition)
/// * `k` - Decomposition width
/// * `kappa` - Number of commitment rows
/// * `B` - Norm bound parameter
///
/// # Returns
/// Tuple of (input Mlin, folded LinB2, proof)
///
/// # Validation
/// Proof is verified to be valid before returning
pub fn setup_mlin_proof(
    L: usize,
    n: usize,
    k: usize,
    kappa: usize,
    B: usize,
) -> (Mlin<R>, LinB2<R>, latticefold_plus::cm::CmProof<R>) {
    let mlin = setup_mlin_input(L, n, k, kappa, B);
    let mut rng = bench_rng();

    // Create M matrix
    let mut m = SparseMatrix::identity(n);
    m.coeffs[0][0].0 = 2u128.into();
    let M = vec![m];

    // Create another A matrix for folding
    let A = Matrix::<R>::rand(&mut rng, kappa, n);

    let mut ts = PoseidonTranscript::empty::<PC>();

    // Execute multilinear folding
    let (linb2, proof) = mlin.mlin(&A, &M, &mut ts);

    // Cryptographic validation: verify proof is valid
    let mut verify_ts = PoseidonTranscript::empty::<PC>();
    proof
        .verify(&M, &mut verify_ts)
        .expect("Generated multilinear folding proof should be valid");

    (mlin, linb2, proof)
}

/// Generate test input for decomposition benchmarks
///
/// Creates valid Decomp<R> input for Construction 5.3 (Π_decomp,B).
/// This decomposes a LinB2 instance (norm B²) into 2 LinB instances (norm B each).
///
/// # Arguments
/// * `n` - Witness size (with norm B² before decomposition)
/// * `k` - Decomposition width
/// * `kappa` - Number of commitment rows
/// * `B` - Norm bound parameter (output will have norm B, input has B²)
///
/// # Returns
/// `Decomp<R>` structure ready for decomposition
///
/// # Panics
/// Panics if n violates the constraint: n >= kappa * k * d * l * d
pub fn setup_decomp_input(n: usize, k: usize, kappa: usize, _B: usize) -> Decomp<R> {
    use ark_std::One;
    use latticefold::arith::r1cs::R1CS;
    use latticefold_plus::r1cs::r1cs_decomposed_square;

    // Create identity R1CS
    let mut r1cs = R1CS::<R> {
        l: 1,
        A: SparseMatrix::identity(n / k),
        B: SparseMatrix::identity(n / k),
        C: SparseMatrix::identity(n / k),
    };
    r1cs.A.coeffs[0][0].0 = 2u128.into();
    r1cs.C.coeffs[0][0].0 = 2u128.into();

    // Apply r1cs_decomposed_square transformation
    let r1cs = r1cs_decomposed_square(r1cs, n, 2, k);

    // Create witness (all ones)
    let z = vec![R::one(); n / k];

    let mut rng = bench_rng();
    let A = Matrix::<R>::rand(&mut rng, kappa, n);

    // Create ComR1CS
    let cr1cs = ComR1CS::new(r1cs, z, 1, 2, k, &A);

    // Linearize to get random challenges
    let mut ts = PoseidonTranscript::empty::<PC>();
    let (linb, lproof) = cr1cs.linearize(&mut ts);

    // Verify linearization
    let mut ts = PoseidonTranscript::empty::<PC>();
    lproof.verify(&mut ts);

    // Create Decomp with witness f for prover
    // Note: f is the witness vector (needed for decompose), not the commitment
    Decomp {
        f: cr1cs.f,
        r: lproof.r.iter().map(|&r| (r, r)).collect::<Vec<_>>(),
        M: cr1cs.x.matrices(),
    }
}

/// Generate decomposition proof for verification benchmarks
///
/// Creates a complete valid proof by running the decomposition prover.
///
/// # Arguments
/// * `n` - Witness size
/// * `k` - Decomposition width
/// * `kappa` - Number of commitment rows
/// * `B` - Norm bound parameter
///
/// # Returns
/// Tuple of (input Decomp, two output LinB instances, proof)
///
/// # Validation
/// Proof is verified to be valid before returning
pub fn setup_decomp_proof(
    n: usize,
    k: usize,
    kappa: usize,
    B: usize,
) -> (Decomp<R>, (LinB<R>, LinB<R>), DecompProof<R>) {
    use ark_std::One;
    use latticefold::arith::r1cs::R1CS;
    use latticefold_plus::r1cs::r1cs_decomposed_square;

    // Create identity R1CS
    let mut r1cs = R1CS::<R> {
        l: 1,
        A: SparseMatrix::identity(n / k),
        B: SparseMatrix::identity(n / k),
        C: SparseMatrix::identity(n / k),
    };
    r1cs.A.coeffs[0][0].0 = 2u128.into();
    r1cs.C.coeffs[0][0].0 = 2u128.into();
    let r1cs = r1cs_decomposed_square(r1cs, n, 2, k);

    let z = vec![R::one(); n / k];

    let mut rng = bench_rng();
    let A = Matrix::<R>::rand(&mut rng, kappa, n);

    let cr1cs = ComR1CS::new(r1cs, z, 1, 2, k, &A);

    let mut ts = PoseidonTranscript::empty::<PC>();
    let (linb, lproof) = cr1cs.linearize(&mut ts);

    let mut ts = PoseidonTranscript::empty::<PC>();
    lproof.verify(&mut ts);

    // Create Decomp with witness f for calling decompose
    let decomp_prover = Decomp {
        f: cr1cs.f,
        r: lproof.r.iter().map(|&r| (r, r)).collect::<Vec<_>>(),
        M: cr1cs.x.matrices(),
    };

    // Execute decomposition with B
    let ((linb0, linb1), proof) = decomp_prover.decompose(&A, B as u128);

    // Verify proof
    proof.verify(&cr1cs.x.cm_f, &linb.x.v, B as u128);

    // Create Decomp with commitment and v for verifier
    let decomp_verifier = Decomp {
        f: cr1cs.x.cm_f,
        r: linb.x.v,
        M: decomp_prover.M,
    };

    (decomp_verifier, (linb0, linb1), proof)
}
