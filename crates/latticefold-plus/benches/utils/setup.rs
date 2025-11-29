use stark_rings::{
    cyclotomic_ring::models::frog_ring::RqPoly as R,
    PolyRing,
};
use stark_rings_linalg::{Matrix, SparseMatrix};

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