use ark_std::log2;
use lattirust_linear_algebra::SparseMatrix;
use lattirust_ring::Ring;

use super::{r1cs::create_dummy_identity_sparse_matrix, CCS};

pub fn get_test_dummy_degree_three_ccs_non_scalar<
    R: Ring,
    const X_LEN: usize,
    const WIT_LEN: usize,
>(
    rows: usize,
    witness: &[R],
) -> CCS<R> {
    let A = create_dummy_identity_sparse_matrix(rows, X_LEN + WIT_LEN + 1);
    let B = A.clone();
    let C = A.clone();
    let D = create_dummy_cubing_sparse_matrix(rows, X_LEN + WIT_LEN + 1, witness);

    CCS {
        m: rows,
        n: X_LEN + WIT_LEN + 1,
        l: 1,
        t: 4,
        q: 2,
        d: 3,
        s: log2(rows) as usize,
        s_prime: (X_LEN + WIT_LEN + 1),
        M: vec![A, B, C, D],
        S: vec![vec![0, 1, 2], vec![3]],
        c: vec![R::one(), R::one().neg()],
    }
}

// Takes a vector and returns a matrix that will square the vector
pub fn create_dummy_cubing_sparse_matrix<R: Ring>(
    rows: usize,
    columns: usize,
    witness: &[R],
) -> SparseMatrix<R> {
    assert_eq!(
        rows,
        witness.len(),
        "Length of witness vector must be equal to ccs width"
    );
    let mut matrix = SparseMatrix {
        n_rows: rows,
        n_cols: columns,
        coeffs: vec![vec![]; rows],
    };
    for (i, row) in matrix.coeffs.iter_mut().enumerate() {
        row.push((witness[i] * witness[i], i));
    }
    matrix
}
