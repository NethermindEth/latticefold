use ark_std::{log2, vec::Vec};
use lattirust_linear_algebra::SparseMatrix;
use lattirust_ring::Ring;

use super::{
    r1cs::{create_dummy_identity_sparse_matrix, to_F_matrix, to_F_vec},
    CCS,
};

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

pub fn get_test_degree_three_z<R: Ring>(input: usize) -> Vec<R> {
    // z = (io, 1, w)
    to_F_vec(vec![
        input, // io
        1,
        input * input * input,             // x^3
        input * input * input + input,     // x^3 + x
        input * input * input + input + 5, // x^3 +x + 5
    ])
}

pub fn get_test_degree_three_z_split<R: Ring>(input: usize) -> (R, Vec<R>, Vec<R>) {
    let z = get_test_degree_three_z(input);
    (z[1], vec![z[0]], z[2..].to_vec())
}

pub fn get_test_degree_three_ccs<R: Ring>() -> CCS<R> {
    // Degree 3 CCS for: x^3 + x + 5 = y
    let A = to_F_matrix::<R>(vec![
        vec![1, 0, 0, 0, 0],
        vec![1, 0, 1, 0, 0],
        vec![0, 5, 0, 1, 0],
    ]);
    let B = to_F_matrix::<R>(vec![
        vec![1, 0, 0, 0, 0],
        vec![0, 1, 0, 0, 0],
        vec![0, 1, 0, 0, 0],
    ]);

    let C = to_F_matrix::<R>(vec![
        vec![1, 0, 0, 0, 0],
        vec![0, 1, 0, 0, 0],
        vec![0, 1, 0, 0, 0],
    ]);
    let D = to_F_matrix::<R>(vec![
        vec![0, 0, 1, 0, 0],
        vec![0, 0, 0, 1, 0],
        vec![0, 0, 0, 0, 1],
    ]);

    CCS {
        m: 3,
        n: 5,
        l: 1,
        t: 4,
        q: 2,
        d: 3,
        s: log2(3) as usize,
        s_prime: log2(5) as usize,
        M: vec![A, B, C, D],
        S: vec![vec![0, 1, 2], vec![3]],
        c: vec![R::one(), R::one().neg()],
    }
}

pub fn get_test_degree_three_ccs_padded<R: Ring>(W: usize, L: usize) -> CCS<R> {
    let mut ccs = get_test_degree_three_ccs();

    ccs.m = W;
    ccs.s = log2(W) as usize;
    let len = usize::max((ccs.n - ccs.l - 1) * L, ccs.m).next_power_of_two();
    ccs.pad_rows_to(len);
    ccs
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
#[allow(clippy::upper_case_acronyms)]
#[cfg(test)]
mod tests {
    use cyclotomic_rings::rings::GoldilocksRingNTT;

    use crate::arith::{ccs::get_test_degree_three_z, Arith, CCS};

    use super::get_test_degree_three_ccs;
    type NTT = GoldilocksRingNTT;

    #[test]
    fn test_degree_three_ccs() {
        let input = 5;
        let ccs: CCS<NTT> = get_test_degree_three_ccs();
        let z = get_test_degree_three_z(input);
        assert!(ccs.check_relation(&z).is_ok())
    }
}
