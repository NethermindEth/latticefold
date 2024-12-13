//! Provide test and benchmark utility for CCS

use ark_std::{log2, vec::Vec};
use cyclotomic_rings::rings::SuitableRing;
use lattirust_linear_algebra::SparseMatrix;
use lattirust_ring::Ring;

use super::{
    r1cs::{create_dummy_identity_sparse_matrix, to_F_matrix, to_F_vec},
    CCS,
};

/// Given a witness, provides a satisfying degree three CCS of arbitrary size
pub fn get_test_dummy_degree_three_ccs_non_scalar<
    R: Ring,
    const X_LEN: usize,
    const WIT_LEN: usize,
    const W: usize,
>(
    witness: &[R],
    L: usize,
    n_rows: usize,
) -> CCS<R> {
    let A = create_dummy_identity_sparse_matrix(n_rows, X_LEN + WIT_LEN + 1);
    let B = A.clone();
    let C = A.clone();
    let D = create_dummy_cubing_sparse_matrix(n_rows, X_LEN + WIT_LEN + 1, witness);

    let mut ccs = CCS {
        m: W,
        n: X_LEN + WIT_LEN + 1,
        l: 1,
        t: 4,
        q: 2,
        d: 3,
        s: log2(W) as usize,
        s_prime: (X_LEN + WIT_LEN + 1),
        M: vec![A, B, C, D],
        S: vec![vec![0, 1, 2], vec![3]],
        c: vec![R::one(), R::one().neg()],
    };
    let len = usize::max((ccs.n - ccs.l - 1) * L, ccs.m).next_power_of_two();
    ccs.pad_rows_to(len);
    ccs
}

pub(crate) fn get_test_degree_three_z<R: Ring>(input: usize) -> Vec<R> {
    // z = (x, 1, w)
    to_F_vec(vec![
        input, // x
        1,
        input * input * input,             // x^3
        input * input * input + input,     // x^3 + x
        input * input * input + input + 5, // x^3 +x + 5
    ])
}

pub(crate) fn get_test_degree_three_z_non_scalar<R: SuitableRing>() -> Vec<R> {
    let mut res = Vec::new();
    for input in 0..R::dimension() {
        // z = (io, 1, w)
        res.push(to_F_vec::<R::BaseRing>(vec![
            input, // x
            1,
            input * input * input,             // x^3
            input * input * input + input,     // x^3 + x
            input * input * input + input + 5, // x^3 +x + 5
        ]))
    }

    let mut ret: Vec<R> = Vec::new();
    for j in 0..res[0].len() {
        let mut vec = Vec::new();
        for witness in &res {
            vec.push(witness[j]);
        }
        ret.push(R::from(vec));
    }

    ret
}

#[allow(dead_code)]
pub(crate) fn get_test_degree_three_z_split<R: Ring>(input: usize) -> (R, Vec<R>, Vec<R>) {
    let z = get_test_degree_three_z(input);
    (z[1], vec![z[0]], z[2..].to_vec())
}

#[allow(dead_code)]
pub(crate) fn get_test_degree_three_z_non_scalar_split<R: SuitableRing>() -> (R, Vec<R>, Vec<R>) {
    let z = get_test_degree_three_z_non_scalar();
    (z[1], vec![z[0]], z[2..].to_vec())
}

#[allow(dead_code)]
pub(crate) fn get_test_degree_three_ccs<R: Ring>() -> CCS<R> {
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

#[cfg(test)]
pub(crate) fn get_test_degree_three_ccs_padded<R: Ring>(W: usize, L: usize) -> CCS<R> {
    let mut ccs = get_test_degree_three_ccs();

    ccs.m = W;
    ccs.s = log2(W) as usize;
    let len = usize::max((ccs.n - ccs.l - 1) * L, ccs.m).next_power_of_two();
    ccs.pad_rows_to(len);
    ccs
}

// Takes a vector and returns a matrix that will square the vector
pub(crate) fn create_dummy_cubing_sparse_matrix<R: Ring>(
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

    use crate::arith::{
        ccs::{get_test_degree_three_z, get_test_degree_three_z_non_scalar},
        r1cs::get_test_dummy_z_split_ntt,
        Arith, CCS,
    };

    use super::{get_test_degree_three_ccs, get_test_dummy_degree_three_ccs_non_scalar};
    type NTT = GoldilocksRingNTT;

    #[test]
    fn test_degree_three_ccs() {
        let input = 5;
        let ccs: CCS<NTT> = get_test_degree_three_ccs();
        let z = get_test_degree_three_z(input);
        assert!(ccs.check_relation(&z).is_ok())
    }

    #[test]
    fn test_degree_three_ccs_non_scalar() {
        let ccs: CCS<NTT> = get_test_degree_three_ccs();
        let z = get_test_degree_three_z_non_scalar();
        assert!(ccs.check_relation(&z).is_ok())
    }
    #[test]
    fn test_degree_three_dummy_ccs_non_scalar() {
        let (one, x_ccs, w_ccs) = get_test_dummy_z_split_ntt::<NTT, 1, 2048>();
        let mut z = vec![one];
        z.extend(&x_ccs);
        z.extend(&w_ccs);
        let ccs = get_test_dummy_degree_three_ccs_non_scalar::<NTT, 1, 2048, 2050>(&z, 1, 2050);
        assert!(ccs.check_relation(&z).is_ok())
    }
}
