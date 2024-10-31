#![allow(non_snake_case)]
use ark_std::UniformRand;
use cyclotomic_rings::SuitableRing;
use latticefold::arith::{r1cs::R1CS, CCS};
use lattirust_linear_algebra::SparseMatrix;
use std::fmt::Debug;

pub fn get_test_dummy_z_split<
    R: Clone + UniformRand + Debug + SuitableRing + for<'a> std::ops::AddAssign<&'a R>,
    const X_LEN: usize,
    const W: usize,
>() -> (R, Vec<R>, Vec<R>) {
    (R::one(), to_f_vec(vec![1; X_LEN]), to_f_vec(vec![1; W]))
}

pub fn get_test_dummy_ccs<
    R: Clone + UniformRand + Debug + SuitableRing + for<'a> std::ops::AddAssign<&'a R>,
    const X_LEN: usize,
    const WIT_LEN: usize,
    const W: usize,
>(
    r1cs_rows: usize,
) -> CCS<R> {
    let r1cs = get_test_dummy_r1cs::<R, X_LEN, WIT_LEN>(r1cs_rows);
    CCS::<R>::from_r1cs(r1cs, W)
}

pub fn get_test_dummy_r1cs<
    R: Clone + UniformRand + Debug + SuitableRing + for<'a> std::ops::AddAssign<&'a R>,
    const X_LEN: usize,
    const WIT_LEN: usize,
>(
    rows: usize,
) -> R1CS<R> {
    let R1CS_A = to_f_matrix::<R>(create_dummy_identity_matrix(rows, X_LEN + WIT_LEN + 1));
    let R1CS_B = R1CS_A.clone();
    let R1CS_C = R1CS_A.clone();

    R1CS::<R> {
        l: 1,
        A: R1CS_A,
        B: R1CS_B,
        C: R1CS_C,
    }
}

pub fn create_dummy_identity_matrix(rows: usize, columns: usize) -> Vec<Vec<usize>> {
    let mut matrix = vec![vec![0; columns]; rows];
    for (i, item) in matrix.iter_mut().enumerate().take(rows) {
        item[i] = 1;
    }
    matrix
}

pub fn to_f_matrix<
    R: Clone + UniformRand + Debug + SuitableRing + for<'a> std::ops::AddAssign<&'a R>,
>(
    M: Vec<Vec<usize>>,
) -> SparseMatrix<R> {
    to_f_dense_matrix::<R>(M).as_slice().into()
}
pub fn to_f_dense_matrix<
    R: Clone + UniformRand + Debug + SuitableRing + for<'a> std::ops::AddAssign<&'a R>,
>(
    M: Vec<Vec<usize>>,
) -> Vec<Vec<R>> {
    M.iter()
        .map(|m| m.iter().map(|r| R::from(*r as u64)).collect())
        .collect()
}
pub fn to_f_vec<
    R: Clone + UniformRand + Debug + SuitableRing + for<'a> std::ops::AddAssign<&'a R>,
>(
    z: Vec<usize>,
) -> Vec<R> {
    z.iter().map(|c| R::from(*c as u64)).collect()
}
