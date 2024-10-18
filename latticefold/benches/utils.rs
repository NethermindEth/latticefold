use ark_std::UniformRand;
use cyclotomic_rings::SuitableRing;
use latticefold::arith::{r1cs::R1CS, CCS};
use lattirust_linear_algebra::SparseMatrix;
use lattirust_ring::PolyRing;
use std::fmt::Debug;

pub fn get_test_dummy_z_split<
    R: Clone + UniformRand + Debug + SuitableRing + for<'a> std::ops::AddAssign<&'a R>,
    const IO: usize,
    const W: usize,
>() -> (R, Vec<R>, Vec<R>) {
    (R::one(), to_f_vec(vec![1; IO]), to_f_vec(vec![1; W]))
}

pub fn get_test_dummy_ccs<
    R: Clone + UniformRand + Debug + SuitableRing + for<'a> std::ops::AddAssign<&'a R>,
    const IO: usize,
    const W: usize,
>(
    r1cs_rows: usize,
) -> CCS<R> {
    let r1cs = get_test_dummy_r1cs::<R, IO, W>(r1cs_rows);
    CCS::<R>::from_r1cs(r1cs, IO + W + 1)
}

pub fn get_test_dummy_r1cs<
    R: Clone + UniformRand + Debug + SuitableRing + for<'a> std::ops::AddAssign<&'a R>,
    const IO: usize,
    const W: usize,
>(
    rows: usize,
) -> R1CS<R> {
    let R1CS_A = to_f_matrix::<R>(create_dummy_matrix(rows, IO + W + 1));
    let R1CS_B = R1CS_A.clone();
    let R1CS_C = R1CS_A.clone();

    R1CS::<R> {
        l: 1,
        A: R1CS_A,
        B: R1CS_B,
        C: R1CS_C,
    }
}

pub fn create_dummy_matrix(rows: usize, columns: usize) -> Vec<Vec<usize>> {
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
