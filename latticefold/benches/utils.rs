use latticefold::arith::{r1cs::R1CS, CCS};
use lattirust_linear_algebra::SparseMatrix;
use lattirust_ring::PolyRing;

pub fn get_test_dummy_z_split<R: PolyRing, const IO: usize, const W: usize>() -> (R, Vec<R>, Vec<R>)
{
    (R::one(), to_F_vec(vec![1; IO]), to_F_vec(vec![1; W]))
}

pub fn get_test_dummy_ccs<R: PolyRing, const IO: usize, const W: usize>(
    r1cs_rows: usize,
) -> CCS<R> {
    let r1cs = get_test_dummy_r1cs::<R, IO, W>(r1cs_rows);
    CCS::<R>::from_r1cs(r1cs, IO + W + 1)
}

pub fn get_test_dummy_r1cs<R: PolyRing, const IO: usize, const W: usize>(rows: usize) -> R1CS<R> {
    let R1CS_A = to_F_matrix::<R>(create_dummy_matrix(rows, IO + W + 1));
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
    for i in 0..rows {
        matrix[i][i] = 1;
    }
    matrix
}

pub fn to_F_matrix<R: PolyRing>(M: Vec<Vec<usize>>) -> SparseMatrix<R> {
    to_F_dense_matrix::<R>(M).as_slice().into()
}
pub fn to_F_dense_matrix<R: PolyRing>(M: Vec<Vec<usize>>) -> Vec<Vec<R>> {
    M.iter()
        .map(|m| m.iter().map(|r| R::from(*r as u64)).collect())
        .collect()
}
pub fn to_F_vec<R: PolyRing>(z: Vec<usize>) -> Vec<R> {
    z.iter().map(|c| R::from(*c as u64)).collect()
}