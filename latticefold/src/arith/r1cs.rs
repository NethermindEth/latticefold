//! Defines behaviour of R1CS, a degree two constraint system

use cyclotomic_rings::rings::SuitableRing;
use lattirust_linear_algebra::sparse_matrix::dense_matrix_u64_to_sparse;
use lattirust_linear_algebra::SparseMatrix;
use lattirust_ring::Ring;

use super::{
    error::CSError as Error,
    utils::{mat_vec_mul, vec_add, vec_scalar_mul},
};
use crate::arith::hadamard;
use crate::ark_base::*;

/// a representation of a R1CS instance
#[derive(Debug, Clone, PartialEq)]
pub struct R1CS<R: Ring> {
    /// Length of public input
    pub l: usize,
    /// First constraint matrix
    pub A: SparseMatrix<R>,
    /// Second constraint matrix
    pub B: SparseMatrix<R>,
    /// Third constraint matrix
    pub C: SparseMatrix<R>,
}

impl<R: Ring> R1CS<R> {
    /// check that a R1CS structure is satisfied by a z vector. Only for testing.
    #[cfg(test)]
    pub(crate) fn check_relation(&self, z: &[R]) -> Result<(), Error> {
        let Az = mat_vec_mul(&self.A, z)?;
        let Bz = mat_vec_mul(&self.B, z)?;

        let Cz = mat_vec_mul(&self.C, z)?;
        let AzBz = hadamard(&Az, &Bz)?;

        if AzBz != Cz {
            Err(Error::NotSatisfied)
        } else {
            Ok(())
        }
    }
    /// Converts the R1CS instance into a RelaxedR1CS as described in
    /// [Nova](https://eprint.iacr.org/2021/370.pdf#page=14)
    pub fn relax(self) -> RelaxedR1CS<R> {
        RelaxedR1CS::<R> {
            l: self.l,
            E: vec![R::zero(); self.A.nrows()],
            A: self.A,
            B: self.B,
            C: self.C,
            u: R::one(),
        }
    }
}

/// A RelaxedR1CS instance as described in
/// [Nova](https://eprint.iacr.org/2021/370.pdf#page=14).
///
/// A witness $z$ is satisfying if $(A \cdot z) \circ (B \cdot z) = u \cdot (C \cdot z) + E$.
#[derive(Debug, Clone, PartialEq)]
pub struct RelaxedR1CS<R: Ring> {
    /// Public input length
    pub l: usize,
    /// First constraint matrix
    pub A: SparseMatrix<R>,
    /// Second constraint matrix
    pub B: SparseMatrix<R>,
    /// Third constraint matrix
    pub C: SparseMatrix<R>,
    /// Scalar coefficient of $(C \cdot z)$
    pub u: R,
    /// The error matrix
    pub E: Vec<R>,
}

impl<R: Ring> RelaxedR1CS<R> {
    /// check that a RelaxedR1CS structure is satisfied by a z vector.
    pub fn check_relation(&self, z: &[R]) -> Result<(), Error> {
        let Az = mat_vec_mul(&self.A, z)?;
        let Bz = mat_vec_mul(&self.B, z)?;
        let Cz = mat_vec_mul(&self.C, z)?;

        let uCz = vec_scalar_mul(&Cz, &self.u);
        let uCzE = vec_add(&uCz, &self.E)?;
        let AzBz = hadamard(&Az, &Bz)?;
        if AzBz != uCzE {
            Err(Error::NotSatisfied)
        } else {
            Ok(())
        }
    }
}

/// Returns a matrix of ring elements given a matrix of unsigned ints
pub fn to_F_matrix<R: Ring>(M: Vec<Vec<usize>>) -> SparseMatrix<R> {
    // dense_matrix_to_sparse(to_F_dense_matrix::<R>(M))
    let M_u64: Vec<Vec<u64>> = M
        .iter()
        .map(|m| m.iter().map(|r| *r as u64).collect())
        .collect();
    dense_matrix_u64_to_sparse(M_u64)
}

/// Returns a dense matrix of ring elements given a matrix of unsigned ints
pub fn to_F_dense_matrix<R: Ring>(M: Vec<Vec<usize>>) -> Vec<Vec<R>> {
    M.iter()
        .map(|m| m.iter().map(|r| R::from(*r as u64)).collect())
        .collect()
}

/// Returns a vector of ring elements given a vector of unsigned ints
pub fn to_F_vec<R: Ring>(z: Vec<usize>) -> Vec<R> {
    z.iter().map(|c| R::from(*c as u64)).collect()
}

#[cfg(test)]
pub(crate) fn get_test_r1cs<R: Ring>() -> R1CS<R> {
    // R1CS for: x^3 + x + 5 = y (example from article
    // https://www.vitalik.ca/general/2016/12/10/qap.html )
    let A = to_F_matrix::<R>(vec![
        vec![1, 0, 0, 0, 0, 0],
        vec![0, 0, 0, 1, 0, 0],
        vec![1, 0, 0, 0, 1, 0],
        vec![0, 5, 0, 0, 0, 1],
    ]);
    let B = to_F_matrix::<R>(vec![
        vec![1, 0, 0, 0, 0, 0],
        vec![1, 0, 0, 0, 0, 0],
        vec![0, 1, 0, 0, 0, 0],
        vec![0, 1, 0, 0, 0, 0],
    ]);
    let C = to_F_matrix::<R>(vec![
        vec![0, 0, 0, 1, 0, 0],
        vec![0, 0, 0, 0, 1, 0],
        vec![0, 0, 0, 0, 0, 1],
        vec![0, 0, 1, 0, 0, 0],
    ]);

    R1CS::<R> { l: 1, A, B, C }
}

/// Return a R1CS instance of arbitrary size, useful for benching.
/// Only works when z vector consists of multiplicative identities.
pub fn get_test_dummy_r1cs<R: Ring, const X_LEN: usize, const WIT_LEN: usize>(
    rows: usize,
) -> R1CS<R> {
    let R1CS_A = create_dummy_identity_sparse_matrix(rows, X_LEN + WIT_LEN + 1);
    let R1CS_B = R1CS_A.clone();
    let R1CS_C = R1CS_A.clone();

    R1CS::<R> {
        l: 1,
        A: R1CS_A,
        B: R1CS_B,
        C: R1CS_C,
    }
}

/// Return a R1CS instance of arbitrary size, useful for benching.
/// Works for arbitrary z vector.
pub fn get_test_dummy_r1cs_non_scalar<R: Ring, const X_LEN: usize, const WIT_LEN: usize>(
    rows: usize,
    witness: &[R],
) -> R1CS<R> {
    let R1CS_A = create_dummy_identity_sparse_matrix(rows, X_LEN + WIT_LEN + 1);
    let R1CS_B = R1CS_A.clone();
    let R1CS_C = create_dummy_squaring_sparse_matrix(rows, X_LEN + WIT_LEN + 1, witness);

    R1CS::<R> {
        l: 1,
        A: R1CS_A,
        B: R1CS_B,
        C: R1CS_C,
    }
}

pub(crate) fn create_dummy_identity_sparse_matrix<R: Ring>(
    rows: usize,
    columns: usize,
) -> SparseMatrix<R> {
    let mut matrix = SparseMatrix {
        n_rows: rows,
        n_cols: columns,
        coeffs: vec![vec![]; rows],
    };
    for (i, row) in matrix.coeffs.iter_mut().enumerate() {
        row.push((R::one(), i));
    }
    matrix
}

// Takes a vector and returns a matrix that will square the vector
pub(crate) fn create_dummy_squaring_sparse_matrix<R: Ring>(
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
        row.push((witness[i], i));
    }
    matrix
}

pub(crate) fn get_test_z<R: Ring>(input: usize) -> Vec<R> {
    // z = (io, 1, w)
    to_F_vec(vec![
        input, // io
        1,
        input * input * input + input + 5, // x^3 + x + 5
        input * input,                     // x^2
        input * input * input,             // x^2 * x
        input * input * input + input,     // x^3 + x
    ])
}

pub(crate) fn get_test_z_ntt<R: SuitableRing>() -> Vec<R> {
    let mut res = Vec::new();
    for input in 0..R::dimension() {
        // z = (io, 1, w)
        res.push(to_F_vec::<R::BaseRing>(vec![
            input, // io
            1,
            input * input * input + input + 5, // x^3 + x + 5
            input * input,                     // x^2
            input * input * input,             // x^2 * x
            input * input * input + input,     // x^3 + x
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

/// Return scalar z vector for Vitalik's [R1CS example](https://medium.com/@VitalikButerin/quadratic-arithmetic-programs-from-zero-to-hero-f6d558cea649#81e4),
/// split into statement, constant, and witness.
pub fn get_test_z_split<R: Ring>(input: usize) -> (R, Vec<R>, Vec<R>) {
    let z = get_test_z(input);
    (z[1], vec![z[0]], z[2..].to_vec())
}

/// Return non-scalar z vector for Vitalik's [R1CS example](https://medium.com/@VitalikButerin/quadratic-arithmetic-programs-from-zero-to-hero-f6d558cea649#81e4),
/// split into statement, constant, and witness.
pub fn get_test_z_ntt_split<R: SuitableRing>() -> (R, Vec<R>, Vec<R>) {
    let z = get_test_z_ntt();
    (z[1], vec![z[0]], z[2..].to_vec())
}

/// Return z vector consisting only of multiplicative identities,
/// split into statement, constant, and witness.
pub fn get_test_dummy_z_split<R: Ring, const X_LEN: usize, const WIT_LEN: usize>(
) -> (R, Vec<R>, Vec<R>) {
    (
        R::one(),
        to_F_vec(vec![1; X_LEN]),
        to_F_vec(vec![1; WIT_LEN]),
    )
}

/// Return z vector consisting of non scalar ring elements,
/// split into statement, constant, and witness.
pub fn get_test_dummy_z_split_ntt<R: SuitableRing, const X_LEN: usize, const WIT_LEN: usize>(
) -> (R, Vec<R>, Vec<R>) {
    let statement_vec = (0..X_LEN).map(|_| R::one()).collect();

    let witness_vec = (0..WIT_LEN)
        .map(|_| {
            R::from(
                (0..R::dimension())
                    .map(|i| R::BaseRing::from(i as u128))
                    .collect::<Vec<R::BaseRing>>(),
            )
        })
        .collect();

    (R::one(), statement_vec, witness_vec)
}

#[cfg(test)]
mod tests {
    use cyclotomic_rings::rings::FrogRingNTT;

    use super::*;

    #[test]
    fn test_check_relation() {
        let r1cs = get_test_r1cs::<FrogRingNTT>();
        let z = get_test_z(5);

        r1cs.check_relation(&z).unwrap();
        r1cs.relax().check_relation(&z).unwrap();
    }
}
