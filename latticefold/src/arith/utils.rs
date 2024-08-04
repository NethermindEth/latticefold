use lattirust_arithmetic::linear_algebra::SparseMatrix;
use lattirust_arithmetic::ring::Ring;

use super::error::CSError as Error;

//  Computes the hadamard product of two ring
pub fn hadamard_vec<R: Ring>(lhs: &[R], rhs: &[R]) -> Vec<R> {
    lhs.iter().zip(rhs).map(|(lhs, rhs)| *lhs * rhs).collect()
}

// Multiplies Vector of rings by another ring
pub fn vec_value_mul<R: Ring>(lhs: &[R], rhs: &R) -> Vec<R> {
    lhs.iter().map(|lhs_i| *lhs_i * rhs).collect()
}

// Adds two ring vectors
pub fn vec_add<R: Ring>(a: &[R], b: &[R]) -> Result<Vec<R>, Error> {
    if a.len() != b.len() {
        return Err(Error::LengthsNotEqual(
            "a".to_string(),
            "b".to_string(),
            a.len(),
            b.len(),
        ));
    }
    Ok(a.iter().zip(b.iter()).map(|(x, y)| *x + y).collect())
}

pub fn vec_scalar_mul<R: Ring>(vec: &[R], c: &R) -> Vec<R> {
    vec.iter().map(|a| *a * c).collect()
}

pub fn hadamard<R: Ring>(a: &[R], b: &[R]) -> Result<Vec<R>, Error> {
    if a.len() != b.len() {
        return Err(Error::LengthsNotEqual(
            "a".to_string(),
            "b".to_string(),
            a.len(),
            b.len(),
        ));
    }
    Ok(a.iter().zip(b).map(|(a, b)| *a * b).collect())
}

pub fn mat_vec_mul<R: Ring>(M: &SparseMatrix<R>, z: &[R]) -> Result<Vec<R>, Error> {
    if M.ncols() != z.len() {
        return Err(Error::LengthsNotEqual(
            "M".to_string(),
            "z".to_string(),
            M.ncols(),
            z.len(),
        ));
    }
    let mut res = vec![R::zero(); M.nrows()];

    for (col, row, val) in M.triplet_iter() {
        res[row] += *val * z[col];
    }

    Ok(res)
}
