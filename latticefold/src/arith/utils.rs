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
        res[col] += *val * z[row];
    }

    Ok(res)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::Zero;
    use lattirust_arithmetic::linear_algebra::SparseMatrix;
    use lattirust_arithmetic::ring::Z2_64;

    #[test]
    fn test_hadamard_vec() {
        let a = [Z2_64::from(2u64), Z2_64::from(3u64), Z2_64::from(4u64)];
        let b = [Z2_64::from(5u64), Z2_64::from(6u64), Z2_64::from(7u64)];
        let result = hadamard_vec(&a, &b);
        let expected = vec![Z2_64::from(10u64), Z2_64::from(18u64), Z2_64::from(28u64)];
        assert_eq!(result, expected);
    }

    // Add similar tests for other functions here...

    #[test]
    fn test_vec_value_mul() {
        let a = [Z2_64::from(2u64), Z2_64::from(3u64), Z2_64::from(4u64)];
        let scalar = Z2_64::from(2u64);
        let result = vec_value_mul(&a, &scalar);
        let expected = vec![Z2_64::from(4u64), Z2_64::from(6u64), Z2_64::from(8u64)];
        assert_eq!(result, expected);
    }

    #[test]
    fn test_vec_add() {
        let a = [Z2_64::from(1u64), Z2_64::from(2u64), Z2_64::from(3u64)];
        let b = [Z2_64::from(4u64), Z2_64::from(5u64), Z2_64::from(6u64)];
        let result = vec_add(&a, &b);
        let expected = vec![Z2_64::from(5u64), Z2_64::from(7u64), Z2_64::from(9u64)];
        assert_eq!(result.unwrap(), expected);

        // Test error case
        let a = [Z2_64::from(1u64), Z2_64::from(2u64)];
        let b = [Z2_64::from(3u64), Z2_64::from(4u64), Z2_64::from(5u64)];
        let result = vec_add(&a, &b);
        assert!(result.is_err());
    }

    #[test]
    fn test_vec_scalar_mul() {
        let vec = [Z2_64::from(1u64), Z2_64::from(2u64), Z2_64::from(3u64)];
        let c = Z2_64::from(3u64);
        let result = vec_scalar_mul(&vec, &c);
        let expected = vec![Z2_64::from(3u64), Z2_64::from(6u64), Z2_64::from(9u64)];
        assert_eq!(result, expected);
    }

    #[test]
    fn test_hadamard() {
        let a = [Z2_64::from(2u64), Z2_64::from(3u64), Z2_64::from(4u64)];
        let b = [Z2_64::from(5u64), Z2_64::from(6u64), Z2_64::from(7u64)];
        let result = hadamard(&a, &b);
        let expected = vec![Z2_64::from(10u64), Z2_64::from(18u64), Z2_64::from(28u64)];
        assert_eq!(result.unwrap(), expected);

        // Test error case
        let a = [Z2_64::from(2u64), Z2_64::from(3u64)];
        let b = [Z2_64::from(5u64), Z2_64::from(6u64), Z2_64::from(7u64)];
        let result = hadamard(&a, &b);
        assert!(result.is_err());
    }

    #[test]
    fn test_mat_vec_mul() {
        // Construct a sparse matrix M
        let dense_matrix = vec![
            vec![Z2_64::from(1u64), Z2_64::zero(), Z2_64::zero()], // Row 0
            vec![Z2_64::zero(), Z2_64::from(2u64), Z2_64::from(1u64)], // Row 1
            vec![Z2_64::zero(), Z2_64::zero(), Z2_64::from(3u64)], // Row 2
        ];

        let M = SparseMatrix::from(dense_matrix.as_slice());

        let z = [Z2_64::from(1u64), Z2_64::from(1u64), Z2_64::from(1u64)];
        let result = mat_vec_mul(&M, &z);
        let expected = vec![Z2_64::from(1u64), Z2_64::from(3u64), Z2_64::from(3u64)];
        assert_eq!(result.unwrap(), expected);

        // Test error case
        let z = [Z2_64::from(1u64), Z2_64::from(1u64)]; // Wrong size vector
        let result = mat_vec_mul(&M, &z);
        assert!(result.is_err());
    }
}
