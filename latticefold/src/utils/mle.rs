/// Some basic MLE utilities
use ark_std::log2;
use lattirust_linear_algebra::SparseMatrix;
use lattirust_poly::mle::{DenseMultilinearExtension, SparseMultilinearExtension};
use lattirust_ring::Ring;

/// Pad matrix so that its columns and rows are powers of two
pub fn pad_matrix<R: Ring>(m: &SparseMatrix<R>) -> SparseMatrix<R> {
    let nrows_new = m.nrows().next_power_of_two();

    let ncols_new = m.ncols().next_power_of_two();

    if m.nrows() == nrows_new && m.ncols() == ncols_new {
        m.clone()
    } else {
        let offsets = pad_vector(m.col_offsets().to_vec(), ncols_new + 1);
        SparseMatrix::try_from_csc_data(
            nrows_new,
            ncols_new,
            offsets,
            m.row_indices().to_vec(),
            m.values().to_vec(),
        )
        .expect("this shouldn't have happened since we're just enlarging the matrix")
    }
}

/// Returns the dense multilinear extension from the given matrix, without modifying the original
/// matrix.
pub fn matrix_to_dense_mle<R: Ring>(matrix: SparseMatrix<R>) -> DenseMultilinearExtension<R> {
    let n_vars: usize = (log2(matrix.nrows()) + log2(matrix.ncols())) as usize; // n_vars = s + s'

    // Matrices might need to get padded before turned into an MLE
    let padded_matrix = pad_matrix(&matrix);

    // build dense vector representing the sparse padded matrix
    let mut v: Vec<R> = vec![R::zero(); padded_matrix.nrows() * padded_matrix.ncols()];

    for (col_i, row_i, val) in matrix.triplet_iter() {
        v[(padded_matrix.ncols() * row_i) + col_i] = *val;
    }

    // convert the dense vector into a mle
    vec_to_dense_mle(n_vars, &v)
}

/// Takes the n_vars and a dense vector and returns its dense MLE.
pub fn vec_to_dense_mle<R: Ring>(n_vars: usize, v: &[R]) -> DenseMultilinearExtension<R> {
    let v_padded: Vec<R> = if v.len() != (1 << n_vars) {
        // pad to 2^n_vars
        [
            v.to_owned(),
            ark_std::iter::repeat(R::zero())
                .take((1 << n_vars) - v.len())
                .collect(),
        ]
        .concat()
    } else {
        v.to_owned()
    };
    DenseMultilinearExtension::<R>::from_evaluations_vec(n_vars, v_padded)
}

/// Returns the sparse multilinear extension from the given matrix, without modifying the original
/// matrix.
pub fn matrix_to_mle<R: Ring>(m: SparseMatrix<R>) -> SparseMultilinearExtension<R> {
    let n_rows = m.nrows().next_power_of_two();
    let n_cols = m.ncols().next_power_of_two();
    let n_vars: usize = (log2(n_rows * n_cols)) as usize; // n_vars = s + s'

    // build the sparse vec representing the sparse matrix
    let mut v: Vec<(usize, R)> = Vec::with_capacity(m.nnz());

    for (col, row, val) in m.triplet_iter().filter(|(_, _, val)| !val.is_zero()) {
        v.push((row * n_cols + col, *val));
    }

    // convert the dense vector into a mle
    vec_to_mle(n_vars, &v)
}

/// Takes the n_vars and a sparse vector and returns its sparse MLE.
pub fn vec_to_mle<R: Ring>(n_vars: usize, v: &[(usize, R)]) -> SparseMultilinearExtension<R> {
    SparseMultilinearExtension::<R>::from_evaluations(n_vars, v)
}

/// Takes the n_vars and a dense vector and returns its dense MLE.
pub fn dense_vec_to_dense_mle<R: Ring>(n_vars: usize, v: &[R]) -> DenseMultilinearExtension<R> {
    // Pad to 2^n_vars
    let v_padded: Vec<R> = [
        v.to_owned(),
        ark_std::iter::repeat(R::zero())
            .take((1 << n_vars) - v.len())
            .collect(),
    ]
    .concat();
    DenseMultilinearExtension::<R>::from_evaluations_vec(n_vars, v_padded)
}

/// Takes the n_vars and a dense vector and returns its sparse MLE.
pub fn dense_vec_to_mle<R: Ring>(n_vars: usize, v: &[R]) -> SparseMultilinearExtension<R> {
    let v_sparse = v
        .iter()
        .enumerate()
        .map(|(i, v_i)| (i, *v_i))
        .collect::<Vec<(usize, R)>>();
    SparseMultilinearExtension::<R>::from_evaluations(n_vars, &v_sparse)
}

// Pad a vector by repeating last element
fn pad_vector(mut vec: Vec<usize>, target_length: usize) -> Vec<usize> {
    if vec.len() < target_length {
        if let Some(&last_element) = vec.last() {
            vec.resize(target_length, last_element);
        }
    }
    vec
}

#[cfg(test)]
#[allow(non_snake_case)]
mod tests {
    use super::*;
    use crate::arith::{r1cs::tests::to_F_matrix, tests::get_test_z};

    use ark_ff::Zero;
    use lattirust_poly::mle::MultilinearExtension;
    use lattirust_ring::Z2_64;

    // Function to convert usize to a binary vector of Ring elements.
    fn usize_to_binary_vector<R: Ring>(n: usize, dimensions: usize) -> Vec<R> {
        let mut bits = Vec::with_capacity(dimensions);
        let mut current = n;

        for _ in 0..dimensions {
            if (current & 1) == 1 {
                bits.push(R::one());
            } else {
                bits.push(R::zero());
            }
            current >>= 1;
        }
        bits
    }

    // Wrapper function to generate a boolean hypercube.
    fn boolean_hypercube<R: Ring>(dimensions: usize) -> Vec<Vec<R>> {
        let max_val = 1 << dimensions; // 2^dimensions
        (0..max_val)
            .map(|i| usize_to_binary_vector::<R>(i, dimensions))
            .collect()
    }
    #[test]
    fn test_matrix_to_mle() {
        type R = Z2_64;
        let A = to_F_matrix::<R>(vec![
            vec![2, 3, 4, 4],
            vec![4, 11, 14, 14],
            vec![2, 8, 17, 17],
            vec![420, 4, 2, 0],
        ]);

        let A_mle = matrix_to_mle(A);
        assert_eq!(A_mle.evaluations.len(), 15); // 15 non-zero elements
        assert_eq!(A_mle.num_vars, 4); // 4x4 matrix, thus 2bit x 2bit, thus 2^4=16 evals

        let A = to_F_matrix::<R>(vec![
            vec![2, 3, 4, 4, 1],
            vec![4, 11, 14, 14, 2],
            vec![2, 8, 17, 17, 3],
            vec![420, 4, 2, 0, 4],
            vec![420, 4, 2, 0, 5],
        ]);
        let A_mle = matrix_to_mle(A.clone());
        assert_eq!(A_mle.evaluations.len(), 23); // 23 non-zero elements
        assert_eq!(A_mle.num_vars, 6); // 5x5 matrix, thus 3bit x 3bit, thus 2^6=64 evals

        // check that the A_mle evaluated over the boolean hypercube equals the matrix A_i_j values
        let _bhc = boolean_hypercube::<R>(2);

        let _A_dense = matrix_to_dense_mle(A);
    }

    #[test]
    fn test_vec_to_mle() {
        type R = Z2_64;
        let z = get_test_z::<R>(3);
        let n_vars = 3;
        let z_mle = dense_vec_to_mle(n_vars, &z);

        // check that the z_mle evaluated over the boolean hypercube equals the vec z_i values
        let bhc = boolean_hypercube(z_mle.num_vars);

        for (i, z_i) in z.iter().enumerate() {
            let s_i = &bhc[i];
            assert_eq!(z_mle.evaluate(s_i), z_i.clone());
        }
        // for the rest of elements of the boolean hypercube, expect it to evaluate to zero
        for s_i in bhc.iter().take(1 << z_mle.num_vars).skip(z.len()) {
            assert_eq!(z_mle.fix_variables(s_i)[0], R::zero());
        }
    }
}
