/// Some basic MLE utilities
use ark_std::log2;
use lattirust_arithmetic::linear_algebra::SparseMatrix;
use lattirust_arithmetic::mle::DenseMultilinearExtension;
use lattirust_arithmetic::mle::SparseMultilinearExtension;
use lattirust_arithmetic::ring::Ring;

/// Pad matrix so that its columns and rows are powers of two
pub fn pad_matrix<R: Ring>(m: &SparseMatrix<R>) -> SparseMatrix<R> {
    let nrows_new = m.nrows().next_power_of_two();
    let ncols_new = m.ncols().next_power_of_two();

    if m.nrows() == nrows_new && m.ncols() == ncols_new {
        m.clone()
    } else {
        SparseMatrix::try_from_csc_data(
            nrows_new,
            ncols_new,
            m.col_offsets().to_vec(),
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
            std::iter::repeat(R::zero())
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
    let n_vars: usize = (log2(n_rows) + log2(n_cols)) as usize; // n_vars = s + s'

    // build the sparse vec representing the sparse matrix
    let mut v: Vec<(usize, R)> = Vec::with_capacity(m.nnz());

    for (col, row, val) in m.triplet_iter() {
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
        std::iter::repeat(R::zero())
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

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use crate::{
//         arith::ccs::tests::get_test_z,
//         utils::multilinear_polynomial::fix_variables,
//         utils::multilinear_polynomial::tests::fix_last_variables,
//         utils::{hypercube::BooleanHypercube, vec::tests::to_F_matrix},
//     };
//     use ark_poly::MultilinearExtension;
//     use ark_std::Zero;

//     use lattirust_arithmetic::ring::{Pow2CyclotomicPolyRing, Zq};

//     #[test]
//     fn test_matrix_to_mle() {
//         let A = to_F_matrix::<Pow2CyclotomicPolyRing<Zq<17912658149728070801u64, 64>>>(vec![
//             vec![2, 3, 4, 4],
//             vec![4, 11, 14, 14],
//             vec![2, 8, 17, 17],
//             vec![420, 4, 2, 0],
//         ]);

//         let A_mle = matrix_to_mle(A);
//         assert_eq!(A_mle.evaluations.len(), 15); // 15 non-zero elements
//         assert_eq!(A_mle.num_vars, 4); // 4x4 matrix, thus 2bit x 2bit, thus 2^4=16 evals

//         let A = to_F_matrix::<Fr>(vec![
//             vec![2, 3, 4, 4, 1],
//             vec![4, 11, 14, 14, 2],
//             vec![2, 8, 17, 17, 3],
//             vec![420, 4, 2, 0, 4],
//             vec![420, 4, 2, 0, 5],
//         ]);
//         let A_mle = matrix_to_mle(A.clone());
//         assert_eq!(A_mle.evaluations.len(), 23); // 23 non-zero elements
//         assert_eq!(A_mle.num_vars, 6); // 5x5 matrix, thus 3bit x 3bit, thus 2^6=64 evals

//         // check that the A_mle evaluated over the boolean hypercube equals the matrix A_i_j values
//         let bhc = BooleanHypercube::new(A_mle.num_vars);
//         let A_padded = pad_matrix(&A);
//         let A_padded_dense = A_padded.to_dense();
//         for (i, A_row) in A_padded_dense.iter().enumerate() {
//             for (j, _) in A_row.iter().enumerate() {
//                 let s_i_j = bhc.at_i(i * A_row.len() + j);
//                 assert_eq!(A_mle.evaluate(&s_i_j).unwrap(), A_padded_dense[i][j]);
//             }
//         }
//     }

//     #[test]
//     fn test_vec_to_mle() {
//         let z = get_test_z::<Fr>(3);
//         let n_vars = 3;
//         let z_mle = dense_vec_to_mle(n_vars, &z);

//         // check that the z_mle evaluated over the boolean hypercube equals the vec z_i values
//         let bhc = BooleanHypercube::new(z_mle.num_vars);
//         for (i, z_i) in z.iter().enumerate() {
//             let s_i = bhc.at_i(i);
//             assert_eq!(z_mle.evaluate(&s_i).unwrap(), z_i.clone());
//         }
//         // for the rest of elements of the boolean hypercube, expect it to evaluate to zero
//         for i in (z.len())..(1 << z_mle.num_vars) {
//             let s_i = bhc.at_i(i);
//             assert_eq!(z_mle.evaluate(&s_i).unwrap(), Fr::zero());
//         }
//     }

//     #[test]
//     fn test_fix_variables() {
//         let A = to_F_matrix(vec![
//             vec![2, 3, 4, 4],
//             vec![4, 11, 14, 14],
//             vec![2, 8, 17, 17],
//             vec![420, 4, 2, 0],
//         ]);

//         let A_mle = matrix_to_dense_mle(A.clone());
//         let A = A.to_dense();
//         let bhc = BooleanHypercube::new(2);
//         for (i, y) in bhc.enumerate() {
//             // First check that the arkworks and espresso funcs match
//             let expected_fix_left = A_mle.fix_variables(&y); // try arkworks fix_variables
//             let fix_left = fix_variables(&A_mle, &y); // try espresso fix_variables
//             assert_eq!(fix_left, expected_fix_left);

//             // Check that fixing first variables pins down a column
//             // i.e. fixing x to 0 will return the first column
//             //      fixing x to 1 will return the second column etc.
//             let column_i: Vec<Fr> = A.clone().iter().map(|x| x[i]).collect();
//             assert_eq!(fix_left.evaluations, column_i);

//             // Now check that fixing last variables pins down a row
//             // i.e. fixing y to 0 will return the first row
//             //      fixing y to 1 will return the second row etc.
//             let row_i: Vec<Fr> = A[i].clone();
//             let fix_right = fix_last_variables(&A_mle, &y);
//             assert_eq!(fix_right.evaluations, row_i);
//         }
//     }
// }
