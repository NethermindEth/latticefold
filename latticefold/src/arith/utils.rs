use lattirust_arithmetic::ring::Ring;

// Takes a matrix in the form of vec by vec and a separate vec, and multiplies them to result in a column vector
pub fn mat_by_vec<R: Ring>(lhs: &Vec<Vec<R>>, rhs: &Vec<R>) -> Vec<R> {
    let rows = lhs.len();
    let cols = lhs[0].len();
    // Initialize the result vector with zeros
    let mut result = vec![R::zero(); rows];

    // Perform matrix-vector multiplication
    for i in 0..lhs.len() {
        for j in 0..cols {
            result[i] += lhs[i][j] * rhs[j];
        }
    }
    result
}

//  Computes the hadamard product of two ring
pub fn hadamard_vec<R: Ring>(lhs: &Vec<R>, rhs: &Vec<R>) -> Vec<R> {
    lhs.iter()
        .zip(rhs.iter())
        .map(|(lhs, rhs)| *lhs * rhs)
        .collect()
}

// Multiplies Vector of rings by another ring
pub fn vec_value_mul<R: Ring>(lhs: &Vec<R>, rhs: &R) -> Vec<R> {
    lhs.iter().map(|lhs_i| *lhs_i * rhs).collect()
}

// Adds two ring vectors
pub fn vec_add<R: Ring>(lhs: &Vec<R>, rhs: &Vec<R>) -> Vec<R> {
    lhs.iter().zip(rhs.iter()).map(|(a, b)| *a + b).collect()
}
