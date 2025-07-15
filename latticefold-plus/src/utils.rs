use ark_std::{
    iter::once,
    log2,
    ops::{Mul, Sub},
    One, Zero,
};
use latticefold::{
    transcript::Transcript,
    utils::sumcheck::{
        utils::{build_eq_x_r, eq_eval},
        MLSumcheck, Proof, SumCheckError,
    },
};
// Import the concrete ring type for testing
#[cfg(test)]
use stark_rings::cyclotomic_ring::models::goldilocks::RqPoly as TestRing;
use stark_rings::{
    balanced_decomposition::{
        convertible_ring::ConvertibleRing, Decompose, DecomposeToVec, GadgetDecompose,
    },
    exp, psi, psi_range_check, CoeffRing, OverField, PolyRing, Ring, Zq,
};
use stark_rings_linalg::{ops::Transpose, Matrix, SparseMatrix};
use stark_rings_poly::mle::{DenseMultilinearExtension, SparseMultilinearExtension};

pub fn split<R: PolyRing>(com: &Matrix<R>, n: usize, b: u128, k: usize) -> Vec<R::BaseRing>
where
    R: Decompose,
{
    let M_prime = com.gadget_decompose(b, k);
    let M_dprime = M_prime.vals.into_iter().fold(vec![], |mut acc, row| {
        // TODO pre-alloc
        acc.extend(row);
        acc
    });
    let mut tau = M_dprime
        .iter()
        .map(|r| r.coeffs().to_vec())
        .into_iter()
        .fold(vec![], |mut acc, row| {
            // TODO pre-alloc
            acc.extend(row);
            acc
        });
    if tau.len() < n {
        // TODO handle when opposite
        tau.resize(n, R::BaseRing::zero());
    } else {
        panic!(
            "small n {} unsupported, must be >= tau unpadded {}",
            n,
            tau.len()
        );
    }
    tau
}

/// Computes the tensor product of two flat vectors.
///
/// If `a` has length `m` and `b` has length `n`, the result is a new vector
/// of length `m * n` containing the element-wise products.
pub fn tensor_product<T>(a: &[T], b: &[T]) -> Vec<T>
where
    T: Clone + Mul<Output = T>,
{
    if a.is_empty() {
        return b.to_vec();
    }
    if b.is_empty() {
        return a.to_vec();
    }

    let mut result = Vec::with_capacity(a.len() * b.len());
    for a_val in a {
        for b_val in b {
            result.push(a_val.clone() * b_val.clone());
        }
    }
    result
}

/// Computes the tensor operation on a vector `r`.
///
/// This corresponds to the `tensor(r)` function, defined as the sequential
/// tensor product of `(1 - r_i, r_i)` for each element `r_i` in the input vector.
pub fn tensor<T>(r: &[T]) -> Vec<T>
where
    T: Clone + One + Sub<Output = T> + Mul<Output = T>,
{
    let mut result = vec![T::one()];

    for r_i in r {
        let one = T::one();
        let term = [one - r_i.clone(), r_i.clone()];
        result = tensor_product(&result, &term);
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tensor_product() {
        let a = vec![1, 2];
        let b = vec![10, 20, 30];
        let expected = vec![10, 20, 30, 20, 40, 60];
        assert_eq!(tensor_product(&a, &b), expected);
    }

    #[test]
    fn test_tensor() {
        let r = vec![10, 2];
        let expected = vec![-9 * -1, -9 * 2, 10 * -1, 10 * 2];
        assert_eq!(tensor(&r), expected);
    }
}
