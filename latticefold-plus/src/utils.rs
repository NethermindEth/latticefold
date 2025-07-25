use ark_std::{
    ops::{Mul, Sub},
    One, Zero,
};
use latticefold::transcript::Transcript;
use stark_rings::{
    balanced_decomposition::{Decompose, GadgetDecompose},
    OverField, PolyRing,
};
use stark_rings_linalg::Matrix;

pub fn split<R: Decompose + PolyRing>(
    com: &Matrix<R>,
    n: usize,
    b: u128,
    k: usize,
) -> Vec<R::BaseRing> {
    let M_prime = com.gadget_decompose(b, k);
    let M_dprime = M_prime.vals.into_iter().fold(vec![], |mut acc, row| {
        // TODO pre-alloc
        acc.extend(row);
        acc
    });
    let mut tau = M_dprime
        .iter()
        .map(|r| r.coeffs().to_vec())
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

pub fn short_challenge<R: OverField>(lambda: usize, transcript: &mut impl Transcript<R>) -> R {
    let u = 2usize.pow(lambda as u32 / R::dimension() as u32);
    let bytes = transcript.squeeze_bytes(R::dimension());

    let coeffs = bytes
        .iter()
        .map(|b| {
            // rough centering
            R::BaseRing::from((*b as usize % u) as u128) - R::BaseRing::from((u / 2) as u128)
        })
        .collect::<Vec<R::BaseRing>>();

    coeffs.into()
}

pub fn estimate_bound(sop: usize, L: usize, d: usize, k: usize) -> u128 {
    let a = sop * L;
    let c = d / 2 + d * k + 1;

    let discriminant = (a * a + 4 * a * c) as f64;
    let sqrt_discriminant = discriminant.sqrt();

    let b = (a as f64 + sqrt_discriminant) / 2.0;

    (b.ceil()) as u128
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
