use std::ops::{ AddAssign, Mul };

use lattirust_arithmetic::mle::DenseMultilinearExtension;
use lattirust_arithmetic::polynomials::VirtualPolynomial;
use lattirust_arithmetic::ring::Ring;
#[derive(Debug, Clone, PartialEq)]

// Represents a univariate polynomial
// Coefficients represented in ascending order
pub struct UnivPoly<R: Ring> {
    pub coeffs: Vec<R>,
}

impl<R: Ring> UnivPoly<R> {
    pub fn new() -> Self {
        Self { coeffs: Vec::new() }
    }

    pub fn from_virtual_polynomial(poly: VirtualPolynomial<R>) -> Self {
        let flattened_ml_extensions: Vec<DenseMultilinearExtension<R>> =
            poly.flattened_ml_extensions
                .iter()
                .map(|x| x.as_ref().clone())
                .collect();
        // Start with an empty polynomial
        let mut result_poly = UnivPoly::new();

        // Iterate over the products in the virtual polynomial
        for (coeff, list) in poly.products.iter() {
            // Start with the polynomial from the first MLE in the list
            let mut unipoly = UnivPoly::from_mle(&flattened_ml_extensions[list[0]]);

            // Multiply by subsequent MLEs
            for &index in &list[1..] {
                unipoly = unipoly * &flattened_ml_extensions[list[index]];
            }

            // Scale the polynomial by the coefficient
            unipoly = unipoly * coeff;

            // Accumulate the result
            result_poly += &unipoly;
        }
        result_poly
    }
    pub fn from_mle(mle: &DenseMultilinearExtension<R>) -> Self {
        assert!(mle.num_vars == 1, "Multilinear extension must be univariate!");
        let coeffs = vec![mle.evaluations[0], mle.evaluations[1] - mle.evaluations[0]];
        Self { coeffs }
    }

    pub fn evaluate(&self, x: R) -> R {
        let mut result = R::zero();

        // Horner's method to evaluate the polynomial
        for &coeff in self.coeffs.iter().rev() {
            result = result * x + coeff;
        }

        result
    }
    pub fn degree(&self) -> usize {
        self.coeffs
            .iter()
            .enumerate()
            .filter_map(|(i, coeff)| (!coeff.is_zero()).then(|| i))
            .last()
            .unwrap_or(0)
    }
}
impl<R: Ring> Mul<&DenseMultilinearExtension<R>> for UnivPoly<R> {
    type Output = Self;

    fn mul(self, mle: &DenseMultilinearExtension<R>) -> Self {
        assert!(mle.num_vars == 1, "Multilinear extension must be univariate!");
        let mut new_coeffs = vec![R::zero(); self.coeffs.len() + 1];
        for i in 0..self.coeffs.len() {
            new_coeffs[i] += self.coeffs[i] * mle.evaluations[0];
            new_coeffs[i + 1] += self.coeffs[i] * (mle.evaluations[1] - mle.evaluations[0]);
        }
        Self { coeffs: new_coeffs }
    }
}

impl<R: Ring> Mul<&R> for UnivPoly<R> {
    type Output = Self;

    fn mul(self, scalar: &R) -> Self {
        let new_coeffs: Vec<R> = self.coeffs
            .iter()
            .map(|&coeff| coeff * scalar)
            .collect();
        Self { coeffs: new_coeffs }
    }
}
impl<R: Ring> AddAssign<&UnivPoly<R>> for UnivPoly<R> {
    fn add_assign(&mut self, other: &UnivPoly<R>) {
        // Ensure that both polynomials have the same degree by resizing the coefficients vectors
        let max_len = std::cmp::max(self.coeffs.len(), other.coeffs.len());
        self.coeffs.resize(max_len, R::zero());
        let mut other_coeffs = other.coeffs.clone();
        other_coeffs.resize(max_len, R::zero());

        for (self_coeff, other_coeff) in self.coeffs.iter_mut().zip(other_coeffs.iter()) {
            *self_coeff += *other_coeff;
        }
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use crate::utils::sumcheck::univ_poly;

    use super::*;
    use lattirust_arithmetic::mle::DenseMultilinearExtension;
    use lattirust_arithmetic::polynomials::VirtualPolynomial;
    use lattirust_arithmetic::ring::Z2_128;

    // Define some sample DenseMultilinearExtension for testing
    fn sample_mle() -> DenseMultilinearExtension<Z2_128> {
        DenseMultilinearExtension {
            num_vars: 1,
            evaluations: vec![Z2_128::from(2u128), Z2_128::from(3u128)],
        }
    }

    // Define a sample VirtualPolynomial for testing
    fn sample_virtual_polynomial() -> VirtualPolynomial<Z2_128> {
        let mut polynomial = VirtualPolynomial::new(1);
        polynomial.flattened_ml_extensions = vec![Arc::new(sample_mle().clone()); 2];
        polynomial.products = vec![(Z2_128::from(1u128), vec![0, 1])];
        return polynomial;
    }

    #[test]
    fn test_univ_poly_from_mle() {
        let mle = sample_mle();
        let poly = UnivPoly::from_mle(&mle);
        assert_eq!(poly.coeffs, vec![Z2_128::from(2u128), Z2_128::from(1u128)]);
    }

    #[test]
    fn test_univ_poly_multiply_by_mle() {
        let mle = sample_mle();
        let poly = UnivPoly {
            coeffs: vec![Z2_128::from(1u128), Z2_128::from(1u128)],
        };
        let result = poly * &mle;
        assert_eq!(
            result.coeffs,
            vec![Z2_128::from(2u128), Z2_128::from(3u128), Z2_128::from(1u128)]
        );
    }

    #[test]
    fn test_univ_poly_multiply_by_scalar() {
        let poly = UnivPoly {
            coeffs: vec![Z2_128::from(1u128), Z2_128::from(2u128)],
        };
        let scalar = Z2_128::from(3u128);
        let result = poly * &scalar;
        assert_eq!(result.coeffs, vec![Z2_128::from(3u128), Z2_128::from(6u128)]);
    }

    #[test]
    fn test_univ_poly_add_assign() {
        let mut poly1 = UnivPoly {
            coeffs: vec![Z2_128::from(1u128), Z2_128::from(2u128)],
        };
        let poly2 = UnivPoly {
            coeffs: vec![Z2_128::from(3u128), Z2_128::from(4u128)],
        };
        poly1 += &poly2;
        assert_eq!(poly1.coeffs, vec![Z2_128::from(4u128), Z2_128::from(6u128)]);
    }

    #[test]
    fn test_univ_poly_from_virtual_polynomial() {
        let virtual_poly = sample_virtual_polynomial();
        let result = UnivPoly::from_virtual_polynomial(virtual_poly);
        assert_eq!(
            result.coeffs,
            vec![Z2_128::from(4u128), Z2_128::from(4u128), Z2_128::from(1u128)]
        );
    }

    #[test]
    fn test_univ_poly_evaluation() {
        let virtual_poly = sample_virtual_polynomial();
        let unipoly = UnivPoly::from_virtual_polynomial(virtual_poly);
        assert_eq!(unipoly.evaluate(Z2_128::from(2u128)), Z2_128::from(16u128));
        assert_eq!(unipoly.evaluate(Z2_128::from(5u128)), Z2_128::from(49u128));
    }

    #[test]
    fn test_degree() {
        let virtual_poly = sample_virtual_polynomial();
        let unipoly = UnivPoly::from_virtual_polynomial(virtual_poly);
        assert_eq!(unipoly.degree(), 2);
        let unipoly = unipoly * &sample_mle();
        assert_eq!(unipoly.degree(), 3);
    }
}
