use ark_std::ops::{AddAssign, Mul};
use lattirust_poly::{
    mle::DenseMultilinearExtension,
    polynomials::{ArithErrors, VirtualPolynomial},
};
use lattirust_ring::Ring;

// Represents a univariate polynomial
// Coefficients represented in ascending order
#[derive(Debug, Clone, PartialEq)]
pub struct UVPolynomial<R: Ring> {
    pub coeffs: Vec<R>,
}

impl<R: Ring> Default for UVPolynomial<R> {
    fn default() -> Self {
        Self::new()
    }
}

impl<R: Ring> UVPolynomial<R> {
    pub fn new() -> Self {
        Self { coeffs: Vec::new() }
    }

    pub fn evaluate(&self, x: R) -> R {
        self.coeffs
            .iter()
            .rev()
            .fold(R::zero(), |result, coeff| result * x + coeff)
    }

    pub fn degree(&self) -> usize {
        self.coeffs
            .iter()
            .enumerate()
            .rev()
            .filter_map(|(i, coeff)| (!coeff.is_zero()).then_some(i))
            .next()
            .unwrap_or(0)
    }
}

impl<R: Ring> TryFrom<&DenseMultilinearExtension<R>> for UVPolynomial<R> {
    type Error = ArithErrors;
    fn try_from(mle: &DenseMultilinearExtension<R>) -> Result<Self, ArithErrors> {
        assert!(
            mle.num_vars == 1,
            "Multilinear extension must be univariate!"
        );
        let coeffs = vec![mle.evaluations[0], mle.evaluations[1] - mle.evaluations[0]];
        Ok(Self { coeffs })
    }
}

impl<R: Ring> TryFrom<VirtualPolynomial<R>> for UVPolynomial<R> {
    type Error = ArithErrors;
    fn try_from(poly: VirtualPolynomial<R>) -> Result<Self, ArithErrors> {
        let flattened_ml_extensions: Vec<DenseMultilinearExtension<R>> = poly
            .flattened_ml_extensions
            .iter()
            .map(|x| x.as_ref().clone())
            .collect();
        // Start with an empty polynomial
        let mut result_poly = UVPolynomial::new();

        // Iterate over the products in the virtual polynomial
        for (coeff, list) in poly.products.iter() {
            // Start with the polynomial from the first MLE in the list
            let mut unipoly = UVPolynomial::try_from(&flattened_ml_extensions[list[0]])?;

            for &index in &list[1..] {
                unipoly = unipoly * &flattened_ml_extensions[index];
            }

            // Scale the polynomial by the coefficient
            unipoly = unipoly * coeff;

            // Accumulate the result
            result_poly += &unipoly;
        }
        Ok(result_poly)
    }
}

impl<R: Ring> Mul<&DenseMultilinearExtension<R>> for UVPolynomial<R> {
    type Output = Self;

    fn mul(self, mle: &DenseMultilinearExtension<R>) -> Self {
        assert!(
            mle.num_vars == 1,
            "Multilinear extension must be univariate!"
        );

        Self {
            coeffs: self.coeffs.iter().enumerate().fold(
                vec![R::zero(); self.coeffs.len() + 1],
                |mut new_coeffs, (i, coeff)| {
                    new_coeffs[i] += *coeff * mle.evaluations[0];
                    new_coeffs[i + 1] += *coeff * (mle.evaluations[1] - mle.evaluations[0]);
                    new_coeffs
                },
            ),
        }
    }
}

impl<R: Ring> Mul<&R> for UVPolynomial<R> {
    type Output = Self;

    fn mul(self, scalar: &R) -> Self {
        let new_coeffs: Vec<R> = self.coeffs.iter().map(|&coeff| coeff * scalar).collect();
        Self { coeffs: new_coeffs }
    }
}

impl<R: Ring> AddAssign<&UVPolynomial<R>> for UVPolynomial<R> {
    fn add_assign(&mut self, other: &UVPolynomial<R>) {
        // Ensure that both polynomials have the same degree by resizing the coefficients vectors
        let max_len = ark_std::cmp::max(self.coeffs.len(), other.coeffs.len());
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
    use ark_std::sync::Arc;

    use super::*;
    use lattirust_poly::{mle::DenseMultilinearExtension, polynomials::VirtualPolynomial};
    use lattirust_ring::cyclotomic_ring::models::goldilocks::Fq;

    // Define some sample DenseMultilinearExtension for testing
    fn sample_mle() -> DenseMultilinearExtension<Fq> {
        DenseMultilinearExtension {
            num_vars: 1,
            evaluations: vec![Fq::from(2u128), Fq::from(3u128)],
        }
    }

    // Define a sample VirtualPolynomial for testing
    fn sample_virtual_polynomial() -> VirtualPolynomial<Fq> {
        let mut polynomial = VirtualPolynomial::new(1);
        polynomial.flattened_ml_extensions = (0..2).map(|_| Arc::new(sample_mle())).collect();
        polynomial.products = vec![(Fq::from(1u128), vec![0, 1])];
        polynomial
    }

    #[test]
    fn test_univ_poly_from_mle() {
        let mle = sample_mle();
        let poly = UVPolynomial::try_from(&mle);
        assert_eq!(
            poly.unwrap().coeffs,
            vec![Fq::from(2u128), Fq::from(1u128)]
        );
    }

    #[test]
    fn test_univ_poly_multiply_by_mle() {
        let mle = sample_mle();
        let poly = UVPolynomial {
            coeffs: vec![Fq::from(1u128), Fq::from(1u128)],
        };
        let result = poly * &mle;
        assert_eq!(
            result.coeffs,
            vec![
                Fq::from(2u128),
                Fq::from(3u128),
                Fq::from(1u128)
            ]
        );
    }

    #[test]
    fn test_univ_poly_multiply_by_scalar() {
        let poly = UVPolynomial {
            coeffs: vec![Fq::from(1u128), Fq::from(2u128)],
        };
        let scalar = Fq::from(3u128);
        let result = poly * &scalar;
        assert_eq!(
            result.coeffs,
            vec![Fq::from(3u128), Fq::from(6u128)]
        );
    }

    #[test]
    fn test_univ_poly_add_assign() {
        let mut poly1 = UVPolynomial {
            coeffs: vec![Fq::from(1u128), Fq::from(2u128)],
        };
        let poly2 = UVPolynomial {
            coeffs: vec![Fq::from(3u128), Fq::from(4u128)],
        };
        poly1 += &poly2;
        assert_eq!(poly1.coeffs, vec![Fq::from(4u128), Fq::from(6u128)]);
    }

    #[test]
    fn test_univ_poly_from_virtual_polynomial() {
        let virtual_poly = sample_virtual_polynomial();
        let result = UVPolynomial::try_from(virtual_poly);
        assert_eq!(
            result.unwrap().coeffs,
            vec![
                Fq::from(4u128),
                Fq::from(4u128),
                Fq::from(1u128)
            ]
        );
    }

    #[test]
    fn test_univ_poly_evaluation() {
        let virtual_poly = sample_virtual_polynomial();
        let unipoly = UVPolynomial::try_from(virtual_poly);
        assert_eq!(
            unipoly.unwrap().evaluate(Fq::from(2u128)),
            Fq::from(16u128)
        );
    }

    #[test]
    fn test_degree() {
        let virtual_poly = sample_virtual_polynomial();
        let unipoly = UVPolynomial::try_from(virtual_poly);
        assert_eq!(&unipoly.unwrap().degree(), &2);
    }
}
