use std::marker::PhantomData;
use std::ops::AddAssign;

use ark_ff::PrimeField;
use lattirust_arithmetic::challenge_set::latticefold_challenge_set::OverField;
use lattirust_arithmetic::polynomials::VirtualPolynomial;
use lattirust_arithmetic::mle::DenseMultilinearExtension;
pub struct UnivPoly<F: PrimeField, R: OverField<F>> {
    coeffs: Vec<R>,
    _marker: PhantomData<F>,
}

impl<F: PrimeField, R: OverField<F>> UnivPoly<F, R> {
    pub fn new() -> Self {
        Self {
            coeffs: Vec::new(),
            _marker: PhantomData::default(),
        }
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
                unipoly = unipoly.multiply_by_mle(&flattened_ml_extensions[list[index]]);
            }

            // Scale the polynomial by the coefficient
            unipoly = unipoly.multiply_by_scalar(*coeff);

            // Accumulate the result
            result_poly += &unipoly;
        }

        result_poly
    }
    pub fn from_mle(mle: &DenseMultilinearExtension<R>) -> Self {
        assert!(mle.num_vars == 1, "Multilinear extension must be univariate!");
        let coeffs = vec![mle.evaluations[0], mle.evaluations[1] - mle.evaluations[0]];
        Self {
            coeffs,
            _marker: PhantomData::default(),
        }
    }

    pub fn multiply_by_mle(self, mle: &DenseMultilinearExtension<R>) -> Self {
        assert!(mle.num_vars == 1, "Multilinear extension must be univariate!");
        let mut new_coeffs = vec![R::zero(); self.coeffs.len() + 1];
        for i in 0..self.coeffs.len() {
            new_coeffs[i] += self.coeffs[i] * mle.evaluations[0];
            new_coeffs[i + 1] += self.coeffs[i] * (mle.evaluations[1] - mle.evaluations[0]);
        }
        Self {
            coeffs: new_coeffs,
            _marker: PhantomData::default(),
        }
    }

    pub fn multiply_by_scalar(self, scalar: R) -> Self {
        let new_coeffs: Vec<R> = self.coeffs
            .iter()
            .map(|&coeff| coeff * scalar)
            .collect();
        Self {
            coeffs: new_coeffs,
            _marker: PhantomData::default(),
        }
    }
}

impl<F: PrimeField, R: OverField<F>> AddAssign<&UnivPoly<F, R>> for UnivPoly<F, R> {
    fn add_assign(&mut self, other: &UnivPoly<F, R>) {
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
