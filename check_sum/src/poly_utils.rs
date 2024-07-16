// ! Utilities for reasoning with virtual polynomials

use lattirust_arithmetic::polynomials::VirtualPolynomial;
use lattirust_arithmetic::ring::Ring;

pub fn degree_by_index<R: Ring>(polynomial: &VirtualPolynomial<R>, i: &usize) {
    todo!()
}

/// Compute the partial evaluation of `poly` by the specified field element and variable.
pub fn partial_eval<R: Ring>(
    polynomial: &VirtualPolynomial<R>,
    value: R,
    variable: usize
) -> VirtualPolynomial<R> {
    todo!()
}
/// Restrict a multivariate polynomial to the univariate polynomial obtained by
/// evaluating at zero and one for all variables except the one specified, then summing
/// This is computed at every round of the sum check protocol
pub fn to_univariate_sum<R: Ring>(
    polynomial: &VirtualPolynomial<R>,
    variable: usize
) -> VirtualPolynomial<R> {
    todo!()
}

#[cfg(test)]
mod tests {}
