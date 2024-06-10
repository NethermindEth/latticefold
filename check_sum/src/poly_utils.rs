// ! Utilities for reasoning with multi and univariate polynomials

use lattirust_arithmetic::ring::Ring;

// Represents a univariate Polynomial
// Represented by a list of coefficients of ascending powers
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct UnivPoly<R: Ring> {
    coeffs: Vec<R>,
}

// Represents a multivariate polynomial
// Represented by a vector of terms
// Each term has a coefficient and a vector corresponding
// To the degree of each variable in that term
//  See tests for examples
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MultiPoly<R: Ring> {
    pub terms: Vec<(R, Vec<(usize, usize)>)>,
}

impl<R: Ring> UnivPoly<R> {
    // Compute the degree of the polynomial
    pub fn degree(&self) -> usize {
        self.coeffs.len() - 1
    }

    // Takes a value and evaluates the polynomial
    pub fn eval(&self, r: &R) -> R {
        let mut result = R::zero();

        for (i, coeff) in self.coeffs.iter().enumerate() {
            result += r.pow(&[i as u64]) * coeff;
        }

        result
    }
    //  Pretty prints a univariate polynomial
    pub fn format_univ_poly(&self, varname: &str) -> String {
        if self.coeffs.is_empty() || !self.coeffs.iter().any(|coeff| *coeff != R::zero()) {
            String::from("0")
        } else {
            self.coeffs
                .iter()
                .enumerate()
                .rev()
                .filter(|(_, coeff)| **coeff != R::zero())
                .map(|(exp, coeff)| {
                    match (exp, coeff) {
                        (0, _) => format!("{}", coeff),
                        (1, coeff) if *coeff == R::one() => format!("{}", varname),
                        (1, _) => format!("{}*{}", coeff, varname),
                        (_, coeff) if *coeff == R::one() => format!("{}^{}", varname, exp),
                        (_, _) => format!("{}*{}^{}", coeff, varname, exp),
                    }
                })
                .collect::<Vec<_>>()
                .join(" + ")
        }
    }
}

/// Compute the sums of k-powers of the list of summands for k up to `max_exponent`.
fn power_sums<R: Ring>(max_exponent: usize, summands: &[R]) -> Vec<R> {
    let mut powers: Vec<R> = vec![R::one(); summands.len()];
    let mut power_sums: Vec<R> = Vec::with_capacity(max_exponent + 1);

    for _ in 0..=max_exponent {
        power_sums.push(powers.iter().sum());

        for (power, term) in powers.iter_mut().zip(summands.iter()) {
            *power = *power * *term;
        }
    }
    power_sums
}

impl<R: Ring> MultiPoly<R> {
    // Takes a list of values and evaluates the multivariate polynomial
    pub fn eval_poly(&self, values: &Vec<R>) -> R {
        self.terms
            .iter()
            .map(|(coeff, term)| {
                let term_product = term
                    .iter()
                    .map(|(value_index, power)| {
                        if *value_index < values.len() {
                            values[*value_index].clone().pow(&[*power as u64])
                        } else {
                            panic!(
                                "value_index {} is out of bounds for values vector of length {}",
                                value_index,
                                values.len()
                            );
                        }
                    })
                    .fold(R::one(), |acc, x| acc * x);
                *coeff * term_product
            })
            .fold(R::zero(), |acc, x| acc + x)
    }

    // The degree of a multivariate polynomaial
    pub fn degree(&self) -> usize {
        self.terms.iter().fold(0, |max_degree, (_, term)| {
            let term_degree = term
                .iter()
                .map(|(_, exp)| exp)
                .sum();
            max_degree.max(term_degree)
        })
    }

    // The degree of a multivariate polynomial
    pub fn multi_degree(&self) -> Vec<usize> {
        // Find the maximum index of the variables used in the polynomial
        let num_variables =
            self.terms
                .iter()
                .flat_map(|(_, term)| term.iter().map(|(var, _)| *var))
                .max()
                .unwrap_or(0) + 1;

        // Initialize a vector to hold the maximum degree for each variable
        let mut degrees = vec![0; num_variables];

        // Iterate over each term in the polynomial
        for (_, term) in &self.terms {
            // Update the degree for each variable in the term
            for (var, exp) in term {
                if *exp > degrees[*var] {
                    degrees[*var] = *exp;
                }
            }
        }

        degrees
    }

    // Simplify the polynomial by summing terms with the same monomials
    pub fn simplify(&self) -> MultiPoly<R> {
        let mut simplified_terms = Vec::new();

        for &(coeff, ref term) in &self.terms {
            // Search for a term with the same monomials in the simplified_terms vector
            if
                let Some((existing_coeff, _)) = simplified_terms
                    .iter_mut()
                    .find(|(_, t)| *t == *term)
            {
                *existing_coeff += coeff;
            } else {
                // Otherwise add the term to the polynomial
                simplified_terms.push((coeff, term.clone()));
            }
        }

        MultiPoly { terms: simplified_terms }
    }

    // Calculates the number of variables in the polynomial
    pub fn num_vars(&self) -> usize {
        self.terms
            .iter()
            .flat_map(|(_, term)| term.iter().map(|(idx, _)| *idx))
            .max()
            .map_or(0, |max_idx| max_idx + 1)
    }

    // Computes partial evaluation/summation of a multivariate polynomial. Each
    // variable either remains unevaluated or is summed over a fixed set of values.
    // Using a single value for the summation set is equivalent to evaluating the
    // corresponding variable at this value.
    pub fn partial_summation(&self, vals: &Vec<Option<Vec<R>>>) -> MultiPoly<R> {
        let monom_terms: Vec<Option<Vec<R>>> = vals
            .iter()
            .map(|val| {
                match val {
                    Some(summands) => Some(power_sums(self.degree(), summands)),
                    None => None,
                }
            })
            .collect();

        let mut dense_monomial_exponents = vec![0usize; self.num_vars()];

        let new_terms = self.terms
            .iter()
            .map(|(coeff, term)| {
                // Load exponents of all variables (including degree 0) into Vec
                dense_monomial_exponents.fill(0usize);
                for (idx, exp) in term.iter() {
                    dense_monomial_exponents[*idx] = *exp;
                }

                // Compute the new coefficient based on the monomial terms
                let new_coeff = dense_monomial_exponents
                    .iter()
                    .enumerate()
                    .map(|(idx, exp)| {
                        match &monom_terms[idx] {
                            Some(power_sums) => power_sums[*exp],
                            None => R::one(),
                        }
                    })
                    .fold(*coeff, |acc, val| acc * val);

                // Filter out the terms with None values and create a new SparseTerm
                let new_term = term
                    .iter()
                    .filter(|(idx, _)| monom_terms[*idx].is_none())
                    .map(|&(idx, exp)| (idx, exp))
                    .collect();

                (new_coeff, new_term)
            })
            .collect();

        MultiPoly { terms: new_terms }
    }
    /// Compute the partial evaluation of `poly` by the specified field element and variable.
    pub fn partial_eval(&self, value: R, variable: usize) -> MultiPoly<R> {
        let vals = (0..self.num_vars())
            .map(|n| if n == variable { Some(vec![value]) } else { None })
            .collect();
        self.partial_summation(&vals)
    }
    /// Restrict a multivariate polynomial to the univariate polynomial obtained by
    /// evaluating at zero for all variables except the one specified.
    pub fn to_univariate(&self, variable: usize) -> UnivPoly<R> {
        let mut coeffs: Vec<R> = vec![R::zero(); self.degree() + 1];
        self.terms
            .iter()
            .filter(|(_, monom)| { !monom.iter().any(|(idx, exp)| *idx != variable && *exp > 0) })
            .fold(&mut coeffs, |acc, (coeff, monom)| {
                let exp = if monom.is_empty() { 0usize } else { monom[0].1 };
                (*acc)[exp] += coeff;
                acc
            });
        UnivPoly { coeffs }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use lattirust_arithmetic::ring::Z2_64;

    #[test]
    fn test_univ_poly_eval() {
        // Define a polynomial 3 + 2x + x^2 in Z2_64
        let poly = UnivPoly {
            coeffs: vec![Z2_64::from(3 as u8), Z2_64::from(2 as u8), Z2_64::from(1 as u8)],
        };

        // Evaluate the polynomial at x = 1
        let x = Z2_64::from(1 as u8);
        assert_eq!(poly.eval(&x), Z2_64::from(6 as u8));

        // Evaluate the polynomial at x = 2
        let x = Z2_64::from(2 as u8);
        assert_eq!(poly.eval(&x), Z2_64::from(11 as u8));
    }
    #[test]
    fn test_multi_poly_eval() {
        // Define a polynomial 3 + 2x1^2 + x1*x2 + 5x2 in Z2_64
        let poly = MultiPoly {
            terms: vec![
                (Z2_64::from(3 as u8), vec![]),
                (Z2_64::from(2 as u8), vec![(0, 2)]),
                (Z2_64::from(1 as u8), vec![(0, 1), (1, 1)]),
                (Z2_64::from(5 as u8), vec![(1, 1)])
            ],
        };

        // Evaluate the polynomial at x1 = 1, x2 = 1
        let values = vec![Z2_64::from(1 as u8), Z2_64::from(1 as u8)];

        assert_eq!(poly.eval_poly(&values), Z2_64::from(11 as u8));

        // Evaluate the polynomial at x1 = 2, x2 = 3
        let values = vec![Z2_64::from(2 as u8), Z2_64::from(3 as u8)];
        assert_eq!(poly.eval_poly(&values), Z2_64::from(32 as u8));
    }
    #[test]
    fn test_partial_summation() {
        // Define a polynomial 3 + 2x1^2 + x1*x2 + 5x2 in Z2_64
        let poly = MultiPoly {
            terms: vec![
                (Z2_64::from(3 as u8), vec![]),
                (Z2_64::from(2 as u8), vec![(0, 2)]),
                (Z2_64::from(1 as u8), vec![(0, 1), (1, 1)]),
                (Z2_64::from(5 as u8), vec![(1, 1)])
            ],
        };

        // Partial summation with x1 = {1, 2} and x2 = 3
        let vals: Vec<Option<Vec<Z2_64>>> = vec![
            Some(vec![Z2_64::from(1 as u8), Z2_64::from(2 as u8)]),
            Some(vec![Z2_64::from(3 as u8)])
        ];

        let result = poly.partial_summation(&vals).simplify();

        // Expected polynomial: 3 + 2*(1^2 + 2^2) + (1*3 + 2*3) + 5*3
        // = 3 + 2*(1 + 4) + (3 + 6) + 15
        // = 3 + 10 + 9 + 15 = 37
        let expected_result = MultiPoly {
            terms: vec![(Z2_64::from(55 as u8), vec![])],
        };

        assert_eq!(result, expected_result);
    }

    // Test cases for UnivPoly
    #[test]
    fn test_univ_poly_degree() {
        // Create a univariate polynomial 3 + 2x + 5x^2
        let poly = UnivPoly {
            coeffs: vec![Z2_64::from(3 as u8), Z2_64::from(2 as u8), Z2_64::from(5 as u8)],
        };
        // The degree of this polynomial is 2
        assert_eq!(poly.degree(), 2);
    }

    // Test cases for MultiPoly
    #[test]
    fn test_multi_poly_degree() {
        // Create a multivariate polynomial
        let poly = MultiPoly {
            terms: vec![
                (Z2_64::from(3 as u8), vec![(0, 2)]), // 3x1^2
                (Z2_64::from(2 as u8), vec![(1, 1)]), // 2x2
                (Z2_64::from(1 as u8), vec![(0, 1), (1, 2)]) // x1x2^2
            ],
        };
        // The highest total degree of this polynomial is 3 (x1x2^2)
        assert_eq!(poly.degree(), 3);
    }

    #[test]
    fn test_num_vars() {
        // Create a multivariate polynomial with 3 variables
        let poly = MultiPoly {
            terms: vec![
                (Z2_64::from(3 as u8), vec![(0, 2), (1, 1)]), // 3x1^2 * x2
                (Z2_64::from(2 as u8), vec![(1, 1), (2, 2)]), // 2x2 * x3^2
                (Z2_64::from(1 as u8), vec![(0, 1), (2, 1)]) // x1 * x3
            ],
        };

        // The number of variables in this polynomial is 3
        assert_eq!(poly.num_vars(), 3);
    }

    #[test]
    fn test_to_univariate() {
        // Define a multivariate polynomial: 2 + 3x1 + 4x1^2 + 5x2 + 6x1x2
        let poly = MultiPoly {
            terms: vec![
                (Z2_64::from(2 as u8), vec![]),
                (Z2_64::from(3 as u8), vec![(0, 1)]),
                (Z2_64::from(4 as u8), vec![(0, 2)]),
                (Z2_64::from(5 as u8), vec![(1, 1)]),
                (Z2_64::from(6 as u8), vec![(0, 1), (1, 1)])
            ],
        };

        // Convert the multivariate polynomial to a univariate polynomial
        // by evaluating all variables except for x1 at zero
        let univariate_poly = poly.to_univariate(0);

        // Expected univariate polynomial: 2 + 4x + 6x^2
        let expected_univariate_poly = UnivPoly {
            coeffs: vec![Z2_64::from(2 as u8), Z2_64::from(3 as u8), Z2_64::from(4 as u8)],
        };

        assert_eq!(univariate_poly, expected_univariate_poly);
    }
}
