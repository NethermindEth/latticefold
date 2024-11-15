//! Prover

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{cfg_into_iter, cfg_iter_mut, vec::Vec};
use lattirust_poly::{
    mle::MultilinearExtension,
    polynomials::{DenseMultilinearExtension, VirtualPolynomial},
};
use lattirust_ring::{OverField, Ring};

use super::{verifier::VerifierMsg, IPForMLSumcheck};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// Prover Message
#[derive(Clone, Debug, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct ProverMsg<R1: Ring> {
    /// evaluations on P(0), P(1), P(2), ...
    pub(crate) evaluations: Vec<R1>,
}

/// Prover State
pub struct ProverState<R: OverField> {
    /// sampled randomness given by the verifier
    pub randomness: Vec<R::BaseRing>,
    /// Stores the list of products that is meant to be added together. Each multiplicand is represented by
    /// the index in flattened_ml_extensions
    pub list_of_products: Vec<(R, Vec<usize>)>,
    /// Stores a list of multilinear extensions in which `self.list_of_products` points to
    pub flattened_ml_extensions: Vec<DenseMultilinearExtension<R>>,
    /// Number of variables
    pub num_vars: usize,
    /// Max number of multiplicands in a product
    pub max_multiplicands: usize,
    /// The current round number
    pub round: usize,
}

impl<R: OverField, T> IPForMLSumcheck<R, T> {
    /// initialize the prover to argue for the sum of polynomial over {0,1}^`num_vars`
    ///
    /// The polynomial is represented by a list of products of polynomials along with its coefficient that is meant to be added together.
    ///
    /// This data structure of the polynomial is a list of list of `(coefficient, DenseMultilinearExtension)`.
    /// * Number of products n = `polynomial.products.len()`,
    /// * Number of multiplicands of ith product m_i = `polynomial.products[i].1.len()`,
    /// * Coefficient of ith product c_i = `polynomial.products[i].0`
    ///
    /// The resulting polynomial is
    ///
    /// $$\sum_{i=0}^{n}C_i\cdot\prod_{j=0}^{m_i}P_{ij}$$
    ///
    pub fn prover_init(polynomial: &VirtualPolynomial<R>) -> ProverState<R> {
        if polynomial.aux_info.num_variables == 0 {
            panic!("Attempt to prove a constant.")
        }

        // create a deep copy of all unique MLExtensions
        let flattened_ml_extensions = polynomial
            .flattened_ml_extensions
            .iter()
            .map(|x| x.as_ref().clone())
            .collect();

        ProverState {
            randomness: Vec::with_capacity(polynomial.aux_info.num_variables),
            list_of_products: polynomial.products.clone(),
            flattened_ml_extensions,
            num_vars: polynomial.aux_info.num_variables,
            max_multiplicands: polynomial.aux_info.max_degree,
            round: 0,
        }
    }

    pub fn prove_round(
        prover_state: &mut ProverState<R>,
        v_msg: &Option<VerifierMsg<R>>,
    ) -> ProverMsg<R> {
        if let Some(msg) = v_msg {
            if prover_state.round == 0 {
                panic!("first round should be prover first.");
            }
            prover_state.randomness.push(msg.randomness);

            // fix argument
            let i = prover_state.round;
            let r = prover_state.randomness[i - 1];
            cfg_iter_mut!(prover_state.flattened_ml_extensions).for_each(|multiplicand| {
                *multiplicand = multiplicand.fix_variables(&[r.into()]);
            });
        } else if prover_state.round > 0 {
            panic!("verifier message is empty");
        }

        prover_state.round += 1;

        if prover_state.round > prover_state.num_vars {
            panic!("Prover is not active");
        }

        let i = prover_state.round;
        let nv = prover_state.num_vars;
        let degree = prover_state.max_multiplicands;

        // Pre-allocate result array since degree is small (4)
        let mut products_sum = vec![R::zero(); degree + 1];

        // Pre-allocate scratch space for intermediate products
        let product_scratch = vec![R::zero(); degree + 1];

        // Calculate number of iterations
        let num_iterations = 1 << (nv - i);

        #[cfg(feature = "parallel")]
        {
            // Process chunks of iterations in parallel
            let chunk_size = (num_iterations + 7) / 8; // Divide work into 8 chunks or less
            products_sum = (0..num_iterations)
                .into_par_iter()
                .chunks(chunk_size)
                .map(|chunk_iter| {
                    let mut chunk_sum = vec![R::zero(); degree + 1];
                    let mut product = vec![R::zero(); degree + 1];
                    for b in chunk_iter {
                        // Process each product term
                        for (coefficient, products) in &prover_state.list_of_products {
                            for p in &mut product {
                                *p = *coefficient;
                            }

                            // Multiply by each factor
                            for &jth_product in products {
                                let table = &prover_state.flattened_ml_extensions[jth_product];
                                let start = table[b << 1];
                                let step = table[(b << 1) + 1] - start;

                                // Unrolled multiplication since degree is 4
                                let mut val = start;
                                for p in product.iter_mut() {
                                    *p *= val;
                                    val += step;
                                }
                            }

                            // Add to chunk sum
                            for (sum, term) in chunk_sum.iter_mut().zip(product.iter()) {
                                *sum += term;
                            }
                        }
                    }
                    chunk_sum
                })
                .reduce(
                    || vec![R::zero(); degree + 1],
                    |mut a, b| {
                        for (x, y) in a.iter_mut().zip(b.iter()) {
                            *x += y;
                        }
                        a
                    },
                );
        }

        #[cfg(not(feature = "parallel"))]
        {
            let mut product = vec![R::zero(); degree + 1];
            for b in 0..num_iterations {
                for (coefficient, products) in &prover_state.list_of_products {
                    // Initialize product array with coefficient
                    product.copy_from_slice(&product_scratch);
                    for p in &mut product {
                        *p = *coefficient;
                    }

                    // Multiply by each factor
                    for &jth_product in products {
                        let table = &prover_state.flattened_ml_extensions[jth_product];
                        let start = table[b << 1];
                        let step = table[(b << 1) + 1] - start;

                        // Unrolled multiplication since degree is 4
                        let mut val = start;
                        for p in product.iter_mut() {
                            *p *= val;
                            val += step;
                        }
                    }

                    // Add to total sum
                    for (sum, term) in products_sum.iter_mut().zip(product.iter()) {
                        *sum += term;
                    }
                }
            }
        }

        ProverMsg {
            evaluations: products_sum,
        }
    }
}
