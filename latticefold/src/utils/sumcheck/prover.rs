//! Prover
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{cfg_into_iter, cfg_iter_mut, vec::Vec};
use lattirust_poly::{mle::MultilinearExtension, polynomials::DenseMultilinearExtension};
use lattirust_ring::{OverField, Ring};

use super::{verifier::VerifierMsg, virtual_polynomial::VirtualPolynomial, IPForMLSumcheck};

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

#[cfg(feature = "jolt-sumcheck")]
impl<R: OverField> ProverState<R> {
    pub fn combine_product(state: &ProverState<R>, vals: &[R]) -> R {
        let mut sum = R::zero();
        for (coefficient, products) in &state.list_of_products {
            let mut prod = *coefficient;
            for j in products {
                prod *= vals[*j];
            }
            sum += prod;
        }
        sum
    }
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
        let flattened_ml_extensions = ark_std::cfg_iter!(polynomial.flattened_ml_extensions)
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

    #[cfg(not(feature = "jolt-sumcheck"))]
    /// receive message from verifier, generate prover message, and proceed to next round
    ///
    /// Main algorithm used is from section 3.2 of [XZZPS19](https://eprint.iacr.org/2019/317.pdf#subsection.3.2).
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
        let degree = prover_state.max_multiplicands; // the degree of univariate polynomial sent by prover at this round

        #[cfg(not(feature = "parallel"))]
        let zeros = (vec![R::zero(); degree + 1], vec![R::zero(); degree + 1]);
        #[cfg(feature = "parallel")]
        let zeros = || (vec![R::zero(); degree + 1], vec![R::zero(); degree + 1]);

        // generate sum
        let (products_with_zero, products_without_zero): (Vec<_>, Vec<_>) = prover_state
            .list_of_products
            .iter()
            .partition(|(_, products)| products.contains(&0));
        let fold_result =
            cfg_into_iter!(0..1 << (nv - i)).fold(zeros, |(mut products_sum, mut product), b| {
                // Instead of creating a new vector, fill the existing one with zeros
                let mut products_inner_sum = vec![R::zero(); degree + 1];

                // In effect, this fold is essentially doing simply:
                // for b in 0..1 << (nv - i) {
                let index = b << 1;

                // Process products with zero
                for (coefficient, products) in &products_with_zero {
                    product.fill(*coefficient);

                    for &jth_product in products {
                        if jth_product == 0 {
                            continue;
                        }
                        let table = &prover_state.flattened_ml_extensions[jth_product];
                        let mut start = table[index];
                        let step = table[index + 1] - start;
                        for p in product.iter_mut() {
                            *p *= start;
                            start += step;
                        }
                    }
                    for (sum, prod) in products_inner_sum.iter_mut().zip(product.iter()) {
                        *sum += prod;
                    }
                }

                // Factor out eq_beta_x multiplication
                if !products_with_zero.is_empty() {
                    let eq_beta_x = &prover_state.flattened_ml_extensions[0];
                    let mut start = eq_beta_x[index];
                    let step = eq_beta_x[index + 1] - start;
                    for (p, sum) in products_inner_sum.iter_mut().zip(products_sum.iter_mut()) {
                        *p *= start;
                        *sum += *p;
                        start += step;
                    }
                }

                // Process products without zero
                for (coefficient, products) in &products_without_zero {
                    product.fill(*coefficient);
                    for &jth_product in products {
                        let table = &prover_state.flattened_ml_extensions[jth_product];
                        let mut start = table[index];
                        let step = table[index + 1] - start;
                        for p in product.iter_mut() {
                            *p *= start;
                            start += step;
                        }
                    }
                    for (sum, prod) in products_sum.iter_mut().zip(product.iter()) {
                        *sum += prod;
                    }
                }
                (products_sum, product)
            });

        #[cfg(not(feature = "parallel"))]
        let products_sum = fold_result.0;

        // When rayon is used, the `fold` operation results in a iterator of `Vec<F>` rather than a single `Vec<F>`. In this case, we simply need to sum them.
        #[cfg(feature = "parallel")]
        let products_sum = fold_result.map(|scratch| scratch.0).reduce(
            || vec![R::zero(); degree + 1],
            |mut overall_products_sum, sublist_sum| {
                overall_products_sum
                    .iter_mut()
                    .zip(sublist_sum.iter())
                    .for_each(|(f, s)| *f += s);
                overall_products_sum
            },
        );

        ProverMsg {
            evaluations: products_sum,
        }
    }

    #[cfg(feature = "jolt-sumcheck")]
    pub fn prove_round(
        prover_state: &mut ProverState<R>,
        v_msg: &Option<VerifierMsg<R>>,
        comb_fn: impl Fn(&ProverState<R>, &[R]) -> R + Sync + Send,
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

        let polys = &prover_state.flattened_ml_extensions;

        let iter = cfg_into_iter!(0..1 << (nv - i)).map(|b| {
            let index = b << 1;
            let mut eval_points = vec![R::zero(); degree + 1];

            let params_zero: Vec<R> = polys.iter().map(|poly| poly[index]).collect();
            eval_points[0] += comb_fn(prover_state, &params_zero);

            let params_one: Vec<R> = polys.iter().map(|poly| poly[index + 1]).collect();
            eval_points[1] += comb_fn(prover_state, &params_one);

            let steps: Vec<R> = params_one
                .iter()
                .zip(params_zero)
                .map(|(p1, p0)| *p1 - p0)
                .collect();

            let mut poly_evals = vec![R::zero(); polys.len()];
            let mut current = params_one;
            for eval_point in eval_points.iter_mut().take(degree + 1).skip(2) {
                for poly_i in 0..polys.len() {
                    poly_evals[poly_i] = current[poly_i] + steps[poly_i];
                }

                *eval_point += comb_fn(prover_state, &poly_evals);
                ark_std::mem::swap(&mut current, &mut poly_evals);
            }

            eval_points
        });

        // Rayon's reduce interface is different from standard's
        #[cfg(feature = "parallel")]
        let products_sum = iter.reduce(
            || vec![R::zero(); degree + 1],
            |mut products_sum, eval_points| {
                products_sum
                    .iter_mut()
                    .zip(eval_points)
                    .for_each(|(s, e)| *s += e);
                products_sum
            },
        );

        #[cfg(not(feature = "parallel"))]
        let products_sum = {
            let mut products_sum = vec![R::zero(); degree + 1];
            iter.for_each(|eval_points| {
                products_sum
                    .iter_mut()
                    .zip(eval_points)
                    .for_each(|(s, e)| *s += e);
            });
            products_sum
        };

        ProverMsg {
            evaluations: products_sum,
        }
    }
}
