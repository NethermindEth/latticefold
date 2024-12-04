//! Prover

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{cfg_into_iter, cfg_iter_mut, vec::Vec};
use lattirust_poly::{mle::MultilinearExtension, polynomials::DenseMultilinearExtension};
use lattirust_ring::{OverField, Ring};

use super::{dense_polynomial::DensePolynomial, verifier::VerifierMsg, IPForMLSumcheck};

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
    /// Stores a list of multilinear extensions in which `self.list_of_products` points to
    pub flattened_ml_extensions: Vec<DenseMultilinearExtension<R>>,
    /// Number of variables
    pub num_vars: usize,
    /// Max number of multiplicands in a product
    pub max_multiplicands: usize,
    /// The current round number
    pub round: usize,
}

//impl<R: OverField> ProverState<R> {
//    pub fn combine_product(&self, vals: &[R]) -> R {
//        let mut sum = R::zero();
//        for (coefficient, products) in &self.list_of_products {
//            let mut prod = *coefficient;
//            for j in products {
//                prod *= vals[*j];
//            }
//            sum += prod;
//        }
//        sum
//    }
//}

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
    pub fn prover_init(polynomial: &DensePolynomial<R>) -> ProverState<R> {
        if polynomial.aux_info.num_variables == 0 {
            panic!("Attempt to prove a constant.")
        }

        // create a deep copy of all unique MLExtensions
        let flattened_ml_extensions = ark_std::cfg_iter!(polynomial.flattened_ml_extensions)
            .map(|x| x.as_ref().clone())
            .collect();

        ProverState {
            randomness: Vec::with_capacity(polynomial.aux_info.num_variables),
            flattened_ml_extensions,
            num_vars: polynomial.aux_info.num_variables,
            max_multiplicands: polynomial.aux_info.max_degree,
            round: 0,
        }
    }

    /// receive message from verifier, generate prover message, and proceed to next round
    ///
    /// Adapted Jolt's sumcheck implementation
    pub fn prove_round(
        prover_state: &mut ProverState<R>,
        v_msg: &Option<VerifierMsg<R>>,
        comb_fn: impl Fn(&[R]) -> R + Sync + Send,
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
            eval_points[0] += comb_fn(&params_zero);

            let params_one: Vec<R> = polys.iter().map(|poly| poly[index + 1]).collect();
            eval_points[1] += comb_fn(&params_one);

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

                *eval_point += comb_fn(&poly_evals);
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
