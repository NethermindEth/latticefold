//! Prover
use std::sync::Arc;

use ark_std::{cfg_iter_mut, vec::Vec};
use lattirust_poly::{
    mle::MultilinearExtension,
    polynomials::{DenseMultilinearExtension, VirtualPolynomial},
};
use lattirust_ring::{OverField, Ring};

use super::{verifier::VerifierMsg, IPForMLSumcheck};

/// Prover Message
#[derive(Clone)]
pub struct ProverMsg<R: Ring> {
    /// evaluations on P(0), P(1), P(2), ...
    pub(crate) evaluations: Vec<R>,
}

/// Prover State
pub struct ProverState<R: OverField> {
    /// sampled randomness given by the verifier
    pub randomness: Vec<R::BaseRing>,
    /// flattened multilinear extensions
    pub products: Vec<(R, DenseMultilinearExtension<R>)>,
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

        ProverState {
            randomness: Vec::with_capacity(polynomial.aux_info.num_variables),
            products: polynomial.products.clone(),
            num_vars: polynomial.aux_info.num_variables,
            max_multiplicands: polynomial.aux_info.max_degree,
            round: 0,
        }
    }

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
            cfg_iter_mut!(prover_state.products).for_each(|(_, mle)| {
                *mle = mle.fix_variables(&[r.into()]);
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

        // #[cfg(not(feature = "parallel"))]
        let zeros= vec![R::zero(); degree + 1];
        // #[cfg(feature = "parallel")]
        // let zeros = || (vec![F::zero(); degree + 1], vec![F::zero(); degree + 1]);

        // generate sum
        let fold_result = ark_std::cfg_into_iter!(0..1 << (nv - i), 1 << 10).fold(
            zeros,
            |mut products_sum, b| {
                // In effect, this fold is essentially doing simply:
                // for b in 0..1 << (nv - i) {
                for (coefficient, mle) in &prover_state.products {
                    let table = &mle.evaluations;
                    let start = table[b << 1];
                    let step = table[(b << 1) + 1] - start;
                    for t in 0..degree + 1 {
                        products_sum[t] += *coefficient * (start + step*R::from(t as u64));
                    }
                }
                products_sum
            },
        );

        // #[cfg(not(feature = "parallel"))]
        let products_sum = fold_result;

        // When rayon is used, the `fold` operation results in a iterator of `Vec<F>` rather than a single `Vec<F>`. In this case, we simply need to sum them.
        // #[cfg(feature = "parallel")]
        // let products_sum = fold_result.map(|scratch| scratch.0).reduce(
        //     || vec![F::zero(); degree + 1],
        //     |mut overall_products_sum, sublist_sum| {
        //         overall_products_sum
        //             .iter_mut()
        //             .zip(sublist_sum.iter())
        //             .for_each(|(f, s)| *f += s);
        //         overall_products_sum
        //     },
        // );

        ProverMsg {
            evaluations: products_sum,
        }
    }
}
