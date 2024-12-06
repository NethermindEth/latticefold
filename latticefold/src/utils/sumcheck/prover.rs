//! Prover

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{cfg_into_iter, cfg_iter_mut, vec::Vec};
use lattirust_poly::{
    mle::MultilinearExtension,
    polynomials::{DenseMultilinearExtension, RefCounter},
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
    /// Stores a list of multilinear extensions
    pub mles: Vec<DenseMultilinearExtension<R>>,
    /// Number of variables
    pub num_vars: usize,
    /// Max degree
    pub max_degree: usize,
    /// The current round number
    pub round: usize,
}

impl<R: OverField, T> IPForMLSumcheck<R, T> {
    /// initialize the prover to argue for the sum of polynomial over {0,1}^`num_vars`
    pub fn prover_init(
        mles: &[RefCounter<DenseMultilinearExtension<R>>],
        nvars: usize,
        degree: usize,
    ) -> ProverState<R> {
        if nvars == 0 {
            panic!("Attempt to prove a constant.")
        }

        // create a deep copy of all unique MLExtensions
        let mles = ark_std::cfg_iter!(mles)
            .map(|x| x.as_ref().clone())
            .collect();

        ProverState {
            randomness: Vec::with_capacity(nvars),
            mles,
            num_vars: nvars,
            max_degree: degree,
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
            cfg_iter_mut!(prover_state.mles).for_each(|multiplicand| {
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
        let degree = prover_state.max_degree;

        let polys = &prover_state.mles;

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
