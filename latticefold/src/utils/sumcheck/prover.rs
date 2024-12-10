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

        struct Scratch<R> {
            evals: Vec<R>,
            steps: Vec<R>,
            vals0: Vec<R>,
            vals1: Vec<R>,
            vals: Vec<R>,
            levals: Vec<R>,
        }
        let scratch = || Scratch {
            evals: vec![R::zero(); degree + 1],
            steps: vec![R::zero(); polys.len()],
            vals0: vec![R::zero(); polys.len()],
            vals1: vec![R::zero(); polys.len()],
            vals: vec![R::zero(); polys.len()],
            levals: vec![R::zero(); degree + 1],
        };

        #[cfg(not(feature = "parallel"))]
        let zeros = scratch();
        #[cfg(feature = "parallel")]
        let zeros = scratch;

        let summer = cfg_into_iter!(0..1 << (nv - i)).fold(zeros, |mut s, b| {
            let index = b << 1;

            s.vals0
                .iter_mut()
                .zip(polys.iter())
                .for_each(|(v0, poly)| *v0 = poly[index]);
            s.levals[0] = comb_fn(&s.vals0);

            s.vals1
                .iter_mut()
                .zip(polys.iter())
                .for_each(|(v1, poly)| *v1 = poly[index + 1]);
            s.levals[1] = comb_fn(&s.vals1);

            for (i, (v1, v0)) in s.vals1.iter().zip(s.vals0.iter()).enumerate() {
                s.steps[i] = *v1 - v0;
                s.vals[i] = *v1;
            }

            for eval_point in s.levals.iter_mut().take(degree + 1).skip(2) {
                for poly_i in 0..polys.len() {
                    s.vals[poly_i] += s.steps[poly_i];
                }
                *eval_point = comb_fn(&s.vals);
            }

            s.evals
                .iter_mut()
                .zip(s.levals.iter())
                .for_each(|(e, l)| *e += l);
            s
        });

        // Rayon's fold outputs an iter which still needs to be summed over
        #[cfg(feature = "parallel")]
        let evaluations = summer.map(|s| s.evals).reduce(
            || vec![R::zero(); degree + 1],
            |mut evaluations, levals| {
                evaluations
                    .iter_mut()
                    .zip(levals)
                    .for_each(|(e, l)| *e += l);
                evaluations
            },
        );

        #[cfg(not(feature = "parallel"))]
        let evaluations = summer.evals;

        ProverMsg { evaluations }
    }
}
