pub mod prover;
pub mod verifier;
pub mod univ_poly;

use std::fmt::Display;
use std::marker::PhantomData;

use crate::transcript::Transcript;
use ark_crypto_primitives::sponge::Absorb;
use ark_ff::PrimeField;
use lattirust_arithmetic::challenge_set::latticefold_challenge_set::OverField;
use lattirust_arithmetic::polynomials::VPAuxInfo;
use lattirust_arithmetic::polynomials::{ ArithErrors, VirtualPolynomial };
use lattirust_arithmetic::ring::Ring;
use thiserror_no_std::Error;

pub struct SumCheckIP<F: PrimeField, R: OverField<F>> where F: Absorb {
    pub _f: PhantomData<F>,
    pub claimed_sum: R,
    pub poly_info: VPAuxInfo<R>,
}

impl<F: PrimeField, R: OverField<F>> SumCheckIP<F, R> where F: Absorb {
    pub fn new(claimed_sum: R, poly_info: VPAuxInfo<R>) -> Self {
        SumCheckIP {
            _f: Default::default(),
            claimed_sum,
            poly_info,
        }
    }
}

pub struct SumCheckProof<F: PrimeField, R: OverField<F>> where F: Absorb {
    pub rounds: Vec<SumCheckRound<F, R>>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SumCheckRound<F: PrimeField, R: OverField<F>> {
    _marker: std::marker::PhantomData<F>,
    pub unipoly: VirtualPolynomial<R>,
}

#[derive(Error, Debug)]
pub enum SumCheckError<R: Ring + Display> {
    #[error("univariate polynomial evaluation error")] EvaluationError(#[from] ArithErrors),
    #[error("incorrect sumcheck sum. Expected `{0}`. Received `{1}`")] SumCheckFailed(R, R),
    #[error("max degree exceeded")]
    MaxDegreeExceeded,
}

impl<F: PrimeField, R: OverField<F>> SumCheckProof<F, R> where F: Absorb {
    pub fn new(num_rounds: usize) -> SumCheckProof<F, R> {
        SumCheckProof {
            rounds: Vec::with_capacity(num_rounds),
        }
    }

    pub fn add_round(
        &mut self,
        transcript: &mut impl Transcript<F, R>,
        unipoly: VirtualPolynomial<R>
    ) {
        let coeffs: Vec<R> = unipoly.products
            .clone()
            .into_iter()
            .map(|(x, _)| x)
            .collect();
        transcript.absorb_ring_vec(&coeffs);
        let round = SumCheckRound {
            unipoly,
            _marker: std::marker::PhantomData,
        };

        self.rounds.push(round);
    }
}

impl<F: PrimeField, R: OverField<F>> SumCheckRound<F, R> {
    pub fn new(unipoly: VirtualPolynomial<R>) -> SumCheckRound<F, R> {
        SumCheckRound {
            unipoly,
            _marker: std::marker::PhantomData,
        }
    }
}
