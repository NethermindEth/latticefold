pub mod prover;
pub mod univ_poly;
pub mod verifier;

use std::fmt::Display;
use std::marker::PhantomData;

use crate::transcript::Transcript;
use ark_crypto_primitives::sponge::Absorb;
use ark_ff::PrimeField;
use lattirust_arithmetic::challenge_set::latticefold_challenge_set::OverField;
use lattirust_arithmetic::polynomials::ArithErrors;
use lattirust_arithmetic::polynomials::VPAuxInfo;
use lattirust_arithmetic::ring::Ring;
use thiserror_no_std::Error;
use univ_poly::UnivPoly;

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
    pub unipoly: UnivPoly<R>,
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

    pub fn add_round(&mut self, transcript: &mut impl Transcript<F, R>, unipoly: UnivPoly<R>) {
        transcript.absorb_ring_vec(&unipoly.coeffs);
        let round = SumCheckRound {
            unipoly,
            _marker: std::marker::PhantomData,
        };

        self.rounds.push(round);
    }
}

impl<F: PrimeField, R: OverField<F>> SumCheckRound<F, R> {
    pub fn new(unipoly: UnivPoly<R>) -> SumCheckRound<F, R> {
        SumCheckRound {
            unipoly,
            _marker: std::marker::PhantomData,
        }
    }
}
