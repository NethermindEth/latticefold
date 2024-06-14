use crate::{
    poly_utils::{ MultiPoly, UnivPoly },
    transcript::{ poseidon::PoseidonTranscript, Transcript },
};

use ark_crypto_primitives::sponge::{ poseidon::PoseidonConfig, Absorb };
use ark_ff::PrimeField;
use lattirust_arithmetic::challenge_set::latticefold_challenge_set::{
    LatticefoldChallengeSet,
    OverField,
};

pub struct SumCheckTranscript<F: PrimeField, R: OverField<F>, CS: LatticefoldChallengeSet<F, R>>
    where F: Absorb {
    pub claimed_sum: R,
    pub polynomial: MultiPoly<R>,
    pub rounds: Vec<SumCheckRound<F, R>>,
    pub hasher: PoseidonTranscript<F, R, CS>,
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct SumCheckRound<F: PrimeField, R: OverField<F>> {
    _marker: std::marker::PhantomData<F>,
    pub challenge: R,
    var_index: usize,
    pub unipoly: UnivPoly<R>,
}

impl<F: PrimeField, R: OverField<F>, CS: LatticefoldChallengeSet<F, R>> SumCheckTranscript<F, R, CS>
    where F: Absorb
{
    pub fn new(
        claimed_sum: R,
        polynomial: MultiPoly<R>,
        num_rounds: usize
    ) -> SumCheckTranscript<F, R, CS> {
        let config = PoseidonConfig {
            full_rounds: 8, // Example values, adjust according to your needs
            partial_rounds: 57,
            alpha: 5,
            ark: vec![vec![F::zero(); 3]; 8 + 57], // Adjust to actual ark parameters
            mds: vec![vec![F::zero(); 3]; 3], // Adjust to actual MDS matrix parameters
            rate: 2,
            capacity: 1,
        };

        SumCheckTranscript {
            claimed_sum,
            polynomial,
            rounds: Vec::with_capacity(num_rounds),
            hasher: PoseidonTranscript::new(&config),
        }
    }

    pub fn add_round(&mut self, challenge: R, var_index: usize, unipoly: UnivPoly<R>) {
        let round = SumCheckRound {
            challenge,
            var_index,
            unipoly,
            _marker: std::marker::PhantomData,
        };
        self.rounds.push(round);
    }
}

impl<F: PrimeField, R: OverField<F>> SumCheckRound<F, R> {
    pub fn new(unipoly: UnivPoly<R>, var_index: usize, challenge: R) -> SumCheckRound<F, R> {
        SumCheckRound {
            challenge,
            var_index,
            unipoly,
            _marker: std::marker::PhantomData,
        }
    }
}
