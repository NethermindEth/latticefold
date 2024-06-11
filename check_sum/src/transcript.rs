use std::iter::Sum;

use crate::poly_utils::{ MultiPoly, UnivPoly };
use lattirust_arithmetic::ring::Ring;
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct SumCheckTranscript<R: Ring> {
    pub claimed_sum: R,
    pub polynomial: MultiPoly<R>,
    pub rounds: Vec<SumCheckRound<R>>,
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct SumCheckRound<R: Ring> {
    pub challenge: R,
    var_index: usize,
    pub unipoly: UnivPoly<R>,
    multipoly: MultiPoly<R>,
}

impl<R: Ring> SumCheckTranscript<R> {
    pub fn new(
        claimed_sum: R,
        polynomial: MultiPoly<R>,
        num_rounds: usize
    ) -> SumCheckTranscript<R> {
        SumCheckTranscript { claimed_sum, polynomial, rounds: Vec::with_capacity(num_rounds) }
    }

    pub fn add_round(&mut self, round: SumCheckRound<R>) {
        self.rounds.push(round);
    }

    // pub fn validate(&self) -> bool;
}

impl<R: Ring> SumCheckRound<R> {
    pub fn new(
        multipoly: MultiPoly<R>,
        unipoly: UnivPoly<R>,
        var_index: usize,
        challenge: R
    ) -> SumCheckRound<R> {
        SumCheckRound {
            challenge,
            var_index,
            unipoly,
            multipoly,
        }
    }
}
