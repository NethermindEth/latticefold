use crate::poly_utils::{ MultiPoly, UnivPoly };
use lattirust_arithmetic::ring::Ring;
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct SumCheckTranscript<R: Ring> {
    claimed_sum: R,
    polynomial: MultiPoly<R>,
    rounds: Vec<Challenge<R>>,
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct Challenge<R: Ring> {
    challenge: R,
    response: UnivPoly<R>,
}

impl<R: Ring> SumCheckTranscript<R> {
    pub fn new(
        claimed_sum: R,
        polynomial: MultiPoly<R>,
        num_rounds: usize
    ) -> SumCheckTranscript<R> {
        SumCheckTranscript { claimed_sum, polynomial, rounds: Vec::with_capacity(num_rounds) }
    }
}
