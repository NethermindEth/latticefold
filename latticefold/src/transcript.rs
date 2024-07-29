use lattirust_arithmetic::challenge_set::latticefold_challenge_set::*;
use std::fmt::Debug;

pub mod poseidon;

pub trait Transcript<R: OverField> {
    type ChallengeSet: LatticefoldChallengeSet<R>;
    type TranscriptConfig: Debug;

    fn new(config: &Self::TranscriptConfig) -> Self;
    fn absorb(&mut self, v: &R::F);
    fn absorb_vec(&mut self, v: &[R::F]);
    fn absorb_ring(&mut self, v: &R);
    fn absorb_ring_vec(&mut self, v: &[R]);
    fn get_big_challenge(&mut self) -> R::BaseRing;
    fn get_small_challenge(&mut self) -> R;

    // /// get_challenge_nbits returns a field element of size nbits
    // fn get_challenge_nbits(&mut self, nbits: usize) -> Vec<bool>;
    // fn get_challenges(&mut self, n: usize) -> Vec<C::ScalarField>;
}
