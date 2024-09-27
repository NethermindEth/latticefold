use ark_std::fmt::Debug;
use ark_std::One;
use lattirust_ring::OverField;

use cyclotomic_rings::challenge_set::LatticefoldChallengeSet;

pub mod poseidon;

pub trait Transcript<R: OverField> {
    type ChallengeSet: LatticefoldChallengeSet<R>;
    type TranscriptConfig: Debug;

    fn new(config: &Self::TranscriptConfig) -> Self;
    fn absorb(&mut self, v: &R);
    fn absorb_field_element(&mut self, e: &R) {
        self.absorb(&From::from(e.from_field_element(*e)))
    }
    fn absorb_slice(&mut self, v: &[R]);
    fn get_big_challenge(&mut self) -> R::BaseRing;
    fn get_small_challenge(&mut self) -> R;
    fn get_small_challenges(&mut self, n: usize) -> Vec<R> {
        let mut challenges = Vec::with_capacity(n);
        challenges.extend((0..n).map(|_| self.get_small_challenge()));
        challenges
    }

    fn get_big_challenges(&mut self, n: usize) -> Vec<R::BaseRing> {
        let mut challenges = Vec::with_capacity(n);
        challenges.extend((0..n).map(|_| self.get_big_challenge()));
        challenges
    }
    // /// get_challenge_nbits returns a field element of size nbits
    // fn get_challenge_nbits(&mut self, nbits: usize) -> Vec<bool>;
}
