// PGold = 2^64 − 2^32 + 1
use crate::challenge_set::LatticefoldChallengeSet;
use ark_ff::Field;
use lattirust_ring::{
    cyclotomic_ring::models::goldilocks::{Fq3, RqNTT},
    PolyRing,
};

pub type PGoldCyclotomicRing<const N: usize> = RqNTT;

#[allow(dead_code)]
pub struct PGoldChallengeSet<const N: usize>;

impl<const N: usize> LatticefoldChallengeSet<PGoldCyclotomicRing<N>> for PGoldChallengeSet<N> {
    fn small_challenge_coefficient_from_random_bytes(
        _i: usize,
        bs: &[u8],
    ) -> <PGoldCyclotomicRing<N> as PolyRing>::BaseRing {
        if bs[0] == 0 {
            Fq3::ZERO
        } else {
            Fq3::ONE
        }
    }
}
