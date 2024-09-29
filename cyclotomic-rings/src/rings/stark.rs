// PGold = 2^64 − 2^32 + 1
use crate::challenge_set::LatticefoldChallengeSet;
use ark_ff::Field;
use lattirust_ring::{
    cyclotomic_ring::models::stark_prime::{RqNTT, RqPoly},
    PolyRing,
};

use super::SuitableRing;

pub type StarkRingNTT = RqNTT;
pub type StarkRingPoly = RqPoly;

impl SuitableRing for StarkRingNTT {
    type CoefficientRepresentation = StarkRingPoly;
}

#[allow(dead_code)]
pub struct StarkChallengeSet;

// impl LatticefoldChallengeSet<GoldilocksRingNTT> for PGoldChallengeSet {
//     fn small_challenge_coefficient_from_random_bytes(
//         _i: usize,
//         bs: &[u8],
//     ) -> <GoldilocksRingPoly as PolyRing>::BaseRing {
//         if bs[0] == 0 {
//             Fq3::ZERO
//         } else {
//             Fq3::ONE
//         }
//     }
// }
