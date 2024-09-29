// PGold = 2^64 − 2^32 + 1
use lattirust_ring::cyclotomic_ring::models::stark_prime::{RqNTT, RqPoly};

use super::SuitableRing;

pub type StarkRingNTT = RqNTT;
pub type StarkRingPoly = RqPoly;

impl SuitableRing for StarkRingNTT {
    type CoefficientRepresentation = StarkRingPoly;

    type PoseidonParams = StarkPoseidonConfig;
}

pub struct StarkPoseidonConfig;

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
