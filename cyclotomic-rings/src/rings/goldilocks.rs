// PGold = 2^64 − 2^32 + 1

use lattirust_ring::cyclotomic_ring::models::goldilocks::{RqNTT, RqPoly};

use super::SuitableRing;

pub type GoldilocksRingNTT = RqNTT;
pub type GoldilocksRingPoly = RqPoly;

impl SuitableRing for GoldilocksRingNTT {
    type CoefficientRepresentation = RqPoly;
    type PoseidonParams = GoldilocksPoseidonConfig;
}

pub struct GoldilocksPoseidonConfig;

#[allow(dead_code)]
pub struct GoldilocksChallengeSet;

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
