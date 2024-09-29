// PGold = 2^64 − 2^32 + 1

use lattirust_ring::{cyclotomic_ring::models::goldilocks::{Fq, RqNTT, RqPoly}, Ring};

use super::SuitableRing;
use crate::challenge_set::LatticefoldChallengeSet;

pub type GoldilocksRingNTT = RqNTT;
pub type GoldilocksRingPoly = RqPoly;

impl SuitableRing for GoldilocksRingNTT {
    type CoefficientRepresentation = RqPoly;
    type PoseidonParams = GoldilocksPoseidonConfig;
}

pub struct GoldilocksPoseidonConfig;

#[allow(dead_code)]
pub struct GoldilocksChallengeSet;

const MAX_COEFF: i16 = 32;

impl LatticefoldChallengeSet<GoldilocksRingNTT> for GoldilocksChallengeSet {
    const BYTES_NEEDED: usize = 18;

    fn small_challenge_from_random_bytes(
        bs: &[u8],
    ) -> Result<GoldilocksRingPoly, crate::challenge_set::error::ChallengeSetError> {
        if bs.len() != Self::BYTES_NEEDED {
            return Err(crate::challenge_set::error::ChallengeSetError::TooFewBytes(bs.len(), Self::BYTES_NEEDED));
        }

        let mut coeffs: Vec<Fq> = Vec::with_capacity(24);

        for i in 0..6 {
            let x0: i16 = ((bs[3 * i] & 0b0011_1111)) as i16 - MAX_COEFF;
            let x1: i16 = (((bs[3 * i] & 0b1100_0000) >> 6) | ((bs[3 * i + 1] & 0b0000_1111) << 2)) as i16 - MAX_COEFF;
            let x2: i16 = (((bs[3 * i + 1] & 0b1111_0000) >> 4) | ((bs[3 * i + 2] & 0b0000_0011) << 4)) as i16 - MAX_COEFF;
            let x3: i16 = (((bs[3 * i + 2] & 0b1111_1100) >> 2)) as i16 - MAX_COEFF;

            coeffs.extend_from_slice(&[Fq::from(x0), Fq::from(x1), Fq::from(x2), Fq::from(x3)]);
        }
        
        Ok(GoldilocksRingPoly::from(coeffs))
    }
}
