// PGold = 2^64 − 2^32 + 1
use lattirust_ring::cyclotomic_ring::models::frog_ring::{Fq, RqNTT, RqPoly};

use crate::challenge_set::{error, LatticefoldChallengeSet};

use super::SuitableRing;

pub type FrogRingNTT = RqNTT;
pub type FrogRingPoly = RqPoly;

impl SuitableRing for FrogRingNTT {
    type CoefficientRepresentation = RqPoly;
    type PoseidonParams = FrogPoseidonConfig;
}

pub struct FrogPoseidonConfig;

pub struct FrogChallengeSet;

/// For Frog prime the challenge set is the set of all
/// ring elements whose coefficients are in the range [-128, 128[.
impl LatticefoldChallengeSet<FrogRingNTT> for FrogChallengeSet {
    const BYTES_NEEDED: usize = 16;

    fn small_challenge_from_random_bytes(
        bs: &[u8],
    ) -> Result<
        <FrogRingNTT as SuitableRing>::CoefficientRepresentation,
        crate::challenge_set::error::ChallengeSetError,
    > {
        if bs.len() != Self::BYTES_NEEDED {
            return Err(error::ChallengeSetError::TooFewBytes(
                bs.len(),
                Self::BYTES_NEEDED,
            ));
        }

        Ok(FrogRingPoly::from(
            bs.iter()
                .map(|&x| Fq::from(x as i16 - 128))
                .collect::<Vec<Fq>>(),
        ))
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::BigInt;
    use lattirust_ring::cyclotomic_ring::models::frog_ring::Fq;

    use super::*;

    #[test]
    fn test_small_challenge_from_random_bytes() {
        let challenge = FrogChallengeSet::small_challenge_from_random_bytes(&[
            0x7b, 0x4b, 0xe5, 0x8e, 0xe5, 0x11, 0xd2, 0xd0, 0x9c, 0x22, 0xba, 0x2e, 0xeb, 0xa8,
            0xba, 0x35,
        ])
        .unwrap();

        let res_coeffs: Vec<Fq> = vec![
            Fq::new(BigInt([15912092521325583636])),
            Fq::new(BigInt([15912092521325583588])),
            Fq::new(BigInt([101])),
            Fq::new(BigInt([14])),
            Fq::new(BigInt([101])),
            Fq::new(BigInt([15912092521325583530])),
            Fq::new(BigInt([82])),
            Fq::new(BigInt([80])),
            Fq::new(BigInt([28])),
            Fq::new(BigInt([15912092521325583547])),
            Fq::new(BigInt([58])),
            Fq::new(BigInt([15912092521325583559])),
            Fq::new(BigInt([107])),
            Fq::new(BigInt([40])),
            Fq::new(BigInt([58])),
            Fq::new(BigInt([15912092521325583566])),
        ];

        let expected = FrogRingPoly::from(res_coeffs);

        assert_eq!(expected, challenge)
    }
}
