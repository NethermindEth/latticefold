// PGold = 2^64 − 2^32 + 1
use lattirust_ring::cyclotomic_ring::models::stark_prime::{Fq, RqNTT, RqPoly};

use crate::{challenge_set::error, challenge_set::LatticefoldChallengeSet};

use super::SuitableRing;

pub type StarkRingNTT = RqNTT;
pub type StarkRingPoly = RqPoly;

impl SuitableRing for StarkRingNTT {
    type CoefficientRepresentation = StarkRingPoly;

    type PoseidonParams = StarkPoseidonConfig;
}

pub struct StarkPoseidonConfig;

pub struct StarkChallengeSet;

/// Small challenges are the ring elements with coefficients in range [0; 2^8[.
impl LatticefoldChallengeSet<StarkRingNTT> for StarkChallengeSet {
    const BYTES_NEEDED: usize = 16;

    fn small_challenge_from_random_bytes(
        bs: &[u8],
    ) -> Result<<StarkRingNTT as SuitableRing>::CoefficientRepresentation, error::ChallengeSetError>
    {
        if bs.len() != Self::BYTES_NEEDED {
            return Err(error::ChallengeSetError::TooFewBytes(
                bs.len(),
                Self::BYTES_NEEDED,
            ));
        }

        Ok(StarkRingPoly::from(
            bs.iter().map(|&x| Fq::from(x)).collect::<Vec<Fq>>(),
        ))
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::BigInt;

    use super::*;

    #[test]
    fn test_small_challenge_from_random_bytes() {
        let challenge = StarkChallengeSet::small_challenge_from_random_bytes(&[
            0x7b, 0x4b, 0xe5, 0x8e, 0xe5, 0x11, 0xd2, 0xd0, 0x9c, 0x22, 0xba, 0x2e, 0xeb, 0xa8,
            0xba, 0x35,
        ])
        .unwrap();

        let res_coeffs: Vec<Fq> = vec![
            Fq::new(BigInt([123, 0, 0, 0])),
            Fq::new(BigInt([75, 0, 0, 0])),
            Fq::new(BigInt([229, 0, 0, 0])),
            Fq::new(BigInt([142, 0, 0, 0])),
            Fq::new(BigInt([229, 0, 0, 0])),
            Fq::new(BigInt([17, 0, 0, 0])),
            Fq::new(BigInt([210, 0, 0, 0])),
            Fq::new(BigInt([208, 0, 0, 0])),
            Fq::new(BigInt([156, 0, 0, 0])),
            Fq::new(BigInt([34, 0, 0, 0])),
            Fq::new(BigInt([186, 0, 0, 0])),
            Fq::new(BigInt([46, 0, 0, 0])),
            Fq::new(BigInt([235, 0, 0, 0])),
            Fq::new(BigInt([168, 0, 0, 0])),
            Fq::new(BigInt([186, 0, 0, 0])),
            Fq::new(BigInt([53, 0, 0, 0])),
        ];

        let expected = StarkRingPoly::from(res_coeffs);

        assert_eq!(expected, challenge)
    }
}