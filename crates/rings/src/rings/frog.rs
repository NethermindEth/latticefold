use stark_rings::cyclotomic_ring::models::frog_ring::{Fq, RqNTT, RqPoly};

use super::SuitableRing;
use crate::{
    ark_base::*,
    challenge_set::{error, LatticefoldChallengeSet},
};

/// Frog ring in the NTT form.
///
/// The base field of the NTT form is a degree-4
/// extension of the Frog field ($p=15912092521325583641$).
///
/// The NTT norm has 4 components.
pub type FrogRingNTT = RqNTT;

/// Frog ring in the coefficient form.
///
/// The cyclotomic polynomial is $X^16+1$ of degree 16.
pub type FrogRingPoly = RqPoly;

impl SuitableRing for FrogRingNTT {
    type CoefficientRepresentation = RqPoly;
    type PoseidonParams = FrogPoseidonConfig;
    type ChallengeSet = FrogChallengeSet;
}

pub struct FrogPoseidonConfig;

#[derive(Clone)]
pub struct FrogChallengeSet;

/// For Frog prime the challenge set is the set of all
/// ring elements whose coefficients are in the range [-128, 128[.
impl LatticefoldChallengeSet<FrogRingNTT> for FrogChallengeSet {
    const BYTES_NEEDED: usize = 16;

    fn short_challenge_from_random_bytes(
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
    use stark_rings::cyclotomic_ring::models::frog_ring::Fq;

    use super::*;

    #[test]
    fn test_small_challenge_from_random_bytes() {
        let challenge = FrogChallengeSet::short_challenge_from_random_bytes(&[
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
