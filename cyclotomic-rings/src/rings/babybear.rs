use stark_rings::cyclotomic_ring::models::babybear::{Fq, RqNTT, RqPoly};

use super::SuitableRing;
use crate::ark_base::*;
use crate::challenge_set::error;
use crate::challenge_set::LatticefoldChallengeSet;

/// BabyBear ring in the NTT form.
///
/// The base field of the NTT form is a degree-9
/// extension of the BabyBear field.
///
/// The NTT form has 8 components.
pub type BabyBearRingNTT = RqNTT;

/// BabyBear ring in the coefficient form.
///
/// The cyclotomic polynomial is $X^72 - X^36 + 1$ of degree 72.
pub type BabyBearRingPoly = RqPoly;

impl SuitableRing for BabyBearRingNTT {
    type CoefficientRepresentation = RqPoly;
    type PoseidonParams = BabyBearPoseidonConfig;
}

pub struct BabyBearPoseidonConfig;

#[derive(Clone)]
pub struct BabyBearChallengeSet;

const MAX_COEFF: i16 = 32;

/// For Babybear prime the challenge set is the set of all
/// ring elements whose coefficients are in the range [-32, 32[.
impl LatticefoldChallengeSet<BabyBearRingNTT> for BabyBearChallengeSet {
    /// To generate an element in [-32, 32[ it is enough to use 6 bits.
    /// Thus to generate 24 coefficients in that range 18 bytes is enough.
    const BYTES_NEEDED: usize = 18;

    fn short_challenge_from_random_bytes(
        bs: &[u8],
    ) -> Result<BabyBearRingPoly, error::ChallengeSetError> {
        if bs.len() != Self::BYTES_NEEDED {
            return Err(error::ChallengeSetError::TooFewBytes(
                bs.len(),
                Self::BYTES_NEEDED,
            ));
        }

        let mut coeffs: Vec<Fq> = Vec::with_capacity(24);

        for i in 0..6 {
            let x0: i16 = (bs[3 * i] & 0b0011_1111) as i16 - MAX_COEFF;
            let x1: i16 = (((bs[3 * i] & 0b1100_0000) >> 6) | ((bs[3 * i + 1] & 0b0000_1111) << 2))
                as i16
                - MAX_COEFF;
            let x2: i16 = (((bs[3 * i + 1] & 0b1111_0000) >> 4)
                | ((bs[3 * i + 2] & 0b0000_0011) << 4)) as i16
                - MAX_COEFF;
            let x3: i16 = ((bs[3 * i + 2] & 0b1111_1100) >> 2) as i16 - MAX_COEFF;

            coeffs.extend_from_slice(&[Fq::from(x0), Fq::from(x1), Fq::from(x2), Fq::from(x3)]);
        }

        Ok(BabyBearRingPoly::from(coeffs))
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::BigInt;
    use stark_rings::cyclotomic_ring::models::babybear::Fq;

    use super::*;

    #[test]
    fn test_small_challenge_from_random_bytes() {
        let challenge = BabyBearChallengeSet::short_challenge_from_random_bytes(&[
            0x7b, 0x4b, 0xe5, 0x8e, 0xe5, 0x11, 0xd2, 0xd0, 0x9c, 0x22, 0xba, 0x2e, 0xeb, 0xa8,
            0xba, 0x35, 0xf2, 0x18,
        ])
        .unwrap();

        let res_coeffs: Vec<Fq> = vec![
            Fq::new(BigInt([27])),
            Fq::new(BigInt([13])),
            Fq::new(BigInt([2013265909])),
            Fq::new(BigInt([25])),
            Fq::new(BigInt([2013265903])),
            Fq::new(BigInt([2013265911])),
            Fq::new(BigInt([2013265919])),
            Fq::new(BigInt([2013265893])),
            Fq::new(BigInt([2013265907])),
            Fq::new(BigInt([2013265892])),
            Fq::new(BigInt([2013265902])),
            Fq::new(BigInt([7])),
            Fq::new(BigInt([2])),
            Fq::new(BigInt([8])),
            Fq::new(BigInt([11])),
            Fq::new(BigInt([2013265900])),
            Fq::new(BigInt([11])),
            Fq::new(BigInt([3])),
            Fq::new(BigInt([10])),
            Fq::new(BigInt([14])),
            Fq::new(BigInt([21])),
            Fq::new(BigInt([2013265897])),
            Fq::new(BigInt([2013265904])),
            Fq::new(BigInt([2013265895])),
        ];
        let expected = BabyBearRingPoly::from(res_coeffs);

        assert_eq!(expected, challenge)
    }
}
