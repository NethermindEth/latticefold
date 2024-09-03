// PGold = 2^64 − 2^32 + 1

use ark_ff::Field;
use lattirust_arithmetic::{
    challenge_set::latticefold_challenge_set::LatticefoldChallengeSet,
    ring::{CyclotomicPolyRingSplittedNTT, PolynomialRing, Zq},
};

const Q: u64 = 18446744069414584321;
const D: usize = 120;
const Z: usize = 225;
const PHI_Z: usize = 120;
// zth root of unity
const ROU: u64 = 12030049465912850429;

pub type PGoldCyclotomicRing<const N: usize> =
    CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>;

#[allow(dead_code)]
pub struct PGoldChallengeSet<const N: usize>;

impl<const N: usize> LatticefoldChallengeSet<CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>>
    for PGoldChallengeSet<N>
{
    fn small_challenge_coefficient_from_random_bytes(
        _i: usize,
        bs: &[u8],
    ) -> <PGoldCyclotomicRing<N> as PolynomialRing>::BaseRing {
        if bs[0] == 0 {
            <Zq<Q> as Field>::ZERO
        } else {
            <Zq<Q> as Field>::ONE
        }
    }
}
