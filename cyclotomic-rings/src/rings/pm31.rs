use ark_ff::Field;
// PM31 = 2^31 -1
use lattirust_arithmetic::{
    challenge_set::latticefold_challenge_set::LatticefoldChallengeSet,
    ring::{CyclotomicPolyRingSplittedNTT, Zq},
};

const Q: u64 = (1 << 31) - 1;
const D: usize = 84;
const Z: usize = 225;
const PHI_Z: usize = 120;
// zth root of unity
const ROU: u64 = 309107220;

pub type PM31CyclotomicRing<const N: usize> = CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>;

#[allow(dead_code)]
pub struct PM31ChallengeSet<const N: usize>;
impl<const N: usize> LatticefoldChallengeSet<CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>>
    for PM31ChallengeSet<N>
{
    fn small_challenge_coefficient_from_random_bytes(
        _i: usize,
        bs: &[u8],
    ) -> <PM31CyclotomicRing<N> as lattirust_arithmetic::ring::PolyRing>::BaseRing {
        if bs[0] == 0 {
            <Zq<Q> as Field>::ZERO
        } else {
            <Zq<Q> as Field>::ONE
        }
    }
}
