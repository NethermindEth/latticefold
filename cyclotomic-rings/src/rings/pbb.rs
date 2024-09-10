// pBB = 15 · 2^27 + 1

use ark_ff::Field;
use lattirust_ring::{CyclotomicPolyRingSplittedNTT, PolyRing, Zq};

use crate::challenge_set::LatticefoldChallengeSet;

const Q: u64 = 15 * (1 << 27) + 1;
const D: usize = 120;
const Z: usize = 225;
const PHI_Z: usize = 120;
// zth root of unity
const ROU: u64 = 1995471372;

pub type PBBCyclotomicRing<const N: usize> = CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>;
#[allow(dead_code)]
pub struct PBBChallengeSet<const N: usize>;
impl<const N: usize> LatticefoldChallengeSet<CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>>
    for PBBChallengeSet<N>
{
    fn small_challenge_coefficient_from_random_bytes(
        _i: usize,
        bs: &[u8],
    ) -> <PBBCyclotomicRing<N> as PolyRing>::BaseRing {
        if bs[0] == 0 {
            <Zq<Q> as Field>::ZERO
        } else {
            <Zq<Q> as Field>::ONE
        }
    }
}
