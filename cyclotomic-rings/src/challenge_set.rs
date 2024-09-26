use ark_ff::Field;

use lattirust_ring::{OverField, cyclotomic_ring::models::pow2_debug::{Pow2CyclotomicPolyRingNTT, Fp64Pow2}};

pub trait LatticefoldChallengeSet<R: OverField> {
    fn small_challenge_coefficient_from_random_bytes(i: usize, bs: &[u8]) -> R::BaseRing;
    fn small_challenge_from_random_bytes(bs: &[u8]) -> R {
        <R as From<Vec<R::BaseRing>>>::from(
            (0..R::dimension())
                .map(|i| Self::small_challenge_coefficient_from_random_bytes(i, &bs[i..]))
                .collect(),
        )
    }
}

pub struct BinarySmallSet<const Q: u64, const N: usize>;

impl<const Q: u64, const N: usize> LatticefoldChallengeSet<Pow2CyclotomicPolyRingNTT<Q, N>>
    for BinarySmallSet<Q, N>
{
    fn small_challenge_coefficient_from_random_bytes(_i: usize, bs: &[u8]) -> Fp64Pow2<Q, N> {
        if bs[0] == 0 {
            <Fp64Pow2<Q, N> as Field>::ZERO
        } else {
            <Fp64Pow2<Q, N> as Field>::ONE
        }
    }
}
