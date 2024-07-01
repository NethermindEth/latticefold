// pBB = 15 · 2^27 + 1

use super::PrimeCyclotomicRing;
use lattirust_arithmetic::ring::{ Zq, CyclotomicPolyRingSplittedNTT };
use rand::Rng;
use lattirust_arithmetic::partial_ntt::PartialNTT;
use std::ops::{ Add, Mul };
const Q: u64 = 18446744069414584321;
const D: usize = 120;
const Z: usize = 225;
const PHI_Z: usize = 120;

type ZqQ = Zq<Q>;
pub struct PGoldCyclotomicRing<const N: usize>(CyclotomicPolyRingSplittedNTT<Q, N, D, Z, PHI_Z>);

impl<const N: usize> PrimeCyclotomicRing<Q, N> for PGoldCyclotomicRing<N> {
    // Challenge can be any polynomial with degree up to 120
    fn get_challenge_set(&self) -> Vec<ZqQ> {
        let mut rng = rand::thread_rng();
        let mut random_bytes = [0u8; 15];
        rng.fill(&mut random_bytes);

        // Convert the bytes to bits
        let mut bits = Vec::new();
        for byte in random_bytes.iter() {
            for i in 0..8 {
                bits.push(ZqQ::from((byte >> (7 - i)) & 1));
            }
        }
        return bits;
    }

    fn ntt(&self, a: &mut [Zq<Q>; N], rou: Zq<Q>) {
        CyclotomicPolyRingSplittedNTT::<Q, N, D, Z, PHI_Z>::ntt(a, rou);
    }
}

impl<const N: usize> Add for PGoldCyclotomicRing<N> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        PGoldCyclotomicRing(rhs.0 + self.0)
    }
}

impl<const N: usize> Mul for PGoldCyclotomicRing<N> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        PGoldCyclotomicRing(rhs.0 * self.0)
    }
}
