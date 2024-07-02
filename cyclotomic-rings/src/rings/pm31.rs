// pBB = 15 · 2^27 + 1
use super::PrimeCyclotomicRing;
use lattirust_arithmetic::ring::{ Zq, CyclotomicPolyRingSplittedNTT };
use rand::Rng;
use lattirust_arithmetic::partial_ntt::PartialNTT;
use std::ops::{ Add, Mul };
const Q: u64 = (1 << 31) - 1;
const D: usize = 84;
const Z: usize = 225;
const PHI_Z: usize = 120;

type ZqQ = Zq<Q>;
pub struct PM31CyclotomicRing<const N: usize>(CyclotomicPolyRingSplittedNTT<Q, N, D, Z, PHI_Z>);

impl<const N: usize> PrimeCyclotomicRing<Q, N> for PM31CyclotomicRing<N> {
    // Challenge can be any polynomial with degree up to 120
    fn get_challenge(&self) -> Vec<ZqQ> {
        let mut rng = rand::thread_rng();
        let mut random_bytes = [0u8; 7];
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

    fn get_challenge_from_random_bytes(&self, bytes: &[u8]) -> Vec<Zq<Q>> {
        todo!()
    }
}

impl<const N: usize> Add for PM31CyclotomicRing<N> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        PM31CyclotomicRing(rhs.0 + self.0)
    }
}

impl<const N: usize> Mul for PM31CyclotomicRing<N> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        PM31CyclotomicRing(rhs.0 * self.0)
    }
}
