// pBB = 15 · 2^27 + 1

use super::PrimeCyclotomicRing;
use lattirust_arithmetic::ring::{ Zq, CyclotomicPolyRingSplittedNTT };
use rand::Rng;
use std::ops::{ Deref, DerefMut };
const Q: u64 = 18446744069414584321;
const D: usize = 120;
const Z: usize = 225;
const PHI_Z: usize = 120;

type ZqQ = Zq<Q>;
pub struct PGoldCyclotomicRing<const N: usize>(CyclotomicPolyRingSplittedNTT<Q, N, D, Z, PHI_Z>);

impl<const N: usize> PrimeCyclotomicRing<Q, N> for PGoldCyclotomicRing<N> {
    // Challenge can be any polynomial with degree up to 120
    fn get_challenge(&self) -> Vec<ZqQ> {
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

    fn try_challenge_from_random_bytes(&self, bytes: &[u8]) -> Vec<Zq<Q>> {
        assert!(bytes.len() >= 15);
        let mut bits = Vec::new();
        for byte in bytes.iter().take(15) {
            for i in 0..8 {
                bits.push(ZqQ::from((byte >> (7 - i)) & 1));
            }
        }
        return bits;
    }
}
impl<const N: usize> Deref for PGoldCyclotomicRing<N> {
    type Target = CyclotomicPolyRingSplittedNTT<Q, N, D, Z, PHI_Z>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<const N: usize> DerefMut for PGoldCyclotomicRing<N> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}
