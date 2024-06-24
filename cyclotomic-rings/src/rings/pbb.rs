// pBB = 15 · 2^27 + 1

use lattirust_arithmetic::{
    ntt::NTT,
    ring::{ Zq, Pow2CyclotomicPolyRingNTT, Pow2CyclotomicPolyRing },
};
use super::CyclotomicRing;
use rand::Rng;
const Q: u64 = 15 * (1 << 27) + 1;

type ZqQ = Zq<Q>;

pub struct PBBCyclotomicRing<const N: usize>(Pow2CyclotomicPolyRing<ZqQ, N>);

impl<const N: usize> CyclotomicRing<Q> for PBBCyclotomicRing<N> {
    // Challenge can be any polynomial with degree up to 120
    fn get_challenge_set(&self) -> Vec<u8> {
        let mut rng = rand::thread_rng();
        let mut random_bytes = [0u8; 15];
        rng.fill(&mut random_bytes);

        // Convert the bytes to bits
        let mut bits = Vec::new();
        for byte in random_bytes.iter() {
            for i in 0..8 {
                bits.push((byte >> (7 - i)) & 1);
            }
        }
        return bits;
    }

    fn to_ntt(&self) -> Vec<ZqQ> {
        let ntt_ring = Pow2CyclotomicPolyRingNTT::from(self.0);
        ntt_ring.ntt_coeffs()
    }
}
