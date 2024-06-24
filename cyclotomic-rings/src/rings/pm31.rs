// pM31 = 2^31 − 1

use lattirust_arithmetic::{
    ntt::NTT,
    ring::{ Zq, Pow2CyclotomicPolyRingNTT, Pow2CyclotomicPolyRing },
};
use super::CyclotomicRing;
use rand::Rng;

const Q: u64 = (1 << 31) - 1;

type ZqQ = Zq<Q>;

pub struct PM31CyclotomicRing<const N: usize>(Pow2CyclotomicPolyRing<ZqQ, N>);

impl<const N: usize> CyclotomicRing<Q> for PM31CyclotomicRing<N> {
    // Challenge can be any polynomial with degree up to 84
    fn get_challenge_set(&self) -> Vec<u8> {
        let mut rng = rand::thread_rng();
        let mut random_bytes = [0u8; 11];
        rng.fill(&mut random_bytes);

        // Convert the bytes to bits
        let mut bits = Vec::new();
        for byte in random_bytes.iter() {
            for i in 0..8 {
                bits.push((byte >> (7 - i)) & 1);
            }
        }

        // Trim the bits to 84
        bits.truncate(84);
        return bits;
    }

    fn to_ntt(&self) -> Vec<ZqQ> {
        let ntt_ring = Pow2CyclotomicPolyRingNTT::from(self.0);
        ntt_ring.ntt_coeffs()
    }
}
