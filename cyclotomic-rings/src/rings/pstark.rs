// pSTARK = 2^251 + 17 · 2^192 + 1

use lattirust_arithmetic::{
    ntt::NTT,
    ring::{ Pow2CyclotomicPolyRing, Pow2CyclotomicPolyRingNTT, ZPStark, Zq },
};
use super::CyclotomicRing;
use rand::Rng;

pub struct PStarkCyclotomicRing<const N: usize>(Pow2CyclotomicPolyRing<ZPStark, N>);
impl<const N: usize> CyclotomicRing<ZPStark> for PStarkCyclotomicRing<N> {
    fn get_challenge_set(&self) -> Vec<u8> {
        // Challenge can be any polynomial with degree up to 84
        // Coefficients can only be one or zero
        let mut rng = rand::thread_rng();
        let mut random_bytes = [0u8; 3];
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

    fn to_ntt(&self) -> Vec<ZPStark> {
        todo!()
    }
}
