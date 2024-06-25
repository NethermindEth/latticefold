// pgold = 2^64 − 2^32 + 1

use lattirust_arithmetic::{
    ntt::NTT,
    ring::{ Pow2CyclotomicPolyRing, Pow2CyclotomicPolyRingNTT, ZPGold, Zq },
};
use super::CyclotomicRing;
use rand::Rng;

pub struct PGoldCyclotomicRing<const N: usize>(Pow2CyclotomicPolyRing<ZPGold, N>);
impl<const N: usize> CyclotomicRing<ZPGold> for PGoldCyclotomicRing<N> {
    fn get_challenge_set(&self) -> Vec<u8> {
        // Challenge can be any polynomial with degree up to 84
        // Coefficients can only be one or zero
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

    fn to_ntt(&self) -> Vec<ZPGold> {
        todo!()
    }
}
