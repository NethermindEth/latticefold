// pM31 = 2^31 − 1

use lattirust_arithmetic::{
    ntt::NTT,
    ring::{ Zq, Pow2CyclotomicPolyRingNTT, Pow2CyclotomicPolyRing },
};
use super::CyclotomicRing;
use num_bigint::BigInt;

const Q: u64 = (1 << 31) - 1;

type ZqQ = Zq<Q>;

pub struct PM31CyclotomicRing<const N: usize>(Pow2CyclotomicPolyRing<ZqQ, N>);

impl<const N: usize> CyclotomicRing<Q> for PM31CyclotomicRing<N> {
    // Returns integer modulus of the set
    // TODO coverth this to the csmall set in the latticefold params document
    fn get_challenge_set(&self) -> BigInt {
        Q.into()
    }

    fn to_ntt(&self) -> Vec<ZqQ> {
        let ntt_ring = Pow2CyclotomicPolyRingNTT::from(self.0);
        ntt_ring.ntt_coeffs()
    }
}
