use lattirust_arithmetic::ring::Zq;

mod pbb;
mod pgold;
mod pm31;

pub trait PrimeCyclotomicRing<const Q: u64, const N: usize> {
    // Challenge is on the form of polynomial with 0 and 1 coefficients
    // TODO This is currently implemented as random zeroes and ones
    // Change this is so it uses a real hash function
    fn get_challenge(&self) -> Vec<Zq<Q>>;
    fn try_challenge_from_random_bytes(&self, bytes: &[u8]) -> Vec<Zq<Q>>;
    fn ntt(&self, a: &mut [Zq<Q>; N], rou: Zq<Q>) -> ();
}
