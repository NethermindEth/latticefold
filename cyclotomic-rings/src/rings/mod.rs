use lattirust_arithmetic::ring::Zq;
use num_bigint::BigInt;

mod pbb;
mod pstark;
mod pgold;
mod pm31;

pub trait CyclotomicRing<const Q: u64> {
    fn get_challenge_set(&self) -> BigInt;
    fn to_ntt(&self) -> Vec<Zq<Q>>;
}
