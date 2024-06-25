use lattirust_arithmetic::ring::{ ConvertibleRing, Zq };

mod pbb;
mod pstark;
mod pgold;
mod pm31;

pub trait CyclotomicRing<R: ConvertibleRing> {
    // Challenge is on the form of polynomial with 0 and 1 coefficients
    // TODO This is currently implemented as random zeroes and ones
    // Change this is so it uses a real hash function
    fn get_challenge_set(&self) -> Vec<u8>;
    fn to_ntt(&self) -> Vec<R>;
}
