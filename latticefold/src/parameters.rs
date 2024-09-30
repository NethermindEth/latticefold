use ark_std::fmt::Display;
use lattirust_ring::cyclotomic_ring::models::pow2_debug::{
    Pow2CyclotomicPolyRing, Pow2CyclotomicPolyRingNTT,
};

/// Decomposition parameters.
/// Convenient to enforce them compile-time.
/// Contains both gadget matrix data and Latticefold decomposition step data.
pub trait DecompositionParams: Clone {
    /// The MSIS bound.
    const B: u128;
    /// The ring modulus should be < B^L.
    const L: usize;
    /// The small b from the decomposition step of LF.
    const B_SMALL: u128;
    /// K = log_b B.
    const K: usize;
}

// Some classic lattice parameter sets.

pub const DILITHIUM_PRIME: u64 = 0x00000000_007FE001;

pub type DilithiumCR = Pow2CyclotomicPolyRing<DILITHIUM_PRIME, 256>;
pub type DilithiumNTT = Pow2CyclotomicPolyRingNTT<DILITHIUM_PRIME, 256>;

#[derive(Clone, Copy)]
pub struct DilithiumTestParams;

// TODO: Revise this later
impl DecompositionParams for DilithiumTestParams {
    const B: u128 = 1 << 13;
    const L: usize = 2;
    const B_SMALL: u128 = 2;
    const K: usize = 13;
}
#[derive(Clone, Copy)]
pub struct StarkPrimeParams;

impl DecompositionParams for StarkPrimeParams {
    const B: u128 = 1 << 54;
    const L: usize = 2;
    const B_SMALL: u128 = 2;
    const K: usize = 54;
}

#[derive(Clone, Copy)]
pub struct GoldilocksParams;

impl DecompositionParams for GoldilocksParams {
    const B: u128 = 1 << 24;
    const L: usize = 2;
    const B_SMALL: u128 = 2;
    const K: usize = 24;
}

#[derive(Clone, Copy)]
pub struct FrogParams;

impl DecompositionParams for crate::parameters::FrogParams {
    const B: u128 = 1 << 16;
    const L: usize = 2;
    const B_SMALL: u128 = 2;
    const K: usize = 16;
}

impl<P: DecompositionParams> From<P> for DecompositionParamData {
    fn from(_: P) -> Self {
        {
            Self { b: P::B, l: P::L }
        }
    }
}

// Nice representation of parameters for printing out in benchmarks.
#[derive(Clone, Copy)]
pub struct DecompositionParamData {
    // The MSIS bound.
    b: u128,
    // The ring modulus should be < B^L.
    l: usize,
}

impl Display for DecompositionParamData {
    fn fmt(&self, f: &mut ark_std::fmt::Formatter<'_>) -> ark_std::fmt::Result {
        write!(f, "B={}, l={}", self.b, self.l,)
    }
}
