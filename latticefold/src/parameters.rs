use lattirust_arithmetic::ring::{Pow2CyclotomicPolyRing, Pow2CyclotomicPolyRingNTT, Zq};
use std::fmt::Display;

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

pub type DilithiumCR = Pow2CyclotomicPolyRing<Zq<DILITHIUM_PRIME>, 256>;
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

pub const GOLDILOCKS_PRIME: u64 = (1 << 32) * ((1 << 32) - 1) + 1;

pub type GoldilocksCR = Pow2CyclotomicPolyRing<Zq<GOLDILOCKS_PRIME>, 256>;
pub type GoldilocksNTT = Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 256>;

#[derive(Clone, Copy)]
pub struct GoldilocksTestParams;

// TODO: Revise this later
impl DecompositionParams for GoldilocksTestParams {
    const B: u128 = 1 << 63; // log2(GOLDILOCKS) ~ 64
    const L: usize = 2;
    const B_SMALL: u128 = 2;
    const K: usize = 13;
}

pub const BABYBEAR_PRIME: u64 = 15 * (1 << 27) + 1;

pub type BabyBearCR = Pow2CyclotomicPolyRing<Zq<BABYBEAR_PRIME>, 256>;
pub type BabyBearNTT = Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 256>;

#[derive(Clone, Copy)]
pub struct BabyBearTestParams;

// TODO: Revise this later
impl DecompositionParams for BabyBearTestParams {
    const B: u128 = 1 << 30; // log2(BABYBEAR) ~ 31
    const L: usize = 2;
    const B_SMALL: u128 = 2;
    const K: usize = 13;
}

// p = 27*2^59 + 1
pub const POW2_59_PRIME: u64 = 0xd800000000000001;

pub type POW2_59CR = Pow2CyclotomicPolyRing<Zq<POW2_59_PRIME>, 256>;
pub type POW2_59NTT = Pow2CyclotomicPolyRingNTT<POW2_59_PRIME, 256>;
#[derive(Clone, Copy)]
pub struct Pow2_59TestParams;

// TODO: Revise this later
impl DecompositionParams for Pow2_59TestParams {
    const B: u128 = 1 << 13;
    const L: usize = 2;
    const B_SMALL: u128 = 2;
    const K: usize = 13;
}

pub const POW2_57_PRIME: u64 = 0xf600000000000001;

pub type POW2_57CR = Pow2CyclotomicPolyRing<Zq<POW2_57_PRIME>, 256>;
pub type POW2_57NTT = Pow2CyclotomicPolyRingNTT<POW2_57_PRIME, 256>;

#[derive(Clone, Copy)]
pub struct Pow2_57TestParams;

// TODO: Revise this later
impl DecompositionParams for Pow2_57TestParams {
    const B: u128 = 1 << 13;
    const L: usize = 2;
    const B_SMALL: u128 = 2;
    const K: usize = 13;
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
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "B={}, l={}", self.b, self.l,)
    }
}
