use lattirust_ring::cyclotomic_ring::models::pow2_debug::{
    Pow2CyclotomicPolyRing, Pow2CyclotomicPolyRingNTT,
};

// Some classic lattice parameter sets.
pub const DILITHIUM_PRIME: u64 = 0x00000000_007FE001; // = 8380417

pub type DilithiumPoly = Pow2CyclotomicPolyRing<DILITHIUM_PRIME, 256>;
pub type DilithiumNTT = Pow2CyclotomicPolyRingNTT<DILITHIUM_PRIME, 256>;
