use ark_crypto_primitives::sponge::{poseidon::PoseidonConfig, Absorb};
use ark_ff::PrimeField;
use ark_ff::{Field, Zero};
use lattirust_ring::zn::z_q::Zq;
use lattirust_ring::{
    balanced_decomposition::Decompose,
    cyclotomic_ring::models::pow2_debug::{Pow2CyclotomicPolyRing, Pow2CyclotomicPolyRingNTT},
    OverField, PolyRing,
};

mod frog;
mod goldilocks;
mod stark;

pub use frog::*;
pub use goldilocks::*;
pub use stark::*;

/// This trait should be used for rings in Latticefold.
/// It contains all the data needed in the protocol.
/// The type itself is meant to be the NTT-representation of a ring.
/// The associated type `CoefficientRepresentation` is the ring in the coefficient basis.
pub trait SuitableRing:
    OverField + From<Self::CoefficientRepresentation> + Into<Self::CoefficientRepresentation>
where
    <<Self as PolyRing>::BaseRing as Field>::BasePrimeField: Absorb,
{
    /// The coefficient basis version of the ring.
    type CoefficientRepresentation: OverField<BaseRing = <<Self as PolyRing>::BaseRing as Field>::BasePrimeField>
        + Decompose;
    type PoseidonParams: GetPoseidonParams<<<Self as PolyRing>::BaseRing as Field>::BasePrimeField>;
}

pub trait GetPoseidonParams<Fq: PrimeField> {
    fn get_poseidon_config() -> PoseidonConfig<Fq>;
}

pub struct Pow2Params<const Q: u64, const N: usize>;

impl<const Q: u64, const N: usize> GetPoseidonParams<Zq<Q>> for Pow2Params<Q, N> {
    fn get_poseidon_config() -> PoseidonConfig<Zq<Q>> {
        PoseidonConfig {
            full_rounds: 8, // Example values, adjust according to your needs
            partial_rounds: 57,
            alpha: 5,
            ark: vec![vec![Zq::zero(); 3]; 8 + 57], // Adjust to actual ark parameters
            mds: vec![vec![Zq::zero(); 3]; 3],      // Adjust to actual MDS matrix parameters
            rate: 2,
            capacity: 1,
        }
    }
}

impl<const Q: u64, const N: usize> SuitableRing for Pow2CyclotomicPolyRingNTT<Q, N>
where
    Pow2CyclotomicPolyRingNTT<Q, N>: From<Pow2CyclotomicPolyRing<Q, N>>,
    Pow2CyclotomicPolyRing<Q, N>: From<Pow2CyclotomicPolyRingNTT<Q, N>>,
{
    type CoefficientRepresentation = Pow2CyclotomicPolyRing<Q, N>;
    type PoseidonParams = Pow2Params<Q, N>;
}

#[cfg(test)]
mod tests {
    use lattirust_ring::Ring;

    use super::*;

    #[test]
    fn test_trait_implementation() {
        fn takes_suitable_ring<R: SuitableRing>(_x: R) {}

        let x: Pow2CyclotomicPolyRingNTT<17, 4> = Pow2CyclotomicPolyRingNTT::ZERO;

        takes_suitable_ring(x);
    }
}
