use ark_crypto_primitives::sponge::Absorb;
use ark_ff::Field;
use lattirust_ring::{
    balanced_decomposition::Decompose,
    cyclotomic_ring::models::pow2_debug::{Pow2CyclotomicPolyRing, Pow2CyclotomicPolyRingNTT},
    OverField, PolyRing,
};

pub mod pgold;

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
}

impl<const Q: u64, const N: usize> SuitableRing for Pow2CyclotomicPolyRingNTT<Q, N>
where
    Pow2CyclotomicPolyRingNTT<Q, N>: From<Pow2CyclotomicPolyRing<Q, N>>,
    Pow2CyclotomicPolyRing<Q, N>: From<Pow2CyclotomicPolyRingNTT<Q, N>>,
{
    type CoefficientRepresentation = Pow2CyclotomicPolyRing<Q, N>;
}

mod test {
    use lattirust_ring::Ring;

    use super::*;

    #[test]
    fn test_trait_implementation() {
        fn takes_suitable_ring<R: SuitableRing>(_x: R) {}

        let x: Pow2CyclotomicPolyRingNTT<17, 4> = Pow2CyclotomicPolyRingNTT::ZERO;

        takes_suitable_ring(x);

        assert!(true)
    }
}
