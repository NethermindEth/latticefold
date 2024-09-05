use lattirust_ring::{
    balanced_decomposition::Decompose, OverField, PolyRing, Pow2CyclotomicPolyRing,
    Pow2CyclotomicPolyRingNTT, Zq,
};

pub mod pbb;
pub mod pgold;
pub mod pm31;

/// This trait should be used for rings in Latticefold.
/// It contains all the data needed in the protocol.
/// The type itself is meant to be the NTT-representation of a ring.
/// The associated type `CoefficientRepresentation` is the ring in the coefficient basis.
pub trait SuitableRing:
    OverField + From<Self::CoefficientRepresentation> + Into<Self::CoefficientRepresentation>
{
    /// The coefficient basis version of the ring.
    type CoefficientRepresentation: PolyRing<BaseRing = Self::BaseRing> + Decompose;
}

impl<const Q: u64, const N: usize> SuitableRing for Pow2CyclotomicPolyRingNTT<Q, N> {
    type CoefficientRepresentation = Pow2CyclotomicPolyRing<Zq<Q>, N>;
}
