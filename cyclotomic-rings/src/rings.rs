use ark_crypto_primitives::sponge::{poseidon::PoseidonConfig, Absorb};
use ark_ff::PrimeField;
use ark_ff::{Field, Zero};
use lattirust_ring::{
    balanced_decomposition::Decompose,
    cyclotomic_ring::models::pow2_debug::{Pow2CyclotomicPolyRing, Pow2CyclotomicPolyRingNTT},
    zn::z_q::Zq,
    Cyclotomic, OverField, PolyRing,
};

mod babybear;
mod frog;
mod goldilocks;
mod poseidon;
mod stark;

pub use babybear::*;
pub use frog::*;
pub use goldilocks::*;
pub use stark::*;

/// An umbrella trait of a ring suitable to be used in the LatticeFold protocol.
///
/// The ring is assumed to be of the form $$\mathbb{Z}_p\[X\]/(f(X)),$$ for a polynomial
/// $f(X) \in \mathbb{Z}_p\[X\],\: d=\mathrm{deg}\ f(X)$, (typically, this is a cyclotomic polynomial $\Phi_m(X)$) so it has
/// two isomorphic forms:
///   * <i>The coefficient form</i>, i.e. a ring element is represented as the unique polynomial $g$ of the
///     degree $\mathrm{deg}\ g < d$.
///   * <i>The NTT form</i>, i.e. a ring element is represented as its image along the Chinese-remainder isomorphism
///     $$\mathbb{Z}_p\[X\]/(f(X))\cong \prod\limits\_{i=1}^k\mathbb{Z}_p\[X\]/(f\_i(X)),$$
///     where $f\_1(X),\ldots, f\_k(X)$ are irreducible polynomials in $ \mathbb{Z}_p\[X\]$ such that
///     $$f(X) = f\_1(X)\cdot\ldots\cdot f\_k(X).$$
///
/// When $f(X)$ is a cyclotomic polynomial the factors $f\_1(X),\ldots, f\_k(X)$ have equal degrees, thus the fields in the RHS of
/// the Chinese-remainder isomorphism are all isomorphic to the same extension of the field $\mathbb{Z}\_p$, implying the NTT form
/// of the ring is a direct product of $k$ instances of $\mathbb{Z}\_{p^\tau}$ for some $\tau$ with componentwise operations.
///
/// If `R: SuitableRing` then we assume that the type `R` represents the NTT form of the ring as the arithmetic operations
/// in the NTT form are much faster and we intend to use the NTT form as much as possible only occasionally turning to the
/// coefficient form (usually, when Ajtai security aspects are discussed). The associated type `CoefficientRepresentation` is the corresponding
/// coefficient form representation of the ring.
///
/// A type `R: SuitableRing` and its `R::CoefficientRepresentation` has to satisfy the following conditions:
///   * `R` has to be an `OverField` to exhibit an algebra over a field `R::BaseRing` structure.  
///   * `R::CoefficientRepresentation` has to be an algebra over the prime field `R::BaseRing::BasePrimeField` of the field `R::BaseRing`.
///   * `R::BaseRing::BasePrimeField` has to be absorbable by sponge hashes (`R::BaseRing::BasePrimeField: Absorb`).
///   * `R` and `R::CoefficientRepresentation` should be convertible into each other.
///   * `R::CoefficientRepresentation` is radix-$B$ decomposable and exhibits cyclotomic structure (`R::CoefficientRepresentation: Decompose + Cyclotomic`).
///
/// In addition to the data above a suitable ring has to provide Poseidon hash parameters for its base prime field (i.e. $\mathbb{Z}\_p$).
pub trait SuitableRing:
    OverField + From<Self::CoefficientRepresentation> + Into<Self::CoefficientRepresentation>
where
    <<Self as PolyRing>::BaseRing as Field>::BasePrimeField: Absorb,
{
    /// The coefficient form version of the ring.
    type CoefficientRepresentation: OverField<BaseRing = <<Self as PolyRing>::BaseRing as Field>::BasePrimeField>
        + Decompose
        + Cyclotomic;

    /// Poseidon sponge parameters for the base prime field.
    type PoseidonParams: GetPoseidonParams<<<Self as PolyRing>::BaseRing as Field>::BasePrimeField>;
}

/// A trait for types with an associated Poseidon sponge configuration.
pub trait GetPoseidonParams<Fq: PrimeField> {
    /// Returns the associated Poseidon sponge configuration.
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
