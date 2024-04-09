mod sparse;

use fields::Field;

pub trait AdditiveGroupElement {
    fn add(&mut self, z: &mut Self) -> &mut Self;

    fn add_invert(&mut self) -> &mut Self;
}

pub trait MultiplicativeGroupElement {
    fn mul(&mut self, z: &mut Self) -> &mut Self;

    fn mul_invert(&mut self) -> &mut Self;
}

pub trait CyclotomicFieldElement: AdditiveGroupElement + MultiplicativeGroupElement {
    // Add NTT
}

/// This corresponds to cyclotomic fields of prime degree only
pub trait PrimeCyclotomicFieldElement<F>
where
    F: Field,
{
    /// Returns $\zeta_n^k$
    fn e(n: usize, k: usize) -> Self;

    /// Multiplies in place by scalar
    fn scalar_mul(&mut self, scalar: &F) -> &mut Self;

    /// Gives zero expressed as an element of $\mathbb{F}(\zeta_n)$
    fn zero_order(n: usize) -> Self;

    /// Gives one expressed as an element of $\mathbb{F}(\zeta_n)$
    fn one_order(n: usize) -> Self;
}
