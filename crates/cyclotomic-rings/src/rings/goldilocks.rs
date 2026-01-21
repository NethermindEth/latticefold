use ark_ff::Field;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use stark_rings::balanced_decomposition::Decompose;
use stark_rings::cyclotomic_ring::models::goldilocks::{Fq, FqConfig, RqNTT, RqPoly};
use stark_rings::cyclotomic_ring::Flatten;
use stark_rings::traits::FromRandomBytes;
use stark_rings::traits::MulUnchecked;
use stark_rings::{OverField, PolyRing, Ring};

use super::SuitableRing;
use crate::{
    ark_base::*,
    challenge_set::{error, LatticefoldChallengeSet},
};

/// Goldilocks ring in the NTT form.
///
/// The base field of the NTT form is a degree-3
/// extension of the Goldilocks field.
///
/// The NTT form has 8 components.
pub type GoldilocksRingNTT = RqNTT;

/// BabyBear ring in the coefficient form.
///
/// The cyclotomic polynomial is $X^24-X^12+1$ of degree 24.
pub type GoldilocksRingPoly = RqPoly;

impl SuitableRing for GoldilocksRingNTT {
    type CoefficientRepresentation = RqPoly;
    type PoseidonParams = GoldilocksPoseidonConfig;
}

pub struct GoldilocksPoseidonConfig;

#[derive(Clone)]
pub struct GoldilocksChallengeSet;

const MAX_COEFF: i16 = 32;

/// For Goldilocks prime the challenge set is the set of all
/// ring elements whose coefficients are in the range [-32, 32[.
impl LatticefoldChallengeSet<GoldilocksRingNTT> for GoldilocksChallengeSet {
    /// To generate an element in [-32, 32[ it is enough to use 6 bits.
    /// Thus to generate 24 coefficients in that range 18 bytes is enough.
    const BYTES_NEEDED: usize = 18;

    fn short_challenge_from_random_bytes(
        bs: &[u8],
    ) -> Result<GoldilocksRingPoly, error::ChallengeSetError> {
        if bs.len() != Self::BYTES_NEEDED {
            return Err(error::ChallengeSetError::TooFewBytes(
                bs.len(),
                Self::BYTES_NEEDED,
            ));
        }

        let mut coeffs: Vec<Fq> = Vec::with_capacity(24);

        for i in 0..6 {
            let x0: i16 = (bs[3 * i] & 0b0011_1111) as i16 - MAX_COEFF;
            let x1: i16 = (((bs[3 * i] & 0b1100_0000) >> 6) | ((bs[3 * i + 1] & 0b0000_1111) << 2))
                as i16
                - MAX_COEFF;
            let x2: i16 = (((bs[3 * i + 1] & 0b1111_0000) >> 4)
                | ((bs[3 * i + 2] & 0b0000_0011) << 4)) as i16
                - MAX_COEFF;
            let x3: i16 = ((bs[3 * i + 2] & 0b1111_1100) >> 2) as i16 - MAX_COEFF;

            coeffs.extend_from_slice(&[Fq::from(x0), Fq::from(x1), Fq::from(x2), Fq::from(x3)]);
        }

        Ok(GoldilocksRingPoly::from(coeffs))
    }
}

/// Goldilocks ring wrapper with **dimension 64** (coefficient representation):
/// \( \mathbb{F}_p[X]/(X^{64}+1) \) over the Goldilocks base prime field.
///
/// This is specifically intended for LF+/WE experiments where we want:
/// - the same power-of-two negacyclic shape as `FrogRing64`,
/// - but with Goldilocks' high 2-adicity to enable an NTT-based `ring_mul_negacyclic` gadget.
///
/// NOTE: Over Goldilocks, \(X^{64}+1\) *splits completely* (since 128 | (p-1)), so the ring is
/// isomorphic to \( \mathbb{F}_p^{64} \) in NTT form. This has security implications for schemes
/// whose hardness depends on extension-field size in the CRT/NTT representation (see Neo paper).
#[repr(transparent)]
#[derive(
    Copy, Clone, Debug, Default, Eq, PartialEq, Hash, CanonicalSerialize, CanonicalDeserialize,
)]
pub struct GoldilocksRing64(
    pub stark_rings::cyclotomic_ring::CyclotomicPolyRingGeneral<Goldilocks64Config, 1, 64>,
);

/// Parameters for \( \mathbb{F}_p[X]/(X^{64}+1) \) over Goldilocks base prime field.
pub struct Goldilocks64Config;

impl stark_rings::cyclotomic_ring::CyclotomicConfig<1> for Goldilocks64Config {
    type BaseFieldConfig = ark_ff::MontBackend<FqConfig, 1>;
    type BaseCRTField = Fq;
    const CRT_FIELD_EXTENSION_DEGREE: usize = 1;

    fn reduce_in_place(coefficients: &mut Vec<Fq>) {
        // Reduce mod (X^64 + 1): x^{64+i} = -x^i.
        if coefficients.len() > 64 {
            let (left, right) = coefficients.split_at_mut(64);
            for (a, b) in left.iter_mut().zip(right.iter()) {
                *a -= *b;
            }
        }
        coefficients.resize(64, <Fq as Field>::ZERO);
    }

    // For this configuration we treat CRT/ICRT as identity: we use coefficient representation
    // directly as the "NTT form" type for LF+ integration experiments.
    fn crt_in_place(_coefficients: &mut [Fq]) {}
    fn crt(coefficients: Vec<Fq>) -> Vec<Fq> {
        coefficients
    }
    fn icrt(evaluations: Vec<Fq>) -> Vec<Fq> {
        evaluations
    }
    fn icrt_in_place(_evaluations: &mut [Fq]) {}
}

// ---- Forwarding impls to satisfy `stark_rings::Ring` + `stark_rings::PolyRing` ----
macro_rules! fwd_binop {
    ($trait:ident, $method:ident, $op:tt) => {
        impl core::ops::$trait for GoldilocksRing64 {
            type Output = Self;
            #[inline(always)]
            fn $method(self, rhs: Self) -> Self::Output {
                Self(self.0 $op rhs.0)
            }
        }
        impl<'a> core::ops::$trait<&'a Self> for GoldilocksRing64 {
            type Output = Self;
            #[inline(always)]
            fn $method(self, rhs: &'a Self) -> Self::Output {
                Self(self.0 $op rhs.0)
            }
        }
        impl<'a> core::ops::$trait<&'a mut Self> for GoldilocksRing64 {
            type Output = Self;
            #[inline(always)]
            fn $method(self, rhs: &'a mut Self) -> Self::Output {
                Self(self.0 $op rhs.0)
            }
        }
    };
}
macro_rules! fwd_binop_assign {
    ($trait:ident, $method:ident, $op:tt) => {
        impl core::ops::$trait for GoldilocksRing64 {
            #[inline(always)]
            fn $method(&mut self, rhs: Self) {
                self.0 $op rhs.0;
            }
        }
        impl<'a> core::ops::$trait<&'a Self> for GoldilocksRing64 {
            #[inline(always)]
            fn $method(&mut self, rhs: &'a Self) {
                self.0 $op rhs.0;
            }
        }
        impl<'a> core::ops::$trait<&'a mut Self> for GoldilocksRing64 {
            #[inline(always)]
            fn $method(&mut self, rhs: &'a mut Self) {
                self.0 $op rhs.0;
            }
        }
    };
}

fwd_binop!(Add, add, +);
fwd_binop!(Sub, sub, -);
fwd_binop_assign!(AddAssign, add_assign, +=);
fwd_binop_assign!(SubAssign, sub_assign, -=);

impl core::ops::Mul for GoldilocksRing64 {
    type Output = Self;
    #[inline(always)]
    fn mul(self, rhs: Self) -> Self::Output {
        Self(self.0 * rhs.0)
    }
}
impl<'a> core::ops::Mul<&'a Self> for GoldilocksRing64 {
    type Output = Self;
    #[inline(always)]
    fn mul(self, rhs: &'a Self) -> Self::Output {
        Self(self.0 * rhs.0)
    }
}
impl<'a> core::ops::Mul<&'a mut Self> for GoldilocksRing64 {
    type Output = Self;
    #[inline(always)]
    fn mul(self, rhs: &'a mut Self) -> Self::Output {
        Self(self.0 * rhs.0)
    }
}
impl core::ops::MulAssign for GoldilocksRing64 {
    #[inline(always)]
    fn mul_assign(&mut self, rhs: Self) {
        self.0 *= rhs.0;
    }
}
impl<'a> core::ops::MulAssign<&'a Self> for GoldilocksRing64 {
    #[inline(always)]
    fn mul_assign(&mut self, rhs: &'a Self) {
        self.0 *= rhs.0;
    }
}
impl<'a> core::ops::MulAssign<&'a mut Self> for GoldilocksRing64 {
    #[inline(always)]
    fn mul_assign(&mut self, rhs: &'a mut Self) {
        self.0 *= rhs.0;
    }
}

impl core::ops::Neg for GoldilocksRing64 {
    type Output = Self;
    #[inline(always)]
    fn neg(self) -> Self::Output {
        Self(-self.0)
    }
}

impl ark_std::fmt::Display for GoldilocksRing64 {
    fn fmt(&self, f: &mut ark_std::fmt::Formatter<'_>) -> ark_std::fmt::Result {
        ark_std::fmt::Display::fmt(&self.0, f)
    }
}

impl ark_std::iter::Sum for GoldilocksRing64 {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        Self(iter.map(|x| x.0).sum())
    }
}
impl<'a> ark_std::iter::Sum<&'a Self> for GoldilocksRing64 {
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        Self(iter.map(|x| x.0).sum())
    }
}
impl ark_std::iter::Product for GoldilocksRing64 {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        Self(iter.map(|x| x.0).product())
    }
}
impl<'a> ark_std::iter::Product<&'a Self> for GoldilocksRing64 {
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        Self(iter.map(|x| x.0).product())
    }
}

impl ark_std::Zero for GoldilocksRing64 {
    #[inline(always)]
    fn zero() -> Self {
        Self(<stark_rings::cyclotomic_ring::CyclotomicPolyRingGeneral<Goldilocks64Config, 1, 64> as Ring>::ZERO)
    }
    #[inline(always)]
    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

impl ark_std::One for GoldilocksRing64 {
    #[inline(always)]
    fn one() -> Self {
        Self(<stark_rings::cyclotomic_ring::CyclotomicPolyRingGeneral<Goldilocks64Config, 1, 64> as Ring>::ONE)
    }
}

impl ark_std::UniformRand for GoldilocksRing64 {
    fn rand<R: rand::Rng + ?Sized>(rng: &mut R) -> Self {
        Self(stark_rings::cyclotomic_ring::CyclotomicPolyRingGeneral::<Goldilocks64Config, 1, 64>::rand(rng))
    }
}

impl FromRandomBytes<Self> for GoldilocksRing64 {
    fn byte_size() -> usize {
        stark_rings::cyclotomic_ring::CyclotomicPolyRingGeneral::<Goldilocks64Config, 1, 64>::byte_size()
    }

    fn try_from_random_bytes(bytes: &[u8]) -> Option<Self> {
        stark_rings::cyclotomic_ring::CyclotomicPolyRingGeneral::<Goldilocks64Config, 1, 64>::try_from_random_bytes(bytes)
            .map(Self)
    }
}

impl Ring for GoldilocksRing64 {
    const ZERO: Self = Self(<stark_rings::cyclotomic_ring::CyclotomicPolyRingGeneral<Goldilocks64Config, 1, 64> as Ring>::ZERO);
    const ONE: Self = Self(<stark_rings::cyclotomic_ring::CyclotomicPolyRingGeneral<Goldilocks64Config, 1, 64> as Ring>::ONE);
}

impl From<u128> for GoldilocksRing64 {
    fn from(value: u128) -> Self {
        Self(stark_rings::cyclotomic_ring::CyclotomicPolyRingGeneral::<Goldilocks64Config, 1, 64>::from(value))
    }
}
impl From<u64> for GoldilocksRing64 {
    fn from(value: u64) -> Self {
        Self::from(value as u128)
    }
}
impl From<u32> for GoldilocksRing64 {
    fn from(value: u32) -> Self {
        Self::from(value as u128)
    }
}
impl From<u16> for GoldilocksRing64 {
    fn from(value: u16) -> Self {
        Self::from(value as u128)
    }
}
impl From<u8> for GoldilocksRing64 {
    fn from(value: u8) -> Self {
        Self::from(value as u128)
    }
}
impl From<bool> for GoldilocksRing64 {
    fn from(value: bool) -> Self {
        Self::from(value as u128)
    }
}

impl From<Fq> for GoldilocksRing64 {
    fn from(value: Fq) -> Self {
        Self(stark_rings::cyclotomic_ring::CyclotomicPolyRingGeneral::<Goldilocks64Config, 1, 64>::from(value))
    }
}

impl From<Vec<Fq>> for GoldilocksRing64 {
    fn from(value: Vec<Fq>) -> Self {
        Self(stark_rings::cyclotomic_ring::CyclotomicPolyRingGeneral::<Goldilocks64Config, 1, 64>::from(value))
    }
}

impl PolyRing for GoldilocksRing64 {
    type BaseRing = Fq;

    fn coeffs(&self) -> &[Self::BaseRing] {
        self.0.coeffs()
    }
    fn coeffs_mut(&mut self) -> &mut [Self::BaseRing] {
        self.0.coeffs_mut()
    }
    fn into_coeffs(self) -> Vec<Self::BaseRing> {
        self.0.into_coeffs()
    }
    fn dimension() -> usize {
        64
    }
    fn from_scalar(scalar: Self::BaseRing) -> Self {
        Self(stark_rings::cyclotomic_ring::CyclotomicPolyRingGeneral::<Goldilocks64Config, 1, 64>::from_scalar(scalar))
    }
}

impl Flatten for GoldilocksRing64 {}

impl core::ops::Mul<Fq> for GoldilocksRing64 {
    type Output = Self;
    fn mul(self, rhs: Fq) -> Self::Output {
        let mut out = self;
        for c in out.coeffs_mut() {
            *c *= rhs;
        }
        out
    }
}

impl OverField for GoldilocksRing64 {}

impl stark_rings::Cyclotomic for GoldilocksRing64 {
    fn rot(&mut self) {
        let d = <Self as PolyRing>::dimension();
        let mut buf = -self.coeffs()[d - 1];
        for i in 0..d {
            ark_std::mem::swap(&mut buf, &mut self.coeffs_mut()[i]);
        }
    }
}

impl stark_rings::cyclotomic_ring::CRT for GoldilocksRing64 {
    type CRTForm = Self;
    fn crt(self) -> Self::CRTForm {
        self
    }
}
impl stark_rings::cyclotomic_ring::ICRT for GoldilocksRing64 {
    type ICRTForm = Self;
    fn icrt(self) -> Self::ICRTForm {
        self
    }
}

impl SuitableRing for GoldilocksRing64 {
    type CoefficientRepresentation = GoldilocksRing64;
    type PoseidonParams = GoldilocksPoseidonConfig;
}

impl<'a> core::ops::MulAssign<&'a u128> for GoldilocksRing64 {
    fn mul_assign(&mut self, rhs: &'a u128) {
        self.0 *= rhs;
    }
}

impl MulUnchecked for GoldilocksRing64 {
    type Output = Self;

    fn mul_unchecked(self, rhs: Self) -> Self::Output {
        // Keep `mul_unchecked` consistent with `Mul` (we're already in coefficient form).
        self * rhs
    }
}

impl Decompose for GoldilocksRing64 {
    fn decompose_to(&self, b: u128, out: &mut [Self]) {
        type Inner = stark_rings::cyclotomic_ring::CyclotomicPolyRingGeneral<Goldilocks64Config, 1, 64>;
        use std::cell::RefCell;
        thread_local! {
            static SCRATCH: RefCell<Vec<Inner>> = const { RefCell::new(Vec::new()) };
        }
        SCRATCH.with(|cell| {
            let mut buf = cell.borrow_mut();
            if buf.len() != out.len() {
                buf.resize(out.len(), Inner::ZERO);
            }
            self.0.decompose_to(b, &mut buf[..]);
            for (o, t) in out.iter_mut().zip(buf.iter()) {
                *o = GoldilocksRing64(t.clone());
            }
        });
    }
}

#[derive(Clone)]
pub struct Goldilocks64ChallengeSet;

impl LatticefoldChallengeSet<GoldilocksRing64> for Goldilocks64ChallengeSet {
    const BYTES_NEEDED: usize = 64;

    fn short_challenge_from_random_bytes(
        bs: &[u8],
    ) -> Result<
        <GoldilocksRing64 as SuitableRing>::CoefficientRepresentation,
        crate::challenge_set::error::ChallengeSetError,
    > {
        if bs.len() != Self::BYTES_NEEDED {
            return Err(error::ChallengeSetError::TooFewBytes(
                bs.len(),
                Self::BYTES_NEEDED,
            ));
        }

        Ok(GoldilocksRing64::from(
            bs.iter()
                .map(|&x| Fq::from(x as i16 - 128))
                .collect::<Vec<Fq>>(),
        ))
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::BigInt;
    use stark_rings::cyclotomic_ring::models::goldilocks::Fq;

    use super::*;

    #[test]
    fn test_small_challenge_from_random_bytes() {
        let challenge = GoldilocksChallengeSet::short_challenge_from_random_bytes(&[
            0x7b, 0x4b, 0xe5, 0x8e, 0xe5, 0x11, 0xd2, 0xd0, 0x9c, 0x22, 0xba, 0x2e, 0xeb, 0xa8,
            0xba, 0x35, 0xf2, 0x18,
        ])
        .unwrap();

        let res_coeffs: Vec<Fq> = vec![
            Fq::new(BigInt([27])),
            Fq::new(BigInt([13])),
            Fq::new(BigInt([18446744069414584309])),
            Fq::new(BigInt([25])),
            Fq::new(BigInt([18446744069414584303])),
            Fq::new(BigInt([18446744069414584311])),
            Fq::new(BigInt([18446744069414584319])),
            Fq::new(BigInt([18446744069414584293])),
            Fq::new(BigInt([18446744069414584307])),
            Fq::new(BigInt([18446744069414584292])),
            Fq::new(BigInt([18446744069414584302])),
            Fq::new(BigInt([7])),
            Fq::new(BigInt([2])),
            Fq::new(BigInt([8])),
            Fq::new(BigInt([11])),
            Fq::new(BigInt([18446744069414584300])),
            Fq::new(BigInt([11])),
            Fq::new(BigInt([3])),
            Fq::new(BigInt([10])),
            Fq::new(BigInt([14])),
            Fq::new(BigInt([21])),
            Fq::new(BigInt([18446744069414584297])),
            Fq::new(BigInt([18446744069414584304])),
            Fq::new(BigInt([18446744069414584295])),
        ];

        let expected = GoldilocksRingPoly::from(res_coeffs);

        assert_eq!(expected, challenge)
    }
}
