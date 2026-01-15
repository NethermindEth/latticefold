use ark_ff::Field;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use stark_rings::balanced_decomposition::Decompose;
use stark_rings::cyclotomic_ring::models::frog_ring::{Fq, FqConfig, RqNTT, RqPoly};
use stark_rings::cyclotomic_ring::Flatten;
use stark_rings::traits::FromRandomBytes;
use stark_rings::traits::MulUnchecked;
use stark_rings::{OverField, PolyRing, Ring};

use super::SuitableRing;
use crate::{
    ark_base::*,
    challenge_set::{error, LatticefoldChallengeSet},
};

/// Frog ring in the NTT form.
///
/// The base field of the NTT form is a degree-4
/// extension of the Frog field ($p=15912092521325583641$).
///
/// The NTT norm has 4 components.
pub type FrogRingNTT = RqNTT;

/// Frog ring in the coefficient form.
///
/// The cyclotomic polynomial is $X^16+1$ of degree 16.
pub type FrogRingPoly = RqPoly;

impl SuitableRing for FrogRingNTT {
    type CoefficientRepresentation = RqPoly;
    type PoseidonParams = FrogPoseidonConfig;
}

pub struct FrogPoseidonConfig;

#[derive(Clone)]
pub struct FrogChallengeSet;

/// For Frog prime the challenge set is the set of all
/// ring elements whose coefficients are in the range [-128, 128[.
impl LatticefoldChallengeSet<FrogRingNTT> for FrogChallengeSet {
    const BYTES_NEEDED: usize = 16;

    fn short_challenge_from_random_bytes(
        bs: &[u8],
    ) -> Result<
        <FrogRingNTT as SuitableRing>::CoefficientRepresentation,
        crate::challenge_set::error::ChallengeSetError,
    > {
        if bs.len() != Self::BYTES_NEEDED {
            return Err(error::ChallengeSetError::TooFewBytes(
                bs.len(),
                Self::BYTES_NEEDED,
            ));
        }

        Ok(FrogRingPoly::from(
            bs.iter()
                .map(|&x| Fq::from(x as i16 - 128))
                .collect::<Vec<Fq>>(),
        ))
    }
}

/// Frog ring wrapper with **dimension 64** (coefficient representation).
///
/// This is a minimal, correct algebraic ring type that plugs into LF+/WE as a `SuitableRing`.
/// The underlying arithmetic is performed by the existing generic cyclotomic coefficient-form
/// implementation from `stark-rings` (schoolbook multiply + reduction).
#[repr(transparent)]
#[derive(
    Copy, Clone, Debug, Default, Eq, PartialEq, Hash, CanonicalSerialize, CanonicalDeserialize,
)]
pub struct FrogRing64(pub stark_rings::cyclotomic_ring::CyclotomicPolyRingGeneral<Frog64Config, 1, 64>);

/// Parameters for \( \mathbb{F}_p[X]/(X^{64}+1) \) over Frog's base prime field.
pub struct Frog64Config;

impl stark_rings::cyclotomic_ring::CyclotomicConfig<1> for Frog64Config {
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
        impl core::ops::$trait for FrogRing64 {
            type Output = Self;
            #[inline(always)]
            fn $method(self, rhs: Self) -> Self::Output {
                Self(self.0 $op rhs.0)
            }
        }
        impl<'a> core::ops::$trait<&'a Self> for FrogRing64 {
            type Output = Self;
            #[inline(always)]
            fn $method(self, rhs: &'a Self) -> Self::Output {
                Self(self.0 $op rhs.0)
            }
        }
        impl<'a> core::ops::$trait<&'a mut Self> for FrogRing64 {
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
        impl core::ops::$trait for FrogRing64 {
            #[inline(always)]
            fn $method(&mut self, rhs: Self) {
                self.0 $op rhs.0;
            }
        }
        impl<'a> core::ops::$trait<&'a Self> for FrogRing64 {
            #[inline(always)]
            fn $method(&mut self, rhs: &'a Self) {
                self.0 $op rhs.0;
            }
        }
        impl<'a> core::ops::$trait<&'a mut Self> for FrogRing64 {
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

// -----------------------------------------------------------------------------
// FrogRing64 multiplication (performance-critical)
// -----------------------------------------------------------------------------
//
// IMPORTANT:
// The generic `CyclotomicPolyRingGeneral` multiplication path can be allocation-heavy. LF+/WE uses
// *a lot* of ring multiplies in sumcheck combiners, so `FrogRing64` needs an allocation-free
// `Mul` implementation to be remotely competitive.
//
// We implement the negacyclic convolution mod (X^64 + 1):
//   (a * b)[k] = Σ_{i+j=k} a[i]b[j]  -  Σ_{i+j=k+64} a[i]b[j]
//
// This is an O(d^2) schoolbook multiply (4096 muls) but avoids heap churn and is a major
// improvement over the generic Vec-based path.
impl core::ops::Mul for FrogRing64 {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        let a = self.coeffs();
        let b = rhs.coeffs();
        debug_assert_eq!(a.len(), 64);
        debug_assert_eq!(b.len(), 64);

        let mut out = [<Fq as Field>::ZERO; 64];
        for i in 0..64 {
            let ai = a[i];
            if ai == <Fq as Field>::ZERO {
                continue;
            }
            for j in 0..64 {
                let prod = ai * b[j];
                let k = i + j;
                if k < 64 {
                    out[k] += prod;
                } else {
                    out[k - 64] -= prod;
                }
            }
        }

        let mut r = FrogRing64::ZERO;
        r.coeffs_mut().copy_from_slice(&out);
        r
    }
}

impl<'a> core::ops::Mul<&'a Self> for FrogRing64 {
    type Output = Self;
    #[inline(always)]
    fn mul(self, rhs: &'a Self) -> Self::Output {
        self * (*rhs)
    }
}
impl<'a> core::ops::Mul<&'a mut Self> for FrogRing64 {
    type Output = Self;
    #[inline(always)]
    fn mul(self, rhs: &'a mut Self) -> Self::Output {
        self * (*rhs)
    }
}

impl core::ops::MulAssign for FrogRing64 {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl<'a> core::ops::MulAssign<&'a Self> for FrogRing64 {
    #[inline(always)]
    fn mul_assign(&mut self, rhs: &'a Self) {
        *self = *self * (*rhs);
    }
}
impl<'a> core::ops::MulAssign<&'a mut Self> for FrogRing64 {
    #[inline(always)]
    fn mul_assign(&mut self, rhs: &'a mut Self) {
        *self = *self * (*rhs);
    }
}

impl core::ops::Neg for FrogRing64 {
    type Output = Self;
    #[inline(always)]
    fn neg(self) -> Self::Output {
        Self(-self.0)
    }
}

impl ark_std::fmt::Display for FrogRing64 {
    fn fmt(&self, f: &mut ark_std::fmt::Formatter<'_>) -> ark_std::fmt::Result {
        ark_std::fmt::Display::fmt(&self.0, f)
    }
}

impl ark_std::iter::Sum for FrogRing64 {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        Self(iter.map(|x| x.0).sum())
    }
}
impl<'a> ark_std::iter::Sum<&'a Self> for FrogRing64 {
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        Self(iter.map(|x| x.0).sum())
    }
}
impl ark_std::iter::Product for FrogRing64 {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        Self(iter.map(|x| x.0).product())
    }
}
impl<'a> ark_std::iter::Product<&'a Self> for FrogRing64 {
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        Self(iter.map(|x| x.0).product())
    }
}

impl ark_std::Zero for FrogRing64 {
    #[inline(always)]
    fn zero() -> Self {
        Self(<stark_rings::cyclotomic_ring::CyclotomicPolyRingGeneral<Frog64Config, 1, 64> as Ring>::ZERO)
    }
    #[inline(always)]
    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

impl ark_std::One for FrogRing64 {
    #[inline(always)]
    fn one() -> Self {
        Self(<stark_rings::cyclotomic_ring::CyclotomicPolyRingGeneral<Frog64Config, 1, 64> as Ring>::ONE)
    }
}

impl ark_std::UniformRand for FrogRing64 {
    fn rand<R: rand::Rng + ?Sized>(rng: &mut R) -> Self {
        Self(stark_rings::cyclotomic_ring::CyclotomicPolyRingGeneral::<Frog64Config, 1, 64>::rand(rng))
    }
}

impl FromRandomBytes<Self> for FrogRing64 {
    fn byte_size() -> usize {
        stark_rings::cyclotomic_ring::CyclotomicPolyRingGeneral::<Frog64Config, 1, 64>::byte_size()
    }
    fn try_from_random_bytes(bytes: &[u8]) -> Option<Self> {
        stark_rings::cyclotomic_ring::CyclotomicPolyRingGeneral::<Frog64Config, 1, 64>::try_from_random_bytes(bytes)
            .map(Self)
    }
}

impl Ring for FrogRing64 {
    const ZERO: Self = Self(<stark_rings::cyclotomic_ring::CyclotomicPolyRingGeneral<Frog64Config, 1, 64> as Ring>::ZERO);
    const ONE: Self = Self(<stark_rings::cyclotomic_ring::CyclotomicPolyRingGeneral<Frog64Config, 1, 64> as Ring>::ONE);
}

impl From<u128> for FrogRing64 {
    fn from(value: u128) -> Self {
        Self(stark_rings::cyclotomic_ring::CyclotomicPolyRingGeneral::<Frog64Config, 1, 64>::from(value))
    }
}
impl From<u64> for FrogRing64 {
    fn from(value: u64) -> Self {
        Self::from(value as u128)
    }
}
impl From<u32> for FrogRing64 {
    fn from(value: u32) -> Self {
        Self::from(value as u128)
    }
}
impl From<u16> for FrogRing64 {
    fn from(value: u16) -> Self {
        Self::from(value as u128)
    }
}
impl From<u8> for FrogRing64 {
    fn from(value: u8) -> Self {
        Self::from(value as u128)
    }
}
impl From<bool> for FrogRing64 {
    fn from(value: bool) -> Self {
        Self::from(value as u128)
    }
}

impl From<Fq> for FrogRing64 {
    fn from(value: Fq) -> Self {
        Self(stark_rings::cyclotomic_ring::CyclotomicPolyRingGeneral::<Frog64Config, 1, 64>::from(value))
    }
}

impl From<Vec<Fq>> for FrogRing64 {
    fn from(value: Vec<Fq>) -> Self {
        Self(stark_rings::cyclotomic_ring::CyclotomicPolyRingGeneral::<Frog64Config, 1, 64>::from(value))
    }
}

impl PolyRing for FrogRing64 {
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
        Self(stark_rings::cyclotomic_ring::CyclotomicPolyRingGeneral::<Frog64Config, 1, 64>::from_scalar(scalar))
    }
}

impl Flatten for FrogRing64 {}

impl core::ops::Mul<Fq> for FrogRing64 {
    type Output = Self;
    fn mul(self, rhs: Fq) -> Self::Output {
        let mut out = self;
        for c in out.coeffs_mut() {
            *c *= rhs;
        }
        out
    }
}

impl OverField for FrogRing64 {}

impl stark_rings::Cyclotomic for FrogRing64 {
    fn rot(&mut self) {
        let d = <Self as PolyRing>::dimension();
        let mut buf = -self.coeffs()[d - 1];
        for i in 0..d {
            ark_std::mem::swap(&mut buf, &mut self.coeffs_mut()[i]);
        }
    }
}

impl stark_rings::cyclotomic_ring::CRT for FrogRing64 {
    type CRTForm = Self;
    fn crt(self) -> Self::CRTForm {
        self
    }
}
impl stark_rings::cyclotomic_ring::ICRT for FrogRing64 {
    type ICRTForm = Self;
    fn icrt(self) -> Self::ICRTForm {
        self
    }
}

impl SuitableRing for FrogRing64 {
    type CoefficientRepresentation = FrogRing64;
    type PoseidonParams = FrogPoseidonConfig;
}

impl<'a> core::ops::MulAssign<&'a u128> for FrogRing64 {
    fn mul_assign(&mut self, rhs: &'a u128) {
        self.0 *= rhs;
    }
}

impl MulUnchecked for FrogRing64 {
    type Output = Self;

    fn mul_unchecked(self, rhs: Self) -> Self::Output {
        // Keep `mul_unchecked` consistent with `Mul` (we're already in coefficient form).
        self * rhs
    }
}

impl Decompose for FrogRing64 {
    fn decompose_to(&self, b: u128, out: &mut [Self]) {
        // Delegate to the underlying cyclotomic coefficient-form implementation.
        let mut tmp = vec![
            stark_rings::cyclotomic_ring::CyclotomicPolyRingGeneral::<Frog64Config, 1, 64>::ZERO;
            out.len()
        ];
        self.0.decompose_to(b, &mut tmp);
        for (o, t) in out.iter_mut().zip(tmp.into_iter()) {
            *o = Self(t);
        }
    }
}

#[derive(Clone)]
pub struct Frog64ChallengeSet;

impl LatticefoldChallengeSet<FrogRing64> for Frog64ChallengeSet {
    const BYTES_NEEDED: usize = 64;

    fn short_challenge_from_random_bytes(
        bs: &[u8],
    ) -> Result<
        <FrogRing64 as SuitableRing>::CoefficientRepresentation,
        crate::challenge_set::error::ChallengeSetError,
    > {
        if bs.len() != Self::BYTES_NEEDED {
            return Err(error::ChallengeSetError::TooFewBytes(
                bs.len(),
                Self::BYTES_NEEDED,
            ));
        }

        Ok(FrogRing64::from(
            bs.iter()
                .map(|&x| Fq::from(x as i16 - 128))
                .collect::<Vec<Fq>>(),
        ))
    }
}

#[cfg(test)]
mod frog64_tests {
    use super::*;
    use ark_std::{test_rng, UniformRand, Zero};

    type Inner = stark_rings::cyclotomic_ring::CyclotomicPolyRingGeneral<Frog64Config, 1, 64>;

    #[test]
    fn test_frog64_mul_matches_inner_generic() {
        let mut rng = test_rng();
        for _ in 0..256 {
            let a = FrogRing64::rand(&mut rng);
            let b = FrogRing64::rand(&mut rng);

            // Reference: underlying generic ring multiplication.
            let ref_out = Inner::from(a.into_coeffs()) * Inner::from(b.into_coeffs());
            let got = a * b;

            assert_eq!(got.coeffs(), ref_out.coeffs());
        }
    }

    #[test]
    fn test_frog64_mul_by_scalar_in_place() {
        let mut rng = test_rng();
        let a = FrogRing64::rand(&mut rng);
        let s = Fq::rand(&mut rng);
        let got = a * s;

        // Reference via inner generic multiplication by scalar (should match).
        let ref_out = Inner::from(a.into_coeffs()) * s;
        assert_eq!(got.coeffs(), ref_out.coeffs());
    }

    /// Rough timing smoke test (ignored by default).
    ///
    /// This is not a benchmark harness; it’s just a sanity check that `FrogRing64` isn't
    /// catastrophically slower than `FrogRingPoly` for multiplication.
    #[test]
    #[ignore]
    fn test_frog64_mul_timing_smoke() {
        use std::time::Instant;

        let mut rng = test_rng();
        let iters: usize = 50_000;

        let mut acc64 = FrogRing64::ONE;
        let a64 = FrogRing64::rand(&mut rng);
        let b64 = FrogRing64::rand(&mut rng);
        let t64 = Instant::now();
        for _ in 0..iters {
            acc64 *= a64;
            acc64 *= b64;
        }
        let dt64 = t64.elapsed();

        let mut acc16 = FrogRingPoly::ONE;
        let a16 = FrogRingPoly::rand(&mut rng);
        let b16 = FrogRingPoly::rand(&mut rng);
        let t16 = Instant::now();
        for _ in 0..iters {
            acc16 *= a16;
            acc16 *= b16;
        }
        let dt16 = t16.elapsed();

        // Keep the values live.
        assert!(!acc64.is_zero() || !acc16.is_zero());
        eprintln!(
            "[frog64_smoke] iters={} mul64={:?} mul16={:?}",
            iters, dt64, dt16
        );
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::BigInt;
    use stark_rings::cyclotomic_ring::models::frog_ring::Fq;

    use super::*;

    #[test]
    fn test_small_challenge_from_random_bytes() {
        let challenge = FrogChallengeSet::short_challenge_from_random_bytes(&[
            0x7b, 0x4b, 0xe5, 0x8e, 0xe5, 0x11, 0xd2, 0xd0, 0x9c, 0x22, 0xba, 0x2e, 0xeb, 0xa8,
            0xba, 0x35,
        ])
        .unwrap();

        let res_coeffs: Vec<Fq> = vec![
            Fq::new(BigInt([15912092521325583636])),
            Fq::new(BigInt([15912092521325583588])),
            Fq::new(BigInt([101])),
            Fq::new(BigInt([14])),
            Fq::new(BigInt([101])),
            Fq::new(BigInt([15912092521325583530])),
            Fq::new(BigInt([82])),
            Fq::new(BigInt([80])),
            Fq::new(BigInt([28])),
            Fq::new(BigInt([15912092521325583547])),
            Fq::new(BigInt([58])),
            Fq::new(BigInt([15912092521325583559])),
            Fq::new(BigInt([107])),
            Fq::new(BigInt([40])),
            Fq::new(BigInt([58])),
            Fq::new(BigInt([15912092521325583566])),
        ];

        let expected = FrogRingPoly::from(res_coeffs);

        assert_eq!(expected, challenge)
    }

    #[test]
    fn test_small_challenge_from_random_bytes_frog64() {
        let mut bs = [0u8; 64];
        for (i, b) in bs.iter_mut().enumerate() {
            *b = (i as u8).wrapping_mul(13);
        }
        let challenge = Frog64ChallengeSet::short_challenge_from_random_bytes(&bs).unwrap();
        assert_eq!(challenge.coeffs().len(), 64);
    }
}
