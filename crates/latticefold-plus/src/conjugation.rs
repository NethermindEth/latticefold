//! Complex conjugation automorphism for power-of-two cyclotomic rings.
//!
//! For `R_q = Z_q[X]/(X^N + 1)` with `N` a power of two, the automorphism
//! `sigma_{-1}: X -> X^{-1} = -X^{N-1}` is the building block of the exact
//! shortening step of the ell-2 norm-check (ePrint 2026/721, section 4.4).
//!
//! The load-bearing identity is `||w||_2^2 = ct(<w, conj(w)>)`, where `ct(.)`
//! is the constant term: for `a = sum_i a_i X^i`, `ct(a * conj(a)) = sum_i a_i^2`.
//!
//! SOUNDNESS: this identity holds only for power-of-two `N` (rings `X^N + 1`).
//! On general cyclotomics (Goldilocks `X^24 - X^12 + 1`, BabyBear
//! `X^72 - X^36 + 1`) `X^{-1}` is not `-X^{N-1}` and the identity is false. The
//! general-cyclotomic path (trace / trace-dual, per KLNO25) is out of scope, see
//! caveat 3 of 2026/721. We guard with a hard `assert!` so misuse fails loudly.

use ark_std::Zero;
use stark_rings::{CoeffRing, PolyRing, Zq};

/// The complex conjugation automorphism `sigma_{-1}: X -> -X^{N-1}` in the
/// coefficient representation.
///
/// For `a = sum_{i=0}^{N-1} a_i X^i` in `Z_q[X]/(X^N + 1)`, the result has
/// coefficient `a_0` at index 0 and `-a_{N-j}` at index `j` for `1 <= j < N`.
/// This is `a(X^{-1})` reduced modulo `X^N + 1`.
///
/// # Panics
/// Panics if `R::dimension()` is not a power of two (see module SOUNDNESS note).
pub fn conjugate<R: PolyRing>(a: &R) -> R {
    let n = R::dimension();
    assert!(
        n.is_power_of_two(),
        "conjugate is only defined for power-of-two rings (X^N + 1), got N = {n}"
    );

    let c = a.coeffs();
    let mut out = vec![R::BaseRing::zero(); n];
    out[0] = c[0];
    for j in 1..n {
        out[j] = -c[n - j];
    }
    R::from(out)
}

/// Squared ell-2 norm of a witness vector over the coefficient embedding:
/// `||w||_2^2 = sum_i ct(w_i * conj(w_i)) = sum_i sum_j (w_i)_j^2`.
///
/// Returned as the base ring element `ct(<w, conj(w)>)`. For the exact-shortening
/// equality to hold over the integers (not just mod q), the caller must ensure
/// `||w||_2^2 < q` (soundness gate 6 of the plan).
///
/// # Panics
/// Panics if `R::dimension()` is not a power of two.
pub fn l2_norm_sq<R: CoeffRing>(w: &[R]) -> R::BaseRing
where
    R::BaseRing: Zq,
{
    w.iter().fold(R::BaseRing::zero(), |acc, w_i| {
        acc + (*w_i * conjugate(w_i)).ct()
    })
}

#[cfg(test)]
mod tests {
    use ark_ff::UniformRand;
    use ark_std::{One, Zero};
    use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as R;

    use super::*;

    fn rand_ring(rng: &mut impl rand::Rng) -> R {
        R::rand(rng)
    }

    // (i) The load-bearing identity: ct(a * conj(a)) == sum_i a_i^2.
    #[test]
    fn test_ct_a_conj_a_is_sum_of_squares() {
        let mut rng = ark_std::test_rng();
        for _ in 0..32 {
            let a = rand_ring(&mut rng);
            let lhs = (a * conjugate(&a)).ct();
            let rhs = a
                .coeffs()
                .iter()
                .fold(<R as PolyRing>::BaseRing::zero(), |acc, &c| acc + c * c);
            assert_eq!(lhs, rhs);
        }
    }

    // (ii) conj is additive and involutive.
    #[test]
    fn test_conj_additive_and_involutive() {
        let mut rng = ark_std::test_rng();
        for _ in 0..32 {
            let a = rand_ring(&mut rng);
            let b = rand_ring(&mut rng);
            assert_eq!(conjugate(&(a + b)), conjugate(&a) + conjugate(&b));
            assert_eq!(conjugate(&conjugate(&a)), a);
        }
    }

    // (iii) l2_norm_sq on a hand-computed small vector.
    #[test]
    fn test_l2_norm_sq_known_vector() {
        // w = [ 2 + 5X , 4 + X^2 ] ; ||w||_2^2 = 4 + 25 + 16 + 1 = 46
        let mut w0 = R::zero();
        w0.coeffs_mut()[0] = 2u128.into();
        w0.coeffs_mut()[1] = 5u128.into();
        let mut w1 = R::zero();
        w1.coeffs_mut()[0] = 4u128.into();
        w1.coeffs_mut()[2] = 1u128.into();

        let got = l2_norm_sq(&[w0, w1]);
        assert_eq!(got, <R as PolyRing>::BaseRing::from(46u128));
    }

    // (iv) A non-power-of-two ring (Goldilocks, X^24 - X^12 + 1) must be rejected.
    #[test]
    #[should_panic(expected = "power-of-two")]
    fn test_non_power_of_two_rejected() {
        use stark_rings::cyclotomic_ring::models::goldilocks::RqPoly as G;
        let mut rng = ark_std::test_rng();
        let a = G::rand(&mut rng);
        let _ = conjugate(&a);
    }

    // conj fixes the constant polynomial and negates the others as expected.
    #[test]
    fn test_conj_basis() {
        // conj(1) = 1
        assert_eq!(conjugate(&R::one()), R::one());
        // conj(X) = -X^{N-1}
        let mut x = R::zero();
        x.coeffs_mut()[1] = 1u128.into();
        let mut expected = R::zero();
        expected.coeffs_mut()[R::dimension() - 1] = -<R as PolyRing>::BaseRing::one();
        assert_eq!(conjugate(&x), expected);
    }
}
