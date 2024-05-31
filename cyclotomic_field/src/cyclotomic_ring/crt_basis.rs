use std::ops::{Add, Mul};

use super::Field;

/// A ring elelemnt can be represented in the CRT basis where the basis consist of the elements
/// (X - w_m^i) where m is define by the m-th root of unity use as an abstract element to extend
/// the field
/// TODO: Cite paper
pub struct RqCRT<F: Field> {
    pub crt_coeffs: Vec<F>,
    pub prime: usize,
    pub prime_power: usize,
}

impl<F> RqCRT<F>
where
    F: Field,
{
    pub fn new(prime: usize, prime_power: usize, crt_coeffs: Vec<F>) -> Self {
        // Assure is a prime power cyclotomic
        assert_eq!(
            prime.pow(prime_power as u32 - 1) * (prime - 1),
            crt_coeffs.len()
        );

        RqCRT {
            crt_coeffs,
            prime,
            prime_power,
        }
    }
}

impl<F> Add for RqCRT<F>
where
    F: Field,
{
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        let coeffs_1 = self.crt_coeffs;
        let coeffs_2 = rhs.crt_coeffs;

        let crt_coeffs = coeffs_1
            .iter()
            .zip(coeffs_2.iter())
            .map(|(&c1, &c2)| c1 + c2)
            .collect::<Vec<_>>();

        RqCRT {
            crt_coeffs,
            prime: self.prime,
            prime_power: self.prime_power,
        }
    }
}
impl<F> Mul for RqCRT<F>
where
    F: Field,
{
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        let coeffs_1 = self.crt_coeffs;
        let coeffs_2 = rhs.crt_coeffs;

        let crt_coeffs = coeffs_1
            .iter()
            .zip(coeffs_2.iter())
            .map(|(&c1, &c2)| c1 * c2)
            .collect::<Vec<_>>();

        RqCRT {
            crt_coeffs,
            prime: self.prime,
            prime_power: self.prime_power,
        }
    }
}
