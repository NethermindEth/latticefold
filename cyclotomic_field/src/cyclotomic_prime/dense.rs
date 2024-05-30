use super::Field;

/// Ring m-th cyclotomic element where m is a power of a prime
/// The minimal polynomial for a prime is define as Phi_p(X) = Sum_{i=0}^{p-1} X^i
/// then the minimal polynomial for a prime is define as Phi_m(X) = Phi_p(x^{m/p})
pub struct RqDense<F: Field> {
    pub coeffs: Vec<F>,
    prime: usize,
    prime_power: usize,
}

impl<F> RqDense<F>
where
    F: Field,
{
    pub fn new(prime: usize, prime_power: usize, coeffs: Vec<F>) -> Self {
        // Assure is a prime power cyclotomic
        assert_eq!(
            prime.pow(prime_power as u32 - 1) * (prime - 1),
            coeffs.len()
        );

        RqDense {
            coeffs,
            prime,
            prime_power,
        }
    }

    pub fn to_crt_basis(&self, prime: usize, prime_power: usize, rou: F) {
        todo!()
    }
}
