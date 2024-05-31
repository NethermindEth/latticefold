use std::ops::{Add, Mul};

use super::{dense::RqDense, Field};

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

    pub fn to_power_basis(&self, rou: F) -> RqDense<F> {
        todo!()
    }

    pub fn intt(prime: usize, prime_power: usize, omega: F, crt_coeffs: &mut [F]) {
        let varphi_m = crt_coeffs.len();
        let m_prime = varphi_m / (prime - 1);

        // Set up omegas powers
        let mut current_omega_power = F::one();
        let mut omega_powers: Vec<F> = Vec::new();
        for _ in 0..m_prime * prime - 1 {
            omega_powers.push(current_omega_power.clone());
            current_omega_power = current_omega_power * omega;
        }
        omega_powers.push(current_omega_power);
        todo!()
    }

    fn icrt(prime: usize, prime_power: usize, omega: F, crt_coeffs: &mut [F]) {
        todo!()
    }

    fn icrt_prime(prime_omegas: &[F], crt_coeffs: &mut [F]) {
        todo!()
    }

    fn radixp_intt(prime: usize, prime_power: usize, omega_powers: &[F], crt_coeffs: &mut [F]) {
        todo!()
    }

    fn stride_permutation(varphi_p: usize, input: &mut [F]) {
        let varphi_m = input.len();
        assert_eq!(varphi_m % varphi_p, 0); // Assure permutation is well define

        let mut temp = vec![F::zero(); varphi_m];

        for i in 0..varphi_m - 1 {
            let new_index = (i * varphi_p) % (varphi_m - 1);
            temp[new_index] = input[i].clone(); // Remove clone, try using swap to do it in-place
        }
        temp[varphi_m - 1] = input[varphi_m - 1].clone();
        input.clone_from_slice(&temp);
    }

    fn inverse_stride_permutation(varphi_p: usize, input: &mut [F]) {
        let varphi_m = input.len();
        assert_eq!(varphi_m % varphi_p, 0); // Assure permutation is well define

        let mut temp = vec![F::zero(); varphi_m];

        let d = varphi_m / varphi_p;
        for i in 0..varphi_m - 1 {
            let new_index = (i * d) % (varphi_m - 1);
            temp[new_index] = input[i].clone(); // Remove clone, try using swap to do it in-place
        }
        temp[varphi_m - 1] = input[varphi_m - 1].clone();
        input.clone_from_slice(&temp);
    }

    fn twiddle_hat_factors(prime: usize, omega_powers: &[F]) -> Vec<F> {
        let m_prime = omega_powers.len() / prime;
        let mut t_hat = Vec::new();
        for r in 1..prime {
            for j in 0..m_prime {
                t_hat.push(omega_powers[r * j].clone());
            }
        }
        t_hat
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
