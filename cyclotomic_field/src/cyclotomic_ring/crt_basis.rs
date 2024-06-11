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
        let prime = self.prime;
        let prime_power = self.prime_power;
        assert_eq!(rou.pow(prime.pow(prime_power as u32) as u128), F::one());
        let mut crt_coeffs = self.crt_coeffs.clone();
        if prime == 2 {
            RqCRT::intt(prime, prime_power - 1, rou, crt_coeffs.as_mut_slice());
        } else {
            RqCRT::icrt(prime, prime_power, rou, crt_coeffs.as_mut_slice());
        }

        RqDense {
            coeffs: crt_coeffs,
            prime,
            prime_power,
        }
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
        RqCRT::radixp_intt(prime, prime_power, omega_powers.as_slice(), crt_coeffs);
    }

    fn icrt(prime: usize, prime_power: usize, omega: F, crt_coeffs: &mut [F]) {
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
        let prime_omegas = omega_powers
            .iter()
            .step_by(m_prime)
            .map(|w| w.clone()) // Remove clone
            .collect::<Vec<_>>();

        if prime_power == 1 {
            RqCRT::<F>::icrt_prime(prime_omegas.as_slice(), crt_coeffs);
            return;
        }

        RqCRT::inverse_stride_permutation(prime - 1, crt_coeffs);

        let m_prime_omega_powers = omega_powers
            .iter()
            .step_by(prime)
            .map(|w| w.clone())
            .collect::<Vec<_>>();

        for coeffs_chunk in crt_coeffs.chunks_exact_mut(m_prime) {
            RqCRT::radixp_intt(prime, prime_power - 1, &m_prime_omega_powers, coeffs_chunk);
        }

        let t_hat = RqCRT::twiddle_hat_factors(prime, omega_powers.as_slice());
        for (twiddle_factor, coeff) in t_hat.iter().zip(crt_coeffs.iter_mut()) {
            *coeff = *coeff * twiddle_factor.inverse();
        }

        RqCRT::stride_permutation(prime - 1, crt_coeffs);

        for crt_coeffs_chunks in crt_coeffs.chunks_exact_mut(prime - 1) {
            RqCRT::icrt_prime(prime_omegas.as_slice(), crt_coeffs_chunks);
        }

        RqCRT::inverse_stride_permutation(prime - 1, crt_coeffs);
    }

    fn icrt_prime(prime_omegas: &[F], crt_coeffs: &mut [F]) {
        let prime = prime_omegas.len();
        // Assert if prime?
        match prime {
            2 => return,
            3 => {
                RqCRT::icrt_three(prime_omegas, crt_coeffs);
            }
            _ => {
                unimplemented!("No support for powers-of-{}", prime);
            }
        }
        // TODO: Easy way to define the inverse of the matrix?
        todo!()
    }

    fn radixp_intt(prime: usize, prime_power: usize, omega_powers: &[F], crt_coeffs: &mut [F]) {
        // Get the inverse powers
        let mut omega_powers = omega_powers.to_vec();
        omega_powers.reverse();
        omega_powers.rotate_right(1);

        RqDense::radixp_ntt(prime, prime_power, &omega_powers, crt_coeffs);
        let mut inv_n = F::zero();
        let one = F::one();
        for _ in 0..crt_coeffs.len() {
            inv_n = inv_n + one;
        }
        inv_n = inv_n.inverse();
        for coeff in crt_coeffs {
            *coeff = *coeff * inv_n;
        }
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

    fn icrt_three(prime_omegas: &[F], crt_coeffs: &mut [F]) {
        let b_0 = crt_coeffs[0] * prime_omegas[2] - crt_coeffs[1] * prime_omegas[1];
        let b_1 = crt_coeffs[1] - crt_coeffs[0];
        let det = prime_omegas[2] - prime_omegas[1];
        let inv_det = det.inverse();
        crt_coeffs[0] = b_0 * inv_det;
        crt_coeffs[1] = b_1 * inv_det;
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
