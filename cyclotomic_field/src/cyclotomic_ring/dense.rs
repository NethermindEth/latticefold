use super::crt_basis::RqCRT;
use super::Field;

/// Ring m-th cyclotomic element where m is a power of a prime
/// The minimal polynomial for a prime is define as Phi_p(X) = Sum_{i=0}^{p-1} X^i
/// then the minimal polynomial for a prime is define as Phi_m(X) = Phi_p(x^{m/p})
pub struct RqDense<F: Field> {
    pub coeffs: Vec<F>,
    pub prime: usize,
    pub prime_power: usize,
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

    pub fn to_crt_basis(&self, rou: F) -> RqCRT<F> {
        let prime = self.prime;
        let prime_power = self.prime_power;
        // Ensure rou is the m-th root of unity
        assert_eq!(rou.pow(prime.pow(prime_power as u32) as u128), F::one());
        let mut coeffs = self.coeffs.clone();
        if prime == 2 {
            // Note that in the prime = 2 case the twiddle factors for CRT are just a scalar
            // multiplication of each entry and because all operations are element wise it can be
            // rescale at a later time or not rescale and perform ICRT without affecting
            // Basically CRT = NTT for po2
            // Same thing happens with reordering/permutations
            RqDense::<F>::ntt(prime, prime_power - 1, rou, coeffs.as_mut_slice());
        } else {
            RqDense::crt(prime, prime_power, rou, coeffs.as_mut_slice())
        }
        RqCRT {
            crt_coeffs: coeffs,
            prime,
            prime_power,
        }
    }

    pub fn ntt(prime: usize, prime_power: usize, omega: F, coeffs: &mut [F]) {
        let varphi_m = coeffs.len();
        let m_prime = varphi_m / (prime - 1);

        // Set up omegas powers
        let mut current_omega_power = F::one();
        let mut omega_powers: Vec<F> = Vec::new();
        for _ in 0..m_prime * prime - 1 {
            omega_powers.push(current_omega_power.clone());
            current_omega_power = current_omega_power * omega;
        }
        omega_powers.push(current_omega_power);
        // TODO: optimize radix2-ntt
        RqDense::radixp_ntt(prime, prime_power, omega_powers.as_slice(), coeffs);
    }

    // Assure omega is the m-th root of unity
    fn crt(prime: usize, prime_power: usize, omega: F, coeffs: &mut [F]) {
        let varphi_m = coeffs.len();
        let m_prime = varphi_m / (prime - 1);

        // Set up omegas powers
        let mut current_omega_power = F::one();
        let mut omega_powers: Vec<F> = Vec::new();
        // Do we need all the powers of omega or just some? review
        // twiddle factors for CRT
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
            RqDense::<F>::crt_prime(prime_omegas.as_slice(), coeffs);
            return;
        }

        // Handle powers diferent than two
        RqDense::stride_permutation(prime - 1, coeffs);

        // In the prime = 2 case the crt matrix is [[1]]
        for coeffs_chunk in coeffs.chunks_exact_mut(prime - 1) {
            RqDense::crt_prime(prime_omegas.as_slice(), coeffs_chunk);
        }

        RqDense::inverse_stride_permutation(prime - 1, coeffs);

        let t_hat = RqDense::twiddle_hat_factors(prime, omega_powers.as_slice());
        for (&twiddle_factor, coeff) in t_hat.iter().zip(coeffs.iter_mut()) {
            *coeff = *coeff * twiddle_factor;
        }

        let m_prime_omega_powers = omega_powers
            .iter()
            .step_by(prime)
            .map(|w| w.clone())
            .collect::<Vec<_>>();

        for coeffs_chunk in coeffs.chunks_exact_mut(m_prime) {
            RqDense::radixp_ntt(prime, prime_power - 1, &m_prime_omega_powers, coeffs_chunk);
        }

        RqDense::stride_permutation(prime - 1, coeffs);
    }

    fn crt_prime(prime_omegas: &[F], coeffs: &mut [F]) {
        assert_eq!(prime_omegas.len(), coeffs.len());
        let p = coeffs.len();
        let mut changed_coeffs = Vec::new();
        for i in 0..p {
            let mut sum = F::zero();
            for j in 0..p {
                sum = sum + prime_omegas[((i + 1) * j) % p] * coeffs[j];
            }
            changed_coeffs.push(sum);
        }
        coeffs.copy_from_slice(&changed_coeffs);
    }

    pub fn radixp_ntt(prime: usize, prime_power: usize, omega_powers: &[F], coeffs: &mut [F]) {
        let n = coeffs.len();
        if n == 1 {
            return;
        }
        let mut decomposed_coeffs = vec![vec![F::zero(); n / prime]; prime];
        for i in 0..n / prime {
            for j in 0..decomposed_coeffs.len() {
                decomposed_coeffs[j][i] = coeffs[prime * i + j].clone();
            }
        }

        let primed_omegas = omega_powers
            .iter()
            .step_by(prime)
            .map(|w| w.clone())
            .collect::<Vec<_>>();

        for i in 0..prime {
            RqDense::radixp_ntt(
                prime,
                prime_power - 1,
                primed_omegas.as_slice(),
                decomposed_coeffs[i].as_mut_slice(),
            );
        }

        for q in 0..n / prime {
            for s in 0..prime {
                let mut sum = F::zero();
                for l in 0..prime {
                    sum = sum
                        + omega_powers[(l * (q + s * (n / prime))) % n] * decomposed_coeffs[l][q];
                }
                coeffs[q + s * (n / prime)] = sum;
            }
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
}
