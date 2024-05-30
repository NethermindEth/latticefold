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

    pub fn to_crt_basis(&self, rou: F) {
        let prime = self.prime;
        let prime_power = self.prime_power;
        // Ensure rou is the m-th root of unity
        assert_eq!(rou.pow(prime.pow(prime_power as u32) as u128), F::one());
        let mut coeffs = self.coeffs.clone();
        todo!()
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

        if prime != 2 {
            // Handle powers of two
        }
    }

    fn crt_prime(prime_omegas: &[F], coeffs: &mut [F]) {
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
}
