mod crt_basis;
mod dense;
mod sparse;

use fields::Field;

pub struct CRTDomain<F: Field> {
    omega_powers: Vec<F>,
}

impl<F: Field> CRTDomain<F> {
    pub fn new(prime: usize, prime_power: usize, rou: F, q: usize) -> Self {
        let n = prime.pow(prime_power as u32);
        assert_eq!(
            (q - 1) / n % prime,
            0,
            "Field does not support powers of {}",
            prime
        );
        let resize_power = (q - 1) / n;
        let omega = rou.pow(resize_power as u128);
        assert_eq!(omega.pow(n as u128), F::one(), "Omega^N != 1");
        let mut current_omega_power = F::one();
        let mut omega_powers = Vec::new();
        for _ in 0..n {
            omega_powers.push(current_omega_power);
            current_omega_power = current_omega_power * omega;
        }

        CRTDomain { omega_powers }
    }
}
