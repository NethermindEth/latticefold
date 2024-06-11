use std::{collections::HashMap, fmt::Debug};

use super::Field;

use super::crt_basis::RqCRT;

#[derive(Clone)]
pub struct RqSparse<F: Field> {
    pub coeffs: HashMap<usize, F>,
    pub prime: usize,
    pub prime_power: usize,
}

fn print_cyclotomic_element<F: Field>(z: &RqSparse<F>) -> String {
    let mut str_list: Vec<String> = vec![];
    let mut exp = 0;
    let order = z.prime.pow(z.prime_power as u32);
    while &exp != &order {
        let zero = F::zero();
        let coeff = z.coeffs.get(&exp).unwrap_or(&zero);
        if !coeff.is_zero() {
            str_list.push(String::from(
                format!("{}*E({})^{}", coeff, order, exp).as_str(),
            ));
        }
        exp += 1;
    }
    "(".to_string() + &str_list.join(" + ") + ")"
}

impl<F: Field> Debug for RqSparse<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "SparseCyclotomic ({})", print_cyclotomic_element(self))
    }
}

impl<F> RqSparse<F>
where
    F: Field,
{
    pub fn new(prime: usize, prime_power: usize, coeffs: &HashMap<usize, F>) -> RqSparse<F> {
        assert!(prime.pow(prime_power as u32 - 1) * (prime - 1) >= coeffs.len());
        RqSparse {
            coeffs: coeffs.clone(),
            prime,
            prime_power,
        }
    }

    pub fn to_crt_basis(&self, rou: F) -> RqCRT<F> {
        let prime = self.prime;
        let prime_power = self.prime_power;
        let varphi_m = self.coeffs.len();
        let m_prime = varphi_m / (prime - 1);
        // Ensure rou is the m-th root of unity
        assert_eq!(rou.pow(prime.pow(prime_power as u32) as u128), F::one());
        let omega = rou;
        let order = prime.pow(prime_power as u32 - 1) * (prime - 1);
        let mut crt_coeffs = vec![F::zero(); order];

        let mut current_omega_power = F::one();
        let mut omega_powers: Vec<F> = Vec::new();
        for _ in 0..m_prime * prime - 1 {
            omega_powers.push(current_omega_power.clone());
            current_omega_power = current_omega_power * omega;
        }
        omega_powers.push(current_omega_power);
        if prime == 2 {
            // Note that in the dense case we use NTT instead of CRT for po2
            for i in 0..varphi_m {
                let mut sum = F::zero();
                for (j, &coeff) in &self.coeffs {
                    sum = sum + omega_powers[(i * j) % varphi_m] * coeff;
                }
                crt_coeffs[i] = sum;
            }
        } else {
            let relatives_set = euler_totient(prime, prime_power);
            for i in relatives_set {
                let mut sum = F::zero();
                for (j, &coeff) in &self.coeffs {
                    sum = sum + omega_powers[(i * j) % varphi_m] * coeff;
                }
                crt_coeffs[i] = sum;
            }
        }
        RqCRT {
            crt_coeffs,
            prime,
            prime_power,
        }
    }
}

fn euler_totient(prime: usize, prime_power: usize) -> impl Iterator<Item = usize> {
    let m_prime = prime.pow(prime_power as u32);
    let relative_set =
        (1..prime).chain((1..(m_prime)).flat_map(move |i| i * prime + 1..prime * (i + 1)));
    relative_set
}
