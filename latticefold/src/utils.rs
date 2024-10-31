pub mod mle;
pub mod sumcheck;
pub mod uv_polynomial;

use num_bigint::BigUint;
use num_traits::ToPrimitive;
use std::f64;
// Position this in a better place ?
pub fn check_ring_modulus_128_bits_security(
    ring_modulus: &BigUint,
    kappa: usize,
    degree: usize,
    num_cols: usize,
    b: u128,
    l: usize,
) -> bool {
    // Calculate the logarithm of stark_modulus
    let ring_modulus_log2 = ring_modulus.bits() as f64;
    let ring_modulus_half = ring_modulus / 2u32;

    // Calculate the left side of the inequality
    let bound_l2 = 2f64.powf(
        2.0 * (1.0045f64.ln() / 2f64.ln()).sqrt()
            * (degree as f64 * kappa as f64 * ring_modulus_log2).sqrt(),
    );
    let bound_l2_ceil = bound_l2.ceil() as u64; // Ceil and convert to u64
    let bound_l2_bigint = BigUint::from(bound_l2_ceil); // Convert to BigUint
    let bound_l2_check = bound_l2_bigint < ring_modulus_half;
    // Calculate bound_inf
    let bound_inf = bound_l2 / ((degree as f64 * num_cols as f64).sqrt());

    let b_check = b.to_f64().unwrap() < bound_inf;
    // Calculate the right side of the inequality
    // Check if b^l > stark_modulus/2
    let b_bigint = BigUint::from(b);
    let b_pow_l = b_bigint.pow(l as u32);
    let b_pow_l_check = b_pow_l > ring_modulus_half;

    // Return the result of the condition
    bound_l2_check && b_check && b_pow_l_check
}
