use crate::ark_base::*;
use num_bigint::BigUint;
use num_traits::ToPrimitive;

#[cfg(feature = "std")]
use f64 as math;
#[cfg(not(feature = "std"))]
use libm as math;

#[cfg(feature = "std")]
fn log(val: f64) -> f64 {
    val.ln()
}
#[cfg(not(feature = "std"))]
fn log(val: f64) -> f64 {
    math::log(val)
}

#[cfg(feature = "std")]
fn pow(val: f64, n: f64) -> f64 {
    val.powf(n)
}
#[cfg(not(feature = "std"))]
fn pow(val: f64, n: f64) -> f64 {
    math::pow(val, n)
}

fn calculate_bound_l2(degree: usize, kappa: usize, ring_modulus_log2: f64) -> BigUint {
    // The current security parameter use log2(delta)
    let delta = 1.0045_f64;
    // Calculate B_{L_2} as 2^{2 \sqrt{\text{log2}(\delta) \times \text{degree} \times \kappa \times \frac{\text{modulus}}{2}}}
    let bound_l2 = pow(
        2f64,
        2.0 * math::sqrt(log(delta) / log(2f64))
            * math::sqrt(degree as f64 * kappa as f64 * ring_modulus_log2),
    );
    let bound_l2_ceil = math::ceil(bound_l2) as u64; // Ceil and convert to u64
    BigUint::from(bound_l2_ceil) // Convert to BigUint
}

pub fn check_ring_modulus_128_bits_security(
    ring_modulus: &BigUint,
    kappa: usize,
    degree: usize,
    num_cols: usize,
    b: u128,
    l: usize,
    already_under_bound: bool,
) -> bool {
    // Modulus bits and half
    let (ring_modulus_log2, ring_modulus_half) = (ring_modulus.bits() as f64, ring_modulus / 2u32);

    // Calculate the left side of the inequality
    let bound_l2_bigint = calculate_bound_l2(degree, kappa, ring_modulus_log2);
    let bound_l2_check = bound_l2_bigint < ring_modulus_half;
    // Calculate bound_inf B_inf as B_{L_2} / \sqrt{\text{degree} \times \text{num_cols}}
    let bound_inf =
        bound_l2_bigint.to_f64().unwrap() / (math::sqrt(degree as f64 * num_cols as f64));

    let b_check = b.to_f64().unwrap() < bound_inf;
    // Check if we need to decompose and b^l > stark_modulus/2
    let b_pow_l_check = if already_under_bound && l == 1 {
        true
    } else {
        BigUint::from(b).pow(l as u32) > ring_modulus_half
    };

    // Return the result of the condition
    bound_l2_check && b_check && b_pow_l_check
}
