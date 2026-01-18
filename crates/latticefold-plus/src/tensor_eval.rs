//! Optimized tensor product MLE evaluation.
//!
//! This module provides O(sum of factor sizes) evaluation for tensor product MLEs,
//! instead of O(product of factor sizes) for dense evaluation.
//!
//! # Key Insight
//!
//! For a tensor product of MLEs:
//! ```text
//! MLE_{a⊗b}(r_x, r_y) = MLE_a(r_x) * MLE_b(r_y)
//! ```
//!
//! So a 4-way tensor product `t = t1 ⊗ t2 ⊗ t3 ⊗ t4` evaluated at point `r`:
//! ```text
//! t(r) = MLE_{t1}(r[0:n1]) * MLE_{t2}(r[n1:n2]) * MLE_{t3}(r[n2:n3]) * MLE_{t4}(r[n3:n4])
//! ```
//!
//! Complexity: O(|t1| + |t2| + |t3| + |t4|) instead of O(|t1| * |t2| * |t3| * |t4|)

use ark_std::log2;
use stark_rings::OverField;

/// Evaluate a small MLE (given as evaluations) at a point.
///
/// Input: `evals` of length `m`, point `r` of length `k` where `2^k >= m`
/// Output: `MLE(r) = Σ_i evals[i] * eq(r, bits(i))`
///
/// When `m < 2^k`, this is equivalent to evaluating the MLE of `evals` padded with zeros
/// to length `2^k`. The padding contributes 0 to the sum, so we only iterate over actual
/// elements but use ALL coordinates of `r` when computing eq weights.
///
/// Complexity: O(|evals| * |r|)
pub fn eval_small_mle<F: OverField>(evals: &[F], r: &[F]) -> F {
    let n = evals.len();
    if n == 0 {
        return F::zero();
    }
    
    // Use ALL coordinates of r (this handles padding correctly)
    let num_vars = r.len();
    
    let mut result = F::ZERO;
    for (i, &eval_i) in evals.iter().enumerate() {
        if eval_i == F::ZERO {
            continue; // Skip zero terms
        }
        
        // Compute eq(r, bits(i)) using ALL num_vars bits
        // For i < n, the high bits (beyond log2(n)) are all 0
        let mut eq_val = F::one();
        for j in 0..num_vars {
            let bit_j = (i >> j) & 1;
            if bit_j == 1 {
                eq_val *= r[j];
            } else {
                eq_val *= F::one() - r[j];
            }
        }
        
        result += eval_i * eq_val;
    }
    
    result
}

/// Evaluate `tensor(c)` MLE at a point.
///
/// `tensor(c)` for `c = [c_0, c_1, ..., c_{k-1}]` produces the vector:
/// `[Π_i (1-c_i if b_i=0 else c_i) for b in {0,1}^k]`
///
/// The MLE of `tensor(c)` evaluated at point `r` is:
/// `Π_i (1 - r_i) + r_i * (some linear combo)`
///
/// Actually, for tensor(c), we have:
/// `MLE_{tensor(c)}(r) = Π_i ((1-c_i)(1-r_i) + c_i * r_i)`
///
/// Complexity: O(|c|)
pub fn eval_tensor_mle<F: OverField>(c: &[F], r: &[F]) -> F {
    assert!(r.len() >= c.len(), "point must have at least {} coordinates", c.len());
    
    let mut result = F::one();
    for (i, &c_i) in c.iter().enumerate() {
        // tensor(c)[j] at position j (with bits b_0, b_1, ...) is:
        // Π_i (c_i if b_i=1 else 1-c_i)
        //
        // MLE evaluation: eq(r, j) * tensor(c)[j] summed over j
        // = Π_i ((1-r_i)(1-c_i) + r_i * c_i)
        let term = (F::one() - r[i]) * (F::one() - c_i) + r[i] * c_i;
        result *= term;
    }
    
    result
}

/// Configuration for tensor product evaluation.
///
/// Describes how the evaluation point `r` is split across the tensor factors.
#[derive(Clone, Debug)]
pub struct TensorConfig {
    /// Sizes of each factor (in elements, not log)
    pub factor_sizes: Vec<usize>,
}

impl TensorConfig {
    /// Create a new tensor config from factor sizes.
    pub fn new(factor_sizes: Vec<usize>) -> Self {
        Self { factor_sizes }
    }
    
    /// Total number of variables (log of product of sizes, rounded up).
    pub fn total_vars(&self) -> usize {
        let total_size: usize = self.factor_sizes.iter().product();
        log2(total_size.next_power_of_two()) as usize
    }
    
    /// Number of variables for each factor.
    pub fn vars_per_factor(&self) -> Vec<usize> {
        self.factor_sizes
            .iter()
            .map(|&s| log2(s.next_power_of_two()) as usize)
            .collect()
    }
    
    /// Starting index in `r` for each factor.
    pub fn factor_offsets(&self) -> Vec<usize> {
        let vars = self.vars_per_factor();
        let mut offsets = Vec::with_capacity(vars.len());
        let mut offset = 0;
        for v in vars {
            offsets.push(offset);
            offset += v;
        }
        offsets
    }
}

/// Evaluate a 4-way tensor product MLE at a point.
///
/// Given:
/// - `t = t1 ⊗ t2 ⊗ t3 ⊗ t4` (conceptually, using the nested loop ordering from utils::tensor_product)
/// - `t1_evals`, `t2_evals`, `t3_evals`, `t4_evals` (factor evaluations)
/// - `r` (evaluation point)
///
/// The tensor product is computed as:
/// ```text
/// for a in t1:
///   for b in t2:
///     for c in t3:
///       for d in t4:
///         result.push(a * b * c * d)
/// ```
///
/// So the index is: i = ((i1 * |t2| + i2) * |t3| + i3) * |t4| + i4
/// And the bit layout is: [t4 bits (lowest)] [t3 bits] [t2 bits] [t1 bits (highest)]
///
/// Returns: `MLE_t(r) = MLE_{t1}(r4) * MLE_{t2}(r3) * MLE_{t3}(r2) * MLE_{t4}(r1)`
/// where `r = (r1, r2, r3, r4)` are chunks from lowest to highest bits.
///
/// Complexity: O(|t1| + |t2| + |t3| + |t4|) instead of O(|t1| * |t2| * |t3| * |t4|)
pub fn eval_tensor4_mle<F: OverField>(
    t1_evals: &[F],
    t2_evals: &[F],
    t3_evals: &[F],
    t4_evals: &[F],
    r: &[F],
) -> F {
    // Factor sizes in REVERSE order (innermost to outermost)
    let config = TensorConfig::new(vec![
        t4_evals.len(),
        t3_evals.len(),
        t2_evals.len(),
        t1_evals.len(),
    ]);
    
    let offsets = config.factor_offsets();
    let vars = config.vars_per_factor();
    
    // Split r into chunks for each factor (innermost to outermost)
    let r4 = &r[offsets[0]..offsets[0] + vars[0]];
    let r3 = &r[offsets[1]..offsets[1] + vars[1]];
    let r2 = &r[offsets[2]..offsets[2] + vars[2]];
    let r1 = &r[offsets[3]..offsets[3] + vars[3]];
    
    // Evaluate each factor's MLE
    let v1 = eval_small_mle(t1_evals, r1);
    let v2 = eval_small_mle(t2_evals, r2);
    let v3 = eval_small_mle(t3_evals, r3);
    let v4 = eval_small_mle(t4_evals, r4);
    
    v1 * v2 * v3 * v4
}

/// Check if n is a power of 2
fn is_power_of_two(n: usize) -> bool {
    n > 0 && (n & (n - 1)) == 0
}

/// Check if tensor factors would use the optimized path
pub fn would_use_optimized_path(
    c_z_len: usize,
    s_prime_len: usize,
    d_prime_powers_len: usize,
    x_powers_len: usize,
) -> bool {
    let kappa = 1 << c_z_len;
    is_power_of_two(kappa) 
        && is_power_of_two(s_prime_len)
        && is_power_of_two(d_prime_powers_len)
        && is_power_of_two(x_powers_len)
}

/// Print tensor factor sizes and optimization status
pub fn print_tensor_optimization_status(
    c_z_len: usize,
    s_prime_len: usize,
    d_prime_powers_len: usize,
    x_powers_len: usize,
) {
    let kappa = 1 << c_z_len;
    let optimized = would_use_optimized_path(c_z_len, s_prime_len, d_prime_powers_len, x_powers_len);
    
    println!("=== Tensor Evaluation Optimization ===");
    println!("  tensor(c_z) size (κ): {} (pow2: {})", kappa, is_power_of_two(kappa));
    println!("  s_prime size (k×d):   {} (pow2: {})", s_prime_len, is_power_of_two(s_prime_len));
    println!("  d_prime_powers (ℓ):   {} (pow2: {})", d_prime_powers_len, is_power_of_two(d_prime_powers_len));
    println!("  x_powers (d):         {} (pow2: {})", x_powers_len, is_power_of_two(x_powers_len));
    println!("  Using optimized path: {}", if optimized { "YES ✓" } else { "NO (dense fallback)" });
    if optimized {
        let ops = kappa + s_prime_len + d_prime_powers_len + x_powers_len;
        println!("  Estimated ops: ~{}", ops);
    } else {
        let ops = kappa * s_prime_len * d_prime_powers_len * x_powers_len;
        println!("  Estimated ops (dense): ~{}", ops);
    }
    println!("======================================");
}

/// Evaluate `t(z) = tensor(c_z) ⊗ s' ⊗ d_powers ⊗ x_powers` at a point.
///
/// This is the optimized version of the calculation in `cm.rs`.
///
/// The tensor product uses the nested loop ordering from utils::tensor_product:
/// ```text
/// for a in tensor(c_z):        // outermost, highest bits
///   for b in s_prime:
///     for c in d_prime_powers:
///       for d in x_powers:     // innermost, lowest bits
///         result.push(a * b * c * d)
/// ```
///
/// # Important: Factored Evaluation
/// The factored approach `MLE(r) = Π MLE_factor(r_factor)` only works correctly
/// when ALL factor sizes are powers of 2. This is because the bit decomposition
/// of tensor indices only aligns with factor indices when sizes are powers of 2.
///
/// When factor sizes aren't all powers of 2, we fall back to dense evaluation
/// of the full tensor product.
///
/// # Important: Zero Padding
/// The prover pads the tensor product with zeros to length 2^nvars.
/// When evaluating at point r of length nvars, we must account for this:
/// - For indices in the tensor range: contribute tensor[i] * eq(r, bits(i))
/// - For indices in padding range: contribute 0 (don't affect the sum)
///
/// # Arguments
/// - `c_z`: Challenge vector (length log κ)
/// - `s_prime`: Decomposition challenges (length k*d)
/// - `d_prime_powers`: Powers of d' (length ℓ)
/// - `x_powers`: Monomial powers (length d)
/// - `r`: Evaluation point (length nvars, may be longer than tensor vars)
///
/// # Complexity
/// - O(κ + k*d + ℓ + d + nvars) when all factor sizes are powers of 2
/// - O(κ * k*d * ℓ * d) otherwise (falls back to dense)
///
/// Note: We precompute tensor(c_z) (size κ = 2^log_κ, typically small like 8-16)
/// instead of using eval_tensor_mle, because the tensor() function has specific
/// bit ordering that doesn't match the standard MLE eq formula.
pub fn eval_t_z_optimized<F: OverField>(
    c_z: &[F],
    s_prime: &[F],
    d_prime_powers: &[F],
    x_powers: &[F],
    r: &[F],
) -> F {
    // Precompute tensor(c_z) - this is O(κ) where κ = 2^|c_z|
    // For typical κ = 8 or 16, this is negligible.
    let tensor_c_z = tensor(c_z);
    let kappa = tensor_c_z.len();
    
    // Check if all factor sizes are powers of 2 (required for factored evaluation)
    let all_pow2 = is_power_of_two(kappa) 
        && is_power_of_two(s_prime.len())
        && is_power_of_two(d_prime_powers.len())
        && is_power_of_two(x_powers.len());
    
    if !all_pow2 {
        // Fall back to dense evaluation when factors aren't all powers of 2
        return eval_t_z_dense(&tensor_c_z, s_prime, d_prime_powers, x_powers, r);
    }
    
    // Factor sizes in REVERSE order (innermost to outermost)
    let config = TensorConfig::new(vec![
        x_powers.len(),
        d_prime_powers.len(),
        s_prime.len(),
        kappa,
    ]);
    
    let offsets = config.factor_offsets();
    let vars = config.vars_per_factor();
    let tensor_vars = vars.iter().sum::<usize>();
    
    // Split r into chunks (innermost to outermost)
    let r4 = &r[offsets[0]..offsets[0] + vars[0]]; // x_powers
    let r3 = &r[offsets[1]..offsets[1] + vars[1]]; // d_prime_powers
    let r2 = &r[offsets[2]..offsets[2] + vars[2]]; // s_prime
    let r1 = &r[offsets[3]..offsets[3] + vars[3]]; // tensor(c_z)
    
    // Evaluate each factor's MLE using standard small MLE eval
    let v1 = eval_small_mle(&tensor_c_z, r1);
    let v2 = eval_small_mle(s_prime, r2);
    let v3 = eval_small_mle(d_prime_powers, r3);
    let v4 = eval_small_mle(x_powers, r4);
    
    let mut result = v1 * v2 * v3 * v4;
    
    // Account for zero padding: for indices in tensor range, high bits are 0.
    // eq(r, bits(i)) for i < tensor_size includes factor (1 - r[j]) for j >= tensor_vars.
    for j in tensor_vars..r.len() {
        result *= F::one() - r[j];
    }
    
    result
}

/// Dense fallback for eval_t_z when factor sizes aren't powers of 2.
/// 
/// Builds the full tensor product and evaluates the MLE directly.
fn eval_t_z_dense<F: OverField>(
    tensor_c_z: &[F],
    s_prime: &[F],
    d_prime_powers: &[F],
    x_powers: &[F],
    r: &[F],
) -> F {
    // Build full tensor product
    let part1 = tensor_product_pair(tensor_c_z, s_prime);
    let part2 = tensor_product_pair(&part1, d_prime_powers);
    let t_z = tensor_product_pair(&part2, x_powers);
    
    // Evaluate MLE (handles padding implicitly since eval_small_mle only iterates over actual elements)
    eval_small_mle(&t_z, r)
}

/// Compute tensor(c) = fold over tensor_product with [1-c_i, c_i]
/// 
/// This matches the `tensor` function in utils.rs
fn tensor<F: OverField>(c: &[F]) -> Vec<F> {
    c.iter().fold(vec![F::one()], |acc, x| {
        tensor_product_pair(&acc, &[F::one() - *x, *x])
    })
}

/// Simple tensor product of two vectors: a ⊗ b
fn tensor_product_pair<F: OverField>(a: &[F], b: &[F]) -> Vec<F> {
    let mut result = Vec::with_capacity(a.len() * b.len());
    for &ai in a {
        for &bi in b {
            result.push(ai * bi);
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use stark_rings::cyclotomic_ring::models::goldilocks::RqPoly as R;
    use stark_rings::Ring;
    
    #[test]
    fn test_eval_small_mle() {
        // MLE of [1, 2] over 1 variable
        // MLE(0) = 1, MLE(1) = 2
        // MLE(r) = (1-r)*1 + r*2 = 1 + r
        let evals = vec![R::from(1u128), R::from(2u128)];
        let r = vec![R::from(3u128)];
        let result = eval_small_mle(&evals, &r);
        assert_eq!(result, R::from(4u128)); // 1 + 3 = 4
    }
    
    #[test]
    fn test_eval_tensor_mle() {
        // tensor([c]) = [1-c, c]
        // MLE at r: (1-r)(1-c) + r*c
        let c = vec![R::from(2u128)];
        let r = vec![R::from(3u128)];
        let result = eval_tensor_mle(&c, &r);
        // (1-3)(1-2) + 3*2 = (-2)(-1) + 6 = 2 + 6 = 8
        assert_eq!(result, R::from(8u128));
    }
    
    #[test]
    fn test_tensor4_consistency() {
        // Test that eval_tensor4_mle matches dense evaluation
        let t1 = vec![R::from(1u128), R::from(2u128)];
        let t2 = vec![R::from(3u128), R::from(4u128)];
        let t3 = vec![R::from(5u128)];
        let t4 = vec![R::from(6u128), R::from(7u128)];
        
        // Full tensor product (for comparison)
        let mut full = Vec::new();
        for &a in &t1 {
            for &b in &t2 {
                for &c in &t3 {
                    for &d in &t4 {
                        full.push(a * b * c * d);
                    }
                }
            }
        }
        
        // Evaluate both ways at a random point
        let r = vec![
            R::from(10u128), // for t1 (1 bit)
            R::from(11u128), // for t2 (1 bit)
            // t3 has 1 element, needs 0 bits
            R::from(12u128), // for t4 (1 bit)
        ];
        
        let optimized = eval_tensor4_mle(&t1, &t2, &t3, &t4, &r);
        let dense = eval_small_mle(&full, &r);
        
        assert_eq!(optimized, dense);
    }
    
    #[test]
    fn test_eval_t_z_vs_dense() {
        use stark_rings_poly::mle::DenseMultilinearExtension;
        use crate::utils::{tensor as utils_tensor, tensor_product as utils_tensor_product};
        
        // Same parameters as cm.rs uses
        let c_z = vec![R::from(2u128), R::from(3u128)]; // log_kappa = 2
        let s_prime = vec![R::from(5u128), R::from(6u128), R::from(7u128), R::from(8u128)]; // k*d = 4
        let d_prime_powers = vec![R::from(1u128), R::from(9u128)]; // l = 2
        let x_powers = vec![R::from(1u128), R::from(10u128)]; // d = 2
        
        // Build dense t_z using same method as cm.rs calculate_t_z
        let tensor_c_z = utils_tensor(&c_z);
        let part1 = utils_tensor_product(&tensor_c_z, &s_prime);
        let part2 = utils_tensor_product(&part1, &d_prime_powers);
        let t_z_dense = utils_tensor_product(&part2, &x_powers);
        
        // Total vars
        let nvars = ark_std::log2(t_z_dense.len()) as usize;
        
        // Random evaluation point
        let ro: Vec<R> = (0..nvars).map(|i| R::from((11 + i) as u128)).collect();
        
        // Dense MLE evaluation (original method)
        let mle = DenseMultilinearExtension::from_evaluations_vec(nvars, t_z_dense.clone());
        let dense_result = mle.evaluate(&ro).unwrap();
        
        // Optimized evaluation
        let optimized_result = eval_t_z_optimized(&c_z, &s_prime, &d_prime_powers, &x_powers, &ro);
        
        assert_eq!(optimized_result, dense_result);
    }
    
    #[test]
    fn test_eval_t_z_vs_dense_with_padding() {
        use stark_rings_poly::mle::DenseMultilinearExtension;
        use crate::utils::{tensor as utils_tensor, tensor_product as utils_tensor_product};
        
        // Parameters that result in non-power-of-2 tensor product
        let c_z = vec![R::from(2u128), R::from(3u128)]; // log_kappa = 2, kappa = 4
        let s_prime = vec![R::from(5u128), R::from(6u128), R::from(7u128)]; // k*d = 3 (not power of 2)
        let d_prime_powers = vec![R::from(1u128), R::from(9u128), R::from(81u128)]; // l = 3 (not power of 2)
        let x_powers = vec![R::from(1u128), R::from(10u128)]; // d = 2
        
        // Build dense t_z using same method as cm.rs calculate_t_z
        let tensor_c_z = utils_tensor(&c_z);
        let part1 = utils_tensor_product(&tensor_c_z, &s_prime);
        let part2 = utils_tensor_product(&part1, &d_prime_powers);
        let mut t_z_dense = utils_tensor_product(&part2, &x_powers);
        
        // Pad to larger power of 2 (simulating cm.rs behavior)
        let padded_nvars = 10; // Much larger than needed
        let padded_len = 1 << padded_nvars;
        t_z_dense.resize(padded_len, R::ZERO);
        
        // Random evaluation point (now has padded_nvars coordinates)
        let ro: Vec<R> = (0..padded_nvars).map(|i| R::from((11 + i) as u128)).collect();
        
        // Dense MLE evaluation (original method with padding)
        let mle = DenseMultilinearExtension::from_evaluations_vec(padded_nvars, t_z_dense.clone());
        let dense_result = mle.evaluate(&ro).unwrap();
        
        // Optimized evaluation (should handle padding automatically)
        let optimized_result = eval_t_z_optimized(&c_z, &s_prime, &d_prime_powers, &x_powers, &ro);
        
        
        assert_eq!(optimized_result, dense_result);
    }
}
