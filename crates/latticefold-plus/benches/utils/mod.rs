//! Parameter sets and utilities for LatticeFold+ benchmarks.
//!
//! This module defines the shared infrastructure for all protocol benchmarks:
//! - Common helper functions and traits (see [`helpers`])
//! - Standardized parameter sets for each protocol
//!
//! All parameter sets follow consistent naming conventions:
//! - `WITNESS_SCALING`: Varies witness size while keeping other params fixed
//! - `K_SCALING`: Varies decomposition width k
//! - `KAPPA_SCALING`: Varies security parameter κ (kappa)
//! - `FOLDING_ARITY`: Varies number of instances L for batched protocols

pub mod helpers;

/// Range check protocol (Constructions 4.3-4.4) parameter sets.
///
/// Benchmarks the double commitment and gadget decomposition protocols
/// for verifying that witness coefficients lie within [-B, B].
pub mod range_check {
    /// Witness size scaling benchmark.
    ///
    /// Measures performance across witness sizes from 32K to 512K.
    /// Fixed parameters: k=2, κ=2.
    ///
    /// Format: `(witness_size, k, kappa)`
    pub const WITNESS_SCALING: &[(usize, usize, usize)] = &[
        (32768, 2, 2),
        (65536, 2, 2),
        (131072, 2, 2),
        (262144, 2, 2),
        (524288, 2, 2),
    ];

    /// Decomposition width scaling benchmark.
    ///
    /// Measures how performance scales with k ∈ [2,3,4,5].
    /// Fixed parameters: witness_size=65536, κ=2.
    ///
    /// Format: `(witness_size, k, kappa)`
    pub const K_SCALING: &[(usize, usize, usize)] =
        &[(65536, 2, 2), (65536, 3, 2), (65536, 4, 2), (65536, 5, 2)];

    /// Security parameter scaling benchmark.
    ///
    /// Measures how performance scales with κ ∈ [2,3,4,5].
    /// Fixed parameters: witness_size=65536, k=2.
    ///
    /// Format: `(witness_size, k, kappa)`
    pub const KAPPA_SCALING: &[(usize, usize, usize)] =
        &[(65536, 2, 2), (65536, 2, 3), (65536, 2, 4), (65536, 2, 5)];
}

/// Set check protocol (Construction 4.2) parameter sets.
///
/// Benchmarks verification that matrices contain monomials (exactly one
/// non-zero entry per row/column) using sumcheck protocols.
pub mod set_check {
    /// Set size and batching benchmark.
    ///
    /// Tests various combinations of set sizes (256, 512, 1024) and
    /// batch counts (1, 2, 4) to measure batching efficiency.
    ///
    /// Format: `(set_size, num_batches)`
    pub const SET_SIZES: &[(usize, usize)] = &[(256, 1), (512, 1), (1024, 1), (512, 2), (512, 4)];
}

/// Split function (Construction 4.1) parameter sets.
///
/// Benchmarks gadget decomposition on double commitments, converting
/// ring element matrices to base ring scalar representations.
pub mod split {
    /// Witness and commitment parameter scaling benchmark.
    ///
    /// Tests various combinations of witness size, first decomposition
    /// width k_first, and security parameter κ.
    ///
    /// Format: `(witness_size, k_first, kappa)`
    pub const WITNESS_SCALING: &[(usize, usize, usize)] =
        &[(16384, 2, 1), (32768, 4, 2), (65536, 6, 3), (131072, 8, 4)];

    /// First decomposition width scaling benchmark.
    ///
    /// Measures how performance scales with k_first ∈ [2,4,6,8].
    /// Fixed parameters: witness_size=131072, κ=2.
    ///
    /// Format: `(witness_size, k_first, kappa)`
    pub const K_FIRST_SCALING: &[(usize, usize, usize)] = &[
        (131072, 2, 2),
        (131072, 4, 2),
        (131072, 6, 2),
        (131072, 8, 2),
    ];

    /// Security parameter scaling benchmark.
    ///
    /// Measures how performance scales with κ ∈ [1,2,3,4].
    /// Fixed parameters: witness_size=131072, k_first=4.
    ///
    /// Format: `(witness_size, k_first, kappa)`
    pub const KAPPA_SCALING: &[(usize, usize, usize)] = &[
        (131072, 4, 1),
        (131072, 4, 2),
        (131072, 4, 3),
        (131072, 4, 4),
    ];
}

/// Commitment transformation (Construction 4.5) parameter sets.
///
/// Benchmarks conversion of double commitments from range check into
/// folded commitments for the main LatticeFold+ protocol.
pub mod commitment_transform {
    /// Folding arity scaling benchmark.
    ///
    /// Measures batching efficiency across L ∈ [2,3,4,5,6,7,8] instances.
    /// Fixed parameters: witness_size=65536, k=2, κ=2.
    ///
    /// Format: `(L, witness_size, k, kappa)`
    pub const FOLDING_ARITY: &[(usize, usize, usize, usize)] = &[
        (2, 65536, 2, 2),
        (3, 65536, 2, 2),
        (4, 65536, 2, 2),
        (5, 65536, 2, 2),
        (6, 65536, 2, 2),
        (7, 65536, 2, 2),
        (8, 65536, 2, 2),
    ];
}

/// Multilinear folding (Construction 5.2) parameter sets.
///
/// Benchmarks the main folding protocol that batches multiple LinB instances
/// into a single LinB2 instance with improved amortization.
pub mod multilinear_fold {
    /// Folding arity scaling benchmark.
    ///
    /// Measures batching efficiency across L ∈ [2,3,4,5,6,7,8] instances.
    /// Fixed parameters: n=65536, k=2, κ=2, B=50.
    ///
    /// Format: `(L, n, k, kappa, B)`
    pub const FOLDING_ARITY: &[(usize, usize, usize, usize, usize)] = &[
        (2, 65536, 2, 2, 50),
        (3, 65536, 2, 2, 50),
        (4, 65536, 2, 2, 50),
        (5, 65536, 2, 2, 50),
        (6, 65536, 2, 2, 50),
        (7, 65536, 2, 2, 50),
        (8, 65536, 2, 2, 50),
    ];

    /// Decomposition width scaling benchmark.
    ///
    /// Measures how performance scales with k ∈ [2,3,4].
    /// Note: witness size n is adjusted proportionally with k.
    /// Fixed parameters: L=4, κ=2, B=50.
    ///
    /// Format: `(L, n, k, kappa, B)`
    pub const K_SCALING: &[(usize, usize, usize, usize, usize)] = &[
        (4, 65536, 2, 2, 50),
        (4, 98304, 3, 2, 50),
        (4, 131072, 4, 2, 50),
    ];

    /// Large witness scaling benchmark.
    ///
    /// Measures performance on very large witnesses from 128K to 512K.
    /// Fixed parameters: L=4, k=2, κ=2, B=50.
    ///
    /// Format: `(L, n, k, kappa, B)`
    pub const LARGE_WITNESS: &[(usize, usize, usize, usize, usize)] = &[
        (4, 131072, 2, 2, 50),
        (4, 262144, 2, 2, 50),
        (4, 524288, 2, 2, 50),
    ];

    /// Security parameter scaling benchmark.
    ///
    /// Measures how performance scales with κ ∈ [2,3,4,5].
    /// Fixed parameters: L=4, n=65536, k=2, B=50.
    ///
    /// Format: `(L, n, k, kappa, B)`
    pub const KAPPA_SCALING: &[(usize, usize, usize, usize, usize)] = &[
        (4, 65536, 2, 2, 50),
        (4, 65536, 2, 3, 50),
        (4, 65536, 2, 4, 50),
        (4, 65536, 2, 5, 50),
    ];
}

/// Single instance folding (Construction 5.1) parameter sets.
///
/// Benchmarks the baseline folding protocol with L=1, reducing one LinB
/// instance to one LinB2 instance.
pub mod single_instance_fold {
    /// Witness size scaling benchmark.
    ///
    /// Measures performance across witness sizes from 32K to 128K.
    /// Fixed parameters: k=2, κ=2, B=50.
    ///
    /// Format: `(n, k, kappa, B)`
    pub const WITNESS_SCALING: &[(usize, usize, usize, usize)] =
        &[(32768, 2, 2, 50), (65536, 2, 2, 50), (131072, 2, 2, 50)];
}

/// Decomposition protocol (Construction 5.3) parameter sets.
///
/// Benchmarks splitting a LinB2 instance (norm B²) into two LinB instances
/// (norm B each), critical for IVC/PCD to prevent norm explosion.
pub mod decomposition {
    /// Witness size scaling benchmark.
    ///
    /// Measures performance across witness sizes from 32K to 128K.
    /// Fixed parameters: k=2, κ=2, B=50.
    ///
    /// Format: `(n, k, kappa, B)`
    pub const WITNESS_SCALING: &[(usize, usize, usize, usize)] =
        &[(32768, 2, 2, 50), (65536, 2, 2, 50), (131072, 2, 2, 50)];
}
