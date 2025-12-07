//! Benchmark parameter sets for LatticeFold+ protocols.
//!
//! All parameter sets satisfy the fundamental constraint:
//! `witness_size >= kappa * k * d * l * d`
//!
//! For FrogRing (d=16, l≈22): `witness_size >= kappa * k * 5632`

/// Range check parameters.
pub mod range_check {
    /// Witness size scaling with fixed k=2, kappa=2.
    ///
    /// Parameters: (witness_size, k, kappa)
    pub const WITNESS_SCALING: &[(usize, usize, usize)] = &[
        (32768, 2, 2),
        (65536, 2, 2),
        (131072, 2, 2),
        (262144, 2, 2),
        (524288, 2, 2),
    ];

    /// Decomposition width scaling with fixed witness_size=65536, kappa=2.
    ///
    /// Parameters: (witness_size, k, kappa)
    pub const K_SCALING: &[(usize, usize, usize)] =
        &[(65536, 2, 2), (65536, 3, 2), (65536, 4, 2), (65536, 5, 2)];

    /// Security parameter scaling with fixed witness_size=65536, k=2.
    ///
    /// Parameters: (witness_size, k, kappa)
    pub const KAPPA_SCALING: &[(usize, usize, usize)] =
        &[(65536, 2, 2), (65536, 2, 3), (65536, 2, 4), (65536, 2, 5)];
}

/// Commitment transformation parameters.
pub mod commitment_transform {
    /// Folding arity scaling with fixed witness_size=65536, k=2, kappa=2.
    ///
    /// Parameters: (L, witness_size, k, kappa)
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

/// Multilinear folding parameters.
pub mod multilinear_fold {
    /// Folding arity scaling with fixed n=65536, k=2, kappa=2, B=50.
    ///
    /// Parameters: (L, n, k, kappa, B)
    pub const FOLDING_ARITY: &[(usize, usize, usize, usize, usize)] = &[
        (2, 65536, 2, 2, 50),
        (3, 65536, 2, 2, 50),
        (4, 65536, 2, 2, 50),
        (5, 65536, 2, 2, 50),
        (6, 65536, 2, 2, 50),
        (7, 65536, 2, 2, 50),
        (8, 65536, 2, 2, 50),
    ];

    /// Decomposition width scaling with fixed L=4, kappa=2, B=50.
    ///
    /// Parameters: (L, n, k, kappa, B)
    pub const K_SCALING: &[(usize, usize, usize, usize, usize)] = &[
        (4, 65536, 2, 2, 50),
        (4, 90112, 3, 2, 50),
        (4, 131072, 4, 2, 50),
    ];

    /// Large witness scaling with fixed L=4, k=2, kappa=2, B=50.
    ///
    /// Parameters: (L, n, k, kappa, B)
    pub const LARGE_WITNESS: &[(usize, usize, usize, usize, usize)] = &[
        (4, 131072, 2, 2, 50),
        (4, 262144, 2, 2, 50),
        (4, 524288, 2, 2, 50),
    ];

    /// Security parameter scaling with fixed L=4, n=65536, k=2, B=50.
    ///
    /// Parameters: (L, n, k, kappa, B)
    pub const KAPPA_SCALING: &[(usize, usize, usize, usize, usize)] = &[
        (4, 65536, 2, 2, 50),
        (4, 65536, 2, 3, 50),
        (4, 65536, 2, 4, 50),
        (4, 65536, 2, 5, 50),
    ];
}

/// Single instance folding parameters (Construction 5.1).
pub mod single_instance_fold {
    /// Witness scaling with fixed k=4, kappa=2, B=50.
    ///
    /// Parameters: (n, k, kappa, B)
    pub const WITNESS_SCALING: &[(usize, usize, usize, usize)] =
        &[(49152, 4, 2, 50), (65536, 4, 2, 50), (98304, 4, 2, 50)];
}

/// Decomposition parameters (Construction 5.3).
pub mod decomposition {
    /// Witness scaling with fixed k=4, kappa=2, B=50.
    ///
    /// Parameters: (n, k, kappa, B)
    pub const WITNESS_SCALING: &[(usize, usize, usize, usize)] =
        &[(49152, 4, 2, 50), (65536, 4, 2, 50), (98304, 4, 2, 50)];
}

/// Set check parameters (Construction 4.2).
pub mod set_check {
    /// Set size and batching scaling.
    ///
    /// Parameters: (set_size, num_batches)
    pub const SET_SIZES: &[(usize, usize)] = &[(256, 1), (512, 2), (1024, 4)];
}

/// Split function parameters (Construction 4.1).
pub mod split {
    /// Witness size scaling with varying k_first and kappa.
    ///
    /// Parameters: (witness_size, k_first, kappa)
    pub const WITNESS_SCALING: &[(usize, usize, usize)] =
        &[(16384, 2, 2), (65536, 4, 2), (131072, 8, 2)];

    /// K_first scaling with fixed witness_size=131072, kappa=2.
    ///
    /// Parameters: (witness_size, k_first, kappa)
    pub const K_FIRST_SCALING: &[(usize, usize, usize)] = &[
        (131072, 2, 2),
        (131072, 4, 2),
        (131072, 6, 2),
        (131072, 8, 2),
    ];

    /// Kappa scaling with fixed witness_size=131072, k_first=4.
    ///
    /// Parameters: (witness_size, k_first, kappa)
    pub const KAPPA_SCALING: &[(usize, usize, usize)] = &[
        (131072, 4, 1),
        (131072, 4, 2),
        (131072, 4, 3),
        (131072, 4, 4),
    ];
}
