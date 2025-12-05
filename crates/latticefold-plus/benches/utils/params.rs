pub mod quick {

    pub const SPLIT: &[(usize, usize, usize)] = &[
        (16384, 2, 1),    
        (32768, 2, 2),    
        (65536, 4, 2),    
        (131072, 4, 3),  
    ];

    /// Set check parameters: (set_size, num_batches)
    /// - set_size: monomial set cardinality
    /// - num_batches: number of sets to check simultaneously (tests batching efficiency)
    pub const SETCHK: &[(usize, usize)] = &[
        (256, 1),
        (512, 2),
        (1024, 4),
    ];

    /// Range check parameters: (witness_size, k, kappa)
    ///
    /// Constraint: witness_size >= kappa * k * d * l * d
    /// For FrogRing (d=16, l=22): witness_size >= kappa * k * 5632
    ///
    /// Parameters:
    /// - witness_size: length of witness vector
    /// - k: decomposition width (determines range B = (d/2)^k)
    /// - kappa: number of commitment rows (security parameter)
    pub const RGCHK: &[(usize, usize, usize)] = &[
        (16384, 2, 1),
        (32768, 2, 2),
        (65536, 4, 2),
    ];

    /// Commitment transformation parameters: (L, witness_size, k, kappa)
    ///
    /// The commitment transformation builds on range check, so the same constraint applies:
    /// witness_size >= kappa * k * d * l * d
    /// For FrogRing (d=16, l=22): witness_size >= kappa * k * 5632
    ///
    /// Parameters:
    /// - L: number of instances to transform/fold
    /// - witness_size: length of witness vector
    /// - k: decomposition width (determines range B = (d/2)^k)
    /// - kappa: number of commitment rows (security parameter)
    pub const CM: &[(usize, usize, usize, usize)] = &[
        (2, 65536, 2, 2),
        (4, 65536, 2, 2),
        (8, 65536, 2, 2),
    ];

    /// Multilinear folding parameters: (L, n, k, κ, B)
    ///
    /// Parameters:
    /// - L: number of instances to fold (higher L = better amortization)
    /// - n: witness size (length of witness vector after decomposition)
    /// - k: decomposition width (determines range bound)
    /// - κ (kappa): number of commitment rows (security parameter)
    /// - B: norm bound parameter
    ///
    /// These parameters hold n, k, kappa, and B constant to isolate L's impact
    pub const MLIN: &[(usize, usize, usize, usize, usize)] = &[
        (2, 65536, 2, 2, 50),
        (4, 65536, 2, 2, 50),
        (8, 65536, 2, 2, 50),
    ];

    /// Decomposition parameters: (n, k, κ, B)
    ///
    /// This is Construction 5.3 (Π_decomp,B) - decomposes norm B² to 2×B
    ///
    /// Constraint: n >= κ * k * d * l * d
    /// For FrogRing (d=16, l=22): n >= κ * k * 5632
    ///
    /// Parameters:
    /// - n: witness size (before decomposition, with norm B²)
    /// - k: decomposition width
    /// - κ (kappa): number of commitment rows
    /// - B: norm bound for output (input has bound B²)
    ///
    /// These parameters hold n, k, and B constant to isolate kappa scaling
    pub const DECOMP_PARAMS: &[(usize, usize, usize, usize)] = &[
        (131072, 4, 2, 50),
        (131072, 4, 3, 50),
        (131072, 4, 4, 50),
    ];
}