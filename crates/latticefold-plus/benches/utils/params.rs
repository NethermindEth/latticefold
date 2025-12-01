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
}