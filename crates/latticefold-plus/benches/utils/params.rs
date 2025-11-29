pub mod quick {

    pub const SPLIT: &[(usize, usize, usize)] = &[
        (16384, 2, 1),    
        (32768, 2, 2),    
        (65536, 4, 2),    
        (131072, 4, 3),  
    ];
}