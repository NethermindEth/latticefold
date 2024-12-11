#[macro_export]
macro_rules! define_params {
    ($w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
        paste::paste! {
            #[derive(Clone)]
            struct [<DecompParamsWithB $b W $w b $b_small K $k>];

            impl DecompositionParams for [<DecompParamsWithB $b W $w b $b_small K $k>] {
                const B: u128 = $b;
                const L: usize = $l;
                const B_SMALL: usize = $b_small;
                const K: usize = $k;
            }
        }
    };
}

#[macro_export]
macro_rules! run_goldilocks_ajai_benchmarks {
    ($group: ident) => {
//         run_single_goldilocks_benchmark!(&mut $group, 1, 32768);
    };
}

#[macro_export]
macro_rules! run_babybear_ajai_benchmarks {
    ($group: ident) => {
//         run_single_babybear_benchmark!(&mut $group, 1, 32768);
    };
}

#[macro_export]
macro_rules! run_starkprime_ajai_benchmarks {
    ($group: ident) => {
//         run_single_starkprime_benchmark!(&mut $group, 1, 32768);
    };
}

#[macro_export]
macro_rules! run_frog_ajai_benchmarks {
    ($group: ident) => {
//         run_single_frog_benchmark!(&mut $group, 1, 32768);
    };
}

#[macro_export]
macro_rules! run_goldilocks_benchmarks {
    ($group: ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        run_single_goldilocks_benchmark!(&mut $group, 1, 7, 512, 256, 8, 2, 8);
        run_single_goldilocks_benchmark!(&mut $group, 1, 9, 512, 512, 8, 2, 9);
        run_single_goldilocks_benchmark!(&mut $group, 1, 10, 512, 2048, 6, 2, 11);
        run_single_goldilocks_benchmark!(&mut $group, 1, 12, 512, 8192, 5, 2, 13);
        run_single_goldilocks_benchmark!(&mut $group, 1, 15, 512, 65536, 4, 2, 16);
        run_single_goldilocks_benchmark!(&mut $group, 1, 22, 512, 2097152, 4, 2, 21);
//         run_single_goldilocks_benchmark!(&mut $group, 1, 39, 512, 4294967296, 2, 2, 32);
        run_single_goldilocks_benchmark!(&mut $group, 1, 8, 1024, 256, 8, 2, 8);
        run_single_goldilocks_benchmark!(&mut $group, 1, 9, 1024, 512, 8, 2, 9);
        run_single_goldilocks_benchmark!(&mut $group, 1, 10, 1024, 2048, 6, 2, 11);
        run_single_goldilocks_benchmark!(&mut $group, 1, 12, 1024, 8192, 5, 2, 13);
        run_single_goldilocks_benchmark!(&mut $group, 1, 15, 1024, 65536, 4, 2, 16);
        run_single_goldilocks_benchmark!(&mut $group, 1, 23, 1024, 2097152, 4, 2, 21);
//         run_single_goldilocks_benchmark!(&mut $group, 1, 40, 1024, 4294967296, 2, 2, 32);
        run_single_goldilocks_benchmark!(&mut $group, 1, 8, 2048, 256, 8, 2, 8);
        run_single_goldilocks_benchmark!(&mut $group, 1, 10, 2048, 512, 8, 2, 9);
        run_single_goldilocks_benchmark!(&mut $group, 1, 11, 2048, 2048, 6, 2, 11);
        run_single_goldilocks_benchmark!(&mut $group, 1, 13, 2048, 8192, 5, 2, 13);
        run_single_goldilocks_benchmark!(&mut $group, 1, 16, 2048, 65536, 4, 2, 16);
        run_single_goldilocks_benchmark!(&mut $group, 1, 24, 2048, 2097152, 4, 2, 21);
//         run_single_goldilocks_benchmark!(&mut $group, 1, 41, 2048, 4294967296, 2, 2, 32);
        run_single_goldilocks_benchmark!(&mut $group, 1, 8, 4096, 256, 8, 2, 8);
        run_single_goldilocks_benchmark!(&mut $group, 1, 10, 4096, 512, 8, 2, 9);
        run_single_goldilocks_benchmark!(&mut $group, 1, 11, 4096, 2048, 6, 2, 11);
        run_single_goldilocks_benchmark!(&mut $group, 1, 13, 4096, 8192, 5, 2, 13);
        run_single_goldilocks_benchmark!(&mut $group, 1, 17, 4096, 65536, 4, 2, 16);
        run_single_goldilocks_benchmark!(&mut $group, 1, 25, 4096, 2097152, 4, 2, 21);
//         run_single_goldilocks_benchmark!(&mut $group, 1, 42, 4096, 4294967296, 2, 2, 32);
        run_single_goldilocks_benchmark!(&mut $group, 1, 9, 8192, 256, 8, 2, 8);
        run_single_goldilocks_benchmark!(&mut $group, 1, 11, 8192, 512, 8, 2, 9);
        run_single_goldilocks_benchmark!(&mut $group, 1, 12, 8192, 2048, 6, 2, 11);
        run_single_goldilocks_benchmark!(&mut $group, 1, 14, 8192, 8192, 5, 2, 13);
        // run_single_goldilocks_benchmark!(&mut $group, 1, 17, 8192, 65536, 4, 2, 16);
        // run_single_goldilocks_benchmark!(&mut $group, 1, 26, 8192, 2097152, 4, 2, 21);
//         // run_single_goldilocks_benchmark!(&mut $group, 1, 43, 8192, 4294967296, 2, 2, 32);
    };
}

#[macro_export]
macro_rules! run_goldilocks_non_scalar_benchmarks {
    ($group: ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 7, 512, 256, 8, 2, 8);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 9, 512, 512, 8, 2, 9);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 10, 512, 2048, 6, 2, 11);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 12, 512, 8192, 5, 2, 13);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 15, 512, 65536, 4, 2, 16);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 22, 512, 2097152, 4, 2, 21);
//         run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 39, 512, 4294967296, 2, 2, 32);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 8, 1024, 256, 8, 2, 8);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 9, 1024, 512, 8, 2, 9);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 10, 1024, 2048, 6, 2, 11);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 12, 1024, 8192, 5, 2, 13);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 15, 1024, 65536, 4, 2, 16);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 23, 1024, 2097152, 4, 2, 21);
//         run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 40, 1024, 4294967296, 2, 2, 32);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 8, 2048, 256, 8, 2, 8);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 10, 2048, 512, 8, 2, 9);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 11, 2048, 2048, 6, 2, 11);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 13, 2048, 8192, 5, 2, 13);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 16, 2048, 65536, 4, 2, 16);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 24, 2048, 2097152, 4, 2, 21);
//         run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 41, 2048, 4294967296, 2, 2, 32);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 8, 4096, 256, 8, 2, 8);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 10, 4096, 512, 8, 2, 9);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 11, 4096, 2048, 6, 2, 11);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 13, 4096, 8192, 5, 2, 13);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 17, 4096, 65536, 4, 2, 16);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 25, 4096, 2097152, 4, 2, 21);
//         run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 42, 4096, 4294967296, 2, 2, 32);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 9, 8192, 256, 8, 2, 8);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 11, 8192, 512, 8, 2, 9);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 12, 8192, 2048, 6, 2, 11);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 14, 8192, 8192, 5, 2, 13);
        // run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 17, 8192, 65536, 4, 2, 16);
        // run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 26, 8192, 2097152, 4, 2, 21);
//         // run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 43, 8192, 4294967296, 2, 2, 32);
    };
}

#[macro_export]
macro_rules! run_goldilocks_degree_three_non_scalar_benchmarks {
    ($group: ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 7, 512, 256, 8, 2, 8);
        run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 9, 512, 512, 8, 2, 9);
        run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 10, 512, 2048, 6, 2, 11);
        run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 12, 512, 8192, 5, 2, 13);
        run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 15, 512, 65536, 4, 2, 16);
        run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 22, 512, 2097152, 4, 2, 21);
//         run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 39, 512, 4294967296, 2, 2, 32);
        run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 8, 1024, 256, 8, 2, 8);
        run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 9, 1024, 512, 8, 2, 9);
        run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 10, 1024, 2048, 6, 2, 11);
        run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 12, 1024, 8192, 5, 2, 13);
        run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 15, 1024, 65536, 4, 2, 16);
        run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 23, 1024, 2097152, 4, 2, 21);
//         run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 40, 1024, 4294967296, 2, 2, 32);
        run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 8, 2048, 256, 8, 2, 8);
        run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 10, 2048, 512, 8, 2, 9);
        run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 11, 2048, 2048, 6, 2, 11);
        run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 13, 2048, 8192, 5, 2, 13);
        run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 16, 2048, 65536, 4, 2, 16);
        run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 24, 2048, 2097152, 4, 2, 21);
//         run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 41, 2048, 4294967296, 2, 2, 32);
        run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 8, 4096, 256, 8, 2, 8);
        run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 10, 4096, 512, 8, 2, 9);
        run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 11, 4096, 2048, 6, 2, 11);
        run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 13, 4096, 8192, 5, 2, 13);
        run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 17, 4096, 65536, 4, 2, 16);
        run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 25, 4096, 2097152, 4, 2, 21);
//         run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 42, 4096, 4294967296, 2, 2, 32);
        run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 9, 8192, 256, 8, 2, 8);
        run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 11, 8192, 512, 8, 2, 9);
        run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 12, 8192, 2048, 6, 2, 11);
        run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 14, 8192, 8192, 5, 2, 13);
        // run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 17, 8192, 65536, 4, 2, 16);
        // run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 26, 8192, 2097152, 4, 2, 21);
//         // run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 43, 8192, 4294967296, 2, 2, 32);
    };
}

#[macro_export]
macro_rules! run_starkprime_benchmarks {
    ($group: ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
    };
}

#[macro_export]
macro_rules! run_starkprime_non_scalar_benchmarks {
    ($group: ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
    };
}

#[macro_export]
macro_rules! run_starkprime_degree_three_non_scalar_benchmarks {
    ($group: ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
    };
}

#[macro_export]
macro_rules! run_frog_benchmarks {
    ($group:ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
    };
}

#[macro_export]
macro_rules! run_frog_non_scalar_benchmarks {
    ($group:ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
    };
}

#[macro_export]
macro_rules! run_frog_degree_three_non_scalar_benchmarks {
    ($group:ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
    };
}

#[macro_export]
macro_rules! run_babybear_benchmarks {
    ($group:ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
    };
}

#[macro_export]
macro_rules! run_babybear_non_scalar_benchmarks {
    ($group:ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
    };
}

#[macro_export]
macro_rules! run_babybear_degree_three_non_scalar_benchmarks {
    ($group:ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
    };
}

#[macro_export]
macro_rules! run_goldilocks_operations_benchmarks {
    ($group:ident) => {
        run_single_operations_goldilocks_benchmark!(&mut $group, 7, 512, 256, 8, 2, 8);
        run_single_operations_goldilocks_benchmark!(&mut $group, 9, 512, 512, 8, 2, 9);
        run_single_operations_goldilocks_benchmark!(&mut $group, 10, 512, 2048, 6, 2, 11);
        run_single_operations_goldilocks_benchmark!(&mut $group, 12, 512, 8192, 5, 2, 13);
        run_single_operations_goldilocks_benchmark!(&mut $group, 15, 512, 65536, 4, 2, 16);
        run_single_operations_goldilocks_benchmark!(&mut $group, 22, 512, 2097152, 4, 2, 21);
//         run_single_operations_goldilocks_benchmark!(&mut $group, 39, 512, 4294967296, 2, 2, 32);
        run_single_operations_goldilocks_benchmark!(&mut $group, 8, 1024, 256, 8, 2, 8);
        run_single_operations_goldilocks_benchmark!(&mut $group, 9, 1024, 512, 8, 2, 9);
        run_single_operations_goldilocks_benchmark!(&mut $group, 10, 1024, 2048, 6, 2, 11);
        run_single_operations_goldilocks_benchmark!(&mut $group, 12, 1024, 8192, 5, 2, 13);
        run_single_operations_goldilocks_benchmark!(&mut $group, 15, 1024, 65536, 4, 2, 16);
        run_single_operations_goldilocks_benchmark!(&mut $group, 23, 1024, 2097152, 4, 2, 21);
//         run_single_operations_goldilocks_benchmark!(&mut $group, 40, 1024, 4294967296, 2, 2, 32);
        run_single_operations_goldilocks_benchmark!(&mut $group, 8, 2048, 256, 8, 2, 8);
        run_single_operations_goldilocks_benchmark!(&mut $group, 10, 2048, 512, 8, 2, 9);
        run_single_operations_goldilocks_benchmark!(&mut $group, 11, 2048, 2048, 6, 2, 11);
        run_single_operations_goldilocks_benchmark!(&mut $group, 13, 2048, 8192, 5, 2, 13);
        run_single_operations_goldilocks_benchmark!(&mut $group, 16, 2048, 65536, 4, 2, 16);
        run_single_operations_goldilocks_benchmark!(&mut $group, 24, 2048, 2097152, 4, 2, 21);
//         run_single_operations_goldilocks_benchmark!(&mut $group, 41, 2048, 4294967296, 2, 2, 32);
        run_single_operations_goldilocks_benchmark!(&mut $group, 8, 4096, 256, 8, 2, 8);
        run_single_operations_goldilocks_benchmark!(&mut $group, 10, 4096, 512, 8, 2, 9);
        run_single_operations_goldilocks_benchmark!(&mut $group, 11, 4096, 2048, 6, 2, 11);
        run_single_operations_goldilocks_benchmark!(&mut $group, 13, 4096, 8192, 5, 2, 13);
        run_single_operations_goldilocks_benchmark!(&mut $group, 17, 4096, 65536, 4, 2, 16);
        run_single_operations_goldilocks_benchmark!(&mut $group, 25, 4096, 2097152, 4, 2, 21);
//         run_single_operations_goldilocks_benchmark!(&mut $group, 42, 4096, 4294967296, 2, 2, 32);
        run_single_operations_goldilocks_benchmark!(&mut $group, 9, 8192, 256, 8, 2, 8);
        run_single_operations_goldilocks_benchmark!(&mut $group, 11, 8192, 512, 8, 2, 9);
        run_single_operations_goldilocks_benchmark!(&mut $group, 12, 8192, 2048, 6, 2, 11);
        run_single_operations_goldilocks_benchmark!(&mut $group, 14, 8192, 8192, 5, 2, 13);
        // run_single_operations_goldilocks_benchmark!(&mut $group, 17, 8192, 65536, 4, 2, 16);
        // run_single_operations_goldilocks_benchmark!(&mut $group, 26, 8192, 2097152, 4, 2, 21);
//         // run_single_operations_goldilocks_benchmark!(&mut $group, 43, 8192, 4294967296, 2, 2, 32);
    };
}

#[macro_export]
macro_rules! run_goldilocks_operations_non_scalar_benchmarks {
    ($group:ident) => {
        run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 7, 512, 256, 8, 2, 8);
        run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 9, 512, 512, 8, 2, 9);
        run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 10, 512, 2048, 6, 2, 11);
        run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 12, 512, 8192, 5, 2, 13);
        run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 15, 512, 65536, 4, 2, 16);
        run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 22, 512, 2097152, 4, 2, 21);
//         run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 39, 512, 4294967296, 2, 2, 32);
        run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 8, 1024, 256, 8, 2, 8);
        run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 9, 1024, 512, 8, 2, 9);
        run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 10, 1024, 2048, 6, 2, 11);
        run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 12, 1024, 8192, 5, 2, 13);
        run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 15, 1024, 65536, 4, 2, 16);
        run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 23, 1024, 2097152, 4, 2, 21);
//         run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 40, 1024, 4294967296, 2, 2, 32);
        run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 8, 2048, 256, 8, 2, 8);
        run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 10, 2048, 512, 8, 2, 9);
        run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 11, 2048, 2048, 6, 2, 11);
        run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 13, 2048, 8192, 5, 2, 13);
        run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 16, 2048, 65536, 4, 2, 16);
        run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 24, 2048, 2097152, 4, 2, 21);
//         run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 41, 2048, 4294967296, 2, 2, 32);
        run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 8, 4096, 256, 8, 2, 8);
        run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 10, 4096, 512, 8, 2, 9);
        run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 11, 4096, 2048, 6, 2, 11);
        run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 13, 4096, 8192, 5, 2, 13);
        run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 17, 4096, 65536, 4, 2, 16);
        run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 25, 4096, 2097152, 4, 2, 21);
//         run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 42, 4096, 4294967296, 2, 2, 32);
        run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 9, 8192, 256, 8, 2, 8);
        run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 11, 8192, 512, 8, 2, 9);
        run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 12, 8192, 2048, 6, 2, 11);
        run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 14, 8192, 8192, 5, 2, 13);
        // run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 17, 8192, 65536, 4, 2, 16);
        // run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 26, 8192, 2097152, 4, 2, 21);
//         // run_single_operations_non_scalar_goldilocks_benchmark!(&mut $group, 43, 8192, 4294967296, 2, 2, 32);
    };
}

#[macro_export]
macro_rules! run_goldilocks_operations_degree_three_non_scalar_benchmarks {
    ($group:ident) => {
        run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 7, 512, 256, 8, 2, 8);
        run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 9, 512, 512, 8, 2, 9);
        run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 10, 512, 2048, 6, 2, 11);
        run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 12, 512, 8192, 5, 2, 13);
        run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 15, 512, 65536, 4, 2, 16);
        run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 22, 512, 2097152, 4, 2, 21);
//         run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 39, 512, 4294967296, 2, 2, 32);
        run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 8, 1024, 256, 8, 2, 8);
        run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 9, 1024, 512, 8, 2, 9);
        run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 10, 1024, 2048, 6, 2, 11);
        run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 12, 1024, 8192, 5, 2, 13);
        run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 15, 1024, 65536, 4, 2, 16);
        run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 23, 1024, 2097152, 4, 2, 21);
//         run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 40, 1024, 4294967296, 2, 2, 32);
        run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 8, 2048, 256, 8, 2, 8);
        run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 10, 2048, 512, 8, 2, 9);
        run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 11, 2048, 2048, 6, 2, 11);
        run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 13, 2048, 8192, 5, 2, 13);
        run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 16, 2048, 65536, 4, 2, 16);
        run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 24, 2048, 2097152, 4, 2, 21);
//         run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 41, 2048, 4294967296, 2, 2, 32);
        run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 8, 4096, 256, 8, 2, 8);
        run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 10, 4096, 512, 8, 2, 9);
        run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 11, 4096, 2048, 6, 2, 11);
        run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 13, 4096, 8192, 5, 2, 13);
        run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 17, 4096, 65536, 4, 2, 16);
        run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 25, 4096, 2097152, 4, 2, 21);
//         run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 42, 4096, 4294967296, 2, 2, 32);
        run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 9, 8192, 256, 8, 2, 8);
        run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 11, 8192, 512, 8, 2, 9);
        run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 12, 8192, 2048, 6, 2, 11);
        run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 14, 8192, 8192, 5, 2, 13);
        // run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 17, 8192, 65536, 4, 2, 16);
        // run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 26, 8192, 2097152, 4, 2, 21);
//         // run_single_operations_degree_three_non_scalar_goldilocks_benchmark!(&mut $group, 43, 8192, 4294967296, 2, 2, 32);
    };
}