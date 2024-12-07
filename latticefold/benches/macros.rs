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
macro_rules! run_goldilocks_benchmarks {
    ($group: ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        run_single_goldilocks_benchmark!(&mut $group, 1, 6, 512, 128, 10, 2, 7);
        run_single_goldilocks_benchmark!(&mut $group, 1, 7, 512, 256, 8, 2, 8);
        run_single_goldilocks_benchmark!(&mut $group, 1, 8, 512, 512, 8, 2, 9);
        run_single_goldilocks_benchmark!(&mut $group, 1, 8, 1024, 512, 8, 2, 9);
        run_single_goldilocks_benchmark!(&mut $group, 1, 8, 2048, 256, 8, 2, 8);
        run_single_goldilocks_benchmark!(&mut $group, 1, 9, 1024, 1024, 7, 2, 10);
        run_single_goldilocks_benchmark!(&mut $group, 1, 9, 2048, 512, 8, 2, 9);
        run_single_goldilocks_benchmark!(&mut $group, 1, 10, 512, 2048, 6, 2, 11);
        run_single_goldilocks_benchmark!(&mut $group, 1, 10, 1024, 2048, 6, 2, 11);
        run_single_goldilocks_benchmark!(&mut $group, 1, 11, 1024, 4096, 6, 2, 12);
        run_single_goldilocks_benchmark!(&mut $group, 1, 11, 2048, 2048, 6, 2, 11);
        run_single_goldilocks_benchmark!(&mut $group, 1, 12, 1024, 8192, 5, 2, 13);
        run_single_goldilocks_benchmark!(&mut $group, 1, 13, 1024, 16384, 5, 2, 14);
        run_single_goldilocks_benchmark!(&mut $group, 1, 13, 2048, 8192, 5, 2, 13);
        run_single_goldilocks_benchmark!(&mut $group, 1, 14, 1024, 32768, 5, 2, 15);
        run_single_goldilocks_benchmark!(&mut $group, 1, 14, 2048, 16384, 5, 2, 14);
        run_single_goldilocks_benchmark!(&mut $group, 1, 15, 2048, 32768, 5, 2, 15);
        run_single_goldilocks_benchmark!(&mut $group, 1, 16, 2048, 65536, 4, 2, 16);
    };
}

#[macro_export]
macro_rules! run_goldilocks_non_scalar_benchmarks {
    ($group: ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 6, 512, 128, 10, 2, 7);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 7, 512, 256, 8, 2, 8);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 8, 512, 512, 8, 2, 9);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 8, 1024, 512, 8, 2, 9);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 8, 2048, 256, 8, 2, 8);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 9, 1024, 1024, 7, 2, 10);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 9, 2048, 512, 8, 2, 9);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 10, 512, 2048, 6, 2, 11);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 10, 1024, 2048, 6, 2, 11);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 11, 1024, 4096, 6, 2, 12);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 11, 2048, 2048, 6, 2, 11);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 12, 1024, 8192, 5, 2, 13);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 13, 1024, 16384, 5, 2, 14);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 13, 2048, 8192, 5, 2, 13);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 14, 1024, 32768, 5, 2, 15);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 14, 2048, 16384, 5, 2, 14);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 15, 2048, 32768, 5, 2, 15);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 16, 2048, 65536, 4, 2, 16);
    };
}

#[macro_export]
macro_rules! run_starkprime_benchmarks {
    ($group: ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        run_single_starkprime_benchmark!(&mut $group, 1, 15, 1024, 1073741824u128, 9, 2, 30);
        run_single_starkprime_benchmark!(&mut $group, 1, 16, 1024, 4294967296u128, 8, 2, 32);
        run_single_starkprime_benchmark!(&mut $group, 1, 17, 2048, 8589934592u128, 8, 2, 33);
        run_single_starkprime_benchmark!(&mut $group, 1, 18, 2048, 17179869184u128, 8, 2, 34);
        run_single_starkprime_benchmark!(&mut $group, 1, 19, 2048, 34359738368u128, 8, 2, 35);
    };
}

#[macro_export]
macro_rules! run_starkprime_non_scalar_benchmarks {
    ($group: ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            15,
            1024,
            1073741824u128,
            9,
            2,
            30
        );
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            16,
            1024,
            4294967296u128,
            8,
            2,
            32
        );
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            17,
            2048,
            8589934592u128,
            8,
            2,
            33
        );
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            18,
            2048,
            17179869184u128,
            8,
            2,
            34
        );
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            19,
            2048,
            34359738368u128,
            8,
            2,
            35
        );
    };
}

#[macro_export]
macro_rules! run_frog_benchmarks {
    ($group:ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        run_single_frog_benchmark!(&mut $group, 1, 5, 512, 8, 23, 2, 3);
        run_single_frog_benchmark!(&mut $group, 1, 9, 1024, 128, 10, 2, 7);
        run_single_frog_benchmark!(&mut $group, 1, 10, 1024, 256, 9, 2, 8);
        run_single_frog_benchmark!(&mut $group, 1, 12, 512, 1024, 7, 2, 10);
        run_single_frog_benchmark!(&mut $group, 1, 15, 1024, 4096, 6, 2, 12);
    };
}

#[macro_export]
macro_rules! run_frog_non_scalar_benchmarks {
    ($group:ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 5, 512, 8, 23, 2, 3);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 9, 1024, 128, 10, 2, 7);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 10, 1024, 256, 9, 2, 8);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 12, 512, 1024, 7, 2, 10);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 15, 1024, 4096, 6, 2, 12);
    };
}

#[macro_export]
macro_rules! run_babybear_benchmarks {
    ($group:ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        run_single_babybear_benchmark!(&mut $group, 1, 6, 1024, 512, 4, 2, 9);
        run_single_babybear_benchmark!(&mut $group, 1, 7, 1024, 2048, 3, 2, 11);
        run_single_babybear_benchmark!(&mut $group, 1, 8, 4096, 2048, 3, 2, 11);
        run_single_babybear_benchmark!(&mut $group, 1, 9, 2048, 8192, 3, 2, 13);
        run_single_babybear_benchmark!(&mut $group, 1, 10, 4096, 16384, 3, 2, 14);
    };
}
#[macro_export]
macro_rules! run_babybear_non_scalar_benchmarks {
    ($group:ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 6, 1024, 512, 4, 2, 9);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 7, 1024, 2048, 3, 2, 11);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 8, 4096, 2048, 3, 2, 11);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 9, 2048, 8192, 3, 2, 13);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 10, 4096, 16384, 3, 2, 14);
    };
}

#[macro_export]
macro_rules! run_goldilocks_linearization_benchmarks {
    ($group:ident) => {
        run_single_linearization_goldilocks_benchmark!(&mut $group, 9, 512, 512, 7, 2, 9);
        run_single_linearization_goldilocks_benchmark!(&mut $group, 12, 512, 8192, 5, 2, 13);
        run_single_linearization_goldilocks_benchmark!(&mut $group, 15, 512, 65536, 4, 2, 16);
        run_single_linearization_goldilocks_benchmark!(&mut $group, 22, 512, 2097152, 3, 2, 21);
        run_single_linearization_goldilocks_benchmark!(&mut $group, 39, 512, 4294967296, 2, 2, 32);
        run_single_linearization_goldilocks_benchmark!(&mut $group, 8, 1024, 256, 8, 2, 8);
        run_single_linearization_goldilocks_benchmark!(&mut $group, 9, 1024, 512, 7, 2, 9);
        run_single_linearization_goldilocks_benchmark!(&mut $group, 15, 1024, 65536, 4, 2, 16);
        run_single_linearization_goldilocks_benchmark!(&mut $group, 23, 1024, 2097152, 3, 2, 21);
        run_single_linearization_goldilocks_benchmark!(&mut $group, 40, 1024, 4294967296, 2, 2, 32);
        run_single_linearization_goldilocks_benchmark!(&mut $group, 8, 2048, 256, 8, 2, 8);
        run_single_linearization_goldilocks_benchmark!(&mut $group, 10, 2048, 512, 7, 2, 9);
        run_single_linearization_goldilocks_benchmark!(&mut $group, 11, 2048, 2048, 6, 2, 11);
        run_single_linearization_goldilocks_benchmark!(&mut $group, 24, 2048, 2097152, 3, 2, 21);
        run_single_linearization_goldilocks_benchmark!(&mut $group, 41, 2048, 4294967296, 2, 2, 32);
        run_single_linearization_goldilocks_benchmark!(&mut $group, 8, 4096, 256, 8, 2, 8);
        run_single_linearization_goldilocks_benchmark!(&mut $group, 10, 4096, 512, 7, 2, 9);
        run_single_linearization_goldilocks_benchmark!(&mut $group, 11, 4096, 2048, 6, 2, 11);
        run_single_linearization_goldilocks_benchmark!(&mut $group, 13, 4096, 8192, 5, 2, 13);
        run_single_linearization_goldilocks_benchmark!(&mut $group, 17, 4096, 65536, 4, 2, 16);
        run_single_linearization_goldilocks_benchmark!(&mut $group, 25, 4096, 2097152, 3, 2, 21);
        run_single_linearization_goldilocks_benchmark!(&mut $group, 42, 4096, 4294967296, 2, 2, 32);
        run_single_linearization_goldilocks_benchmark!(&mut $group, 9, 8192, 256, 8, 2, 8);
        run_single_linearization_goldilocks_benchmark!(&mut $group, 11, 8192, 512, 7, 2, 9);
        run_single_linearization_goldilocks_benchmark!(&mut $group, 12, 8192, 2048, 6, 2, 11);
        run_single_linearization_goldilocks_benchmark!(&mut $group, 14, 8192, 8192, 5, 2, 13);
        run_single_linearization_goldilocks_benchmark!(&mut $group, 17, 8192, 65536, 4, 2, 16);
        run_single_linearization_goldilocks_benchmark!(&mut $group, 26, 8192, 2097152, 3, 2, 21);
        run_single_linearization_goldilocks_benchmark!(&mut $group, 43, 8192, 4294967296, 2, 2, 32);
        run_single_linearization_goldilocks_benchmark!(&mut $group, 9, 16384, 256, 8, 2, 8);
        run_single_linearization_goldilocks_benchmark!(&mut $group, 11, 16384, 512, 7, 2, 9);
    };
}

#[macro_export]
macro_rules! run_goldilocks_decomposition_benchmarks {
    ($group:ident) => {
        run_single_decomposition_goldilocks_benchmark!(&mut $group, 9, 512, 512, 7, 2, 9);
        run_single_decomposition_goldilocks_benchmark!(&mut $group, 12, 512, 8192, 5, 2, 13);
        run_single_decomposition_goldilocks_benchmark!(&mut $group, 15, 512, 65536, 4, 2, 16);
        run_single_decomposition_goldilocks_benchmark!(&mut $group, 22, 512, 2097152, 3, 2, 21);
        run_single_decomposition_goldilocks_benchmark!(&mut $group, 39, 512, 4294967296, 2, 2, 32);
        run_single_decomposition_goldilocks_benchmark!(&mut $group, 8, 1024, 256, 8, 2, 8);
        run_single_decomposition_goldilocks_benchmark!(&mut $group, 9, 1024, 512, 7, 2, 9);
        run_single_decomposition_goldilocks_benchmark!(&mut $group, 15, 1024, 65536, 4, 2, 16);
        run_single_decomposition_goldilocks_benchmark!(&mut $group, 23, 1024, 2097152, 3, 2, 21);
        run_single_decomposition_goldilocks_benchmark!(&mut $group, 40, 1024, 4294967296, 2, 2, 32);
        run_single_decomposition_goldilocks_benchmark!(&mut $group, 8, 2048, 256, 8, 2, 8);
        run_single_decomposition_goldilocks_benchmark!(&mut $group, 10, 2048, 512, 7, 2, 9);
        run_single_decomposition_goldilocks_benchmark!(&mut $group, 11, 2048, 2048, 6, 2, 11);
        run_single_decomposition_goldilocks_benchmark!(&mut $group, 24, 2048, 2097152, 3, 2, 21);
        run_single_decomposition_goldilocks_benchmark!(&mut $group, 41, 2048, 4294967296, 2, 2, 32);
        run_single_decomposition_goldilocks_benchmark!(&mut $group, 8, 4096, 256, 8, 2, 8);
        run_single_decomposition_goldilocks_benchmark!(&mut $group, 10, 4096, 512, 7, 2, 9);
        run_single_decomposition_goldilocks_benchmark!(&mut $group, 11, 4096, 2048, 6, 2, 11);
        run_single_decomposition_goldilocks_benchmark!(&mut $group, 13, 4096, 8192, 5, 2, 13);
        run_single_decomposition_goldilocks_benchmark!(&mut $group, 17, 4096, 65536, 4, 2, 16);
        run_single_decomposition_goldilocks_benchmark!(&mut $group, 25, 4096, 2097152, 3, 2, 21);
        run_single_decomposition_goldilocks_benchmark!(&mut $group, 42, 4096, 4294967296, 2, 2, 32);
        run_single_decomposition_goldilocks_benchmark!(&mut $group, 9, 8192, 256, 8, 2, 8);
        run_single_decomposition_goldilocks_benchmark!(&mut $group, 11, 8192, 512, 7, 2, 9);
        run_single_decomposition_goldilocks_benchmark!(&mut $group, 12, 8192, 2048, 6, 2, 11);
        run_single_decomposition_goldilocks_benchmark!(&mut $group, 14, 8192, 8192, 5, 2, 13);
        run_single_decomposition_goldilocks_benchmark!(&mut $group, 17, 8192, 65536, 4, 2, 16);
        run_single_decomposition_goldilocks_benchmark!(&mut $group, 26, 8192, 2097152, 3, 2, 21);
        run_single_decomposition_goldilocks_benchmark!(&mut $group, 43, 8192, 4294967296, 2, 2, 32);
        run_single_decomposition_goldilocks_benchmark!(&mut $group, 9, 16384, 256, 8, 2, 8);
        run_single_decomposition_goldilocks_benchmark!(&mut $group, 11, 16384, 512, 7, 2, 9);
    };
}

#[macro_export]
macro_rules! run_goldilocks_folding_benchmarks {
    ($group:ident) => {
        run_single_folding_goldilocks_benchmark!(&mut $group, 9, 512, 512, 7, 2, 9);
        run_single_folding_goldilocks_benchmark!(&mut $group, 12, 512, 8192, 5, 2, 13);
        run_single_folding_goldilocks_benchmark!(&mut $group, 15, 512, 65536, 4, 2, 16);
        run_single_folding_goldilocks_benchmark!(&mut $group, 22, 512, 2097152, 3, 2, 21);
        run_single_folding_goldilocks_benchmark!(&mut $group, 39, 512, 4294967296, 2, 2, 32);
        run_single_folding_goldilocks_benchmark!(&mut $group, 8, 1024, 256, 8, 2, 8);
        run_single_folding_goldilocks_benchmark!(&mut $group, 9, 1024, 512, 7, 2, 9);
        run_single_folding_goldilocks_benchmark!(&mut $group, 15, 1024, 65536, 4, 2, 16);
        run_single_folding_goldilocks_benchmark!(&mut $group, 23, 1024, 2097152, 3, 2, 21);
        run_single_folding_goldilocks_benchmark!(&mut $group, 40, 1024, 4294967296, 2, 2, 32);
        run_single_folding_goldilocks_benchmark!(&mut $group, 8, 2048, 256, 8, 2, 8);
        run_single_folding_goldilocks_benchmark!(&mut $group, 10, 2048, 512, 7, 2, 9);
        run_single_folding_goldilocks_benchmark!(&mut $group, 11, 2048, 2048, 6, 2, 11);
        run_single_folding_goldilocks_benchmark!(&mut $group, 24, 2048, 2097152, 3, 2, 21);
        run_single_folding_goldilocks_benchmark!(&mut $group, 41, 2048, 4294967296, 2, 2, 32);
        run_single_folding_goldilocks_benchmark!(&mut $group, 8, 4096, 256, 8, 2, 8);
        run_single_folding_goldilocks_benchmark!(&mut $group, 10, 4096, 512, 7, 2, 9);
        run_single_folding_goldilocks_benchmark!(&mut $group, 11, 4096, 2048, 6, 2, 11);
        run_single_folding_goldilocks_benchmark!(&mut $group, 13, 4096, 8192, 5, 2, 13);
        run_single_folding_goldilocks_benchmark!(&mut $group, 17, 4096, 65536, 4, 2, 16);
        run_single_folding_goldilocks_benchmark!(&mut $group, 25, 4096, 2097152, 3, 2, 21);
        run_single_folding_goldilocks_benchmark!(&mut $group, 42, 4096, 4294967296, 2, 2, 32);
        run_single_folding_goldilocks_benchmark!(&mut $group, 9, 8192, 256, 8, 2, 8);
        run_single_folding_goldilocks_benchmark!(&mut $group, 11, 8192, 512, 7, 2, 9);
        run_single_folding_goldilocks_benchmark!(&mut $group, 12, 8192, 2048, 6, 2, 11);
        run_single_folding_goldilocks_benchmark!(&mut $group, 14, 8192, 8192, 5, 2, 13);
        run_single_folding_goldilocks_benchmark!(&mut $group, 17, 8192, 65536, 4, 2, 16);
        run_single_folding_goldilocks_benchmark!(&mut $group, 26, 8192, 2097152, 3, 2, 21);
        run_single_folding_goldilocks_benchmark!(&mut $group, 43, 8192, 4294967296, 2, 2, 32);
        run_single_folding_goldilocks_benchmark!(&mut $group, 9, 16384, 256, 8, 2, 8);
        run_single_folding_goldilocks_benchmark!(&mut $group, 11, 16384, 512, 7, 2, 9);
    };
}
