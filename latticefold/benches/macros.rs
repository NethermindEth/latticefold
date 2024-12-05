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
        run_single_goldilocks_benchmark!(&mut $group, 1, 5, 512, 16, 16, 2, 4);
        run_single_goldilocks_benchmark!(&mut $group, 1, 7, 512, 256, 8, 2, 8);
        run_single_goldilocks_benchmark!(&mut $group, 1, 15, 512, 65536, 4, 2, 16);
        run_single_goldilocks_benchmark!(&mut $group, 1, 39, 512, 4294967296, 2, 2, 32);
        run_single_goldilocks_benchmark!(&mut $group, 1, 5, 1024, 16, 16, 2, 4);
        run_single_goldilocks_benchmark!(&mut $group, 1, 8, 1024, 256, 8, 2, 8);
        run_single_goldilocks_benchmark!(&mut $group, 1, 15, 1024, 65536, 4, 2, 16);
        run_single_goldilocks_benchmark!(&mut $group, 1, 40, 1024, 4294967296, 2, 2, 32);
        run_single_goldilocks_benchmark!(&mut $group, 1, 5, 2048, 16, 16, 2, 4);
        run_single_goldilocks_benchmark!(&mut $group, 1, 8, 2048, 256, 8, 2, 8);
        run_single_goldilocks_benchmark!(&mut $group, 1, 16, 2048, 65536, 4, 2, 16);
        run_single_goldilocks_benchmark!(&mut $group, 1, 41, 2048, 4294967296, 2, 2, 32);
        run_single_goldilocks_benchmark!(&mut $group, 1, 6, 4096, 16, 16, 2, 4);
        run_single_goldilocks_benchmark!(&mut $group, 1, 8, 4096, 256, 8, 2, 8);
        run_single_goldilocks_benchmark!(&mut $group, 1, 17, 4096, 65536, 4, 2, 16);
        run_single_goldilocks_benchmark!(&mut $group, 1, 42, 4096, 4294967296, 2, 2, 32);
        run_single_goldilocks_benchmark!(&mut $group, 1, 6, 8192, 16, 16, 2, 4);
        run_single_goldilocks_benchmark!(&mut $group, 1, 9, 8192, 256, 8, 2, 8);
        run_single_goldilocks_benchmark!(&mut $group, 1, 17, 8192, 65536, 4, 2, 16);
        run_single_goldilocks_benchmark!(&mut $group, 1, 43, 8192, 4294967296, 2, 2, 32);
        run_single_goldilocks_benchmark!(&mut $group, 1, 6, 16384, 16, 16, 2, 4);
        run_single_goldilocks_benchmark!(&mut $group, 1, 9, 16384, 256, 8, 2, 8);
        run_single_goldilocks_benchmark!(&mut $group, 1, 18, 16384, 65536, 4, 2, 16);
        run_single_goldilocks_benchmark!(&mut $group, 1, 44, 16384, 4294967296, 2, 2, 32);
    };
}

#[macro_export]
macro_rules! run_goldilocks_non_scalar_benchmarks {
    ($group: ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 5, 512, 16, 16, 2, 4);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 7, 512, 256, 8, 2, 8);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 15, 512, 65536, 4, 2, 16);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 39, 512, 4294967296, 2, 2, 32);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 5, 1024, 16, 16, 2, 4);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 8, 1024, 256, 8, 2, 8);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 15, 1024, 65536, 4, 2, 16);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 40, 1024, 4294967296, 2, 2, 32);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 5, 2048, 16, 16, 2, 4);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 8, 2048, 256, 8, 2, 8);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 16, 2048, 65536, 4, 2, 16);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 41, 2048, 4294967296, 2, 2, 32);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 6, 4096, 16, 16, 2, 4);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 8, 4096, 256, 8, 2, 8);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 17, 4096, 65536, 4, 2, 16);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 42, 4096, 4294967296, 2, 2, 32);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 6, 8192, 16, 16, 2, 4);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 9, 8192, 256, 8, 2, 8);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 17, 8192, 65536, 4, 2, 16);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 43, 8192, 4294967296, 2, 2, 32);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 6, 16384, 16, 16, 2, 4);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 9, 16384, 256, 8, 2, 8);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 18, 16384, 65536, 4, 2, 16);
        run_single_goldilocks_non_scalar_benchmark!(
            &mut $group,
            1,
            44,
            16384,
            4294967296,
            2,
            2,
            32
        );
    };
}

#[macro_export]
macro_rules! run_starkprime_benchmarks {
    ($group: ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        run_single_starkprime_benchmark!(&mut $group, 1, 6, 512, 65536, 16, 2, 16);
        run_single_starkprime_benchmark!(&mut $group, 1, 16, 512, 4294967296, 8, 2, 32);
        run_single_starkprime_benchmark!(&mut $group, 1, 48, 512, 9223372036854775808, 4, 2, 63);
        run_single_starkprime_benchmark!(
            &mut $group,
            1,
            170,
            512,
            85070591730234615865843651857942052864,
            2,
            2,
            126
        );
        run_single_starkprime_benchmark!(&mut $group, 1, 7, 1024, 65536, 16, 2, 16);
        run_single_starkprime_benchmark!(&mut $group, 1, 16, 1024, 4294967296, 8, 2, 32);
        run_single_starkprime_benchmark!(&mut $group, 1, 49, 1024, 9223372036854775808, 4, 2, 63);
        run_single_starkprime_benchmark!(
            &mut $group,
            1,
            172,
            1024,
            85070591730234615865843651857942052864,
            2,
            2,
            126
        );
        run_single_starkprime_benchmark!(&mut $group, 1, 7, 2048, 65536, 16, 2, 16);
        run_single_starkprime_benchmark!(&mut $group, 1, 17, 2048, 4294967296, 8, 2, 32);
        run_single_starkprime_benchmark!(&mut $group, 1, 50, 2048, 9223372036854775808, 4, 2, 63);
        run_single_starkprime_benchmark!(
            &mut $group,
            1,
            173,
            2048,
            85070591730234615865843651857942052864,
            2,
            2,
            126
        );
        run_single_starkprime_benchmark!(&mut $group, 1, 7, 4096, 65536, 16, 2, 16);
        run_single_starkprime_benchmark!(&mut $group, 1, 17, 4096, 4294967296, 8, 2, 32);
        run_single_starkprime_benchmark!(&mut $group, 1, 50, 4096, 9223372036854775808, 4, 2, 63);
        run_single_starkprime_benchmark!(
            &mut $group,
            1,
            174,
            4096,
            85070591730234615865843651857942052864,
            2,
            2,
            126
        );
        run_single_starkprime_benchmark!(&mut $group, 1, 7, 8192, 65536, 16, 2, 16);
        run_single_starkprime_benchmark!(&mut $group, 1, 17, 8192, 4294967296, 8, 2, 32);
        run_single_starkprime_benchmark!(&mut $group, 1, 51, 8192, 9223372036854775808, 4, 2, 63);
        run_single_starkprime_benchmark!(
            &mut $group,
            1,
            176,
            8192,
            85070591730234615865843651857942052864,
            2,
            2,
            126
        );
        run_single_starkprime_benchmark!(&mut $group, 1, 8, 16384, 65536, 16, 2, 16);
        run_single_starkprime_benchmark!(&mut $group, 1, 18, 16384, 4294967296, 8, 2, 32);
        run_single_starkprime_benchmark!(&mut $group, 1, 52, 16384, 9223372036854775808, 4, 2, 63);
        run_single_starkprime_benchmark!(
            &mut $group,
            1,
            177,
            16384,
            85070591730234615865843651857942052864,
            2,
            2,
            126
        );
    };
}

#[macro_export]
macro_rules! run_starkprime_non_scalar_benchmarks {
    ($group: ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 6, 512, 65536, 16, 2, 16);
        run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 16, 512, 4294967296, 8, 2, 32);
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            48,
            512,
            9223372036854775808,
            4,
            2,
            63
        );
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            170,
            512,
            85070591730234615865843651857942052864,
            2,
            2,
            126
        );
        run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 7, 1024, 65536, 16, 2, 16);
        run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 16, 1024, 4294967296, 8, 2, 32);
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            49,
            1024,
            9223372036854775808,
            4,
            2,
            63
        );
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            172,
            1024,
            85070591730234615865843651857942052864,
            2,
            2,
            126
        );
        run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 7, 2048, 65536, 16, 2, 16);
        run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 17, 2048, 4294967296, 8, 2, 32);
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            50,
            2048,
            9223372036854775808,
            4,
            2,
            63
        );
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            173,
            2048,
            85070591730234615865843651857942052864,
            2,
            2,
            126
        );
        run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 7, 4096, 65536, 16, 2, 16);
        run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 17, 4096, 4294967296, 8, 2, 32);
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            50,
            4096,
            9223372036854775808,
            4,
            2,
            63
        );
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            174,
            4096,
            85070591730234615865843651857942052864,
            2,
            2,
            126
        );
        run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 7, 8192, 65536, 16, 2, 16);
        run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 17, 8192, 4294967296, 8, 2, 32);
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            51,
            8192,
            9223372036854775808,
            4,
            2,
            63
        );
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            176,
            8192,
            85070591730234615865843651857942052864,
            2,
            2,
            126
        );
        run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 8, 16384, 65536, 16, 2, 16);
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            18,
            16384,
            4294967296,
            8,
            2,
            32
        );
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            52,
            16384,
            9223372036854775808,
            4,
            2,
            63
        );
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            177,
            16384,
            85070591730234615865843651857942052864,
            2,
            2,
            126
        );
    };
}

#[macro_export]
macro_rules! run_frog_benchmarks {
    ($group:ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        run_single_frog_benchmark!(&mut $group, 1, 7, 512, 32, 14, 2, 5);
        run_single_frog_benchmark!(&mut $group, 1, 11, 512, 512, 8, 2, 9);
        run_single_frog_benchmark!(&mut $group, 1, 22, 512, 131072, 4, 2, 17);
        run_single_frog_benchmark!(&mut $group, 1, 61, 512, 17179869184, 2, 2, 34);
        run_single_frog_benchmark!(&mut $group, 1, 8, 1024, 32, 14, 2, 5);
        run_single_frog_benchmark!(&mut $group, 1, 12, 1024, 512, 8, 2, 9);
        run_single_frog_benchmark!(&mut $group, 1, 23, 1024, 131072, 4, 2, 17);
        run_single_frog_benchmark!(&mut $group, 1, 62, 1024, 17179869184, 2, 2, 34);
        run_single_frog_benchmark!(&mut $group, 1, 8, 2048, 32, 14, 2, 5);
        run_single_frog_benchmark!(&mut $group, 1, 12, 2048, 512, 8, 2, 9);
        run_single_frog_benchmark!(&mut $group, 1, 24, 2048, 131072, 4, 2, 17);
        run_single_frog_benchmark!(&mut $group, 1, 64, 2048, 17179869184, 2, 2, 34);
        run_single_frog_benchmark!(&mut $group, 1, 9, 4096, 32, 14, 2, 5);
        run_single_frog_benchmark!(&mut $group, 1, 13, 4096, 512, 8, 2, 9);
        run_single_frog_benchmark!(&mut $group, 1, 25, 4096, 131072, 4, 2, 17);
        run_single_frog_benchmark!(&mut $group, 1, 65, 4096, 17179869184, 2, 2, 34);
        run_single_frog_benchmark!(&mut $group, 1, 9, 8192, 32, 14, 2, 5);
        run_single_frog_benchmark!(&mut $group, 1, 13, 8192, 512, 8, 2, 9);
        run_single_frog_benchmark!(&mut $group, 1, 26, 8192, 131072, 4, 2, 17);
        run_single_frog_benchmark!(&mut $group, 1, 67, 8192, 17179869184, 2, 2, 34);
        run_single_frog_benchmark!(&mut $group, 1, 10, 16384, 32, 14, 2, 5);
        run_single_frog_benchmark!(&mut $group, 1, 14, 16384, 512, 8, 2, 9);
        run_single_frog_benchmark!(&mut $group, 1, 27, 16384, 131072, 4, 2, 17);
        run_single_frog_benchmark!(&mut $group, 1, 69, 16384, 17179869184, 2, 2, 34);
    };
}

#[macro_export]
macro_rules! run_frog_non_scalar_benchmarks {
    ($group:ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 7, 512, 32, 14, 2, 5);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 11, 512, 512, 8, 2, 9);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 22, 512, 131072, 4, 2, 17);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 61, 512, 17179869184, 2, 2, 34);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 8, 1024, 32, 14, 2, 5);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 12, 1024, 512, 8, 2, 9);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 23, 1024, 131072, 4, 2, 17);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 62, 1024, 17179869184, 2, 2, 34);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 8, 2048, 32, 14, 2, 5);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 12, 2048, 512, 8, 2, 9);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 24, 2048, 131072, 4, 2, 17);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 64, 2048, 17179869184, 2, 2, 34);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 9, 4096, 32, 14, 2, 5);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 13, 4096, 512, 8, 2, 9);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 25, 4096, 131072, 4, 2, 17);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 65, 4096, 17179869184, 2, 2, 34);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 9, 8192, 32, 14, 2, 5);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 13, 8192, 512, 8, 2, 9);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 26, 8192, 131072, 4, 2, 17);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 67, 8192, 17179869184, 2, 2, 34);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 10, 16384, 32, 14, 2, 5);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 14, 16384, 512, 8, 2, 9);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 27, 16384, 131072, 4, 2, 17);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 69, 16384, 17179869184, 2, 2, 34);
    };
}

#[macro_export]
macro_rules! run_babybear_benchmarks {
    ($group:ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        run_single_babybear_benchmark!(&mut $group, 1, 3, 512, 16, 8, 2, 4);
        run_single_babybear_benchmark!(&mut $group, 1, 5, 512, 256, 4, 2, 8);
        run_single_babybear_benchmark!(&mut $group, 1, 10, 512, 32768, 2, 2, 15);
        run_single_babybear_benchmark!(&mut $group, 1, 3, 1024, 4, 15, 2, 2);
        run_single_babybear_benchmark!(&mut $group, 1, 4, 1024, 16, 8, 2, 4);
        run_single_babybear_benchmark!(&mut $group, 1, 6, 1024, 256, 4, 2, 8);
        run_single_babybear_benchmark!(&mut $group, 1, 10, 1024, 32768, 2, 2, 15);
        run_single_babybear_benchmark!(&mut $group, 1, 3, 2048, 4, 15, 2, 2);
        run_single_babybear_benchmark!(&mut $group, 1, 4, 2048, 16, 8, 2, 4);
        run_single_babybear_benchmark!(&mut $group, 1, 6, 2048, 256, 4, 2, 8);
        run_single_babybear_benchmark!(&mut $group, 1, 11, 2048, 32768, 2, 2, 15);
        run_single_babybear_benchmark!(&mut $group, 1, 3, 4096, 4, 15, 2, 2);
        run_single_babybear_benchmark!(&mut $group, 1, 4, 4096, 16, 8, 2, 4);
        run_single_babybear_benchmark!(&mut $group, 1, 6, 4096, 256, 4, 2, 8);
        run_single_babybear_benchmark!(&mut $group, 1, 11, 4096, 32768, 2, 2, 15);
        run_single_babybear_benchmark!(&mut $group, 1, 4, 8192, 16, 8, 2, 4);
        run_single_babybear_benchmark!(&mut $group, 1, 6, 8192, 256, 4, 2, 8);
        run_single_babybear_benchmark!(&mut $group, 1, 11, 8192, 32768, 2, 2, 15);
        run_single_babybear_benchmark!(&mut $group, 1, 4, 16384, 4, 15, 2, 2);
        run_single_babybear_benchmark!(&mut $group, 1, 5, 16384, 16, 8, 2, 4);
        run_single_babybear_benchmark!(&mut $group, 1, 7, 16384, 256, 4, 2, 8);
        run_single_babybear_benchmark!(&mut $group, 1, 12, 16384, 32768, 2, 2, 15);
    };
}
#[macro_export]
macro_rules! run_babybear_non_scalar_benchmarks {
    ($group:ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 3, 512, 16, 8, 2, 4);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 5, 512, 256, 4, 2, 8);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 10, 512, 32768, 2, 2, 15);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 3, 1024, 4, 15, 2, 2);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 4, 1024, 16, 8, 2, 4);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 6, 1024, 256, 4, 2, 8);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 10, 1024, 32768, 2, 2, 15);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 3, 2048, 4, 15, 2, 2);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 4, 2048, 16, 8, 2, 4);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 6, 2048, 256, 4, 2, 8);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 11, 2048, 32768, 2, 2, 15);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 3, 4096, 4, 15, 2, 2);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 4, 4096, 16, 8, 2, 4);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 6, 4096, 256, 4, 2, 8);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 11, 4096, 32768, 2, 2, 15);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 4, 8192, 16, 8, 2, 4);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 6, 8192, 256, 4, 2, 8);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 11, 8192, 32768, 2, 2, 15);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 4, 16384, 4, 15, 2, 2);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 5, 16384, 16, 8, 2, 4);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 7, 16384, 256, 4, 2, 8);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 12, 16384, 32768, 2, 2, 15);
    };
}
