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
        run_single_goldilocks_benchmark!(&mut $group, 1, 7, 512, 256, 8, 2, 8);
        run_single_goldilocks_benchmark!(&mut $group, 1, 9, 512, 512, 8, 2, 9);
        run_single_goldilocks_benchmark!(&mut $group, 1, 10, 512, 2048, 6, 2, 11);
        run_single_goldilocks_benchmark!(&mut $group, 1, 12, 512, 8192, 5, 2, 13);
        run_single_goldilocks_benchmark!(&mut $group, 1, 15, 512, 65536, 4, 2, 16);
        run_single_goldilocks_benchmark!(&mut $group, 1, 22, 512, 2097152, 4, 2, 21);
        run_single_goldilocks_benchmark!(&mut $group, 1, 39, 512, 4294967296, 2, 2, 32);
        run_single_goldilocks_benchmark!(&mut $group, 1, 8, 1024, 256, 8, 2, 8);
        run_single_goldilocks_benchmark!(&mut $group, 1, 9, 1024, 512, 8, 2, 9);
        run_single_goldilocks_benchmark!(&mut $group, 1, 10, 1024, 2048, 6, 2, 11);
        run_single_goldilocks_benchmark!(&mut $group, 1, 12, 1024, 8192, 5, 2, 13);
        run_single_goldilocks_benchmark!(&mut $group, 1, 15, 1024, 65536, 4, 2, 16);
        run_single_goldilocks_benchmark!(&mut $group, 1, 23, 1024, 2097152, 4, 2, 21);
        run_single_goldilocks_benchmark!(&mut $group, 1, 40, 1024, 4294967296, 2, 2, 32);
        run_single_goldilocks_benchmark!(&mut $group, 1, 8, 2048, 256, 8, 2, 8);
        run_single_goldilocks_benchmark!(&mut $group, 1, 10, 2048, 512, 8, 2, 9);
        run_single_goldilocks_benchmark!(&mut $group, 1, 11, 2048, 2048, 6, 2, 11);
        run_single_goldilocks_benchmark!(&mut $group, 1, 13, 2048, 8192, 5, 2, 13);
        run_single_goldilocks_benchmark!(&mut $group, 1, 16, 2048, 65536, 4, 2, 16);
        run_single_goldilocks_benchmark!(&mut $group, 1, 24, 2048, 2097152, 4, 2, 21);
        run_single_goldilocks_benchmark!(&mut $group, 1, 41, 2048, 4294967296, 2, 2, 32);
        run_single_goldilocks_benchmark!(&mut $group, 1, 8, 4096, 256, 8, 2, 8);
        run_single_goldilocks_benchmark!(&mut $group, 1, 10, 4096, 512, 8, 2, 9);
        run_single_goldilocks_benchmark!(&mut $group, 1, 11, 4096, 2048, 6, 2, 11);
        run_single_goldilocks_benchmark!(&mut $group, 1, 13, 4096, 8192, 5, 2, 13);
        run_single_goldilocks_benchmark!(&mut $group, 1, 17, 4096, 65536, 4, 2, 16);
        run_single_goldilocks_benchmark!(&mut $group, 1, 25, 4096, 2097152, 4, 2, 21);
        run_single_goldilocks_benchmark!(&mut $group, 1, 42, 4096, 4294967296, 2, 2, 32);
        run_single_goldilocks_benchmark!(&mut $group, 1, 9, 8192, 256, 8, 2, 8);
        run_single_goldilocks_benchmark!(&mut $group, 1, 11, 8192, 512, 8, 2, 9);
        run_single_goldilocks_benchmark!(&mut $group, 1, 12, 8192, 2048, 6, 2, 11);
        run_single_goldilocks_benchmark!(&mut $group, 1, 14, 8192, 8192, 5, 2, 13);
        run_single_goldilocks_benchmark!(&mut $group, 1, 17, 8192, 65536, 4, 2, 16);
        run_single_goldilocks_benchmark!(&mut $group, 1, 26, 8192, 2097152, 4, 2, 21);
        run_single_goldilocks_benchmark!(&mut $group, 1, 43, 8192, 4294967296, 2, 2, 32);
        run_single_goldilocks_benchmark!(&mut $group, 1, 9, 16384, 256, 8, 2, 8);
        run_single_goldilocks_benchmark!(&mut $group, 1, 11, 16384, 512, 8, 2, 9);
        run_single_goldilocks_benchmark!(&mut $group, 1, 12, 16384, 2048, 6, 2, 11);
        run_single_goldilocks_benchmark!(&mut $group, 1, 14, 16384, 8192, 5, 2, 13);
        run_single_goldilocks_benchmark!(&mut $group, 1, 18, 16384, 65536, 4, 2, 16);
        run_single_goldilocks_benchmark!(&mut $group, 1, 26, 16384, 2097152, 4, 2, 21);
        run_single_goldilocks_benchmark!(&mut $group, 1, 44, 16384, 4294967296, 2, 2, 32);
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
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 39, 512, 4294967296, 2, 2, 32);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 8, 1024, 256, 8, 2, 8);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 9, 1024, 512, 8, 2, 9);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 10, 1024, 2048, 6, 2, 11);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 12, 1024, 8192, 5, 2, 13);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 15, 1024, 65536, 4, 2, 16);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 23, 1024, 2097152, 4, 2, 21);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 40, 1024, 4294967296, 2, 2, 32);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 8, 2048, 256, 8, 2, 8);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 10, 2048, 512, 8, 2, 9);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 11, 2048, 2048, 6, 2, 11);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 13, 2048, 8192, 5, 2, 13);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 16, 2048, 65536, 4, 2, 16);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 24, 2048, 2097152, 4, 2, 21);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 41, 2048, 4294967296, 2, 2, 32);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 8, 4096, 256, 8, 2, 8);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 10, 4096, 512, 8, 2, 9);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 11, 4096, 2048, 6, 2, 11);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 13, 4096, 8192, 5, 2, 13);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 17, 4096, 65536, 4, 2, 16);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 25, 4096, 2097152, 4, 2, 21);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 42, 4096, 4294967296, 2, 2, 32);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 9, 8192, 256, 8, 2, 8);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 11, 8192, 512, 8, 2, 9);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 12, 8192, 2048, 6, 2, 11);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 14, 8192, 8192, 5, 2, 13);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 17, 8192, 65536, 4, 2, 16);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 26, 8192, 2097152, 4, 2, 21);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 43, 8192, 4294967296, 2, 2, 32);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 9, 16384, 256, 8, 2, 8);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 11, 16384, 512, 8, 2, 9);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 12, 16384, 2048, 6, 2, 11);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 14, 16384, 8192, 5, 2, 13);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 18, 16384, 65536, 4, 2, 16);
        run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 26, 16384, 2097152, 4, 2, 21);
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
macro_rules! run_goldilocks_degree_three_non_scalar_benchmarks {
    ($group: ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        run_single_goldilocks_degree_three_non_scalar_benchmark!(
            &mut $group,
            1,
            6,
            512,
            128,
            10,
            2,
            7
        );
        run_single_goldilocks_degree_three_non_scalar_benchmark!(
            &mut $group,
            1,
            7,
            512,
            256,
            8,
            2,
            8
        );
        run_single_goldilocks_degree_three_non_scalar_benchmark!(
            &mut $group,
            1,
            8,
            512,
            512,
            8,
            2,
            9
        );
        run_single_goldilocks_degree_three_non_scalar_benchmark!(
            &mut $group,
            1,
            8,
            1024,
            512,
            8,
            2,
            9
        );
        run_single_goldilocks_degree_three_non_scalar_benchmark!(
            &mut $group,
            1,
            8,
            2048,
            256,
            8,
            2,
            8
        );
        run_single_goldilocks_degree_three_non_scalar_benchmark!(
            &mut $group,
            1,
            9,
            1024,
            1024,
            7,
            2,
            10
        );
        run_single_goldilocks_degree_three_non_scalar_benchmark!(
            &mut $group,
            1,
            9,
            2048,
            512,
            8,
            2,
            9
        );
        run_single_goldilocks_degree_three_non_scalar_benchmark!(
            &mut $group,
            1,
            10,
            512,
            2048,
            6,
            2,
            11
        );
        run_single_goldilocks_degree_three_non_scalar_benchmark!(
            &mut $group,
            1,
            10,
            1024,
            2048,
            6,
            2,
            11
        );
        run_single_goldilocks_degree_three_non_scalar_benchmark!(
            &mut $group,
            1,
            11,
            1024,
            4096,
            6,
            2,
            12
        );
        run_single_goldilocks_degree_three_non_scalar_benchmark!(
            &mut $group,
            1,
            11,
            2048,
            2048,
            6,
            2,
            11
        );
        run_single_goldilocks_degree_three_non_scalar_benchmark!(
            &mut $group,
            1,
            12,
            1024,
            8192,
            5,
            2,
            13
        );
        run_single_goldilocks_degree_three_non_scalar_benchmark!(
            &mut $group,
            1,
            13,
            1024,
            16384,
            5,
            2,
            14
        );
        run_single_goldilocks_degree_three_non_scalar_benchmark!(
            &mut $group,
            1,
            13,
            2048,
            8192,
            5,
            2,
            13
        );
        run_single_goldilocks_degree_three_non_scalar_benchmark!(
            &mut $group,
            1,
            14,
            1024,
            32768,
            5,
            2,
            15
        );
        run_single_goldilocks_degree_three_non_scalar_benchmark!(
            &mut $group,
            1,
            14,
            2048,
            16384,
            5,
            2,
            14
        );
        run_single_goldilocks_degree_three_non_scalar_benchmark!(
            &mut $group,
            1,
            15,
            2048,
            32768,
            5,
            2,
            15
        );
        run_single_goldilocks_degree_three_non_scalar_benchmark!(
            &mut $group,
            1,
            16,
            2048,
            65536,
            4,
            2,
            16
        );
    };
}

#[macro_export]
macro_rules! run_starkprime_benchmarks {
    ($group: ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        run_single_starkprime_benchmark!(&mut $group, 1, 13, 512, 268435456, 9, 2, 28);
        run_single_starkprime_benchmark!(&mut $group, 1, 16, 512, 4294967296, 8, 2, 32);
        run_single_starkprime_benchmark!(&mut $group, 1, 19, 512, 68719476736, 7, 2, 36);
        run_single_starkprime_benchmark!(&mut $group, 1, 24, 512, 4398046511104, 6, 2, 42);
        run_single_starkprime_benchmark!(&mut $group, 1, 34, 512, 2251799813685248, 5, 2, 51);
        run_single_starkprime_benchmark!(&mut $group, 1, 48, 512, 9223372036854775808, 4, 2, 63);
        run_single_starkprime_benchmark!(
            &mut $group,
            1,
            81,
            512,
            19342813113834066795298816,
            3,
            2,
            84
        );
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
        run_single_starkprime_benchmark!(&mut $group, 1, 13, 1024, 268435456, 9, 2, 28);
        run_single_starkprime_benchmark!(&mut $group, 1, 16, 1024, 4294967296, 8, 2, 32);
        run_single_starkprime_benchmark!(&mut $group, 1, 19, 1024, 68719476736, 7, 2, 36);
        run_single_starkprime_benchmark!(&mut $group, 1, 25, 1024, 4398046511104, 6, 2, 42);
        run_single_starkprime_benchmark!(&mut $group, 1, 34, 1024, 2251799813685248, 5, 2, 51);
        run_single_starkprime_benchmark!(&mut $group, 1, 49, 1024, 9223372036854775808, 4, 2, 63);
        run_single_starkprime_benchmark!(
            &mut $group,
            1,
            81,
            1024,
            19342813113834066795298816,
            3,
            2,
            84
        );
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
        run_single_starkprime_benchmark!(&mut $group, 1, 14, 2048, 268435456, 9, 2, 28);
        run_single_starkprime_benchmark!(&mut $group, 1, 17, 2048, 4294967296, 8, 2, 32);
        run_single_starkprime_benchmark!(&mut $group, 1, 20, 2048, 68719476736, 7, 2, 36);
        run_single_starkprime_benchmark!(&mut $group, 1, 25, 2048, 4398046511104, 6, 2, 42);
        run_single_starkprime_benchmark!(&mut $group, 1, 35, 2048, 2251799813685248, 5, 2, 51);
        run_single_starkprime_benchmark!(&mut $group, 1, 50, 2048, 9223372036854775808, 4, 2, 63);
        run_single_starkprime_benchmark!(
            &mut $group,
            1,
            82,
            2048,
            19342813113834066795298816,
            3,
            2,
            84
        );
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
        run_single_starkprime_benchmark!(&mut $group, 1, 14, 4096, 268435456, 9, 2, 28);
        run_single_starkprime_benchmark!(&mut $group, 1, 17, 4096, 4294967296, 8, 2, 32);
        run_single_starkprime_benchmark!(&mut $group, 1, 20, 4096, 68719476736, 7, 2, 36);
        run_single_starkprime_benchmark!(&mut $group, 1, 26, 4096, 4398046511104, 6, 2, 42);
        run_single_starkprime_benchmark!(&mut $group, 1, 35, 4096, 2251799813685248, 5, 2, 51);
        run_single_starkprime_benchmark!(&mut $group, 1, 50, 4096, 9223372036854775808, 4, 2, 63);
        run_single_starkprime_benchmark!(
            &mut $group,
            1,
            83,
            4096,
            19342813113834066795298816,
            3,
            2,
            84
        );
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
        run_single_starkprime_benchmark!(&mut $group, 1, 14, 8192, 268435456, 9, 2, 28);
        run_single_starkprime_benchmark!(&mut $group, 1, 17, 8192, 4294967296, 8, 2, 32);
        run_single_starkprime_benchmark!(&mut $group, 1, 21, 8192, 68719476736, 7, 2, 36);
        run_single_starkprime_benchmark!(&mut $group, 1, 26, 8192, 4398046511104, 6, 2, 42);
        run_single_starkprime_benchmark!(&mut $group, 1, 36, 8192, 2251799813685248, 5, 2, 51);
        run_single_starkprime_benchmark!(&mut $group, 1, 51, 8192, 9223372036854775808, 4, 2, 63);
        run_single_starkprime_benchmark!(
            &mut $group,
            1,
            84,
            8192,
            19342813113834066795298816,
            3,
            2,
            84
        );
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
        run_single_starkprime_benchmark!(&mut $group, 1, 15, 16384, 268435456, 9, 2, 28);
        run_single_starkprime_benchmark!(&mut $group, 1, 18, 16384, 4294967296, 8, 2, 32);
        run_single_starkprime_benchmark!(&mut $group, 1, 21, 16384, 68719476736, 7, 2, 36);
        run_single_starkprime_benchmark!(&mut $group, 1, 27, 16384, 4398046511104, 6, 2, 42);
        run_single_starkprime_benchmark!(&mut $group, 1, 36, 16384, 2251799813685248, 5, 2, 51);
        run_single_starkprime_benchmark!(&mut $group, 1, 52, 16384, 9223372036854775808, 4, 2, 63);
        run_single_starkprime_benchmark!(
            &mut $group,
            1,
            85,
            16384,
            19342813113834066795298816,
            3,
            2,
            84
        );
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
        run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 13, 512, 268435456, 9, 2, 28);
        run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 16, 512, 4294967296, 8, 2, 32);
        run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 19, 512, 68719476736, 7, 2, 36);
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            24,
            512,
            4398046511104,
            6,
            2,
            42
        );
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            34,
            512,
            2251799813685248,
            5,
            2,
            51
        );
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
            81,
            512,
            19342813113834066795298816,
            3,
            2,
            84
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
        run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 13, 1024, 268435456, 9, 2, 28);
        run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 16, 1024, 4294967296, 8, 2, 32);
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            19,
            1024,
            68719476736,
            7,
            2,
            36
        );
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            25,
            1024,
            4398046511104,
            6,
            2,
            42
        );
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            34,
            1024,
            2251799813685248,
            5,
            2,
            51
        );
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
            81,
            1024,
            19342813113834066795298816,
            3,
            2,
            84
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
        run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 14, 2048, 268435456, 9, 2, 28);
        run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 17, 2048, 4294967296, 8, 2, 32);
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            20,
            2048,
            68719476736,
            7,
            2,
            36
        );
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            25,
            2048,
            4398046511104,
            6,
            2,
            42
        );
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            35,
            2048,
            2251799813685248,
            5,
            2,
            51
        );
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
            82,
            2048,
            19342813113834066795298816,
            3,
            2,
            84
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
        run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 14, 4096, 268435456, 9, 2, 28);
        run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 17, 4096, 4294967296, 8, 2, 32);
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            20,
            4096,
            68719476736,
            7,
            2,
            36
        );
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            26,
            4096,
            4398046511104,
            6,
            2,
            42
        );
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            35,
            4096,
            2251799813685248,
            5,
            2,
            51
        );
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
            83,
            4096,
            19342813113834066795298816,
            3,
            2,
            84
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
        run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 14, 8192, 268435456, 9, 2, 28);
        run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 17, 8192, 4294967296, 8, 2, 32);
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            21,
            8192,
            68719476736,
            7,
            2,
            36
        );
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            26,
            8192,
            4398046511104,
            6,
            2,
            42
        );
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            36,
            8192,
            2251799813685248,
            5,
            2,
            51
        );
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
            84,
            8192,
            19342813113834066795298816,
            3,
            2,
            84
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
        run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 15, 16384, 268435456, 9, 2, 28);
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
            21,
            16384,
            68719476736,
            7,
            2,
            36
        );
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            27,
            16384,
            4398046511104,
            6,
            2,
            42
        );
        run_single_starkprime_non_scalar_benchmark!(
            &mut $group,
            1,
            36,
            16384,
            2251799813685248,
            5,
            2,
            51
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
            85,
            16384,
            19342813113834066795298816,
            3,
            2,
            84
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
macro_rules! run_starkprime_degree_three_non_scalar_benchmarks {
    ($group: ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        run_single_starkprime_degree_three_non_scalar_benchmark!(
            &mut $group,
            1,
            15,
            1024,
            1073741824u128,
            9,
            2,
            30
        );
        run_single_starkprime_degree_three_non_scalar_benchmark!(
            &mut $group,
            1,
            16,
            1024,
            4294967296u128,
            8,
            2,
            32
        );
        run_single_starkprime_degree_three_non_scalar_benchmark!(
            &mut $group,
            1,
            17,
            2048,
            8589934592u128,
            8,
            2,
            33
        );
        run_single_starkprime_degree_three_non_scalar_benchmark!(
            &mut $group,
            1,
            18,
            2048,
            17179869184u128,
            8,
            2,
            34
        );
        run_single_starkprime_degree_three_non_scalar_benchmark!(
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
        run_single_frog_benchmark!(&mut $group, 1, 10, 512, 256, 9, 2, 8);
        run_single_frog_benchmark!(&mut $group, 1, 11, 512, 512, 8, 2, 9);
        run_single_frog_benchmark!(&mut $group, 1, 12, 512, 1024, 7, 2, 10);
        run_single_frog_benchmark!(&mut $group, 1, 15, 512, 4096, 6, 2, 12);
        run_single_frog_benchmark!(&mut $group, 1, 17, 512, 16384, 5, 2, 14);
        run_single_frog_benchmark!(&mut $group, 1, 22, 512, 131072, 4, 2, 17);
        run_single_frog_benchmark!(&mut $group, 1, 33, 512, 8388608, 3, 2, 23);
        run_single_frog_benchmark!(&mut $group, 1, 61, 512, 17179869184, 2, 2, 34);
        run_single_frog_benchmark!(&mut $group, 1, 10, 1024, 256, 9, 2, 8);
        run_single_frog_benchmark!(&mut $group, 1, 12, 1024, 512, 8, 2, 9);
        run_single_frog_benchmark!(&mut $group, 1, 13, 1024, 1024, 7, 2, 10);
        run_single_frog_benchmark!(&mut $group, 1, 15, 1024, 4096, 6, 2, 12);
        run_single_frog_benchmark!(&mut $group, 1, 18, 1024, 16384, 5, 2, 14);
        run_single_frog_benchmark!(&mut $group, 1, 23, 1024, 131072, 4, 2, 17);
        run_single_frog_benchmark!(&mut $group, 1, 35, 1024, 8388608, 3, 2, 23);
        run_single_frog_benchmark!(&mut $group, 1, 62, 1024, 17179869184, 2, 2, 34);
        run_single_frog_benchmark!(&mut $group, 1, 11, 2048, 256, 9, 2, 8);
        run_single_frog_benchmark!(&mut $group, 1, 12, 2048, 512, 8, 2, 9);
        run_single_frog_benchmark!(&mut $group, 1, 13, 2048, 1024, 7, 2, 10);
        run_single_frog_benchmark!(&mut $group, 1, 16, 2048, 4096, 6, 2, 12);
        run_single_frog_benchmark!(&mut $group, 1, 19, 2048, 16384, 5, 2, 14);
        run_single_frog_benchmark!(&mut $group, 1, 24, 2048, 131072, 4, 2, 17);
        run_single_frog_benchmark!(&mut $group, 1, 36, 2048, 8388608, 3, 2, 23);
        run_single_frog_benchmark!(&mut $group, 1, 64, 2048, 17179869184, 2, 2, 34);
        run_single_frog_benchmark!(&mut $group, 1, 12, 4096, 256, 9, 2, 8);
        run_single_frog_benchmark!(&mut $group, 1, 13, 4096, 512, 8, 2, 9);
        run_single_frog_benchmark!(&mut $group, 1, 14, 4096, 1024, 7, 2, 10);
        run_single_frog_benchmark!(&mut $group, 1, 17, 4096, 4096, 6, 2, 12);
        run_single_frog_benchmark!(&mut $group, 1, 20, 4096, 16384, 5, 2, 14);
        run_single_frog_benchmark!(&mut $group, 1, 25, 4096, 131072, 4, 2, 17);
        run_single_frog_benchmark!(&mut $group, 1, 37, 4096, 8388608, 3, 2, 23);
        run_single_frog_benchmark!(&mut $group, 1, 65, 4096, 17179869184, 2, 2, 34);
        run_single_frog_benchmark!(&mut $group, 1, 12, 8192, 256, 9, 2, 8);
        run_single_frog_benchmark!(&mut $group, 1, 13, 8192, 512, 8, 2, 9);
        run_single_frog_benchmark!(&mut $group, 1, 15, 8192, 1024, 7, 2, 10);
        run_single_frog_benchmark!(&mut $group, 1, 18, 8192, 4096, 6, 2, 12);
        run_single_frog_benchmark!(&mut $group, 1, 21, 8192, 16384, 5, 2, 14);
        run_single_frog_benchmark!(&mut $group, 1, 26, 8192, 131072, 4, 2, 17);
        run_single_frog_benchmark!(&mut $group, 1, 38, 8192, 8388608, 3, 2, 23);
        run_single_frog_benchmark!(&mut $group, 1, 67, 8192, 17179869184, 2, 2, 34);
        run_single_frog_benchmark!(&mut $group, 1, 13, 16384, 256, 9, 2, 8);
        run_single_frog_benchmark!(&mut $group, 1, 14, 16384, 512, 8, 2, 9);
        run_single_frog_benchmark!(&mut $group, 1, 15, 16384, 1024, 7, 2, 10);
        run_single_frog_benchmark!(&mut $group, 1, 18, 16384, 4096, 6, 2, 12);
        run_single_frog_benchmark!(&mut $group, 1, 21, 16384, 16384, 5, 2, 14);
        run_single_frog_benchmark!(&mut $group, 1, 27, 16384, 131072, 4, 2, 17);
        run_single_frog_benchmark!(&mut $group, 1, 39, 16384, 8388608, 3, 2, 23);
        run_single_frog_benchmark!(&mut $group, 1, 69, 16384, 17179869184, 2, 2, 34);
    };
}

#[macro_export]
macro_rules! run_frog_non_scalar_benchmarks {
    ($group:ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 10, 512, 256, 9, 2, 8);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 11, 512, 512, 8, 2, 9);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 12, 512, 1024, 7, 2, 10);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 15, 512, 4096, 6, 2, 12);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 17, 512, 16384, 5, 2, 14);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 22, 512, 131072, 4, 2, 17);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 33, 512, 8388608, 3, 2, 23);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 61, 512, 17179869184, 2, 2, 34);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 10, 1024, 256, 9, 2, 8);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 12, 1024, 512, 8, 2, 9);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 13, 1024, 1024, 7, 2, 10);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 15, 1024, 4096, 6, 2, 12);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 18, 1024, 16384, 5, 2, 14);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 23, 1024, 131072, 4, 2, 17);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 35, 1024, 8388608, 3, 2, 23);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 62, 1024, 17179869184, 2, 2, 34);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 11, 2048, 256, 9, 2, 8);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 12, 2048, 512, 8, 2, 9);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 13, 2048, 1024, 7, 2, 10);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 16, 2048, 4096, 6, 2, 12);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 19, 2048, 16384, 5, 2, 14);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 24, 2048, 131072, 4, 2, 17);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 36, 2048, 8388608, 3, 2, 23);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 64, 2048, 17179869184, 2, 2, 34);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 12, 4096, 256, 9, 2, 8);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 13, 4096, 512, 8, 2, 9);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 14, 4096, 1024, 7, 2, 10);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 17, 4096, 4096, 6, 2, 12);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 20, 4096, 16384, 5, 2, 14);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 25, 4096, 131072, 4, 2, 17);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 37, 4096, 8388608, 3, 2, 23);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 65, 4096, 17179869184, 2, 2, 34);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 12, 8192, 256, 9, 2, 8);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 13, 8192, 512, 8, 2, 9);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 15, 8192, 1024, 7, 2, 10);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 18, 8192, 4096, 6, 2, 12);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 21, 8192, 16384, 5, 2, 14);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 26, 8192, 131072, 4, 2, 17);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 38, 8192, 8388608, 3, 2, 23);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 67, 8192, 17179869184, 2, 2, 34);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 13, 16384, 256, 9, 2, 8);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 14, 16384, 512, 8, 2, 9);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 15, 16384, 1024, 7, 2, 10);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 18, 16384, 4096, 6, 2, 12);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 21, 16384, 16384, 5, 2, 14);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 27, 16384, 131072, 4, 2, 17);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 39, 16384, 8388608, 3, 2, 23);
        run_single_frog_non_scalar_benchmark!(&mut $group, 1, 69, 16384, 17179869184, 2, 2, 34);
    };
}

#[macro_export]
macro_rules! run_frog_degree_three_non_scalar_benchmarks {
    ($group:ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 5, 512, 8, 23, 2, 3);
        run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 9, 1024, 128, 10, 2, 7);
        run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 10, 1024, 256, 9, 2, 8);
        run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 12, 512, 1024, 7, 2, 10);
        run_single_frog_degree_three_non_scalar_benchmark!(
            &mut $group,
            1,
            15,
            1024,
            4096,
            6,
            2,
            12
        );
    };
}

#[macro_export]
macro_rules! run_babybear_benchmarks {
    ($group:ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        run_single_babybear_benchmark!(&mut $group, 1, 3, 512, 16, 8, 2, 4);
        run_single_babybear_benchmark!(&mut $group, 1, 4, 512, 32, 7, 2, 5);
        run_single_babybear_benchmark!(&mut $group, 1, 5, 512, 256, 4, 2, 8);
        run_single_babybear_benchmark!(&mut $group, 1, 7, 512, 1024, 4, 2, 10);
        run_single_babybear_benchmark!(&mut $group, 1, 11, 512, 32768, 3, 2, 15);
        run_single_babybear_benchmark!(&mut $group, 1, 4, 1024, 32, 7, 2, 5);
        run_single_babybear_benchmark!(&mut $group, 1, 5, 1024, 64, 6, 2, 6);
        run_single_babybear_benchmark!(&mut $group, 1, 6, 1024, 256, 4, 2, 8);
        run_single_babybear_benchmark!(&mut $group, 1, 7, 1024, 1024, 4, 2, 10);
        run_single_babybear_benchmark!(&mut $group, 1, 11, 1024, 32768, 3, 2, 15);
        run_single_babybear_benchmark!(&mut $group, 1, 4, 2048, 32, 7, 2, 5);
        run_single_babybear_benchmark!(&mut $group, 1, 5, 2048, 64, 6, 2, 6);
        run_single_babybear_benchmark!(&mut $group, 1, 6, 2048, 256, 4, 2, 8);
        run_single_babybear_benchmark!(&mut $group, 1, 8, 2048, 1024, 4, 2, 10);
        run_single_babybear_benchmark!(&mut $group, 1, 11, 2048, 32768, 3, 2, 15);
        run_single_babybear_benchmark!(&mut $group, 1, 4, 4096, 16, 8, 2, 4);
        run_single_babybear_benchmark!(&mut $group, 1, 5, 4096, 32, 7, 2, 5);
        run_single_babybear_benchmark!(&mut $group, 1, 6, 4096, 256, 4, 2, 8);
        run_single_babybear_benchmark!(&mut $group, 1, 8, 4096, 1024, 4, 2, 10);
        run_single_babybear_benchmark!(&mut $group, 1, 12, 4096, 32768, 3, 2, 15);
        run_single_babybear_benchmark!(&mut $group, 1, 4, 8192, 16, 8, 2, 4);
        run_single_babybear_benchmark!(&mut $group, 1, 5, 8192, 32, 7, 2, 5);
        run_single_babybear_benchmark!(&mut $group, 1, 6, 8192, 256, 4, 2, 8);
        run_single_babybear_benchmark!(&mut $group, 1, 8, 8192, 1024, 4, 2, 10);
        run_single_babybear_benchmark!(&mut $group, 1, 12, 8192, 32768, 3, 2, 15);
        run_single_babybear_benchmark!(&mut $group, 1, 5, 16384, 32, 7, 2, 5);
        run_single_babybear_benchmark!(&mut $group, 1, 6, 16384, 64, 6, 2, 6);
        run_single_babybear_benchmark!(&mut $group, 1, 7, 16384, 256, 4, 2, 8);
        run_single_babybear_benchmark!(&mut $group, 1, 9, 16384, 1024, 4, 2, 10);
        run_single_babybear_benchmark!(&mut $group, 1, 13, 16384, 32768, 3, 2, 15);
    };
}

#[macro_export]
macro_rules! run_babybear_non_scalar_benchmarks {
    ($group:ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 3, 512, 16, 8, 2, 4);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 4, 512, 32, 7, 2, 5);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 5, 512, 256, 4, 2, 8);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 7, 512, 1024, 4, 2, 10);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 11, 512, 32768, 3, 2, 15);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 4, 1024, 32, 7, 2, 5);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 5, 1024, 64, 6, 2, 6);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 6, 1024, 256, 4, 2, 8);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 7, 1024, 1024, 4, 2, 10);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 11, 1024, 32768, 3, 2, 15);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 4, 2048, 32, 7, 2, 5);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 5, 2048, 64, 6, 2, 6);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 6, 2048, 256, 4, 2, 8);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 8, 2048, 1024, 4, 2, 10);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 11, 2048, 32768, 3, 2, 15);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 4, 4096, 16, 8, 2, 4);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 5, 4096, 32, 7, 2, 5);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 6, 4096, 256, 4, 2, 8);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 8, 4096, 1024, 4, 2, 10);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 12, 4096, 32768, 3, 2, 15);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 4, 8192, 16, 8, 2, 4);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 5, 8192, 32, 7, 2, 5);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 6, 8192, 256, 4, 2, 8);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 8, 8192, 1024, 4, 2, 10);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 12, 8192, 32768, 3, 2, 15);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 5, 16384, 32, 7, 2, 5);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 6, 16384, 64, 6, 2, 6);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 7, 16384, 256, 4, 2, 8);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 9, 16384, 1024, 4, 2, 10);
        run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 13, 16384, 32768, 3, 2, 15);
    };
}

#[macro_export]
macro_rules! run_babybear_degree_three_non_scalar_benchmarks {
    ($group:ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
        run_single_babybear_degree_three_non_scalar_benchmark!(
            &mut $group,
            1,
            6,
            1024,
            512,
            4,
            2,
            9
        );
        run_single_babybear_degree_three_non_scalar_benchmark!(
            &mut $group,
            1,
            7,
            1024,
            2048,
            3,
            2,
            11
        );
        run_single_babybear_degree_three_non_scalar_benchmark!(
            &mut $group,
            1,
            8,
            4096,
            2048,
            3,
            2,
            11
        );
        run_single_babybear_degree_three_non_scalar_benchmark!(
            &mut $group,
            1,
            9,
            2048,
            8192,
            3,
            2,
            13
        );
        run_single_babybear_degree_three_non_scalar_benchmark!(
            &mut $group,
            1,
            10,
            4096,
            16384,
            3,
            2,
            14
        );
    };
}
