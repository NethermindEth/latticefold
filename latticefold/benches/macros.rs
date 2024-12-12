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
         run_single_goldilocks_benchmark!(&mut $group, 1, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 1, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 1, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 1, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 1, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 1, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 1, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 1, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 1, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 1, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 1, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 1, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 1, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 1, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 1, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 1, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 1, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 1, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 1, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 1, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 1, {1 << 20});
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
		run_single_goldilocks_benchmark!(&mut $group, 1, 8, 512, 256, 8, 2, 8);
		run_single_goldilocks_benchmark!(&mut $group, 1, 9, 512, 1024, 7, 2, 10);
		run_single_goldilocks_benchmark!(&mut $group, 1, 11, 512, 2048, 6, 2, 11);
		run_single_goldilocks_benchmark!(&mut $group, 1, 14, 512, 8192, 5, 2, 13);
		run_single_goldilocks_benchmark!(&mut $group, 1, 21, 512, 65536, 4, 2, 16);
		run_single_goldilocks_benchmark!(&mut $group, 1, 38, 512, 4194304, 3, 2, 22);
		run_single_goldilocks_benchmark!(&mut $group, 1, 99, 512, 4294967296, 2, 2, 32);
		run_single_goldilocks_benchmark!(&mut $group, 1, 8, 1024, 256, 8, 2, 8);
		run_single_goldilocks_benchmark!(&mut $group, 1, 9, 1024, 1024, 7, 2, 10);
		run_single_goldilocks_benchmark!(&mut $group, 1, 11, 1024, 2048, 6, 2, 11);
		run_single_goldilocks_benchmark!(&mut $group, 1, 14, 1024, 8192, 5, 2, 13);
		run_single_goldilocks_benchmark!(&mut $group, 1, 22, 1024, 65536, 4, 2, 16);
		run_single_goldilocks_benchmark!(&mut $group, 1, 39, 1024, 4194304, 3, 2, 22);
		run_single_goldilocks_benchmark!(&mut $group, 1, 99, 1024, 4294967296, 2, 2, 32);
		run_single_goldilocks_benchmark!(&mut $group, 1, 9, 2048, 256, 8, 2, 8);
		run_single_goldilocks_benchmark!(&mut $group, 1, 10, 2048, 1024, 7, 2, 10);
		run_single_goldilocks_benchmark!(&mut $group, 1, 12, 2048, 2048, 6, 2, 11);
		run_single_goldilocks_benchmark!(&mut $group, 1, 15, 2048, 8192, 5, 2, 13);
		run_single_goldilocks_benchmark!(&mut $group, 1, 23, 2048, 65536, 4, 2, 16);
		run_single_goldilocks_benchmark!(&mut $group, 1, 40, 2048, 4194304, 3, 2, 22);
		run_single_goldilocks_benchmark!(&mut $group, 1, 99, 2048, 4294967296, 2, 2, 32);
		run_single_goldilocks_benchmark!(&mut $group, 1, 9, 4096, 256, 8, 2, 8);
		run_single_goldilocks_benchmark!(&mut $group, 1, 10, 4096, 1024, 7, 2, 10);
		run_single_goldilocks_benchmark!(&mut $group, 1, 12, 4096, 2048, 6, 2, 11);
		run_single_goldilocks_benchmark!(&mut $group, 1, 16, 4096, 8192, 5, 2, 13);
		run_single_goldilocks_benchmark!(&mut $group, 1, 24, 4096, 65536, 4, 2, 16);
		run_single_goldilocks_benchmark!(&mut $group, 1, 41, 4096, 4194304, 3, 2, 22);
		run_single_goldilocks_benchmark!(&mut $group, 1, 99, 4096, 4294967296, 2, 2, 32);
		run_single_goldilocks_benchmark!(&mut $group, 1, 10, 8192, 256, 8, 2, 8);
		run_single_goldilocks_benchmark!(&mut $group, 1, 11, 8192, 1024, 7, 2, 10);
		run_single_goldilocks_benchmark!(&mut $group, 1, 13, 8192, 2048, 6, 2, 11);
		run_single_goldilocks_benchmark!(&mut $group, 1, 16, 8192, 8192, 5, 2, 13);
		run_single_goldilocks_benchmark!(&mut $group, 1, 25, 8192, 65536, 4, 2, 16);
		run_single_goldilocks_benchmark!(&mut $group, 1, 42, 8192, 4194304, 3, 2, 22);
		run_single_goldilocks_benchmark!(&mut $group, 1, 99, 8192, 4294967296, 2, 2, 32);
		run_single_goldilocks_benchmark!(&mut $group, 1, 10, 16384, 256, 8, 2, 8);
		run_single_goldilocks_benchmark!(&mut $group, 1, 11, 16384, 1024, 7, 2, 10);
		run_single_goldilocks_benchmark!(&mut $group, 1, 13, 16384, 2048, 6, 2, 11);
		run_single_goldilocks_benchmark!(&mut $group, 1, 17, 16384, 8192, 5, 2, 13);
		run_single_goldilocks_benchmark!(&mut $group, 1, 25, 16384, 65536, 4, 2, 16);
		run_single_goldilocks_benchmark!(&mut $group, 1, 43, 16384, 4194304, 3, 2, 22);
		run_single_goldilocks_benchmark!(&mut $group, 1, 99, 16384, 4294967296, 2, 2, 32);
		run_single_goldilocks_benchmark!(&mut $group, 1, 11, 32768, 256, 8, 2, 8);
		run_single_goldilocks_benchmark!(&mut $group, 1, 12, 32768, 1024, 7, 2, 10);
		run_single_goldilocks_benchmark!(&mut $group, 1, 14, 32768, 2048, 6, 2, 11);
		run_single_goldilocks_benchmark!(&mut $group, 1, 18, 32768, 8192, 5, 2, 13);
		run_single_goldilocks_benchmark!(&mut $group, 1, 26, 32768, 65536, 4, 2, 16);
		run_single_goldilocks_benchmark!(&mut $group, 1, 44, 32768, 4194304, 3, 2, 22);
		run_single_goldilocks_benchmark!(&mut $group, 1, 99, 32768, 4294967296, 2, 2, 32);
		run_single_goldilocks_benchmark!(&mut $group, 1, 11, 65536, 256, 8, 2, 8);
		run_single_goldilocks_benchmark!(&mut $group, 1, 12, 65536, 1024, 7, 2, 10);
		run_single_goldilocks_benchmark!(&mut $group, 1, 15, 65536, 2048, 6, 2, 11);
		run_single_goldilocks_benchmark!(&mut $group, 1, 18, 65536, 8192, 5, 2, 13);
		run_single_goldilocks_benchmark!(&mut $group, 1, 27, 65536, 65536, 4, 2, 16);
		run_single_goldilocks_benchmark!(&mut $group, 1, 46, 65536, 4194304, 3, 2, 22);
		run_single_goldilocks_benchmark!(&mut $group, 1, 99, 65536, 4294967296, 2, 2, 32);
		run_single_goldilocks_benchmark!(&mut $group, 1, 12, 131072, 256, 8, 2, 8);
		run_single_goldilocks_benchmark!(&mut $group, 1, 13, 131072, 1024, 7, 2, 10);
		run_single_goldilocks_benchmark!(&mut $group, 1, 15, 131072, 2048, 6, 2, 11);
		run_single_goldilocks_benchmark!(&mut $group, 1, 19, 131072, 8192, 5, 2, 13);
		run_single_goldilocks_benchmark!(&mut $group, 1, 28, 131072, 65536, 4, 2, 16);
		run_single_goldilocks_benchmark!(&mut $group, 1, 47, 131072, 4194304, 3, 2, 22);
		run_single_goldilocks_benchmark!(&mut $group, 1, 99, 131072, 4294967296, 2, 2, 32);
		run_single_goldilocks_benchmark!(&mut $group, 1, 12, 262144, 256, 8, 2, 8);
		run_single_goldilocks_benchmark!(&mut $group, 1, 13, 262144, 1024, 7, 2, 10);
		run_single_goldilocks_benchmark!(&mut $group, 1, 16, 262144, 2048, 6, 2, 11);
		run_single_goldilocks_benchmark!(&mut $group, 1, 20, 262144, 8192, 5, 2, 13);
		run_single_goldilocks_benchmark!(&mut $group, 1, 29, 262144, 65536, 4, 2, 16);
		run_single_goldilocks_benchmark!(&mut $group, 1, 48, 262144, 4194304, 3, 2, 22);
		run_single_goldilocks_benchmark!(&mut $group, 1, 99, 262144, 4294967296, 2, 2, 32);
		run_single_goldilocks_benchmark!(&mut $group, 1, 13, 524288, 256, 8, 2, 8);
		run_single_goldilocks_benchmark!(&mut $group, 1, 14, 524288, 1024, 7, 2, 10);
		run_single_goldilocks_benchmark!(&mut $group, 1, 16, 524288, 2048, 6, 2, 11);
		run_single_goldilocks_benchmark!(&mut $group, 1, 20, 524288, 8192, 5, 2, 13);
		run_single_goldilocks_benchmark!(&mut $group, 1, 30, 524288, 65536, 4, 2, 16);
		run_single_goldilocks_benchmark!(&mut $group, 1, 49, 524288, 4194304, 3, 2, 22);
		run_single_goldilocks_benchmark!(&mut $group, 1, 99, 524288, 4294967296, 2, 2, 32);
		run_single_goldilocks_benchmark!(&mut $group, 1, 14, 1048576, 256, 8, 2, 8);
		run_single_goldilocks_benchmark!(&mut $group, 1, 15, 1048576, 1024, 7, 2, 10);
		run_single_goldilocks_benchmark!(&mut $group, 1, 17, 1048576, 2048, 6, 2, 11);
		run_single_goldilocks_benchmark!(&mut $group, 1, 21, 1048576, 8192, 5, 2, 13);
		run_single_goldilocks_benchmark!(&mut $group, 1, 30, 1048576, 65536, 4, 2, 16);
		run_single_goldilocks_benchmark!(&mut $group, 1, 50, 1048576, 4194304, 3, 2, 22);
		run_single_goldilocks_benchmark!(&mut $group, 1, 99, 1048576, 4294967296, 2, 2, 32);
    };
}

#[macro_export]
macro_rules! run_goldilocks_non_scalar_benchmarks {
    ($group: ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 8, 512, 256, 8, 2, 8);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 9, 512, 1024, 7, 2, 10);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 11, 512, 2048, 6, 2, 11);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 14, 512, 8192, 5, 2, 13);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 21, 512, 65536, 4, 2, 16);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 38, 512, 4194304, 3, 2, 22);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 99, 512, 4294967296, 2, 2, 32);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 8, 1024, 256, 8, 2, 8);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 9, 1024, 1024, 7, 2, 10);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 11, 1024, 2048, 6, 2, 11);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 14, 1024, 8192, 5, 2, 13);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 22, 1024, 65536, 4, 2, 16);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 39, 1024, 4194304, 3, 2, 22);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 99, 1024, 4294967296, 2, 2, 32);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 9, 2048, 256, 8, 2, 8);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 10, 2048, 1024, 7, 2, 10);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 12, 2048, 2048, 6, 2, 11);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 15, 2048, 8192, 5, 2, 13);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 23, 2048, 65536, 4, 2, 16);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 40, 2048, 4194304, 3, 2, 22);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 99, 2048, 4294967296, 2, 2, 32);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 9, 4096, 256, 8, 2, 8);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 10, 4096, 1024, 7, 2, 10);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 12, 4096, 2048, 6, 2, 11);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 16, 4096, 8192, 5, 2, 13);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 24, 4096, 65536, 4, 2, 16);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 41, 4096, 4194304, 3, 2, 22);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 99, 4096, 4294967296, 2, 2, 32);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 10, 8192, 256, 8, 2, 8);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 11, 8192, 1024, 7, 2, 10);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 13, 8192, 2048, 6, 2, 11);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 16, 8192, 8192, 5, 2, 13);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 25, 8192, 65536, 4, 2, 16);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 42, 8192, 4194304, 3, 2, 22);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 99, 8192, 4294967296, 2, 2, 32);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 10, 16384, 256, 8, 2, 8);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 11, 16384, 1024, 7, 2, 10);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 13, 16384, 2048, 6, 2, 11);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 17, 16384, 8192, 5, 2, 13);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 25, 16384, 65536, 4, 2, 16);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 43, 16384, 4194304, 3, 2, 22);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 99, 16384, 4294967296, 2, 2, 32);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 11, 32768, 256, 8, 2, 8);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 12, 32768, 1024, 7, 2, 10);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 14, 32768, 2048, 6, 2, 11);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 18, 32768, 8192, 5, 2, 13);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 26, 32768, 65536, 4, 2, 16);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 44, 32768, 4194304, 3, 2, 22);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 99, 32768, 4294967296, 2, 2, 32);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 11, 65536, 256, 8, 2, 8);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 12, 65536, 1024, 7, 2, 10);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 15, 65536, 2048, 6, 2, 11);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 18, 65536, 8192, 5, 2, 13);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 27, 65536, 65536, 4, 2, 16);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 46, 65536, 4194304, 3, 2, 22);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 99, 65536, 4294967296, 2, 2, 32);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 12, 131072, 256, 8, 2, 8);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 13, 131072, 1024, 7, 2, 10);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 15, 131072, 2048, 6, 2, 11);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 19, 131072, 8192, 5, 2, 13);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 28, 131072, 65536, 4, 2, 16);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 47, 131072, 4194304, 3, 2, 22);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 99, 131072, 4294967296, 2, 2, 32);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 12, 262144, 256, 8, 2, 8);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 13, 262144, 1024, 7, 2, 10);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 16, 262144, 2048, 6, 2, 11);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 20, 262144, 8192, 5, 2, 13);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 29, 262144, 65536, 4, 2, 16);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 48, 262144, 4194304, 3, 2, 22);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 99, 262144, 4294967296, 2, 2, 32);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 13, 524288, 256, 8, 2, 8);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 14, 524288, 1024, 7, 2, 10);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 16, 524288, 2048, 6, 2, 11);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 20, 524288, 8192, 5, 2, 13);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 30, 524288, 65536, 4, 2, 16);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 49, 524288, 4194304, 3, 2, 22);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 99, 524288, 4294967296, 2, 2, 32);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 14, 1048576, 256, 8, 2, 8);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 15, 1048576, 1024, 7, 2, 10);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 17, 1048576, 2048, 6, 2, 11);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 21, 1048576, 8192, 5, 2, 13);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 30, 1048576, 65536, 4, 2, 16);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 50, 1048576, 4194304, 3, 2, 22);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 99, 1048576, 4294967296, 2, 2, 32);
    };
}

#[macro_export]
macro_rules! run_goldilocks_degree_three_non_scalar_benchmarks {
    ($group: ident) => {
        // Parameters Criterion, X_LEN, C, W, B, L, B_small, K
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