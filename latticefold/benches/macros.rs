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

//--- BabyBear cyclotomic ring (modulus p = 2013265921, degree = 72) ---
//	Maximum kappa for which bound_{l_2} < p/2: 15
#[macro_export]
macro_rules! run_babybear_benchmarks {
   ($group:ident) => {
        // Parameters: Criterion group, X_LEN, Kappa, W_CCS, B, L, b, k
		run_single_babybear_benchmark!(&mut $group, 1, 3, 512, 16, 8, 2, 4);
		run_single_babybear_benchmark!(&mut $group, 1, 4, 512, 32, 7, 2, 5);
		run_single_babybear_benchmark!(&mut $group, 1, 6, 512, 256, 4, 2, 8);
		run_single_babybear_benchmark!(&mut $group, 1, 10, 512, 1024, 4, 2, 10);
		run_single_babybear_benchmark!(&mut $group, 1, 15, 512, 32768, 3, 2, 15);
		run_single_babybear_benchmark!(&mut $group, 1, 4, 1024, 32, 7, 2, 5);
		run_single_babybear_benchmark!(&mut $group, 1, 5, 1024, 64, 6, 2, 6);
		run_single_babybear_benchmark!(&mut $group, 1, 6, 1024, 256, 4, 2, 8);
		run_single_babybear_benchmark!(&mut $group, 1, 10, 1024, 1024, 4, 2, 10);
		run_single_babybear_benchmark!(&mut $group, 1, 15, 1024, 32768, 3, 2, 15);
		run_single_babybear_benchmark!(&mut $group, 1, 4, 2048, 32, 7, 2, 5);
		run_single_babybear_benchmark!(&mut $group, 1, 5, 2048, 64, 6, 2, 6);
		run_single_babybear_benchmark!(&mut $group, 1, 7, 2048, 256, 4, 2, 8);
		run_single_babybear_benchmark!(&mut $group, 1, 10, 2048, 1024, 4, 2, 10);
		run_single_babybear_benchmark!(&mut $group, 1, 15, 2048, 32768, 3, 2, 15);
		run_single_babybear_benchmark!(&mut $group, 1, 4, 4096, 16, 8, 2, 4);
		run_single_babybear_benchmark!(&mut $group, 1, 5, 4096, 32, 7, 2, 5);
		run_single_babybear_benchmark!(&mut $group, 1, 7, 4096, 256, 4, 2, 8);
		run_single_babybear_benchmark!(&mut $group, 1, 11, 4096, 1024, 4, 2, 10);
		run_single_babybear_benchmark!(&mut $group, 1, 15, 4096, 32768, 3, 2, 15);
		run_single_babybear_benchmark!(&mut $group, 1, 4, 8192, 16, 8, 2, 4);
		run_single_babybear_benchmark!(&mut $group, 1, 5, 8192, 32, 7, 2, 5);
		run_single_babybear_benchmark!(&mut $group, 1, 7, 8192, 256, 4, 2, 8);
		run_single_babybear_benchmark!(&mut $group, 1, 11, 8192, 1024, 4, 2, 10);
		run_single_babybear_benchmark!(&mut $group, 1, 15, 8192, 32768, 3, 2, 15);
   };
}

#[macro_export]
macro_rules! run_babybear_non_scalar_benchmarks {
   ($group:ident) => {
        // Parameters: Criterion group, X_LEN, Kappa, W_CCS, B, L, b, k
		run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 3, 512, 16, 8, 2, 4);
		run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 4, 512, 32, 7, 2, 5);
		run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 6, 512, 256, 4, 2, 8);
		run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 10, 512, 1024, 4, 2, 10);
		run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 15, 512, 32768, 3, 2, 15);
		run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 4, 1024, 32, 7, 2, 5);
		run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 5, 1024, 64, 6, 2, 6);
		run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 6, 1024, 256, 4, 2, 8);
		run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 10, 1024, 1024, 4, 2, 10);
		run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 15, 1024, 32768, 3, 2, 15);
		run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 4, 2048, 32, 7, 2, 5);
		run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 5, 2048, 64, 6, 2, 6);
		run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 7, 2048, 256, 4, 2, 8);
		run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 10, 2048, 1024, 4, 2, 10);
		run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 15, 2048, 32768, 3, 2, 15);
		run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 4, 4096, 16, 8, 2, 4);
		run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 5, 4096, 32, 7, 2, 5);
		run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 7, 4096, 256, 4, 2, 8);
		run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 11, 4096, 1024, 4, 2, 10);
		run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 15, 4096, 32768, 3, 2, 15);
		run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 4, 8192, 16, 8, 2, 4);
		run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 5, 8192, 32, 7, 2, 5);
		run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 7, 8192, 256, 4, 2, 8);
		run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 11, 8192, 1024, 4, 2, 10);
		run_single_babybear_non_scalar_benchmark!(&mut $group, 1, 15, 8192, 32768, 3, 2, 15);
   };
}

#[macro_export]
macro_rules! run_babybear_degree_three_non_scalar_benchmarks {
   ($group:ident) => {
        // Parameters: Criterion group, X_LEN, Kappa, W_CCS, B, L, b, k
		run_single_babybear_degree_three_non_scalar_benchmark!(&mut $group, 1, 3, 512, 16, 8, 2, 4);
		run_single_babybear_degree_three_non_scalar_benchmark!(&mut $group, 1, 4, 512, 32, 7, 2, 5);
		run_single_babybear_degree_three_non_scalar_benchmark!(&mut $group, 1, 6, 512, 256, 4, 2, 8);
		run_single_babybear_degree_three_non_scalar_benchmark!(&mut $group, 1, 10, 512, 1024, 4, 2, 10);
		run_single_babybear_degree_three_non_scalar_benchmark!(&mut $group, 1, 15, 512, 32768, 3, 2, 15);
		run_single_babybear_degree_three_non_scalar_benchmark!(&mut $group, 1, 4, 1024, 32, 7, 2, 5);
		run_single_babybear_degree_three_non_scalar_benchmark!(&mut $group, 1, 5, 1024, 64, 6, 2, 6);
		run_single_babybear_degree_three_non_scalar_benchmark!(&mut $group, 1, 6, 1024, 256, 4, 2, 8);
		run_single_babybear_degree_three_non_scalar_benchmark!(&mut $group, 1, 10, 1024, 1024, 4, 2, 10);
		run_single_babybear_degree_three_non_scalar_benchmark!(&mut $group, 1, 15, 1024, 32768, 3, 2, 15);
		run_single_babybear_degree_three_non_scalar_benchmark!(&mut $group, 1, 4, 2048, 32, 7, 2, 5);
		run_single_babybear_degree_three_non_scalar_benchmark!(&mut $group, 1, 5, 2048, 64, 6, 2, 6);
		run_single_babybear_degree_three_non_scalar_benchmark!(&mut $group, 1, 7, 2048, 256, 4, 2, 8);
		run_single_babybear_degree_three_non_scalar_benchmark!(&mut $group, 1, 10, 2048, 1024, 4, 2, 10);
		run_single_babybear_degree_three_non_scalar_benchmark!(&mut $group, 1, 15, 2048, 32768, 3, 2, 15);
		run_single_babybear_degree_three_non_scalar_benchmark!(&mut $group, 1, 4, 4096, 16, 8, 2, 4);
		run_single_babybear_degree_three_non_scalar_benchmark!(&mut $group, 1, 5, 4096, 32, 7, 2, 5);
		run_single_babybear_degree_three_non_scalar_benchmark!(&mut $group, 1, 7, 4096, 256, 4, 2, 8);
		run_single_babybear_degree_three_non_scalar_benchmark!(&mut $group, 1, 11, 4096, 1024, 4, 2, 10);
		run_single_babybear_degree_three_non_scalar_benchmark!(&mut $group, 1, 15, 4096, 32768, 3, 2, 15);
		run_single_babybear_degree_three_non_scalar_benchmark!(&mut $group, 1, 4, 8192, 16, 8, 2, 4);
		run_single_babybear_degree_three_non_scalar_benchmark!(&mut $group, 1, 5, 8192, 32, 7, 2, 5);
		run_single_babybear_degree_three_non_scalar_benchmark!(&mut $group, 1, 7, 8192, 256, 4, 2, 8);
		run_single_babybear_degree_three_non_scalar_benchmark!(&mut $group, 1, 11, 8192, 1024, 4, 2, 10);
		run_single_babybear_degree_three_non_scalar_benchmark!(&mut $group, 1, 15, 8192, 32768, 3, 2, 15);
   };
}

//--- Goldilocks cyclotomic ring (modulus p = 18446744069414584321, degree = 24) ---
//	Maximum kappa for which bound_{l_2} < p/2: 99
#[macro_export]
macro_rules! run_goldilocks_benchmarks {
   ($group:ident) => {
        // Parameters: Criterion group, X_LEN, Kappa, W_CCS, B, L, b, k
		run_single_goldilocks_benchmark!(&mut $group, 1, 8, 512, 256, 8, 2, 8);
		run_single_goldilocks_benchmark!(&mut $group, 1, 9, 512, 512, 8, 2, 9);
		run_single_goldilocks_benchmark!(&mut $group, 1, 11, 512, 2048, 6, 2, 11);
		run_single_goldilocks_benchmark!(&mut $group, 1, 14, 512, 8192, 5, 2, 13);
		run_single_goldilocks_benchmark!(&mut $group, 1, 21, 512, 65536, 4, 2, 16);
		run_single_goldilocks_benchmark!(&mut $group, 1, 38, 512, 2097152, 4, 2, 21);
		run_single_goldilocks_benchmark!(&mut $group, 1, 99, 512, 4294967296, 2, 2, 32);
		run_single_goldilocks_benchmark!(&mut $group, 1, 8, 1024, 256, 8, 2, 8);
		run_single_goldilocks_benchmark!(&mut $group, 1, 9, 1024, 512, 8, 2, 9);
		run_single_goldilocks_benchmark!(&mut $group, 1, 11, 1024, 2048, 6, 2, 11);
		run_single_goldilocks_benchmark!(&mut $group, 1, 14, 1024, 8192, 5, 2, 13);
		run_single_goldilocks_benchmark!(&mut $group, 1, 22, 1024, 65536, 4, 2, 16);
		run_single_goldilocks_benchmark!(&mut $group, 1, 39, 1024, 2097152, 4, 2, 21);
		run_single_goldilocks_benchmark!(&mut $group, 1, 99, 1024, 4294967296, 2, 2, 32);
		run_single_goldilocks_benchmark!(&mut $group, 1, 9, 2048, 256, 8, 2, 8);
		run_single_goldilocks_benchmark!(&mut $group, 1, 10, 2048, 512, 8, 2, 9);
		run_single_goldilocks_benchmark!(&mut $group, 1, 12, 2048, 2048, 6, 2, 11);
		run_single_goldilocks_benchmark!(&mut $group, 1, 15, 2048, 8192, 5, 2, 13);
		run_single_goldilocks_benchmark!(&mut $group, 1, 23, 2048, 65536, 4, 2, 16);
		run_single_goldilocks_benchmark!(&mut $group, 1, 40, 2048, 2097152, 4, 2, 21);
		run_single_goldilocks_benchmark!(&mut $group, 1, 99, 2048, 4294967296, 2, 2, 32);
		run_single_goldilocks_benchmark!(&mut $group, 1, 9, 4096, 256, 8, 2, 8);
		run_single_goldilocks_benchmark!(&mut $group, 1, 10, 4096, 512, 8, 2, 9);
		run_single_goldilocks_benchmark!(&mut $group, 1, 12, 4096, 2048, 6, 2, 11);
		run_single_goldilocks_benchmark!(&mut $group, 1, 16, 4096, 8192, 5, 2, 13);
		run_single_goldilocks_benchmark!(&mut $group, 1, 24, 4096, 65536, 4, 2, 16);
		run_single_goldilocks_benchmark!(&mut $group, 1, 41, 4096, 2097152, 4, 2, 21);
		run_single_goldilocks_benchmark!(&mut $group, 1, 99, 4096, 4294967296, 2, 2, 32);
		run_single_goldilocks_benchmark!(&mut $group, 1, 10, 8192, 256, 8, 2, 8);
		run_single_goldilocks_benchmark!(&mut $group, 1, 11, 8192, 512, 8, 2, 9);
		run_single_goldilocks_benchmark!(&mut $group, 1, 13, 8192, 2048, 6, 2, 11);
		run_single_goldilocks_benchmark!(&mut $group, 1, 16, 8192, 8192, 5, 2, 13);
		run_single_goldilocks_benchmark!(&mut $group, 1, 25, 8192, 65536, 4, 2, 16);
		run_single_goldilocks_benchmark!(&mut $group, 1, 42, 8192, 2097152, 4, 2, 21);
		run_single_goldilocks_benchmark!(&mut $group, 1, 99, 8192, 4294967296, 2, 2, 32);
   };
}

#[macro_export]
macro_rules! run_goldilocks_non_scalar_benchmarks {
   ($group:ident) => {
        // Parameters: Criterion group, X_LEN, Kappa, W_CCS, B, L, b, k
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 8, 512, 256, 8, 2, 8);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 9, 512, 512, 8, 2, 9);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 11, 512, 2048, 6, 2, 11);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 14, 512, 8192, 5, 2, 13);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 21, 512, 65536, 4, 2, 16);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 38, 512, 2097152, 4, 2, 21);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 99, 512, 4294967296, 2, 2, 32);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 8, 1024, 256, 8, 2, 8);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 9, 1024, 512, 8, 2, 9);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 11, 1024, 2048, 6, 2, 11);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 14, 1024, 8192, 5, 2, 13);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 22, 1024, 65536, 4, 2, 16);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 39, 1024, 2097152, 4, 2, 21);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 99, 1024, 4294967296, 2, 2, 32);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 9, 2048, 256, 8, 2, 8);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 10, 2048, 512, 8, 2, 9);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 12, 2048, 2048, 6, 2, 11);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 15, 2048, 8192, 5, 2, 13);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 23, 2048, 65536, 4, 2, 16);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 40, 2048, 2097152, 4, 2, 21);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 99, 2048, 4294967296, 2, 2, 32);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 9, 4096, 256, 8, 2, 8);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 10, 4096, 512, 8, 2, 9);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 12, 4096, 2048, 6, 2, 11);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 16, 4096, 8192, 5, 2, 13);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 24, 4096, 65536, 4, 2, 16);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 41, 4096, 2097152, 4, 2, 21);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 99, 4096, 4294967296, 2, 2, 32);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 10, 8192, 256, 8, 2, 8);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 11, 8192, 512, 8, 2, 9);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 13, 8192, 2048, 6, 2, 11);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 16, 8192, 8192, 5, 2, 13);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 25, 8192, 65536, 4, 2, 16);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 42, 8192, 2097152, 4, 2, 21);
		run_single_goldilocks_non_scalar_benchmark!(&mut $group, 1, 99, 8192, 4294967296, 2, 2, 32);
   };
}

#[macro_export]
macro_rules! run_goldilocks_degree_three_non_scalar_benchmarks {
   ($group:ident) => {
        // Parameters: Criterion group, X_LEN, Kappa, W_CCS, B, L, b, k
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 8, 512, 256, 8, 2, 8);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 9, 512, 512, 8, 2, 9);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 11, 512, 2048, 6, 2, 11);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 14, 512, 8192, 5, 2, 13);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 21, 512, 65536, 4, 2, 16);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 38, 512, 2097152, 4, 2, 21);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 99, 512, 4294967296, 2, 2, 32);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 8, 1024, 256, 8, 2, 8);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 9, 1024, 512, 8, 2, 9);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 11, 1024, 2048, 6, 2, 11);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 14, 1024, 8192, 5, 2, 13);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 22, 1024, 65536, 4, 2, 16);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 39, 1024, 2097152, 4, 2, 21);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 99, 1024, 4294967296, 2, 2, 32);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 9, 2048, 256, 8, 2, 8);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 10, 2048, 512, 8, 2, 9);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 12, 2048, 2048, 6, 2, 11);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 15, 2048, 8192, 5, 2, 13);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 23, 2048, 65536, 4, 2, 16);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 40, 2048, 2097152, 4, 2, 21);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 99, 2048, 4294967296, 2, 2, 32);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 9, 4096, 256, 8, 2, 8);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 10, 4096, 512, 8, 2, 9);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 12, 4096, 2048, 6, 2, 11);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 16, 4096, 8192, 5, 2, 13);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 24, 4096, 65536, 4, 2, 16);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 41, 4096, 2097152, 4, 2, 21);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 99, 4096, 4294967296, 2, 2, 32);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 10, 8192, 256, 8, 2, 8);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 11, 8192, 512, 8, 2, 9);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 13, 8192, 2048, 6, 2, 11);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 16, 8192, 8192, 5, 2, 13);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 25, 8192, 65536, 4, 2, 16);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 42, 8192, 2097152, 4, 2, 21);
		run_single_goldilocks_degree_three_non_scalar_benchmark!(&mut $group, 1, 99, 8192, 4294967296, 2, 2, 32);
   };
}

//--- StarkPrime cyclotomic ring (modulus p = 3618502788666131213697322783095070105623107215331596699973092056135872020481, degree = 16) ---
//	Maximum kappa for which bound_{l_2} < p/2: 600
#[macro_export]
macro_rules! run_starkprime_benchmarks {
   ($group:ident) => {
        // Parameters: Criterion group, X_LEN, Kappa, W_CCS, B, L, b, k
		run_single_starkprime_benchmark!(&mut $group, 1, 15, 512, 268435456, 9, 2, 28);
		run_single_starkprime_benchmark!(&mut $group, 1, 18, 512, 4294967296, 8, 2, 32);
		run_single_starkprime_benchmark!(&mut $group, 1, 23, 512, 68719476736, 7, 2, 36);
		run_single_starkprime_benchmark!(&mut $group, 1, 33, 512, 4398046511104, 6, 2, 42);
		run_single_starkprime_benchmark!(&mut $group, 1, 47, 512, 2251799813685248, 5, 2, 51);
		run_single_starkprime_benchmark!(&mut $group, 1, 80, 512, 9223372036854775808, 4, 2, 63);
		run_single_starkprime_benchmark!(&mut $group, 1, 169, 512, 19342813113834066795298816, 3, 2, 84);
		run_single_starkprime_benchmark!(&mut $group, 1, 600, 512, 85070591730234615865843651857942052864, 2, 2, 126);
		run_single_starkprime_benchmark!(&mut $group, 1, 15, 1024, 268435456, 9, 2, 28);
		run_single_starkprime_benchmark!(&mut $group, 1, 18, 1024, 4294967296, 8, 2, 32);
		run_single_starkprime_benchmark!(&mut $group, 1, 24, 1024, 68719476736, 7, 2, 36);
		run_single_starkprime_benchmark!(&mut $group, 1, 33, 1024, 4398046511104, 6, 2, 42);
		run_single_starkprime_benchmark!(&mut $group, 1, 48, 1024, 2251799813685248, 5, 2, 51);
		run_single_starkprime_benchmark!(&mut $group, 1, 80, 1024, 9223372036854775808, 4, 2, 63);
		run_single_starkprime_benchmark!(&mut $group, 1, 171, 1024, 19342813113834066795298816, 3, 2, 84);
		run_single_starkprime_benchmark!(&mut $group, 1, 600, 1024, 85070591730234615865843651857942052864, 2, 2, 126);
		run_single_starkprime_benchmark!(&mut $group, 1, 16, 2048, 268435456, 9, 2, 28);
		run_single_starkprime_benchmark!(&mut $group, 1, 19, 2048, 4294967296, 8, 2, 32);
		run_single_starkprime_benchmark!(&mut $group, 1, 24, 2048, 68719476736, 7, 2, 36);
		run_single_starkprime_benchmark!(&mut $group, 1, 34, 2048, 4398046511104, 6, 2, 42);
		run_single_starkprime_benchmark!(&mut $group, 1, 49, 2048, 2251799813685248, 5, 2, 51);
		run_single_starkprime_benchmark!(&mut $group, 1, 81, 2048, 9223372036854775808, 4, 2, 63);
		run_single_starkprime_benchmark!(&mut $group, 1, 172, 2048, 19342813113834066795298816, 3, 2, 84);
		run_single_starkprime_benchmark!(&mut $group, 1, 600, 2048, 85070591730234615865843651857942052864, 2, 2, 126);
		run_single_starkprime_benchmark!(&mut $group, 1, 16, 4096, 268435456, 9, 2, 28);
		run_single_starkprime_benchmark!(&mut $group, 1, 19, 4096, 4294967296, 8, 2, 32);
		run_single_starkprime_benchmark!(&mut $group, 1, 25, 4096, 68719476736, 7, 2, 36);
		run_single_starkprime_benchmark!(&mut $group, 1, 34, 4096, 4398046511104, 6, 2, 42);
		run_single_starkprime_benchmark!(&mut $group, 1, 49, 4096, 2251799813685248, 5, 2, 51);
		run_single_starkprime_benchmark!(&mut $group, 1, 82, 4096, 9223372036854775808, 4, 2, 63);
		run_single_starkprime_benchmark!(&mut $group, 1, 173, 4096, 19342813113834066795298816, 3, 2, 84);
		run_single_starkprime_benchmark!(&mut $group, 1, 600, 4096, 85070591730234615865843651857942052864, 2, 2, 126);
		run_single_starkprime_benchmark!(&mut $group, 1, 16, 8192, 268435456, 9, 2, 28);
		run_single_starkprime_benchmark!(&mut $group, 1, 20, 8192, 4294967296, 8, 2, 32);
		run_single_starkprime_benchmark!(&mut $group, 1, 25, 8192, 68719476736, 7, 2, 36);
		run_single_starkprime_benchmark!(&mut $group, 1, 35, 8192, 4398046511104, 6, 2, 42);
		run_single_starkprime_benchmark!(&mut $group, 1, 50, 8192, 2251799813685248, 5, 2, 51);
		run_single_starkprime_benchmark!(&mut $group, 1, 83, 8192, 9223372036854775808, 4, 2, 63);
		run_single_starkprime_benchmark!(&mut $group, 1, 175, 8192, 19342813113834066795298816, 3, 2, 84);
		run_single_starkprime_benchmark!(&mut $group, 1, 600, 8192, 85070591730234615865843651857942052864, 2, 2, 126);
   };
}

#[macro_export]
macro_rules! run_starkprime_non_scalar_benchmarks {
   ($group:ident) => {
        // Parameters: Criterion group, X_LEN, Kappa, W_CCS, B, L, b, k
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 15, 512, 268435456, 9, 2, 28);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 18, 512, 4294967296, 8, 2, 32);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 23, 512, 68719476736, 7, 2, 36);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 33, 512, 4398046511104, 6, 2, 42);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 47, 512, 2251799813685248, 5, 2, 51);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 80, 512, 9223372036854775808, 4, 2, 63);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 169, 512, 19342813113834066795298816, 3, 2, 84);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 600, 512, 85070591730234615865843651857942052864, 2, 2, 126);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 15, 1024, 268435456, 9, 2, 28);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 18, 1024, 4294967296, 8, 2, 32);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 24, 1024, 68719476736, 7, 2, 36);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 33, 1024, 4398046511104, 6, 2, 42);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 48, 1024, 2251799813685248, 5, 2, 51);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 80, 1024, 9223372036854775808, 4, 2, 63);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 171, 1024, 19342813113834066795298816, 3, 2, 84);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 600, 1024, 85070591730234615865843651857942052864, 2, 2, 126);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 16, 2048, 268435456, 9, 2, 28);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 19, 2048, 4294967296, 8, 2, 32);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 24, 2048, 68719476736, 7, 2, 36);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 34, 2048, 4398046511104, 6, 2, 42);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 49, 2048, 2251799813685248, 5, 2, 51);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 81, 2048, 9223372036854775808, 4, 2, 63);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 172, 2048, 19342813113834066795298816, 3, 2, 84);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 600, 2048, 85070591730234615865843651857942052864, 2, 2, 126);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 16, 4096, 268435456, 9, 2, 28);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 19, 4096, 4294967296, 8, 2, 32);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 25, 4096, 68719476736, 7, 2, 36);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 34, 4096, 4398046511104, 6, 2, 42);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 49, 4096, 2251799813685248, 5, 2, 51);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 82, 4096, 9223372036854775808, 4, 2, 63);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 173, 4096, 19342813113834066795298816, 3, 2, 84);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 600, 4096, 85070591730234615865843651857942052864, 2, 2, 126);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 16, 8192, 268435456, 9, 2, 28);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 20, 8192, 4294967296, 8, 2, 32);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 25, 8192, 68719476736, 7, 2, 36);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 35, 8192, 4398046511104, 6, 2, 42);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 50, 8192, 2251799813685248, 5, 2, 51);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 83, 8192, 9223372036854775808, 4, 2, 63);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 175, 8192, 19342813113834066795298816, 3, 2, 84);
		run_single_starkprime_non_scalar_benchmark!(&mut $group, 1, 600, 8192, 85070591730234615865843651857942052864, 2, 2, 126);
   };
}

#[macro_export]
macro_rules! run_starkprime_degree_three_non_scalar_benchmarks {
   ($group:ident) => {
        // Parameters: Criterion group, X_LEN, Kappa, W_CCS, B, L, b, k
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 15, 512, 268435456, 9, 2, 28);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 18, 512, 4294967296, 8, 2, 32);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 23, 512, 68719476736, 7, 2, 36);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 33, 512, 4398046511104, 6, 2, 42);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 47, 512, 2251799813685248, 5, 2, 51);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 80, 512, 9223372036854775808, 4, 2, 63);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 169, 512, 19342813113834066795298816, 3, 2, 84);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 600, 512, 85070591730234615865843651857942052864, 2, 2, 126);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 15, 1024, 268435456, 9, 2, 28);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 18, 1024, 4294967296, 8, 2, 32);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 24, 1024, 68719476736, 7, 2, 36);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 33, 1024, 4398046511104, 6, 2, 42);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 48, 1024, 2251799813685248, 5, 2, 51);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 80, 1024, 9223372036854775808, 4, 2, 63);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 171, 1024, 19342813113834066795298816, 3, 2, 84);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 600, 1024, 85070591730234615865843651857942052864, 2, 2, 126);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 16, 2048, 268435456, 9, 2, 28);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 19, 2048, 4294967296, 8, 2, 32);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 24, 2048, 68719476736, 7, 2, 36);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 34, 2048, 4398046511104, 6, 2, 42);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 49, 2048, 2251799813685248, 5, 2, 51);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 81, 2048, 9223372036854775808, 4, 2, 63);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 172, 2048, 19342813113834066795298816, 3, 2, 84);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 600, 2048, 85070591730234615865843651857942052864, 2, 2, 126);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 16, 4096, 268435456, 9, 2, 28);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 19, 4096, 4294967296, 8, 2, 32);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 25, 4096, 68719476736, 7, 2, 36);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 34, 4096, 4398046511104, 6, 2, 42);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 49, 4096, 2251799813685248, 5, 2, 51);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 82, 4096, 9223372036854775808, 4, 2, 63);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 173, 4096, 19342813113834066795298816, 3, 2, 84);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 600, 4096, 85070591730234615865843651857942052864, 2, 2, 126);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 16, 8192, 268435456, 9, 2, 28);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 20, 8192, 4294967296, 8, 2, 32);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 25, 8192, 68719476736, 7, 2, 36);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 35, 8192, 4398046511104, 6, 2, 42);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 50, 8192, 2251799813685248, 5, 2, 51);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 83, 8192, 9223372036854775808, 4, 2, 63);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 175, 8192, 19342813113834066795298816, 3, 2, 84);
		run_single_starkprime_degree_three_non_scalar_benchmark!(&mut $group, 1, 600, 8192, 85070591730234615865843651857942052864, 2, 2, 126);
   };
}

//--- Frog cyclotomic ring (modulus p = 159120925213255836417, degree = 16) ---
//	Maximum kappa for which bound_{l_2} < p/2: 157
#[macro_export]
macro_rules! run_frog_benchmarks {
   ($group:ident) => {
        // Parameters: Criterion group, X_LEN, Kappa, W_CCS, B, L, b, k
		run_single_frog_benchmark!(&mut $group, 1, 10, 512, 256, 9, 2, 8);
		run_single_frog_benchmark!(&mut $group, 1, 11, 512, 512, 8, 2, 9);
		run_single_frog_benchmark!(&mut $group, 1, 14, 512, 1024, 7, 2, 10);
		run_single_frog_benchmark!(&mut $group, 1, 16, 512, 4096, 6, 2, 12);
		run_single_frog_benchmark!(&mut $group, 1, 21, 512, 16384, 5, 2, 14);
		run_single_frog_benchmark!(&mut $group, 1, 32, 512, 131072, 4, 2, 17);
		run_single_frog_benchmark!(&mut $group, 1, 60, 512, 8388608, 3, 2, 23);
		run_single_frog_benchmark!(&mut $group, 1, 157, 512, 17179869184, 2, 2, 34);
		run_single_frog_benchmark!(&mut $group, 1, 11, 1024, 256, 9, 2, 8);
		run_single_frog_benchmark!(&mut $group, 1, 12, 1024, 512, 8, 2, 9);
		run_single_frog_benchmark!(&mut $group, 1, 14, 1024, 1024, 7, 2, 10);
		run_single_frog_benchmark!(&mut $group, 1, 17, 1024, 4096, 6, 2, 12);
		run_single_frog_benchmark!(&mut $group, 1, 22, 1024, 16384, 5, 2, 14);
		run_single_frog_benchmark!(&mut $group, 1, 34, 1024, 131072, 4, 2, 17);
		run_single_frog_benchmark!(&mut $group, 1, 61, 1024, 8388608, 3, 2, 23);
		run_single_frog_benchmark!(&mut $group, 1, 157, 1024, 17179869184, 2, 2, 34);
		run_single_frog_benchmark!(&mut $group, 1, 11, 2048, 256, 9, 2, 8);
		run_single_frog_benchmark!(&mut $group, 1, 12, 2048, 512, 8, 2, 9);
		run_single_frog_benchmark!(&mut $group, 1, 15, 2048, 1024, 7, 2, 10);
		run_single_frog_benchmark!(&mut $group, 1, 18, 2048, 4096, 6, 2, 12);
		run_single_frog_benchmark!(&mut $group, 1, 23, 2048, 16384, 5, 2, 14);
		run_single_frog_benchmark!(&mut $group, 1, 35, 2048, 131072, 4, 2, 17);
		run_single_frog_benchmark!(&mut $group, 1, 63, 2048, 8388608, 3, 2, 23);
		run_single_frog_benchmark!(&mut $group, 1, 157, 2048, 17179869184, 2, 2, 34);
		run_single_frog_benchmark!(&mut $group, 1, 12, 4096, 256, 9, 2, 8);
		run_single_frog_benchmark!(&mut $group, 1, 13, 4096, 512, 8, 2, 9);
		run_single_frog_benchmark!(&mut $group, 1, 16, 4096, 1024, 7, 2, 10);
		run_single_frog_benchmark!(&mut $group, 1, 19, 4096, 4096, 6, 2, 12);
		run_single_frog_benchmark!(&mut $group, 1, 24, 4096, 16384, 5, 2, 14);
		run_single_frog_benchmark!(&mut $group, 1, 36, 4096, 131072, 4, 2, 17);
		run_single_frog_benchmark!(&mut $group, 1, 64, 4096, 8388608, 3, 2, 23);
		run_single_frog_benchmark!(&mut $group, 1, 157, 4096, 17179869184, 2, 2, 34);
		run_single_frog_benchmark!(&mut $group, 1, 12, 8192, 256, 9, 2, 8);
		run_single_frog_benchmark!(&mut $group, 1, 14, 8192, 512, 8, 2, 9);
		run_single_frog_benchmark!(&mut $group, 1, 17, 8192, 1024, 7, 2, 10);
		run_single_frog_benchmark!(&mut $group, 1, 20, 8192, 4096, 6, 2, 12);
		run_single_frog_benchmark!(&mut $group, 1, 25, 8192, 16384, 5, 2, 14);
		run_single_frog_benchmark!(&mut $group, 1, 37, 8192, 131072, 4, 2, 17);
		run_single_frog_benchmark!(&mut $group, 1, 66, 8192, 8388608, 3, 2, 23);
		run_single_frog_benchmark!(&mut $group, 1, 157, 8192, 17179869184, 2, 2, 34);
   };
}

#[macro_export]
macro_rules! run_frog_non_scalar_benchmarks {
   ($group:ident) => {
        // Parameters: Criterion group, X_LEN, Kappa, W_CCS, B, L, b, k
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 10, 512, 256, 9, 2, 8);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 11, 512, 512, 8, 2, 9);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 14, 512, 1024, 7, 2, 10);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 16, 512, 4096, 6, 2, 12);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 21, 512, 16384, 5, 2, 14);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 32, 512, 131072, 4, 2, 17);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 60, 512, 8388608, 3, 2, 23);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 157, 512, 17179869184, 2, 2, 34);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 11, 1024, 256, 9, 2, 8);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 12, 1024, 512, 8, 2, 9);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 14, 1024, 1024, 7, 2, 10);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 17, 1024, 4096, 6, 2, 12);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 22, 1024, 16384, 5, 2, 14);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 34, 1024, 131072, 4, 2, 17);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 61, 1024, 8388608, 3, 2, 23);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 157, 1024, 17179869184, 2, 2, 34);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 11, 2048, 256, 9, 2, 8);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 12, 2048, 512, 8, 2, 9);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 15, 2048, 1024, 7, 2, 10);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 18, 2048, 4096, 6, 2, 12);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 23, 2048, 16384, 5, 2, 14);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 35, 2048, 131072, 4, 2, 17);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 63, 2048, 8388608, 3, 2, 23);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 157, 2048, 17179869184, 2, 2, 34);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 12, 4096, 256, 9, 2, 8);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 13, 4096, 512, 8, 2, 9);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 16, 4096, 1024, 7, 2, 10);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 19, 4096, 4096, 6, 2, 12);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 24, 4096, 16384, 5, 2, 14);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 36, 4096, 131072, 4, 2, 17);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 64, 4096, 8388608, 3, 2, 23);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 157, 4096, 17179869184, 2, 2, 34);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 12, 8192, 256, 9, 2, 8);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 14, 8192, 512, 8, 2, 9);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 17, 8192, 1024, 7, 2, 10);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 20, 8192, 4096, 6, 2, 12);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 25, 8192, 16384, 5, 2, 14);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 37, 8192, 131072, 4, 2, 17);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 66, 8192, 8388608, 3, 2, 23);
		run_single_frog_non_scalar_benchmark!(&mut $group, 1, 157, 8192, 17179869184, 2, 2, 34);
   };
}

#[macro_export]
macro_rules! run_frog_degree_three_non_scalar_benchmarks {
   ($group:ident) => {
        // Parameters: Criterion group, X_LEN, Kappa, W_CCS, B, L, b, k
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 10, 512, 256, 9, 2, 8);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 11, 512, 512, 8, 2, 9);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 14, 512, 1024, 7, 2, 10);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 16, 512, 4096, 6, 2, 12);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 21, 512, 16384, 5, 2, 14);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 32, 512, 131072, 4, 2, 17);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 60, 512, 8388608, 3, 2, 23);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 157, 512, 17179869184, 2, 2, 34);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 11, 1024, 256, 9, 2, 8);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 12, 1024, 512, 8, 2, 9);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 14, 1024, 1024, 7, 2, 10);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 17, 1024, 4096, 6, 2, 12);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 22, 1024, 16384, 5, 2, 14);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 34, 1024, 131072, 4, 2, 17);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 61, 1024, 8388608, 3, 2, 23);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 157, 1024, 17179869184, 2, 2, 34);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 11, 2048, 256, 9, 2, 8);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 12, 2048, 512, 8, 2, 9);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 15, 2048, 1024, 7, 2, 10);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 18, 2048, 4096, 6, 2, 12);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 23, 2048, 16384, 5, 2, 14);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 35, 2048, 131072, 4, 2, 17);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 63, 2048, 8388608, 3, 2, 23);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 157, 2048, 17179869184, 2, 2, 34);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 12, 4096, 256, 9, 2, 8);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 13, 4096, 512, 8, 2, 9);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 16, 4096, 1024, 7, 2, 10);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 19, 4096, 4096, 6, 2, 12);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 24, 4096, 16384, 5, 2, 14);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 36, 4096, 131072, 4, 2, 17);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 64, 4096, 8388608, 3, 2, 23);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 157, 4096, 17179869184, 2, 2, 34);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 12, 8192, 256, 9, 2, 8);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 14, 8192, 512, 8, 2, 9);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 17, 8192, 1024, 7, 2, 10);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 20, 8192, 4096, 6, 2, 12);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 25, 8192, 16384, 5, 2, 14);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 37, 8192, 131072, 4, 2, 17);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 66, 8192, 8388608, 3, 2, 23);
		run_single_frog_degree_three_non_scalar_benchmark!(&mut $group, 1, 157, 8192, 17179869184, 2, 2, 34);
   };
}
