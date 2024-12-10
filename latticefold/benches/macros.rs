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
		run_single_babybear_benchmarks(&mut $group, 1, 3, 512, 16, 8, 2, 4);
		run_single_babybear_benchmarks(&mut $group, 1, 4, 512, 32, 7, 2, 5);
		run_single_babybear_benchmarks(&mut $group, 1, 6, 512, 256, 4, 2, 8);
		run_single_babybear_benchmarks(&mut $group, 1, 10, 512, 1024, 4, 2, 10);
		run_single_babybear_benchmarks(&mut $group, 1, 15, 512, 32768, 3, 2, 15);
		run_single_babybear_benchmarks(&mut $group, 1, 4, 1024, 32, 7, 2, 5);
		run_single_babybear_benchmarks(&mut $group, 1, 5, 1024, 64, 6, 2, 6);
		run_single_babybear_benchmarks(&mut $group, 1, 6, 1024, 256, 4, 2, 8);
		run_single_babybear_benchmarks(&mut $group, 1, 10, 1024, 1024, 4, 2, 10);
		run_single_babybear_benchmarks(&mut $group, 1, 15, 1024, 32768, 3, 2, 15);
		run_single_babybear_benchmarks(&mut $group, 1, 4, 2048, 32, 7, 2, 5);
		run_single_babybear_benchmarks(&mut $group, 1, 5, 2048, 64, 6, 2, 6);
		run_single_babybear_benchmarks(&mut $group, 1, 7, 2048, 256, 4, 2, 8);
		run_single_babybear_benchmarks(&mut $group, 1, 10, 2048, 1024, 4, 2, 10);
		run_single_babybear_benchmarks(&mut $group, 1, 15, 2048, 32768, 3, 2, 15);
		run_single_babybear_benchmarks(&mut $group, 1, 4, 4096, 16, 8, 2, 4);
		run_single_babybear_benchmarks(&mut $group, 1, 5, 4096, 32, 7, 2, 5);
		run_single_babybear_benchmarks(&mut $group, 1, 7, 4096, 256, 4, 2, 8);
		run_single_babybear_benchmarks(&mut $group, 1, 11, 4096, 1024, 4, 2, 10);
		run_single_babybear_benchmarks(&mut $group, 1, 15, 4096, 32768, 3, 2, 15);
		run_single_babybear_benchmarks(&mut $group, 1, 4, 8192, 16, 8, 2, 4);
		run_single_babybear_benchmarks(&mut $group, 1, 5, 8192, 32, 7, 2, 5);
		run_single_babybear_benchmarks(&mut $group, 1, 7, 8192, 256, 4, 2, 8);
		run_single_babybear_benchmarks(&mut $group, 1, 11, 8192, 1024, 4, 2, 10);
		run_single_babybear_benchmarks(&mut $group, 1, 15, 8192, 32768, 3, 2, 15);
   };
}

#[macro_export]
macro_rules! run_babybear_non_scalar_benchmarks {
   ($group:ident) => {
        // Parameters: Criterion group, X_LEN, Kappa, W_CCS, B, L, b, k
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 3, 512, 16, 8, 2, 4);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 4, 512, 32, 7, 2, 5);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 6, 512, 256, 4, 2, 8);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 10, 512, 1024, 4, 2, 10);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 15, 512, 32768, 3, 2, 15);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 4, 1024, 32, 7, 2, 5);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 5, 1024, 64, 6, 2, 6);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 6, 1024, 256, 4, 2, 8);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 10, 1024, 1024, 4, 2, 10);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 15, 1024, 32768, 3, 2, 15);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 4, 2048, 32, 7, 2, 5);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 5, 2048, 64, 6, 2, 6);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 7, 2048, 256, 4, 2, 8);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 10, 2048, 1024, 4, 2, 10);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 15, 2048, 32768, 3, 2, 15);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 4, 4096, 16, 8, 2, 4);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 5, 4096, 32, 7, 2, 5);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 7, 4096, 256, 4, 2, 8);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 11, 4096, 1024, 4, 2, 10);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 15, 4096, 32768, 3, 2, 15);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 4, 8192, 16, 8, 2, 4);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 5, 8192, 32, 7, 2, 5);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 7, 8192, 256, 4, 2, 8);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 11, 8192, 1024, 4, 2, 10);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 15, 8192, 32768, 3, 2, 15);
   };
}

#[macro_export]
macro_rules! run_babybear_degree_three_non_scalar_benchmarks {
   ($group:ident) => {
        // Parameters: Criterion group, X_LEN, Kappa, W_CCS, B, L, b, k
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 3, 512, 16, 8, 2, 4);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 4, 512, 32, 7, 2, 5);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 6, 512, 256, 4, 2, 8);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 10, 512, 1024, 4, 2, 10);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 15, 512, 32768, 3, 2, 15);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 4, 1024, 32, 7, 2, 5);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 5, 1024, 64, 6, 2, 6);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 6, 1024, 256, 4, 2, 8);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 10, 1024, 1024, 4, 2, 10);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 15, 1024, 32768, 3, 2, 15);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 4, 2048, 32, 7, 2, 5);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 5, 2048, 64, 6, 2, 6);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 7, 2048, 256, 4, 2, 8);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 10, 2048, 1024, 4, 2, 10);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 15, 2048, 32768, 3, 2, 15);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 4, 4096, 16, 8, 2, 4);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 5, 4096, 32, 7, 2, 5);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 7, 4096, 256, 4, 2, 8);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 11, 4096, 1024, 4, 2, 10);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 15, 4096, 32768, 3, 2, 15);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 4, 8192, 16, 8, 2, 4);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 5, 8192, 32, 7, 2, 5);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 7, 8192, 256, 4, 2, 8);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 11, 8192, 1024, 4, 2, 10);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 15, 8192, 32768, 3, 2, 15);
   };
}

//--- Frog cyclotomic ring (modulus p = 159120925213255836417, degree = 16) ---
//	Maximum kappa for which bound_{l_2} < p/2: 157
#[macro_export]
macro_rules! run_frog_benchmarks {
   ($group:ident) => {
        // Parameters: Criterion group, X_LEN, Kappa, W_CCS, B, L, b, k
		run_single_frog_benchmarks(&mut $group, 1, 10, 512, 256, 9, 2, 8);
		run_single_frog_benchmarks(&mut $group, 1, 11, 512, 512, 8, 2, 9);
		run_single_frog_benchmarks(&mut $group, 1, 14, 512, 1024, 7, 2, 10);
		run_single_frog_benchmarks(&mut $group, 1, 16, 512, 4096, 6, 2, 12);
		run_single_frog_benchmarks(&mut $group, 1, 21, 512, 16384, 5, 2, 14);
		run_single_frog_benchmarks(&mut $group, 1, 32, 512, 131072, 4, 2, 17);
		run_single_frog_benchmarks(&mut $group, 1, 60, 512, 8388608, 3, 2, 23);
		run_single_frog_benchmarks(&mut $group, 1, 157, 512, 17179869184, 2, 2, 34);
		run_single_frog_benchmarks(&mut $group, 1, 11, 1024, 256, 9, 2, 8);
		run_single_frog_benchmarks(&mut $group, 1, 12, 1024, 512, 8, 2, 9);
		run_single_frog_benchmarks(&mut $group, 1, 14, 1024, 1024, 7, 2, 10);
		run_single_frog_benchmarks(&mut $group, 1, 17, 1024, 4096, 6, 2, 12);
		run_single_frog_benchmarks(&mut $group, 1, 22, 1024, 16384, 5, 2, 14);
		run_single_frog_benchmarks(&mut $group, 1, 34, 1024, 131072, 4, 2, 17);
		run_single_frog_benchmarks(&mut $group, 1, 61, 1024, 8388608, 3, 2, 23);
		run_single_frog_benchmarks(&mut $group, 1, 157, 1024, 17179869184, 2, 2, 34);
		run_single_frog_benchmarks(&mut $group, 1, 11, 2048, 256, 9, 2, 8);
		run_single_frog_benchmarks(&mut $group, 1, 12, 2048, 512, 8, 2, 9);
		run_single_frog_benchmarks(&mut $group, 1, 15, 2048, 1024, 7, 2, 10);
		run_single_frog_benchmarks(&mut $group, 1, 18, 2048, 4096, 6, 2, 12);
		run_single_frog_benchmarks(&mut $group, 1, 23, 2048, 16384, 5, 2, 14);
		run_single_frog_benchmarks(&mut $group, 1, 35, 2048, 131072, 4, 2, 17);
		run_single_frog_benchmarks(&mut $group, 1, 63, 2048, 8388608, 3, 2, 23);
		run_single_frog_benchmarks(&mut $group, 1, 157, 2048, 17179869184, 2, 2, 34);
		run_single_frog_benchmarks(&mut $group, 1, 12, 4096, 256, 9, 2, 8);
		run_single_frog_benchmarks(&mut $group, 1, 13, 4096, 512, 8, 2, 9);
		run_single_frog_benchmarks(&mut $group, 1, 16, 4096, 1024, 7, 2, 10);
		run_single_frog_benchmarks(&mut $group, 1, 19, 4096, 4096, 6, 2, 12);
		run_single_frog_benchmarks(&mut $group, 1, 24, 4096, 16384, 5, 2, 14);
		run_single_frog_benchmarks(&mut $group, 1, 36, 4096, 131072, 4, 2, 17);
		run_single_frog_benchmarks(&mut $group, 1, 64, 4096, 8388608, 3, 2, 23);
		run_single_frog_benchmarks(&mut $group, 1, 157, 4096, 17179869184, 2, 2, 34);
		run_single_frog_benchmarks(&mut $group, 1, 12, 8192, 256, 9, 2, 8);
		run_single_frog_benchmarks(&mut $group, 1, 14, 8192, 512, 8, 2, 9);
		run_single_frog_benchmarks(&mut $group, 1, 17, 8192, 1024, 7, 2, 10);
		run_single_frog_benchmarks(&mut $group, 1, 20, 8192, 4096, 6, 2, 12);
		run_single_frog_benchmarks(&mut $group, 1, 25, 8192, 16384, 5, 2, 14);
		run_single_frog_benchmarks(&mut $group, 1, 37, 8192, 131072, 4, 2, 17);
		run_single_frog_benchmarks(&mut $group, 1, 66, 8192, 8388608, 3, 2, 23);
		run_single_frog_benchmarks(&mut $group, 1, 157, 8192, 17179869184, 2, 2, 34);
   };
}

#[macro_export]
macro_rules! run_frog_non_scalar_benchmarks {
   ($group:ident) => {
        // Parameters: Criterion group, X_LEN, Kappa, W_CCS, B, L, b, k
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 10, 512, 256, 9, 2, 8);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 11, 512, 512, 8, 2, 9);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 14, 512, 1024, 7, 2, 10);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 16, 512, 4096, 6, 2, 12);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 21, 512, 16384, 5, 2, 14);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 32, 512, 131072, 4, 2, 17);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 60, 512, 8388608, 3, 2, 23);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 157, 512, 17179869184, 2, 2, 34);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 11, 1024, 256, 9, 2, 8);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 12, 1024, 512, 8, 2, 9);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 14, 1024, 1024, 7, 2, 10);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 17, 1024, 4096, 6, 2, 12);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 22, 1024, 16384, 5, 2, 14);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 34, 1024, 131072, 4, 2, 17);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 61, 1024, 8388608, 3, 2, 23);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 157, 1024, 17179869184, 2, 2, 34);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 11, 2048, 256, 9, 2, 8);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 12, 2048, 512, 8, 2, 9);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 15, 2048, 1024, 7, 2, 10);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 18, 2048, 4096, 6, 2, 12);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 23, 2048, 16384, 5, 2, 14);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 35, 2048, 131072, 4, 2, 17);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 63, 2048, 8388608, 3, 2, 23);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 157, 2048, 17179869184, 2, 2, 34);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 12, 4096, 256, 9, 2, 8);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 13, 4096, 512, 8, 2, 9);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 16, 4096, 1024, 7, 2, 10);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 19, 4096, 4096, 6, 2, 12);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 24, 4096, 16384, 5, 2, 14);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 36, 4096, 131072, 4, 2, 17);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 64, 4096, 8388608, 3, 2, 23);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 157, 4096, 17179869184, 2, 2, 34);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 12, 8192, 256, 9, 2, 8);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 14, 8192, 512, 8, 2, 9);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 17, 8192, 1024, 7, 2, 10);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 20, 8192, 4096, 6, 2, 12);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 25, 8192, 16384, 5, 2, 14);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 37, 8192, 131072, 4, 2, 17);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 66, 8192, 8388608, 3, 2, 23);
		run_single_frog_non_scalar_benchmarks(&mut $group, 1, 157, 8192, 17179869184, 2, 2, 34);
   };
}

#[macro_export]
macro_rules! run_frog_degree_three_non_scalar_benchmarks {
   ($group:ident) => {
        // Parameters: Criterion group, X_LEN, Kappa, W_CCS, B, L, b, k
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 10, 512, 256, 9, 2, 8);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 11, 512, 512, 8, 2, 9);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 14, 512, 1024, 7, 2, 10);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 16, 512, 4096, 6, 2, 12);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 21, 512, 16384, 5, 2, 14);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 32, 512, 131072, 4, 2, 17);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 60, 512, 8388608, 3, 2, 23);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 157, 512, 17179869184, 2, 2, 34);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 11, 1024, 256, 9, 2, 8);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 12, 1024, 512, 8, 2, 9);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 14, 1024, 1024, 7, 2, 10);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 17, 1024, 4096, 6, 2, 12);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 22, 1024, 16384, 5, 2, 14);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 34, 1024, 131072, 4, 2, 17);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 61, 1024, 8388608, 3, 2, 23);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 157, 1024, 17179869184, 2, 2, 34);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 11, 2048, 256, 9, 2, 8);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 12, 2048, 512, 8, 2, 9);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 15, 2048, 1024, 7, 2, 10);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 18, 2048, 4096, 6, 2, 12);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 23, 2048, 16384, 5, 2, 14);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 35, 2048, 131072, 4, 2, 17);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 63, 2048, 8388608, 3, 2, 23);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 157, 2048, 17179869184, 2, 2, 34);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 12, 4096, 256, 9, 2, 8);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 13, 4096, 512, 8, 2, 9);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 16, 4096, 1024, 7, 2, 10);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 19, 4096, 4096, 6, 2, 12);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 24, 4096, 16384, 5, 2, 14);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 36, 4096, 131072, 4, 2, 17);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 64, 4096, 8388608, 3, 2, 23);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 157, 4096, 17179869184, 2, 2, 34);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 12, 8192, 256, 9, 2, 8);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 14, 8192, 512, 8, 2, 9);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 17, 8192, 1024, 7, 2, 10);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 20, 8192, 4096, 6, 2, 12);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 25, 8192, 16384, 5, 2, 14);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 37, 8192, 131072, 4, 2, 17);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 66, 8192, 8388608, 3, 2, 23);
		run_single_frog_degree_three_non_scalar_benchmarks(&mut $group, 1, 157, 8192, 17179869184, 2, 2, 34);
   };
}
