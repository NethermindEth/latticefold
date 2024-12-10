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
        // Parameters: Criteion group, X_LEN, Kappa, W_CCS, B, L, b, k
		run_single_babybear_benchmarks(&mut $group, 1, 3, 512, 16, 8, 2, 4);
		run_single_babybear_benchmarks(&mut $group, 1, 4, 512, 32, 7, 2, 5);
		run_single_babybear_benchmarks(&mut $group, 1, 5, 512, 256, 4, 2, 8);
		run_single_babybear_benchmarks(&mut $group, 1, 7, 512, 1024, 4, 2, 10);
		run_single_babybear_benchmarks(&mut $group, 1, 11, 512, 32768, 3, 2, 15);
		run_single_babybear_benchmarks(&mut $group, 1, 4, 1024, 32, 7, 2, 5);
		run_single_babybear_benchmarks(&mut $group, 1, 5, 1024, 64, 6, 2, 6);
		run_single_babybear_benchmarks(&mut $group, 1, 6, 1024, 256, 4, 2, 8);
		run_single_babybear_benchmarks(&mut $group, 1, 7, 1024, 1024, 4, 2, 10);
		run_single_babybear_benchmarks(&mut $group, 1, 11, 1024, 32768, 3, 2, 15);
		run_single_babybear_benchmarks(&mut $group, 1, 4, 2048, 32, 7, 2, 5);
		run_single_babybear_benchmarks(&mut $group, 1, 5, 2048, 64, 6, 2, 6);
		run_single_babybear_benchmarks(&mut $group, 1, 6, 2048, 256, 4, 2, 8);
		run_single_babybear_benchmarks(&mut $group, 1, 8, 2048, 1024, 4, 2, 10);
		run_single_babybear_benchmarks(&mut $group, 1, 11, 2048, 32768, 3, 2, 15);
		run_single_babybear_benchmarks(&mut $group, 1, 4, 4096, 16, 8, 2, 4);
		run_single_babybear_benchmarks(&mut $group, 1, 5, 4096, 32, 7, 2, 5);
		run_single_babybear_benchmarks(&mut $group, 1, 6, 4096, 256, 4, 2, 8);
		run_single_babybear_benchmarks(&mut $group, 1, 8, 4096, 1024, 4, 2, 10);
		run_single_babybear_benchmarks(&mut $group, 1, 12, 4096, 32768, 3, 2, 15);
		run_single_babybear_benchmarks(&mut $group, 1, 4, 8192, 16, 8, 2, 4);
		run_single_babybear_benchmarks(&mut $group, 1, 5, 8192, 32, 7, 2, 5);
		run_single_babybear_benchmarks(&mut $group, 1, 6, 8192, 256, 4, 2, 8);
		run_single_babybear_benchmarks(&mut $group, 1, 8, 8192, 1024, 4, 2, 10);
		run_single_babybear_benchmarks(&mut $group, 1, 12, 8192, 32768, 3, 2, 15);
   };
}


#[macro_export]
macro_rules! run_babybear_non_scalar_benchmarks {
   ($group:ident) => {
        // Parameters: Criteion group, X_LEN, Kappa, W_CCS, B, L, b, k
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 3, 512, 16, 8, 2, 4);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 4, 512, 32, 7, 2, 5);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 5, 512, 256, 4, 2, 8);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 7, 512, 1024, 4, 2, 10);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 11, 512, 32768, 3, 2, 15);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 4, 1024, 32, 7, 2, 5);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 5, 1024, 64, 6, 2, 6);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 6, 1024, 256, 4, 2, 8);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 7, 1024, 1024, 4, 2, 10);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 11, 1024, 32768, 3, 2, 15);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 4, 2048, 32, 7, 2, 5);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 5, 2048, 64, 6, 2, 6);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 6, 2048, 256, 4, 2, 8);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 8, 2048, 1024, 4, 2, 10);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 11, 2048, 32768, 3, 2, 15);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 4, 4096, 16, 8, 2, 4);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 5, 4096, 32, 7, 2, 5);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 6, 4096, 256, 4, 2, 8);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 8, 4096, 1024, 4, 2, 10);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 12, 4096, 32768, 3, 2, 15);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 4, 8192, 16, 8, 2, 4);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 5, 8192, 32, 7, 2, 5);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 6, 8192, 256, 4, 2, 8);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 8, 8192, 1024, 4, 2, 10);
		run_single_babybear_non_scalar_benchmarks(&mut $group, 1, 12, 8192, 32768, 3, 2, 15);
   };
}


#[macro_export]
macro_rules! run_babybear_degree_three_non_scalar_benchmarks {
   ($group:ident) => {
        // Parameters: Criteion group, X_LEN, Kappa, W_CCS, B, L, b, k
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 3, 512, 16, 8, 2, 4);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 4, 512, 32, 7, 2, 5);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 5, 512, 256, 4, 2, 8);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 7, 512, 1024, 4, 2, 10);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 11, 512, 32768, 3, 2, 15);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 4, 1024, 32, 7, 2, 5);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 5, 1024, 64, 6, 2, 6);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 6, 1024, 256, 4, 2, 8);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 7, 1024, 1024, 4, 2, 10);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 11, 1024, 32768, 3, 2, 15);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 4, 2048, 32, 7, 2, 5);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 5, 2048, 64, 6, 2, 6);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 6, 2048, 256, 4, 2, 8);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 8, 2048, 1024, 4, 2, 10);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 11, 2048, 32768, 3, 2, 15);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 4, 4096, 16, 8, 2, 4);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 5, 4096, 32, 7, 2, 5);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 6, 4096, 256, 4, 2, 8);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 8, 4096, 1024, 4, 2, 10);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 12, 4096, 32768, 3, 2, 15);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 4, 8192, 16, 8, 2, 4);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 5, 8192, 32, 7, 2, 5);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 6, 8192, 256, 4, 2, 8);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 8, 8192, 1024, 4, 2, 10);
		run_single_babybear_degree_three_non_scalar_benchmarks(&mut $group, 1, 12, 8192, 32768, 3, 2, 15);
   };
}

