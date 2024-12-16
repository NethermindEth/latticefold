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

         run_single_goldilocks_benchmark!(&mut $group, 2, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 2, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 2, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 2, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 2, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 2, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 2, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 2, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 2, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 2, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 2, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 2, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 2, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 2, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 2, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 2, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 2, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 2, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 2, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 2, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 2, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 3, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 3, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 3, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 3, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 3, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 3, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 3, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 3, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 3, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 3, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 3, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 3, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 3, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 3, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 3, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 3, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 3, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 3, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 3, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 3, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 3, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 4, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 5, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 5, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 5, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 5, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 5, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 5, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 5, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 5, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 5, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 5, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 5, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 5, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 5, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 5, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 5, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 5, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 5, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 5, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 5, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 5, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 5, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 6, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 6, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 6, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 6, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 6, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 6, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 6, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 6, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 6, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 6, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 6, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 6, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 6, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 6, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 6, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 6, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 6, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 6, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 6, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 6, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 6, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 7, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 7, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 7, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 7, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 7, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 7, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 7, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 7, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 7, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 7, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 7, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 7, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 7, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 7, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 7, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 7, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 7, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 7, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 7, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 7, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 7, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 8, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 8, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 8, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 8, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 8, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 8, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 8, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 8, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 8, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 8, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 8, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 8, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 8, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 8, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 8, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 8, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 8, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 8, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 8, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 8, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 8, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 9, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 9, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 9, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 9, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 9, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 9, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 9, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 9, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 9, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 9, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 9, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 9, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 9, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 9, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 9, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 9, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 9, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 9, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 9, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 9, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 9, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 10, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 10, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 10, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 10, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 10, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 10, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 10, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 10, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 10, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 10, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 10, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 10, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 10, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 10, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 10, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 10, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 10, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 10, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 10, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 10, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 10, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 11, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 11, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 11, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 11, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 11, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 11, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 11, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 11, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 11, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 11, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 11, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 11, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 11, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 11, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 11, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 11, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 11, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 11, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 11, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 11, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 11, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 12, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 12, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 12, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 12, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 12, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 12, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 12, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 12, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 12, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 12, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 12, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 12, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 12, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 12, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 12, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 12, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 12, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 12, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 12, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 12, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 12, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 13, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 13, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 13, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 13, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 13, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 13, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 13, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 13, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 13, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 13, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 13, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 13, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 13, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 13, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 13, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 13, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 13, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 13, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 13, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 13, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 13, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 14, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 14, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 14, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 14, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 14, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 14, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 14, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 14, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 14, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 14, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 14, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 14, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 14, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 14, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 14, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 14, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 14, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 14, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 14, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 14, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 14, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 15, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 15, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 15, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 15, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 15, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 15, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 15, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 15, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 15, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 15, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 15, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 15, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 15, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 15, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 15, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 15, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 15, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 15, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 15, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 15, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 15, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 16, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 16, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 16, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 16, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 16, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 16, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 16, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 16, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 16, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 16, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 16, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 16, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 16, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 16, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 16, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 16, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 16, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 16, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 16, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 16, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 16, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 17, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 17, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 17, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 17, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 17, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 17, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 17, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 17, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 17, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 17, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 17, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 17, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 17, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 17, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 17, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 17, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 17, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 17, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 17, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 17, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 17, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 18, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 18, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 18, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 18, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 18, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 18, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 18, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 18, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 18, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 18, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 18, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 18, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 18, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 18, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 18, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 18, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 18, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 18, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 18, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 18, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 18, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 19, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 19, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 19, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 19, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 19, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 19, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 19, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 19, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 19, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 19, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 19, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 19, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 19, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 19, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 19, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 19, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 19, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 19, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 19, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 19, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 19, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 20, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 20, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 20, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 20, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 20, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 20, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 20, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 20, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 20, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 20, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 20, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 20, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 20, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 20, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 20, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 20, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 20, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 20, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 20, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 20, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 20, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 21, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 21, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 21, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 21, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 21, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 21, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 21, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 21, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 21, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 21, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 21, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 21, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 21, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 21, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 21, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 21, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 21, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 21, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 21, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 21, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 21, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 22, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 22, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 22, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 22, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 22, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 22, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 22, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 22, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 22, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 22, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 22, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 22, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 22, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 22, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 22, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 22, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 22, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 22, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 22, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 22, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 22, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 23, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 23, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 23, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 23, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 23, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 23, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 23, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 23, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 23, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 23, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 23, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 23, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 23, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 23, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 23, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 23, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 23, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 23, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 23, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 23, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 23, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 24, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 24, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 24, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 24, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 24, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 24, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 24, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 24, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 24, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 24, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 24, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 24, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 24, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 24, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 24, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 24, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 24, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 24, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 24, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 24, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 24, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 25, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 25, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 25, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 25, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 25, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 25, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 25, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 25, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 25, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 25, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 25, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 25, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 25, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 25, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 25, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 25, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 25, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 25, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 25, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 25, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 25, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 26, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 26, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 26, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 26, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 26, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 26, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 26, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 26, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 26, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 26, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 26, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 26, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 26, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 26, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 26, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 26, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 26, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 26, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 26, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 26, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 26, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 27, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 27, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 27, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 27, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 27, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 27, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 27, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 27, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 27, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 27, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 27, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 27, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 27, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 27, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 27, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 27, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 27, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 27, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 27, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 27, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 27, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 28, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 28, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 28, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 28, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 28, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 28, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 28, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 28, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 28, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 28, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 28, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 28, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 28, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 28, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 28, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 28, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 28, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 28, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 28, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 28, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 28, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 29, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 29, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 29, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 29, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 29, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 29, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 29, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 29, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 29, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 29, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 29, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 29, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 29, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 29, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 29, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 29, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 29, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 29, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 29, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 29, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 29, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 30, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 30, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 30, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 30, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 30, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 30, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 30, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 30, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 30, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 30, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 30, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 30, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 30, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 30, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 30, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 30, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 30, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 30, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 30, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 30, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 30, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 31, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 31, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 31, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 31, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 31, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 31, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 31, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 31, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 31, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 31, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 31, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 31, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 31, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 31, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 31, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 31, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 31, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 31, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 31, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 31, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 31, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 32, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 32, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 32, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 32, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 32, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 32, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 32, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 32, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 32, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 32, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 32, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 32, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 32, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 32, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 32, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 32, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 32, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 32, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 32, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 32, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 32, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 33, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 33, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 33, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 33, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 33, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 33, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 33, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 33, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 33, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 33, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 33, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 33, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 33, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 33, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 33, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 33, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 33, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 33, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 33, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 33, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 33, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 34, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 34, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 34, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 34, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 34, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 34, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 34, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 34, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 34, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 34, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 34, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 34, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 34, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 34, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 34, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 34, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 34, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 34, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 34, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 34, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 34, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 35, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 35, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 35, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 35, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 35, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 35, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 35, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 35, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 35, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 35, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 35, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 35, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 35, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 35, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 35, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 35, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 35, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 35, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 35, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 35, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 35, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 36, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 36, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 36, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 36, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 36, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 36, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 36, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 36, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 36, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 36, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 36, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 36, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 36, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 36, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 36, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 36, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 36, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 36, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 36, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 36, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 36, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 37, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 37, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 37, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 37, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 37, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 37, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 37, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 37, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 37, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 37, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 37, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 37, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 37, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 37, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 37, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 37, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 37, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 37, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 37, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 37, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 37, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 38, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 38, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 38, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 38, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 38, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 38, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 38, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 38, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 38, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 38, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 38, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 38, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 38, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 38, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 38, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 38, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 38, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 38, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 38, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 38, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 38, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 39, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 39, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 39, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 39, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 39, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 39, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 39, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 39, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 39, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 39, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 39, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 39, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 39, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 39, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 39, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 39, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 39, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 39, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 39, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 39, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 39, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 40, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 40, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 40, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 40, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 40, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 40, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 40, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 40, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 40, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 40, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 40, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 40, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 40, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 40, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 40, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 40, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 40, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 40, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 40, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 40, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 40, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 41, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 41, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 41, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 41, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 41, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 41, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 41, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 41, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 41, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 41, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 41, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 41, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 41, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 41, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 41, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 41, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 41, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 41, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 41, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 41, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 41, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 42, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 42, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 42, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 42, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 42, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 42, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 42, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 42, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 42, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 42, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 42, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 42, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 42, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 42, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 42, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 42, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 42, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 42, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 42, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 42, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 42, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 43, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 43, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 43, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 43, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 43, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 43, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 43, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 43, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 43, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 43, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 43, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 43, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 43, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 43, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 43, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 43, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 43, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 43, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 43, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 43, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 43, {1 << 20});

         run_single_goldilocks_benchmark!(&mut $group, 4, {44 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 4, {44 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 4, {44 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 4, {44 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 4, {44 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 4, {44 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 4, {44 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 4, {44 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 4, {44 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 4, {44 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 4, {44 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 4, {44 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 4, {44 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 4, {44 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 4, {44 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 4, {44 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 4, {44 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 4, {44 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 4, {44 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 4, {44 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 4, {44 << 20});


         run_single_goldilocks_benchmark!(&mut $group, 45, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 45, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 45, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 45, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 45, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 45, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 45, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 45, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 45, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 45, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 45, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 45, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 45, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 45, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 45, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 45, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 45, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 45, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 45, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 45, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 45, {1 << 20});


         run_single_goldilocks_benchmark!(&mut $group, 46, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 46, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 46, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 46, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 46, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 46, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 46, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 46, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 46, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 46, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 46, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 46, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 46, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 46, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 46, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 46, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 46, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 46, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 46, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 46, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 46, {1 << 20});


         run_single_goldilocks_benchmark!(&mut $group, 47, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 47, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 47, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 47, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 47, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 47, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 47, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 47, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 47, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 47, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 47, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 47, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 47, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 47, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 47, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 47, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 47, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 47, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 47, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 47, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 47, {1 << 20});


         run_single_goldilocks_benchmark!(&mut $group, 48, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 48, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 48, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 48, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 48, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 48, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 48, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 48, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 48, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 48, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 48, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 48, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 48, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 48, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 48, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 48, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 48, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 48, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 48, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 48, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 48, {1 << 20});


         run_single_goldilocks_benchmark!(&mut $group, 49, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 49, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 49, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 49, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 49, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 49, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 49, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 49, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 49, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 49, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 49, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 49, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 49, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 49, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 49, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 49, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 49, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 49, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 49, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 49, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 49, {1 << 20});


         run_single_goldilocks_benchmark!(&mut $group, 50, {1 << 0});
         run_single_goldilocks_benchmark!(&mut $group, 50, {1 << 1});
         run_single_goldilocks_benchmark!(&mut $group, 50, {1 << 2});
         run_single_goldilocks_benchmark!(&mut $group, 50, {1 << 3});
         run_single_goldilocks_benchmark!(&mut $group, 50, {1 << 4});
         run_single_goldilocks_benchmark!(&mut $group, 50, {1 << 5});
         run_single_goldilocks_benchmark!(&mut $group, 50, {1 << 6});
         run_single_goldilocks_benchmark!(&mut $group, 50, {1 << 7});
         run_single_goldilocks_benchmark!(&mut $group, 50, {1 << 8});
         run_single_goldilocks_benchmark!(&mut $group, 50, {1 << 9});
         run_single_goldilocks_benchmark!(&mut $group, 50, {1 << 10});
         run_single_goldilocks_benchmark!(&mut $group, 50, {1 << 11});
         run_single_goldilocks_benchmark!(&mut $group, 50, {1 << 12});
         run_single_goldilocks_benchmark!(&mut $group, 50, {1 << 13});
         run_single_goldilocks_benchmark!(&mut $group, 50, {1 << 14});
         run_single_goldilocks_benchmark!(&mut $group, 50, {1 << 15});
         run_single_goldilocks_benchmark!(&mut $group, 50, {1 << 16});
         run_single_goldilocks_benchmark!(&mut $group, 50, {1 << 17});
         run_single_goldilocks_benchmark!(&mut $group, 50, {1 << 18});
         run_single_goldilocks_benchmark!(&mut $group, 50, {1 << 19});
         run_single_goldilocks_benchmark!(&mut $group, 50, {1 << 20});


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