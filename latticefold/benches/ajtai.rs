#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

use ark_ff::Field;
use ark_ff::PrimeField;
use ark_std::{time::Duration, UniformRand};
use criterion::{
    criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion, PlotConfiguration,
};
use cyclotomic_rings::{BabyBearRingNTT, GoldilocksRingNTT, StarkRingNTT, FrogRingNTT, SuitableRing};
use lattirust_ring::cyclotomic_ring::models::pow2_debug::Pow2CyclotomicPolyRingNTT;
use latticefold::{commitment::AjtaiCommitmentScheme, parameters::{DecompositionParams, DILITHIUM_PRIME}};
use lattirust_ring::PolyRing;
use rand::{thread_rng, RngCore};
use std::fmt::Debug;

fn draw_bellow_bound<R: SuitableRing, Rng>(rng: &mut Rng, bound: u128, degree: usize) -> R
where
    Rng: rand::Rng + ?Sized,
{
    let bound_as_field = <<R as PolyRing>::BaseRing as Field>::BasePrimeField::from(bound - 1);
    let coeffs = vec![bound_as_field; degree];
    let mut poly = R::CoefficientRepresentation::from(coeffs);
    while !all_elements_bellow_bound::<R>(&poly, bound) {
        poly = R::CoefficientRepresentation::rand(rng);
    }
    R::from(poly)
}

fn all_elements_bellow_bound<R: SuitableRing>(
    poly: &R::CoefficientRepresentation,
    bound: u128,
) -> bool {
    let coeffs = poly.coeffs();
    let big_int = coeffs
        .iter()
        .map(|coeff: &<<R as PolyRing>::BaseRing as Field>::BasePrimeField| coeff.into_bigint())
        .max()
        .unwrap();
    let bound = <<R as PolyRing>::BaseRing as Field>::BasePrimeField::from(bound).into_bigint();
    big_int < bound
}

fn ajtai_benchmark<
    const C: usize, // rows
    const W: usize, // columns
    R: Clone + UniformRand + Debug + SuitableRing + for<'a> std::ops::AddAssign<&'a R>,
>(
    group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>,
) {
    let mut rng = thread_rng();

    let witness: Vec<R> = (0..W).map(|_| R::rand(&mut rng)).collect();
    let witness_2 = witness.clone();
    let witness_3 = witness.clone();
    let ajtai_data: AjtaiCommitmentScheme<C, W, R> = AjtaiCommitmentScheme::rand(&mut rng);

    group.bench_with_input(
        BenchmarkId::new("CommitNTT", format!("C={}, W={}", C, W)),
        &(ajtai_data.clone(), witness),
        |b, (ajtai_data, witness)| {
            b.iter(|| {
                let _ = ajtai_data.commit_ntt(witness);
            })
        },
    );

    // NTT -> INTT (coefficients)
    group.bench_with_input(
        BenchmarkId::new("NTT->INTT", format!("C={}, W={}", C, W)),
        &(witness_2),
        |b, witness| {
            b.iter(|| {
                let _ = witness
                    .iter()
                    .map(|&x| x.into())
                    .collect::<Vec<R::CoefficientRepresentation>>();
            })
        },
    );

    // INTT -> NTT
    let coeff = witness_3
        .iter()
        .map(|&x| x.into())
        .collect::<Vec<R::CoefficientRepresentation>>();
    group.bench_with_input(
        BenchmarkId::new("INTT->NTT", format!("C={}, W={}", C, W)),
        &(coeff),
        |b, coeff| {
            b.iter(|| {
                let _: Vec<R> = coeff.iter().map(|&x| R::from(x)).collect();
            })
        },
    );
}

macro_rules! define_starkprime_params {
    ($w:expr, $b:expr, $l:expr) => {
        paste::paste! {
            #[derive(Clone)]
            struct [<StarkPrimeParamsWithB $b W $w>];

            impl DecompositionParams for [<StarkPrimeParamsWithB $b W $w>] {
                const B: u128 = $b;
                const L: usize = $l;
                const B_SMALL: usize = 2; // This is not use in decompose or commit
                const K: usize = 28;// This is not use in decompose or commit
            }
        }
    };
}

macro_rules! run_single_starkprime_benchmark {
    ($crit:expr, $cw:expr, $w:expr, $b:expr, $l:expr) => {
        define_starkprime_params!($w, $b, $l);
        paste::paste! {
            //const [<W $w B $b L $l>]: usize = $w * [<StarkPrimeParamsWithB $b W $w>]::L; // Define the padded witness
            //ajtai_benchmark::<$cw, $w, StarkRingNTT, [<StarkPrimeParamsWithB $b W $w>]>($crit);
            ajtai_benchmark::<$cw, $w, StarkRingNTT>($crit);            
        }
    };
}

macro_rules! define_goldilocks_params {
    ($w:expr, $b:expr, $l:expr) => {
        paste::paste! {
            #[derive(Clone)]
            struct [<GoldilocksParamsWithB $b W $w>];

            impl DecompositionParams for [<GoldilocksParamsWithB $b W $w>] {
                const B: u128 = $b;
                const L: usize = $l;
                const B_SMALL: usize = 2; // This is not use in decompose or commit
                const K: usize = 28;// This is not use in decompose or commit
            }
        }
    };
}

macro_rules! run_single_goldilocks_benchmark {
    ($crit:expr, $cw:expr, $w:expr, $b:expr, $l:expr) => {
        define_goldilocks_params!($w, $b, $l);
        paste::paste! {
            //const [<W $w B $b L $l>]: usize = $w * [<GoldilocksParamsWithB $b W $w>]::L; // Define the padded witness
            //ajtai_benchmark::<$cw, $w, GoldilocksRingNTT, [<GoldilocksParamsWithB $b W $w>]>($crit);
            ajtai_benchmark::<$cw, $w, GoldilocksRingNTT>($crit);
            
        }
    };
}

macro_rules! define_babybear_params {
    ($w:expr, $b:expr, $l:expr) => {
        paste::paste! {
            #[derive(Clone)]
            struct [<BabyBearParamsWithB $b W $w>];

            impl DecompositionParams for [<BabyBearParamsWithB $b W $w>] {
                const B: u128 = $b;
                const L: usize = $l;
                const B_SMALL: usize = 2; // This is not use in decompose or commit
                const K: usize = 28;// This is not use in decompose or commit
            }
        }
    };
}

macro_rules! run_single_babybear_benchmark {
    ($crit:expr, $cw:expr, $w:expr, $b:expr, $l:expr) => {
        define_babybear_params!($w, $b, $l);
        paste::paste! {
            //const [<W $w B $b L $l>]: usize = $w * [<BabyBearParamsWithB $b W $w>]::L; // Define the padded witness
            //ajtai_benchmark::<$cw, $w, BabyBearRingNTT, [<BabyBearParamsWithB $b W $w>]>($crit);
            ajtai_benchmark::<$cw, $w, BabyBearRingNTT>($crit);
            
        }
    };
}

macro_rules! define_frog_params {
    ($w:expr, $b:expr, $l:expr) => {
        paste::paste! {
            #[derive(Clone)]
            struct [<FrogParamsWithB $b W $w>];

            impl DecompositionParams for [<FrogParamsWithB $b W $w>] {
                const B: u128 = $b;
                const L: usize = $l;
                const B_SMALL: usize = 2; // This is not use in decompose or commit
                const K: usize = 28;// This is not use in decompose or commit
            }
        }
    };
}

macro_rules! run_single_frog_benchmark {
    ($crit:expr, $cw:expr, $w:expr, $b:expr, $l:expr) => {
        define_frog_params!($w, $b, $l);
        paste::paste! {
            //const [<W $w B $b L $l>]: usize = $w * [<FrogParamsWithB $b W $w>]::L; // Define the padded witness
            //ajtai_benchmark::<$cw, $w, FrogRingNTT, [<FrogParamsWithB $b W $w>]>($crit);
            ajtai_benchmark::<$cw, $w, FrogRingNTT>($crit);
        }
    };
}

macro_rules! define_dilithium_params {
    ($w:expr, $b:expr, $l:expr) => {
        paste::paste! {
            #[derive(Clone)]
            struct [<DilithiumParamsWithB $b W $w>];

            impl DecompositionParams for [<DilithiumParamsWithB $b W $w>] {
                const B: u128 = $b;
                const L: usize = $l;
                const B_SMALL: usize = 2; // This is not use in decompose or commit
                const K: usize = 28;// This is not use in decompose or commit
            }
        }
    };
}

macro_rules! run_single_dilithium_benchmark {
    ($crit:expr, $cw:expr, $w:expr, $b:expr, $l:expr) => {
        define_dilithium_params!($w, $b, $l);
        paste::paste! {
            //const [<W $w B $b L $l>]: usize = $w * [<DilithiumParamsWithB $b W $w>]::L; // Define the padded witness
            //ajtai_benchmark::<$cw, $w, Pow2CyclotomicPolyRingNTT<DILITHIUM_PRIME, 256> , [<DilithiumParamsWithB $b W $w>]>($crit);
            ajtai_benchmark::<$cw, $w, Pow2CyclotomicPolyRingNTT<DILITHIUM_PRIME, 256> >($crit);            
        }
    };
}

fn ajtai_benchmarks(c: &mut Criterion) {
    // Parameters are describe in the order C, W, B, L
    // Where:
    //  p: prime modulus
    //  C: number of columns
    //  W: witness length
    //  B: biggest even number less than B_infty from 128 bits of security ` 2^{ 2 * sqrt{ log(1.01) * D * C * log(p) } }/sqrt{ D * W }`
    //  D: Ring degree
    //  L: smallest int such that B^L > p This must be even as well?
    
    // Babybear
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Ajtai BabyBear");
        group.plot_config(plot_config.clone());
        
        run_single_babybear_benchmark!(&mut group, 1, 32768, 2, 31);
        run_single_babybear_benchmark!(&mut group, 2, 32768, 42, 6);
        run_single_babybear_benchmark!(&mut group, 3, 32768, 508, 4);
        run_single_babybear_benchmark!(&mut group, 4, 32768, 4160, 3);
        run_single_babybear_benchmark!(&mut group, 5, 32768, 26458, 3);
        run_single_babybear_benchmark!(&mut group, 6, 32768, 140862, 2);

        run_single_babybear_benchmark!(&mut group, 1, 65536, 2, 31);
        run_single_babybear_benchmark!(&mut group, 2, 65536, 28, 7);
        run_single_babybear_benchmark!(&mut group, 3, 65536, 360, 4);
        run_single_babybear_benchmark!(&mut group, 4, 65536, 2942, 3);
        run_single_babybear_benchmark!(&mut group, 5, 65536, 18708, 3);
        run_single_babybear_benchmark!(&mut group, 6, 65536, 99604, 2);

        run_single_babybear_benchmark!(&mut group, 1, 131072, 0, 31); // Bound not practical
        run_single_babybear_benchmark!(&mut group, 2, 131072, 20, 8);
        run_single_babybear_benchmark!(&mut group, 3, 131072, 254, 4);
        run_single_babybear_benchmark!(&mut group, 4, 131072, 2080, 3);
        run_single_babybear_benchmark!(&mut group, 5, 131072, 13228, 3);
        run_single_babybear_benchmark!(&mut group, 6, 131072, 70430, 2);

        run_single_babybear_benchmark!(&mut group, 1, 262144, 0, 31); // Bound not practical
        run_single_babybear_benchmark!(&mut group, 2, 262144, 14, 9);
        run_single_babybear_benchmark!(&mut group, 3, 262144, 180, 5);
        run_single_babybear_benchmark!(&mut group, 4, 262144, 1470, 3);
        run_single_babybear_benchmark!(&mut group, 5, 262144, 9354, 3);
        run_single_babybear_benchmark!(&mut group, 6, 262144, 49802, 2);

        run_single_babybear_benchmark!(&mut group, 1, 524288, 0, 31); // Bound not practical
        run_single_babybear_benchmark!(&mut group, 2, 524288, 10, 10);
        run_single_babybear_benchmark!(&mut group, 3, 524288, 126, 5);
        run_single_babybear_benchmark!(&mut group, 4, 524288, 1040, 4);
        run_single_babybear_benchmark!(&mut group, 5, 524288, 6614, 3);
        run_single_babybear_benchmark!(&mut group, 6, 524288, 35214, 3);

        run_single_babybear_benchmark!(&mut group, 1, 1048576, 0, 31); // Bound not practical
        run_single_babybear_benchmark!(&mut group, 2, 1048576, 6, 12);
        run_single_babybear_benchmark!(&mut group, 3, 1048576, 90, 5);
        run_single_babybear_benchmark!(&mut group, 4, 1048576, 734, 4);
        run_single_babybear_benchmark!(&mut group, 5, 1048576, 4676, 3);
        run_single_babybear_benchmark!(&mut group, 6, 1048576, 24900, 3);

        group.finish();
    }

    // Goldilocks
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Ajtai Goldilocks");
        group.plot_config(plot_config.clone());

        run_single_goldilocks_benchmark!(&mut group, 1, 32768, 0, 0); // Impractical bound
        run_single_goldilocks_benchmark!(&mut group, 2, 32768, 12, 18);
        run_single_goldilocks_benchmark!(&mut group, 3, 32768, 88, 10);
        run_single_goldilocks_benchmark!(&mut group, 4, 32768, 508, 8);
        run_single_goldilocks_benchmark!(&mut group, 5, 32768, 2364, 6);
        run_single_goldilocks_benchmark!(&mut group, 6, 32768, 9486, 5);

        run_single_goldilocks_benchmark!(&mut group, 1, 65536, 0, 0); // Impractical bound
        run_single_goldilocks_benchmark!(&mut group, 2, 65536, 8, 22);
        run_single_goldilocks_benchmark!(&mut group, 3, 65536, 62, 11);
        run_single_goldilocks_benchmark!(&mut group, 4, 65536, 360, 8);
        run_single_goldilocks_benchmark!(&mut group, 5, 65536, 1672, 6);
        run_single_goldilocks_benchmark!(&mut group, 6, 65536, 6708, 6);

        run_single_goldilocks_benchmark!(&mut group, 1, 131072, 0, 31); // Bound not practical
        run_single_goldilocks_benchmark!(&mut group, 2, 131072, 6, 25);
        run_single_goldilocks_benchmark!(&mut group, 3, 131072, 44, 12);
        run_single_goldilocks_benchmark!(&mut group, 4, 131072, 254, 9);
        run_single_goldilocks_benchmark!(&mut group, 5, 131072, 1182, 7);
        run_single_goldilocks_benchmark!(&mut group, 6, 131072, 4744, 6);

        run_single_goldilocks_benchmark!(&mut group, 1, 262144, 0, 31); // Bound not practical
        run_single_goldilocks_benchmark!(&mut group, 2, 262144, 4, 32);
        run_single_goldilocks_benchmark!(&mut group, 3, 262144, 32, 13);
        run_single_goldilocks_benchmark!(&mut group, 4, 262144, 180, 9);
        run_single_goldilocks_benchmark!(&mut group, 5, 262144, 836, 7);
        run_single_goldilocks_benchmark!(&mut group, 6, 262144, 3354, 6);

        run_single_goldilocks_benchmark!(&mut group, 1, 524288, 0, 31); // Bound not practical
        run_single_goldilocks_benchmark!(&mut group, 2, 524288, 2, 64);
        run_single_goldilocks_benchmark!(&mut group, 3, 524288, 22, 15);
        run_single_goldilocks_benchmark!(&mut group, 4, 524288, 128, 10);
        run_single_goldilocks_benchmark!(&mut group, 5, 524288, 592, 7);
        run_single_goldilocks_benchmark!(&mut group, 6, 524288, 2372, 6);

        run_single_goldilocks_benchmark!(&mut group, 1, 1048576, 0, 31); // Bound not practical
        run_single_goldilocks_benchmark!(&mut group, 2, 1048576, 2, 64);
        run_single_goldilocks_benchmark!(&mut group, 3, 1048576, 16, 16);
        run_single_goldilocks_benchmark!(&mut group, 4, 1048576, 90, 10);
        run_single_goldilocks_benchmark!(&mut group, 5, 1048576, 418, 8);
        run_single_goldilocks_benchmark!(&mut group, 6, 1048576, 1678, 6);
    }

    // StarkPrime
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Ajtai StarkPrime");
        group.plot_config(plot_config.clone());

        run_single_starkprime_benchmark!(&mut group, 1, 32768, 52, 45);
        run_single_starkprime_benchmark!(&mut group, 2, 32768, 4028, 21);
        run_single_starkprime_benchmark!(&mut group, 3, 32768, 114286, 15);
        run_single_starkprime_benchmark!(&mut group, 4, 32768, 1918124, 13);
        run_single_starkprime_benchmark!(&mut group, 5, 32768, 23015556, 11);
        run_single_starkprime_benchmark!(&mut group, 6, 32768, 217592018, 10);

        run_single_starkprime_benchmark!(&mut group, 1, 65536, 36, 49);
        run_single_starkprime_benchmark!(&mut group, 2, 65536, 2848, 22);
        run_single_starkprime_benchmark!(&mut group, 3, 65536, 80812, 16);
        run_single_starkprime_benchmark!(&mut group, 4, 65536, 1356318, 13);
        run_single_starkprime_benchmark!(&mut group, 5, 65536, 16274456, 11);
        run_single_starkprime_benchmark!(&mut group, 6, 65536, 153860792, 10);

        run_single_starkprime_benchmark!(&mut group, 1, 131072, 26, 54);
        run_single_starkprime_benchmark!(&mut group, 2, 131072, 2014, 23);
        run_single_starkprime_benchmark!(&mut group, 3, 131072, 57142, 16);
        run_single_starkprime_benchmark!(&mut group, 4, 131072, 959062, 13);
        run_single_starkprime_benchmark!(&mut group, 5, 131072, 11507778, 11);
        run_single_starkprime_benchmark!(&mut group, 6, 131072, 108796010, 10);

        run_single_starkprime_benchmark!(&mut group, 1, 262144, 18, 61);
        run_single_starkprime_benchmark!(&mut group, 2, 262144, 1424, 24);
        run_single_starkprime_benchmark!(&mut group, 3, 262144, 40406, 17);
        run_single_starkprime_benchmark!(&mut group, 4, 262144, 678160, 13);
        run_single_starkprime_benchmark!(&mut group, 5, 262144, 8137228, 11);
        run_single_starkprime_benchmark!(&mut group, 6, 262144, 76930396, 10);

        run_single_starkprime_benchmark!(&mut group, 1, 524288, 12, 71);
        run_single_starkprime_benchmark!(&mut group, 2, 524288, 1006, 26);
        run_single_starkprime_benchmark!(&mut group, 3, 524288, 28572, 17);
        run_single_starkprime_benchmark!(&mut group, 4, 524288, 479530, 14);
        run_single_starkprime_benchmark!(&mut group, 5, 524288, 5753890, 12);
        run_single_starkprime_benchmark!(&mut group, 6, 524288, 54398004, 10);

        run_single_starkprime_benchmark!(&mut group, 1, 1048576, 10, 76);
        run_single_starkprime_benchmark!(&mut group, 2, 1048576, 712, 27);
        run_single_starkprime_benchmark!(&mut group, 3, 1048576, 20204, 18);
        run_single_starkprime_benchmark!(&mut group, 4, 1048576, 339080, 14);
        run_single_starkprime_benchmark!(&mut group, 5, 1048576, 4068614, 12);
        run_single_starkprime_benchmark!(&mut group, 6, 1048576, 38465198, 10);
    }

    // Frog
    {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Ajtai Frog");
        group.plot_config(plot_config.clone());

        run_single_frog_benchmark!(&mut group, 1, 32768, 0, 1);
        run_single_frog_benchmark!(&mut group, 2, 32768, 2, 64);
        run_single_frog_benchmark!(&mut group, 3, 32768, 14, 17);
        run_single_frog_benchmark!(&mut group, 4, 32768, 56, 11);
        run_single_frog_benchmark!(&mut group, 5, 32768, 196, 9);
        run_single_frog_benchmark!(&mut group, 6, 32768, 610, 7);

        run_single_frog_benchmark!(&mut group, 1, 65536, 0, 1);
        run_single_frog_benchmark!(&mut group, 2, 65536, 2, 64);
        run_single_frog_benchmark!(&mut group, 3, 65536, 10, 20);
        run_single_frog_benchmark!(&mut group, 4, 65536, 40, 12);
        run_single_frog_benchmark!(&mut group, 5, 65536, 138, 9);
        run_single_frog_benchmark!(&mut group, 6, 65536, 430, 8);

        //run_single_frog_benchmark!(&mut group, 1, 131072, 0, 1); // Bound not practical
        run_single_frog_benchmark!(&mut group, 2, 131072, 2, 64);
        run_single_frog_benchmark!(&mut group, 3, 131072, 6, 25);
        run_single_frog_benchmark!(&mut group, 4, 131072, 28, 14);
        run_single_frog_benchmark!(&mut group, 5, 131072, 98, 10);
        run_single_frog_benchmark!(&mut group, 6, 131072, 304, 8);

        run_single_frog_benchmark!(&mut group, 1, 262144, 0, 1); // Bound not practical
        run_single_frog_benchmark!(&mut group, 2, 262144, 1, 1);
        run_single_frog_benchmark!(&mut group, 3, 262144, 4, 32);
        run_single_frog_benchmark!(&mut group, 4, 262144, 20, 15);
        run_single_frog_benchmark!(&mut group, 5, 262144, 70, 11);
        run_single_frog_benchmark!(&mut group, 6, 262144, 216, 9);

        run_single_frog_benchmark!(&mut group, 1, 524288, 0, 1); // Bound not practical
        run_single_frog_benchmark!(&mut group, 2, 524288, 1, 1);
        run_single_frog_benchmark!(&mut group, 3, 524288, 4, 32);
        run_single_frog_benchmark!(&mut group, 4, 524288, 14, 17);
        run_single_frog_benchmark!(&mut group, 5, 524288, 50, 12);
        run_single_frog_benchmark!(&mut group, 6, 524288, 152, 9);

        run_single_frog_benchmark!(&mut group, 1, 1048576, 0, 1); // Bound not practical
        run_single_frog_benchmark!(&mut group, 2, 1048576, 1, 1);
        run_single_frog_benchmark!(&mut group, 3, 1048576, 2, 64);
        run_single_frog_benchmark!(&mut group, 4, 1048576, 10, 20);
        run_single_frog_benchmark!(&mut group, 5, 1048576, 34, 13);
        run_single_frog_benchmark!(&mut group, 6, 1048576, 108, 10);

        group.finish();
    }

    // Dilithium
    // TODO: Implement Suitable Ring for Dilithium
    // TODO: Compute proper bounds
    /*{
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("Ajtai BabyBear");
        group.plot_config(plot_config.clone());

        run_single_dilithium_benchmark!(&mut group, 1, 32768, 2, 31);
        run_single_dilithium_benchmark!(&mut group, 2, 32768, 42, 6);
        run_single_dilithium_benchmark!(&mut group, 3, 32768, 508, 4);
        run_single_dilithium_benchmark!(&mut group, 4, 32768, 4160, 3);
        run_single_dilithium_benchmark!(&mut group, 5, 32768, 26458, 3);
        run_single_dilithium_benchmark!(&mut group, 6, 32768, 140862, 2);

        run_single_dilithium_benchmark!(&mut group, 1, 65536, 2, 31);
        run_single_dilithium_benchmark!(&mut group, 2, 65536, 28, 7);
        run_single_dilithium_benchmark!(&mut group, 3, 65536, 360, 4);
        run_single_dilithium_benchmark!(&mut group, 4, 65536, 2942, 3);
        run_single_dilithium_benchmark!(&mut group, 5, 65536, 18708, 3);
        run_single_dilithium_benchmark!(&mut group, 6, 65536, 99604, 2);

        run_single_dilithium_benchmark!(&mut group, 1, 131072, 0, 31); // Bound not practical
        run_single_dilithium_benchmark!(&mut group, 2, 131072, 20, 8);
        run_single_dilithium_benchmark!(&mut group, 3, 131072, 254, 4);
        run_single_dilithium_benchmark!(&mut group, 4, 131072, 2080, 3);
        run_single_dilithium_benchmark!(&mut group, 5, 131072, 13228, 3);
        run_single_dilithium_benchmark!(&mut group, 6, 131072, 70430, 2);

        run_single_dilithium_benchmark!(&mut group, 1, 262144, 0, 31); // Bound not practical
        run_single_dilithium_benchmark!(&mut group, 2, 262144, 14, 9);
        run_single_dilithium_benchmark!(&mut group, 3, 262144, 180, 5);
        run_single_dilithium_benchmark!(&mut group, 4, 262144, 1470, 3);
        run_single_dilithium_benchmark!(&mut group, 5, 262144, 9354, 3);
        run_single_dilithium_benchmark!(&mut group, 6, 262144, 49802, 2);

        run_single_dilithium_benchmark!(&mut group, 1, 524288, 0, 31); // Bound not practical
        run_single_dilithium_benchmark!(&mut group, 2, 524288, 10, 10);
        run_single_dilithium_benchmark!(&mut group, 3, 524288, 126, 5);
        run_single_dilithium_benchmark!(&mut group, 4, 524288, 1040, 4);
        run_single_dilithium_benchmark!(&mut group, 5, 524288, 6614, 3);
        run_single_dilithium_benchmark!(&mut group, 6, 524288, 35214, 3);

        run_single_dilithium_benchmark!(&mut group, 1, 1048576, 0, 31); // Bound not practical
        run_single_dilithium_benchmark!(&mut group, 2, 1048576, 6, 12);
        run_single_dilithium_benchmark!(&mut group, 3, 1048576, 90, 5);
        run_single_dilithium_benchmark!(&mut group, 4, 1048576, 734, 4);
        run_single_dilithium_benchmark!(&mut group, 5, 1048576, 4676, 3);
        run_single_dilithium_benchmark!(&mut group, 6, 1048576, 24900, 3);

        group.finish();
    }*/
}

pub fn benchmarks_main(c: &mut Criterion) {
    ajtai_benchmarks(c);
}

criterion_group!(
    name=benches;
    config = Criterion::default()
            .sample_size(10)
            .measurement_time(Duration::from_secs(50))
            .warm_up_time(Duration::from_secs(1));
    targets = benchmarks_main
);
criterion_main!(benches);
