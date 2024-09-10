use ark_ff::{Field, UniformRand};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use lattirust_arithmetic::{
    challenge_set::latticefold_challenge_set::OverField,
    ring::{CyclotomicPolyRingSplittedNTT, PolyRing, Pow2CyclotomicPolyRingNTT, Zq},
};
use num_traits::Pow;
use rand::thread_rng;
use std::time::Duration;

use latticefold::{
    commitment::AjtaiCommitmentScheme,
    parameters::{
        BabyBearTestParams, DecompositionParamData, DecompositionParams, DilithiumTestParams,
        GoldilocksCR, GoldilocksTestParams, Pow2_57TestParams, Pow2_59TestParams, BABYBEAR_PRIME,
        DILITHIUM_PRIME, GOLDILOCKS_PRIME, POW2_57_PRIME, POW2_59_PRIME,
    },
};

fn ajtai_benchmark<const C: usize, const W: usize, P: DecompositionParams, R: OverField>(
    c: &mut Criterion,
    p: P,
    prime_name: &str,
) {
    let delta = 1.01_f64;
    let degree = R::dimension() as f64;
    let rows = C as f64;
    let mod_p: f64 = (15 * (1 << 27) + 1) as f64;
    let inner = delta.log2() * degree * rows * mod_p.log2();
    let B = 2_f64.pow(inner.sqrt()).floor();
    println!("Bound B = {}", B);
    let ajtai_data: AjtaiCommitmentScheme<C, W, R> = AjtaiCommitmentScheme::rand(&mut thread_rng());

    let witness: Vec<R> = (0..W)
        .map(|_| {
            let mut value = R::rand(&mut thread_rng());
            let mut inf_norm = value.linf_norm() as f64;
            println!("Todo bien!");
            while inf_norm > B {
                value = R::rand(&mut thread_rng());
                inf_norm = value.linf_norm() as f64;
            }
            value
        })
        .collect();

    let bench_name = format!(
        "Ajtai {}x{} for {} with degree {} and bound {}",
        C,
        W,
        prime_name,
        R::dimension(),
        B,
    );
    c.bench_with_input(
        BenchmarkId::new(bench_name, DecompositionParamData::from(p)),
        &(ajtai_data, witness),
        |b, (ajtai_data, witness)| b.iter(|| ajtai_data.commit_ntt(witness)),
    );
}

fn ajtai_benchmarks(c: &mut Criterion) {
    // ajtai_benchmark::<DILITHIUM_PRIME, 256, 9, { 1 << 15 }, _>(
    //     c,
    //     DilithiumTestParams,
    //     "Dilithium prime",
    // );
    // ajtai_benchmark::<POW2_59_PRIME, 256, 9, { 1 << 15 }, _>(
    //     c,
    //     Pow2_59TestParams,
    //     "p = 27*(1<<59) + 1",
    // );
    // ajtai_benchmark::<POW2_57_PRIME, 256, 9, { 1 << 15 }, _>(
    //     c,
    //     Pow2_57TestParams,
    //     "p = 71*(1<<57) + 1",
    // );

    // // TODO: more benchmarks with different params.
    // // BabyBear benches
    // // Baby bear for 72 degree in NTT splitted
    ajtai_benchmark::<
        1,
        { 1 << 10 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");
    ajtai_benchmark::<
        2,
        { 1 << 10 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");
    ajtai_benchmark::<
        3,
        { 1 << 10 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");
    ajtai_benchmark::<
        4,
        { 1 << 10 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");
    ajtai_benchmark::<
        5,
        { 1 << 10 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");
    ajtai_benchmark::<
        6,
        { 1 << 10 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");

    ajtai_benchmark::<
        1,
        { 1 << 11 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");
    ajtai_benchmark::<
        2,
        { 1 << 11 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");
    ajtai_benchmark::<
        3,
        { 1 << 11 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");
    ajtai_benchmark::<
        4,
        { 1 << 11 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");
    ajtai_benchmark::<
        5,
        { 1 << 11 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");
    ajtai_benchmark::<
        6,
        { 1 << 11 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");

    ajtai_benchmark::<
        1,
        { 1 << 12 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");
    ajtai_benchmark::<
        2,
        { 1 << 12 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");
    ajtai_benchmark::<
        3,
        { 1 << 12 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");
    ajtai_benchmark::<
        4,
        { 1 << 12 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");
    ajtai_benchmark::<
        5,
        { 1 << 12 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");
    ajtai_benchmark::<
        6,
        { 1 << 12 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");

    ajtai_benchmark::<
        1,
        { 1 << 13 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");
    ajtai_benchmark::<
        2,
        { 1 << 13 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");
    ajtai_benchmark::<
        3,
        { 1 << 13 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");
    ajtai_benchmark::<
        4,
        { 1 << 13 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");
    ajtai_benchmark::<
        5,
        { 1 << 13 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");
    ajtai_benchmark::<
        6,
        { 1 << 13 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");

    ajtai_benchmark::<
        1,
        { 1 << 14 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");
    ajtai_benchmark::<
        2,
        { 1 << 14 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");
    ajtai_benchmark::<
        3,
        { 1 << 14 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");
    ajtai_benchmark::<
        4,
        { 1 << 14 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");
    ajtai_benchmark::<
        5,
        { 1 << 14 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");
    ajtai_benchmark::<
        6,
        { 1 << 14 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");

    ajtai_benchmark::<
        1,
        { 1 << 15 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");
    ajtai_benchmark::<
        2,
        { 1 << 15 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");
    ajtai_benchmark::<
        3,
        { 1 << 15 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");
    ajtai_benchmark::<
        4,
        { 1 << 15 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");
    ajtai_benchmark::<
        5,
        { 1 << 15 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");
    ajtai_benchmark::<
        6,
        { 1 << 15 },
        _,
        CyclotomicPolyRingSplittedNTT<BABYBEAR_PRIME, 420899707, 72, 9, 24, 8>,
    >(c, BabyBearTestParams, "BabyBear");

    // Goldilocks for degre 24
    // println!("Run first GL");
    // println!("Goldilocks prime: {}", GOLDILOCKS_PRIME);
    // let g = Zq::<GOLDILOCKS_PRIME>::from(7);
    // println!("Goldilocks generator: {}", g);
    // println!("{}^{} = {}", g, GOLDILOCKS_PRIME-1, g.pow([GOLDILOCKS_PRIME - 1]));

    // let z = 24;
    // let one = Zq::<24>::from(GOLDILOCKS_PRIME);
    // println!("{} ~= {} mod {}", GOLDILOCKS_PRIME, one, m);
    // println!("{} ~= {} mod {}", GOLDILOCKS_PRIME, GOLDILOCKS_PRIME % m, m);

    // ajtai_benchmark::<
    //     1,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // println!("Run second GL");
    // ajtai_benchmark::<
    //     2,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     3,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     4,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     5,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     6,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     7,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     8,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     9,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     10,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     11,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     12,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     13,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     14,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     15,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     16,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     17,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     18,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     19,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     20,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     20,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     21,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     22,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     23,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     24,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     25,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     26,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     27,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     28,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     29,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     30,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     31,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     32,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     33,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     34,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     35,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     36,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     37,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     38,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     39,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     40,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     41,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     42,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     43,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     44,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     45,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     5,
    //     { 1 << 9 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");

    // ajtai_benchmark::<
    //     1,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     2,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     3,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     4,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     5,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     6,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");

    // ajtai_benchmark::<
    //     7,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     8,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     9,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     10,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     11,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     12,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     13,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     14,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     15,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     16,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     17,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     18,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     19,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     20,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     20,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     21,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     22,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     23,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     24,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     25,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     26,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     27,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     28,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     29,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     30,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     31,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     32,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     33,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     34,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     35,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     36,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     37,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     38,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     39,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     40,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     41,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     42,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     43,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     44,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     45,
    //     { 1 << 10 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");

    // ajtai_benchmark::<
    //     1,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     2,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     3,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     4,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     5,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     6,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     7,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     8,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     9,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     10,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     11,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     12,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     13,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     14,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     15,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     16,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     17,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     18,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     19,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     20,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     20,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     21,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     22,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     23,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     24,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     25,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     26,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     27,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     28,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     29,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     30,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     31,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     32,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     33,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     34,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     35,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     36,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     37,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     38,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     39,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     40,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     41,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     42,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     43,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     44,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     45,
    //     { 1 << 11 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");

    // ajtai_benchmark::<
    //     1,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     2,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     3,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     4,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     5,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     6,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     7,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     8,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     9,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     10,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     11,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     12,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     13,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     14,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     15,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     16,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     17,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     18,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     19,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     20,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     20,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     21,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     22,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     23,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     24,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     25,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     26,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     27,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     28,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     29,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     30,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     31,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     32,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     33,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     34,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     35,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     36,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     37,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     38,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     39,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     40,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     41,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     42,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     43,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     44,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     45,
    //     { 1 << 12 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");

    // ajtai_benchmark::<
    //     1,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     2,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     3,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     4,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     5,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     6,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     7,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     8,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     9,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     10,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     11,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     12,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     13,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     14,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     15,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     16,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     17,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     18,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     19,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     20,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     20,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     21,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     22,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     23,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     24,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     25,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     26,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     27,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     28,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     29,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     30,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     31,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     32,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     33,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     34,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     35,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     36,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     37,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     38,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     39,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     40,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     41,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     42,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     43,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     44,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     45,
    //     { 1 << 13 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");

    // ajtai_benchmark::<
    //     1,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     2,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     3,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     4,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     5,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     6,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     7,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     8,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     9,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     10,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     11,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     12,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     13,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     14,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     15,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     16,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     17,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     18,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     19,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     20,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     20,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     21,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     22,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     23,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     24,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     25,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     26,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     27,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     28,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     29,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     30,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     31,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     32,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     33,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     34,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     35,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     36,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     37,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     38,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     39,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     40,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     41,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     42,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     43,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     44,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     45,
    //     { 1 << 14 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");

    // ajtai_benchmark::<
    //     1,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     2,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     3,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     4,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     5,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     6,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     7,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     8,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     9,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     10,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     11,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     12,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     13,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     14,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     15,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     16,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     17,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     18,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     19,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     20,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     20,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     21,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     22,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     23,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     24,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     25,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     26,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     27,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     28,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     29,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     30,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     31,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     32,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     33,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     34,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     35,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     36,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     37,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     38,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     39,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     40,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     41,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     42,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     43,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     44,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<
    //     45,
    //     { 1 << 15 },
    //     _,
    //     CyclotomicPolyRingSplittedNTT<GOLDILOCKS_PRIME, 1099511627776, 24, 3, 24, 8>,
    // >(c, GoldilocksTestParams, "Goldilocks");

    // degree 64
    ajtai_benchmark::<1, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<1, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<1, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<1, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<1, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<1, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<1, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<2, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<2, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<2, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<2, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<2, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<2, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<2, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<3, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<3, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<3, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<3, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<3, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<3, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<3, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<4, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<4, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<4, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<4, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<4, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<4, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<4, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<5, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<5, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<5, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<5, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<5, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<5, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<5, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<6, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<6, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<6, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<6, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<6, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<6, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<6, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<7, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<7, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<7, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<7, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<7, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<7, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<7, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 64>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );

    // degree 128
    ajtai_benchmark::<1, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 128>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<1, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 128>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<1, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 128>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<1, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 128>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<1, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 128>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<1, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 128>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<1, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 128>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<2, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 128>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<2, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 128>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<2, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 128>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<2, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 128>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<2, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 128>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<2, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 128>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<2, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 128>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<3, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 128>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<3, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 128>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<3, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 128>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<3, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 128>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<3, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 128>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<3, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 128>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );
    ajtai_benchmark::<3, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 128>>(
        c,
        BabyBearTestParams,
        "BabyBear",
    );

    // ajtai_benchmark::<10, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<10, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<10, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<10, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<10, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<10, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<10, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<11, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<11, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<11, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<11, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<11, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<11, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<11, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<12, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<12, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<12, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<12, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<12, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<12, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<12, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<13, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<13, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<13, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<13, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<13, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<13, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<13, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<14, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<14, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<14, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<14, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<14, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<14, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<14, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<15, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<15, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<15, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<15, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<15, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<15, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");
    // ajtai_benchmark::<15, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<BABYBEAR_PRIME, 32>>(c, BabyBearTestParams, "BabyBear");

    // // Goldilocks
    // ajtai_benchmark::<1, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<1, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<1, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<1, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<1, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<1, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<1, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<1, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<2, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<2, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<2, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<2, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<2, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<2, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<2, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<3, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<3, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<3, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<3, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<3, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<3, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<3, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<4, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<4, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<4, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<4, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<4, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<4, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<4, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<5, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<5, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<5, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<5, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<5, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<5, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<5, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<6, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<6, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<6, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<6, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<6, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<6, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<6, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<7, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<7, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<7, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<7, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<7, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<7, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<7, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<10, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<10, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<10, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<10, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<10, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<10, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<10, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<11, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<11, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<11, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<11, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<11, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<11, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<11, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<12, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<12, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<12, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<12, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<12, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<12, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<12, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<13, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<13, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<13, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<13, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<13, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<13, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<13, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<14, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<14, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<14, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<14, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<14, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<14, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<14, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<15, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<15, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<15, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<15, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<15, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<15, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<15, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<16, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<16, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<16, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<16, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<16, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<16, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<17, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<17, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<17, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<17, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<17, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<17, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<17, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<18, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<18, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<18, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<18, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<18, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<18, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<18, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<19, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<19, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<19, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<19, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<19, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<19, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<19, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<20, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<20, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<20, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<20, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<20, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<20, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<20, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<21, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<21, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<21, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<21, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<21, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<21, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<21, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<22, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<22, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<22, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<22, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<22, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<22, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<22, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<23, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<23, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<23, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<23, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<23, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<23, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<23, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<24, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<24, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<24, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<24, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<24, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<24, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<24, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<25, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<25, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<25, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<25, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<25, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<25, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<25, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<26, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<26, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<26, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<26, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<26, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<26, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<26, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<27, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<27, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<27, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<27, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<27, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<27, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<27, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<28, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<28, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<28, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<28, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<28, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<28, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<28, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<29, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<29, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<29, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<29, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<29, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<29, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<29, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<30, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<30, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<30, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<30, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<30, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<30, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<30, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<31, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<31, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<31, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<31, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<31, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<31, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<31, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<32, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<32, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<32, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<32, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<32, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<32, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<32, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<33, { 1 << 9 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<33, { 1 << 10 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<33, { 1 << 11 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<33, { 1 << 12 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<33, { 1 << 13 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<33, { 1 << 14 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
    // ajtai_benchmark::<33, { 1 << 15 }, _, Pow2CyclotomicPolyRingNTT<GOLDILOCKS_PRIME, 32>>(c, GoldilocksTestParams, "Goldilocks");
}

pub fn benchmarks_main(c: &mut Criterion) {
    ajtai_benchmarks(c);
}

criterion_group!(
    name=benches;
    config = Criterion::default().sample_size(10).measurement_time(Duration::from_secs(50)).warm_up_time(Duration::from_secs(1));
    targets = benchmarks_main);
criterion_main!(benches);
