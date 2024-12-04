use ark_ff::UniformRand;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use cyclotomic_rings::rings::{GoldilocksRingNTT, SuitableRing};
use lattirust_ring::cyclotomic_ring::{CRT, ICRT};

fn ring_ntt_operations_benchmark<R: SuitableRing>(c: &mut Criterion, ring_name: &str) {
    let mut rng = rand::thread_rng();
    let a = R::rand(&mut rng);
    let b = R::rand(&mut rng);

    c.bench_with_input(
        BenchmarkId::new("NTT Multiplication", ring_name),
        &(a, b),
        |b, (x, y)| {
            b.iter_batched(
                || (x.clone(), y.clone()),
                |input| input.0 * input.1,
                criterion::BatchSize::SmallInput,
            )
        },
    );
    c.bench_with_input(
        BenchmarkId::new("NTT Addition", ring_name),
        &(a, b),
        |b, (x, y)| {
            b.iter_batched(
                || (x.clone(), y.clone()),
                |input| input.0 + input.1,
                criterion::BatchSize::SmallInput,
            )
        },
    );
    c.bench_with_input(
        BenchmarkId::new("NTT Substraction", ring_name),
        &(a, b),
        |b, (x, y)| {
            b.iter_batched(
                || (x.clone(), y.clone()),
                |input| input.0 - input.1,
                criterion::BatchSize::SmallInput,
            )
        },
    );
}

fn mle_evaluation_benchmark<R: SuitableRing>(c: &mut Criterion, ring_name: &str, nv: usize) {
    use lattirust_poly::mle::DenseMultilinearExtension;

    let mut rng = rand::thread_rng();
    let evals: Vec<R> = (0..1 << nv).map(|_| R::rand(&mut rng)).collect();
    let mle = DenseMultilinearExtension::<R>::from_evaluations_vec(nv, evals);
    let point: Vec<R> = (0..nv).map(|_| R::rand(&mut rng)).collect();

    c.bench_with_input(
        BenchmarkId::new(
            "MLE Evaluation",
            format!("{} with 2^{} evals", ring_name, nv),
        ),
        &(mle, point),
        |b, (mle, point)| b.iter(|| mle.evaluate(point)),
    );
}

fn all_ring_mle_evaluations(c: &mut Criterion) {
    for nv in 9..=24 {
        mle_evaluation_benchmark::<GoldilocksRingNTT>(c, &format!("Goldilocks NV={}", nv), nv);
    }
}

fn ring_crt_icrt_benchmark<R: SuitableRing>(c: &mut Criterion, ring_name: &str, nv: usize) {
    let mut rng = rand::thread_rng();
    let vec_ntt_form = (0..(1 << nv))
        .map(|_| R::rand(&mut rng))
        .collect::<Vec<R>>();
    let vec_coeff_form = (0..(1 << nv))
        .map(|_| R::rand(&mut rng).icrt())
        .collect::<Vec<_>>();
    c.bench_with_input(
        BenchmarkId::new(
            "Elementwise CRT",
            format!("{} of size = {}", ring_name, 1 << nv),
        ),
        &vec_coeff_form,
        |b, input| {
            b.iter_batched(
                || input.clone(),
                |input| CRT::elementwise_crt(input),
                criterion::BatchSize::SmallInput,
            )
        },
    );
    c.bench_with_input(
        BenchmarkId::new(
            "Elementwise ICRT",
            format!("{} of size = {}", ring_name, 1 << nv),
        ),
        &vec_ntt_form,
        |b, input| {
            b.iter_batched(
                || input.clone(),
                |input| ICRT::elementwise_icrt(input),
                criterion::BatchSize::SmallInput,
            )
        },
    );
}

fn all_ring_crt_icrt_operations(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let ring_name = "Goldilocks";
    let ntt_form = GoldilocksRingNTT::rand(&mut rng);
    c.bench_with_input(
        BenchmarkId::new(format!("ICRT"), ring_name),
        &ntt_form,
        |b, input| b.iter(|| input.icrt()),
    );
    let coeff_form = ntt_form.icrt();
    c.bench_with_input(
        BenchmarkId::new(format!("CRT"), ring_name),
        &coeff_form,
        |b, input| b.iter(|| input.crt()),
    );
    for nv in 9..=24 {
        ring_crt_icrt_benchmark::<GoldilocksRingNTT>(c, &ring_name, nv);
    }
}

fn all_ring_ntt_operations(c: &mut Criterion) {
    ring_ntt_operations_benchmark::<GoldilocksRingNTT>(c, "Goldilocks");
}

criterion_group!(
    benches,
    all_ring_ntt_operations,
    all_ring_mle_evaluations,
    all_ring_crt_icrt_operations
);
criterion_main!(benches);
