use std::{fmt::Debug, marker::PhantomData};

use ark_std::UniformRand;
use criterion::{measurement::WallTime, BenchmarkGroup, BenchmarkId};
use latticefold::{
    arith::{
        ccs::get_test_dummy_degree_three_ccs_non_scalar,
        r1cs::{
            get_test_dummy_r1cs, get_test_dummy_r1cs_non_scalar, get_test_dummy_z_split,
            get_test_dummy_z_split_ntt,
        },
        Arith, Witness, CCCS, CCS, LCCCS,
    },
    commitment::AjtaiCommitmentScheme,
    decomposition_parameters::DecompositionParams,
    nifs::{
        decomposition::{
            DecompositionProver, DecompositionVerifier, LFDecompositionProver,
            LFDecompositionVerifier,
        },
        folding::{FoldingProver, FoldingVerifier, LFFoldingProver, LFFoldingVerifier},
        linearization::{
            LFLinearizationProver, LFLinearizationVerifier, LinearizationProof,
            LinearizationProver, LinearizationVerifier,
        },
        NIFSProver, NIFSVerifier,
    },
    transcript::poseidon::PoseidonTranscript,
};
use latticefold_rings::rings::SuitableRing;

#[derive(Clone, Copy)]
pub enum R1CS {
    Scalar,
    NonScalar,
    DegreeThreeNonScalar,
}

#[allow(dead_code)]
#[allow(non_snake_case)]
pub fn get_test_dummy_ccs<R: Clone + UniformRand + Debug + SuitableRing>(
    x_len: usize,
    n: usize,
    wit_len: usize,
    r1cs_rows: usize,
    L: usize,
) -> CCS<R> {
    let r1cs = get_test_dummy_r1cs::<R>(x_len, wit_len, r1cs_rows);
    CCS::<R>::from_r1cs_padded(r1cs, n, L)
}

#[allow(dead_code)]
pub fn wit_and_ccs_gen<P: DecompositionParams, R: Clone + UniformRand + Debug + SuitableRing>(
    x_len: usize,
    n: usize,
    wit_len: usize,
    r1cs_rows: usize,
    kappa: usize,
) -> (CCCS<R>, Witness<R>, CCS<R>, AjtaiCommitmentScheme<R>) {
    let mut rng = ark_std::test_rng();

    let new_r1cs_rows = if P::L == 1 && (wit_len > 0 && (wit_len & (wit_len - 1)) == 0) {
        r1cs_rows - 2
    } else {
        r1cs_rows // This makes a square matrix but is too much memory
    };
    let ccs: CCS<R> = get_test_dummy_ccs::<R>(x_len, n, wit_len, new_r1cs_rows, P::L);
    let (one, x_ccs, w_ccs) = get_test_dummy_z_split::<R>(x_len, wit_len);
    let mut z = vec![one];
    z.extend(&x_ccs);
    z.extend(&w_ccs);
    ccs.check_relation(&z).expect("R1CS invalid!");

    let scheme: AjtaiCommitmentScheme<R> = AjtaiCommitmentScheme::rand(kappa, n, &mut rng);
    let wit: Witness<R> = Witness::from_w_ccs::<P>(w_ccs);

    let cm_i: CCCS<R> = CCCS {
        cm: wit.commit::<P>(&scheme).unwrap(),
        x_ccs,
    };

    (cm_i, wit, ccs, scheme)
}

#[allow(dead_code)]
pub fn wit_and_ccs_gen_non_scalar<
    P: DecompositionParams,
    R: Clone + UniformRand + Debug + SuitableRing,
>(
    x_len: usize,
    n: usize,
    wit_len: usize,
    r1cs_rows: usize,
    kappa: usize,
) -> (CCCS<R>, Witness<R>, CCS<R>, AjtaiCommitmentScheme<R>) {
    let mut rng = ark_std::test_rng();

    let new_r1cs_rows = if P::L == 1 && (wit_len > 0 && (wit_len & (wit_len - 1)) == 0) {
        r1cs_rows - 2
    } else {
        r1cs_rows // This makes a square matrix but is too much memory
    };
    let (one, x_ccs, w_ccs) = get_test_dummy_z_split_ntt::<R>(x_len, wit_len);

    let mut z = vec![one];
    z.extend(&x_ccs);
    z.extend(&w_ccs);
    let ccs: CCS<R> = get_test_dummy_ccs_non_scalar::<R>(x_len, n, new_r1cs_rows, P::L, &z);
    ccs.check_relation(&z).expect("R1CS invalid!");

    let scheme: AjtaiCommitmentScheme<R> = AjtaiCommitmentScheme::rand(kappa, n, &mut rng);
    let wit: Witness<R> = Witness::from_w_ccs::<P>(w_ccs);

    let cm_i: CCCS<R> = CCCS {
        cm: wit.commit::<P>(&scheme).unwrap(),
        x_ccs,
    };

    (cm_i, wit, ccs, scheme)
}

#[allow(dead_code)]
#[allow(non_snake_case)]
pub fn get_test_dummy_ccs_non_scalar<R: Clone + UniformRand + Debug + SuitableRing>(
    x_len: usize,
    n: usize,
    r1cs_rows: usize,
    L: usize,
    witness: &[R],
) -> CCS<R> {
    let r1cs = get_test_dummy_r1cs_non_scalar::<R>(x_len, n, r1cs_rows, witness);
    CCS::<R>::from_r1cs_padded(r1cs, n, L)
}

#[allow(dead_code)]
pub fn wit_and_ccs_gen_degree_three_non_scalar<
    P: DecompositionParams,
    R: Clone + UniformRand + Debug + SuitableRing,
>(
    x_len: usize,
    n: usize,
    wit_len: usize,
    r1cs_rows: usize,
    kappa: usize,
) -> (CCCS<R>, Witness<R>, CCS<R>, AjtaiCommitmentScheme<R>) {
    let mut rng = ark_std::test_rng();

    let new_r1cs_rows = if P::L == 1 && (wit_len > 0 && (wit_len & (wit_len - 1)) == 0) {
        r1cs_rows - 2
    } else {
        r1cs_rows // This makes a square matrix but is too much memory
    };
    let (one, x_ccs, w_ccs) = get_test_dummy_z_split_ntt::<R>(x_len, wit_len);

    let mut z = vec![one];
    z.extend(&x_ccs);
    z.extend(&w_ccs);
    let ccs: CCS<R> =
        get_test_dummy_degree_three_ccs_non_scalar::<R>(&z, x_len, n, wit_len, P::L, new_r1cs_rows);
    ccs.check_relation(&z).expect("R1CS invalid!");

    let scheme: AjtaiCommitmentScheme<R> = AjtaiCommitmentScheme::rand(kappa, n, &mut rng);
    let wit: Witness<R> = Witness::from_w_ccs::<P>(w_ccs);

    let cm_i: CCCS<R> = CCCS {
        cm: wit.commit::<P>(&scheme).unwrap(),
        x_ccs,
    };

    (cm_i, wit, ccs, scheme)
}

pub struct Bencher<P: DecompositionParams, R: SuitableRing + Clone> {
    phantom_data: PhantomData<(P, R)>,
}

#[allow(dead_code)]
impl<P: DecompositionParams, R: SuitableRing + Clone> Bencher<P, R> {
    pub fn setup_r1cs(
        x_len: usize,
        n: usize,
        wit_len: usize,
        kappa: usize,
        t: R1CS,
    ) -> (CCCS<R>, Witness<R>, CCS<R>, AjtaiCommitmentScheme<R>) {
        let r1cs_rows = x_len + wit_len + 1;

        match t {
            R1CS::Scalar => wit_and_ccs_gen::<P, R>(x_len, n, wit_len, r1cs_rows, kappa),
            R1CS::NonScalar => {
                wit_and_ccs_gen_non_scalar::<P, R>(x_len, n, wit_len, r1cs_rows, kappa)
            }
            R1CS::DegreeThreeNonScalar => {
                wit_and_ccs_gen_degree_three_non_scalar::<P, R>(x_len, n, wit_len, r1cs_rows, kappa)
            }
        }
    }

    pub fn verify_linearization_proof(
        proof: &LinearizationProof<R>,
        transcript: &mut PoseidonTranscript<R>,
        cm_i: &CCCS<R>,
        ccs: &CCS<R>,
    ) -> LCCCS<R> {
        LFLinearizationVerifier::<_, PoseidonTranscript<R>>::verify(cm_i, proof, transcript, ccs)
            .expect("Failed to verify linearization proof")
    }

    pub fn bench_linearization_prover(
        group: &mut BenchmarkGroup<WallTime>,
        x_len: usize,
        n: usize,
        wit_len: usize,
        kappa: usize,
        t: R1CS,
    ) {
        group.bench_function(
            BenchmarkId::new(
                "Linearization Prover",
                format!(
                    "Param. Kappa={}, Cols={},  B={}, L={}, B_small={}, K={}",
                    kappa,
                    { n / P::L },
                    P::B,
                    P::L,
                    P::B_SMALL,
                    P::K
                ),
            ),
            |b| {
                let transcript = PoseidonTranscript::<R>::default();
                let (cm_i, wit, ccs, _) = Self::setup_r1cs(x_len, n, wit_len, kappa, t);
                b.iter_batched(
                    || transcript.clone(),
                    |mut bench_transcript| {
                        let _ = LFLinearizationProver::<_, PoseidonTranscript<R>>::prove(
                            &cm_i,
                            &wit,
                            &mut bench_transcript,
                            &ccs,
                        )
                        .expect("Failed to generate linearization proof");
                    },
                    criterion::BatchSize::SmallInput,
                );
            },
        );
    }

    pub fn bench_linearization_verifier(
        group: &mut BenchmarkGroup<WallTime>,
        x_len: usize,
        n: usize,
        wit_len: usize,
        kappa: usize,
        t: R1CS,
    ) {
        group.bench_function(
            BenchmarkId::new(
                "Linearization Verifier",
                format!(
                    "Param. Kappa={}, Cols={},  B={}, L={}, B_small={}, K={}",
                    kappa,
                    { n / P::L },
                    P::B,
                    P::L,
                    P::B_SMALL,
                    P::K
                ),
            ),
            |b| {
                let mut transcript = PoseidonTranscript::<R>::default();
                let (cm_i, wit, ccs, _) = Self::setup_r1cs(x_len, n, wit_len, kappa, t);
                let (_, proof) = LFLinearizationProver::<_, PoseidonTranscript<R>>::prove(
                    &cm_i,
                    &wit,
                    &mut transcript,
                    &ccs,
                )
                .expect("Failed to generate linearization proof");

                b.iter(|| {
                    let mut transcript = PoseidonTranscript::<R>::default();
                    let _ = LFLinearizationVerifier::<_, PoseidonTranscript<R>>::verify(
                        &cm_i,
                        &proof,
                        &mut transcript,
                        &ccs,
                    )
                    .expect("Failed to verify linearization proof");
                });
            },
        );
    }

    pub fn bench_decomposition_prover(
        group: &mut BenchmarkGroup<WallTime>,
        x_len: usize,
        n: usize,
        wit_len: usize,
        kappa: usize,
        t: R1CS,
    ) {
        group.bench_function(
            BenchmarkId::new(
                "Decomposition Prover",
                format!(
                    "Param. Kappa={}, Cols={},  B={}, L={}, B_small={}, K={}",
                    kappa,
                    { n / P::L },
                    P::B,
                    P::L,
                    P::B_SMALL,
                    P::K
                ),
            ),
            |b| {
                let mut prover_transcript = PoseidonTranscript::<R>::default();
                let mut verifier_transcript = PoseidonTranscript::<R>::default();
                let (cm_i, wit, ccs, scheme) = Self::setup_r1cs(x_len, n, wit_len, kappa, t);
                let (_, proof) = LFLinearizationProver::<_, PoseidonTranscript<R>>::prove(
                    &cm_i,
                    &wit,
                    &mut prover_transcript,
                    &ccs,
                )
                .expect("Failed to generate linearization proof");
                let lcccs =
                    Self::verify_linearization_proof(&proof, &mut verifier_transcript, &cm_i, &ccs);

                b.iter(|| {
                    let _ = LFDecompositionProver::<_, PoseidonTranscript<R>>::prove::<P>(
                        &lcccs,
                        &wit,
                        &mut prover_transcript,
                        &ccs,
                        &scheme,
                    )
                    .expect("Failed to generate decomposition proof");
                });
            },
        );
    }

    pub fn bench_decomposition_verifier(
        group: &mut BenchmarkGroup<WallTime>,
        x_len: usize,
        n: usize,
        wit_len: usize,
        kappa: usize,
        t: R1CS,
    ) {
        group.bench_function(
            BenchmarkId::new(
                "Decomposition Verifier",
                format!(
                    "Param. Kappa={}, Cols={},  B={}, L={}, B_small={}, K={}",
                    kappa,
                    { n / P::L },
                    P::B,
                    P::L,
                    P::B_SMALL,
                    P::K
                ),
            ),
            |b| {
                let mut prover_transcript = PoseidonTranscript::<R>::default();
                let mut verifier_transcript = PoseidonTranscript::<R>::default();
                let (cm_i, wit, ccs, scheme) = Self::setup_r1cs(x_len, n, wit_len, kappa, t);
                let (_, linearization_proof) =
                    LFLinearizationProver::<_, PoseidonTranscript<R>>::prove(
                        &cm_i,
                        &wit,
                        &mut prover_transcript,
                        &ccs,
                    )
                    .expect("Failed to generate linearization proof");
                let lcccs = LFLinearizationVerifier::<_, PoseidonTranscript<R>>::verify(
                    &cm_i,
                    &linearization_proof,
                    &mut verifier_transcript,
                    &ccs,
                )
                .expect("Failed to verify linearization proof");

                let (_, _, _, proof) =
                    LFDecompositionProver::<_, PoseidonTranscript<R>>::prove::<P>(
                        &lcccs,
                        &wit,
                        &mut prover_transcript,
                        &ccs,
                        &scheme,
                    )
                    .expect("Failed to generate decomposition proof");

                b.iter_batched(
                    || verifier_transcript.clone(),
                    |mut bench_verifier_transcript| {
                        let _ = LFDecompositionVerifier::<_, PoseidonTranscript<R>>::verify::<P>(
                            &lcccs,
                            &proof,
                            &mut bench_verifier_transcript,
                            &ccs,
                        )
                        .expect("Failed to verify decomposition proof");
                    },
                    criterion::BatchSize::SmallInput,
                );
            },
        );
    }

    pub fn bench_folding_prover(
        group: &mut BenchmarkGroup<WallTime>,
        x_len: usize,
        n: usize,
        wit_len: usize,
        kappa: usize,
        t: R1CS,
    ) {
        group.bench_function(
            BenchmarkId::new(
                "Folding Prover",
                format!(
                    "Param. Kappa={}, Cols={},  B={}, L={}, B_small={}, K={}",
                    kappa,
                    { n / P::L },
                    P::B,
                    P::L,
                    P::B_SMALL,
                    P::K
                ),
            ),
            |b| {
                let (cm_i, wit, ccs, scheme) = Self::setup_r1cs(x_len, n, wit_len, kappa, t);
                let mut prover_transcript = PoseidonTranscript::<R>::default();
                let mut verifier_transcript = PoseidonTranscript::<R>::default();

                let (_, linearization_proof) =
                    LFLinearizationProver::<_, PoseidonTranscript<R>>::prove(
                        &cm_i,
                        &wit,
                        &mut prover_transcript,
                        &ccs,
                    )
                    .unwrap();

                let lcccs = LFLinearizationVerifier::<_, PoseidonTranscript<R>>::verify(
                    &cm_i,
                    &linearization_proof,
                    &mut verifier_transcript,
                    &ccs,
                )
                .unwrap();

                let (mz_mles, _, wit_vec, decomposition_proof) =
                    LFDecompositionProver::<_, PoseidonTranscript<R>>::prove::<P>(
                        &lcccs,
                        &wit,
                        &mut prover_transcript,
                        &ccs,
                        &scheme,
                    )
                    .unwrap();

                let lcccs_vec = LFDecompositionVerifier::<_, PoseidonTranscript<R>>::verify::<P>(
                    &lcccs,
                    &decomposition_proof,
                    &mut verifier_transcript,
                    &ccs,
                )
                .unwrap();
                let (lcccs, wit_s, mz_mles) = {
                    let mut lcccs = lcccs_vec.clone();
                    let mut lcccs_r = lcccs_vec.clone();
                    lcccs.append(&mut lcccs_r);

                    let mut wit_s = wit_vec.clone();
                    let mut wit_s_r = wit_vec;
                    wit_s.append(&mut wit_s_r);

                    let mut mz_mles_vec = mz_mles.clone();
                    let mut mz_mles_r = mz_mles;
                    mz_mles_vec.append(&mut mz_mles_r);
                    (lcccs, wit_s, mz_mles_vec)
                };

                b.iter_batched(
                    || (lcccs.clone(), wit_s.clone()),
                    |(lcccs_vec, wit_vec)| {
                        let _ = LFFoldingProver::<R, PoseidonTranscript<R>>::prove::<P>(
                            &lcccs_vec,
                            wit_vec,
                            &mut prover_transcript,
                            &ccs,
                            &mz_mles,
                        )
                        .unwrap();
                    },
                    criterion::BatchSize::SmallInput,
                );
            },
        );
    }

    pub fn bench_folding_verifier(
        group: &mut BenchmarkGroup<WallTime>,
        x_len: usize,
        n: usize,
        wit_len: usize,
        kappa: usize,
        t: R1CS,
    ) {
        group.bench_function(
            BenchmarkId::new(
                "Folding Verifier",
                format!(
                    "Param. Kappa={}, Cols={},  B={}, L={}, B_small={}, K={}",
                    kappa,
                    { n / P::L },
                    P::B,
                    P::L,
                    P::B_SMALL,
                    P::K
                ),
            ),
            |b| {
                let mut prover_transcript = PoseidonTranscript::<R>::default();
                let mut verifier_transcript = PoseidonTranscript::<R>::default();
                let (cm_i, wit, ccs, scheme) = Self::setup_r1cs(x_len, n, wit_len, kappa, t);
                let (_, linearization_proof) =
                    LFLinearizationProver::<_, PoseidonTranscript<R>>::prove(
                        &cm_i,
                        &wit,
                        &mut prover_transcript,
                        &ccs,
                    )
                    .expect("Failed to generate linearization proof");

                let lcccs = LFLinearizationVerifier::<_, PoseidonTranscript<R>>::verify(
                    &cm_i,
                    &linearization_proof,
                    &mut verifier_transcript,
                    &ccs,
                )
                .expect("Failed to verify linearization proof");

                let (mz_mles, _, wit_vec, decomposition_proof) =
                    LFDecompositionProver::<_, PoseidonTranscript<R>>::prove::<P>(
                        &lcccs,
                        &wit,
                        &mut prover_transcript,
                        &ccs,
                        &scheme,
                    )
                    .expect("Failed to generate decomposition proof");

                let lcccs_vec = LFDecompositionVerifier::<_, PoseidonTranscript<R>>::verify::<P>(
                    &lcccs,
                    &decomposition_proof,
                    &mut verifier_transcript,
                    &ccs,
                )
                .expect("Failed to verify decomposition proof");

                let (lcccs, wit_s, mz_mles) = {
                    let mut lcccs = lcccs_vec.clone();
                    let mut lcccs_r = lcccs_vec.clone();
                    lcccs.append(&mut lcccs_r);

                    let mut wit_s = wit_vec.clone();
                    let mut wit_s_r = wit_vec;
                    wit_s.append(&mut wit_s_r);

                    let mut mz_mles_vec = mz_mles.clone();
                    let mut mz_mles_r = mz_mles;
                    mz_mles_vec.append(&mut mz_mles_r);
                    (lcccs, wit_s, mz_mles_vec)
                };

                let (_, _, folding_proof) =
                    LFFoldingProver::<R, PoseidonTranscript<R>>::prove::<P>(
                        &lcccs,
                        wit_s,
                        &mut prover_transcript,
                        &ccs,
                        &mz_mles,
                    )
                    .expect("Failed to generate folding proof");

                b.iter_batched(
                    || verifier_transcript.clone(),
                    |mut bench_verifier_transcript| {
                        let _ = LFFoldingVerifier::<_, PoseidonTranscript<R>>::verify::<P>(
                            &lcccs,
                            &folding_proof,
                            &mut bench_verifier_transcript,
                            &ccs,
                        )
                        .expect("Failed to verify folding proof");
                    },
                    criterion::BatchSize::SmallInput,
                );
            },
        );
    }

    pub fn bench_e2e_prover(
        group: &mut BenchmarkGroup<WallTime>,
        x_len: usize,
        n: usize,
        wit_len: usize,
        kappa: usize,
        t: R1CS,
    ) {
        group.bench_function(
            BenchmarkId::new(
                "E2E Prover",
                format!(
                    "Param. Kappa={}, Cols={}, B={}, L={}, B_small={}, K={}",
                    kappa,
                    { n / P::L },
                    P::B,
                    P::L,
                    P::B_SMALL,
                    P::K
                ),
            ),
            |b| {
                let mut prover_transcript = PoseidonTranscript::<R>::default();
                let mut verifier_transcript = PoseidonTranscript::<R>::default();
                let (cm_i, wit, ccs, scheme) = Self::setup_r1cs(x_len, n, wit_len, kappa, t);

                let (_, linearization_proof) =
                    LFLinearizationProver::<_, PoseidonTranscript<R>>::prove(
                        &cm_i,
                        &wit,
                        &mut prover_transcript,
                        &ccs,
                    )
                    .expect("Failed to generate linearization proof");

                let lcccs = LFLinearizationVerifier::<_, PoseidonTranscript<R>>::verify(
                    &cm_i,
                    &linearization_proof,
                    &mut verifier_transcript,
                    &ccs,
                )
                .expect("Failed to verify linearization");

                b.iter_batched(
                    || prover_transcript.clone(),
                    |mut bench_prover_transcript| {
                        let _ = NIFSProver::<R, P, PoseidonTranscript<R>>::prove(
                            &lcccs,
                            &wit,
                            &cm_i,
                            &wit,
                            &mut bench_prover_transcript,
                            &ccs,
                            &scheme,
                        )
                        .expect("Failed to generate proof");
                    },
                    criterion::BatchSize::SmallInput,
                );
            },
        );
    }

    pub fn bench_e2e_verifier(
        group: &mut BenchmarkGroup<WallTime>,
        x_len: usize,
        n: usize,
        wit_len: usize,
        kappa: usize,
        t: R1CS,
    ) {
        group.bench_function(
            BenchmarkId::new(
                "E2E Verifier",
                format!(
                    "Param. Kappa={}, Cols={}, B={}, L={}, B_small={}, K={}",
                    kappa,
                    { n / P::L },
                    P::B,
                    P::L,
                    P::B_SMALL,
                    P::K
                ),
            ),
            |b| {
                let mut prover_transcript = PoseidonTranscript::<R>::default();
                let mut verifier_transcript = PoseidonTranscript::<R>::default();
                let (cm_i, wit, ccs, scheme) = Self::setup_r1cs(x_len, n, wit_len, kappa, t);

                let (acc_lcccs, linearization_proof) = LFLinearizationProver::<
                    _,
                    PoseidonTranscript<R>,
                >::prove(
                    &cm_i, &wit, &mut prover_transcript, &ccs
                )
                .expect("Failed to generate linearization proof");

                let lcccs = LFLinearizationVerifier::<_, PoseidonTranscript<R>>::verify(
                    &cm_i,
                    &linearization_proof,
                    &mut verifier_transcript,
                    &ccs,
                )
                .expect("Failed to verify linearization");

                let (_, _, proof) = NIFSProver::<R, P, PoseidonTranscript<R>>::prove(
                    &lcccs,
                    &wit,
                    &cm_i,
                    &wit,
                    &mut prover_transcript,
                    &ccs,
                    &scheme,
                )
                .expect("Failed to generate proof");

                b.iter_batched(
                    || verifier_transcript.clone(),
                    |mut bench_verifier_transcript| {
                        let result = NIFSVerifier::<R, P, PoseidonTranscript<R>>::verify(
                            &acc_lcccs,
                            &cm_i,
                            &proof,
                            &mut bench_verifier_transcript,
                            &ccs,
                        );
                        assert!(result.is_ok());
                    },
                    criterion::BatchSize::SmallInput,
                );
            },
        );
    }
}
