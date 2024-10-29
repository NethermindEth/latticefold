#![allow(non_snake_case)]
use ark_std::iter::successors;
use ark_std::iterable::Iterable;
use ark_std::marker::PhantomData;
use cyclotomic_rings::SuitableRing;
use lattirust_ring::OverField;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use std::{f64, time::Instant};
use utils::get_alphas_betas_zetas_mus;

use super::error::FoldingError;
use crate::transcript::TranscriptWithSmallChallenges;
use crate::utils::sumcheck::{MLSumcheck, SumCheckError::SumCheckFailed};
use crate::{
    arith::{utils::mat_vec_mul, Instance, Witness, CCS, LCCCS},
    parameters::DecompositionParams,
    utils::{mle::dense_vec_to_dense_mle, sumcheck},
};

use lattirust_poly::{
    mle::DenseMultilinearExtension,
    polynomials::{eq_eval, VPAuxInfo},
};
use utils::*;

mod utils;

#[derive(Clone)]
pub struct FoldingProof<NTT: OverField> {
    // Step 2.
    pub pointshift_sumcheck_proof: sumcheck::Proof<NTT>,
    // Step 3
    pub theta_s: Vec<NTT>,
    pub eta_s: Vec<Vec<NTT>>,
}

pub trait FoldingProver<NTT: SuitableRing, T: TranscriptWithSmallChallenges<NTT>> {
    fn prove<const C: usize, P: DecompositionParams>(
        cm_i_s: &[LCCCS<C, NTT>],
        w_s: &[Witness<NTT>],
        transcript: &mut impl TranscriptWithSmallChallenges<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(LCCCS<C, NTT>, Witness<NTT>, FoldingProof<NTT>), FoldingError<NTT>>;
}

pub trait FoldingVerifier<NTT: SuitableRing, T: TranscriptWithSmallChallenges<NTT>> {
    fn verify<const C: usize, P: DecompositionParams>(
        cm_i_s: &[LCCCS<C, NTT>],
        proof: &FoldingProof<NTT>,
        transcript: &mut impl TranscriptWithSmallChallenges<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<LCCCS<C, NTT>, FoldingError<NTT>>;
}

pub struct LFFoldingProver<NTT, T> {
    _ntt: PhantomData<NTT>,
    _t: PhantomData<T>,
}

pub struct LFFoldingVerifier<NTT, T> {
    _ntt: PhantomData<NTT>,
    _t: PhantomData<T>,
}

impl<NTT: SuitableRing, T: TranscriptWithSmallChallenges<NTT>> FoldingProver<NTT, T>
    for LFFoldingProver<NTT, T>
{
    fn prove<const C: usize, P: DecompositionParams>(
        cm_i_s: &[LCCCS<C, NTT>],
        w_s: &[Witness<NTT>],
        transcript: &mut impl TranscriptWithSmallChallenges<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(LCCCS<C, NTT>, Witness<NTT>, FoldingProof<NTT>), FoldingError<NTT>> {
        if cm_i_s.len() != 2 * P::K {
            return Err(FoldingError::IncorrectLength);
        }

        let log_m = ccs.s;

        // Step 1: Generate alpha, zeta, mu, beta challenges
        let folding_start = Instant::now();
        let challenges_start = Instant::now();
        let (alpha_s, beta_s, zeta_s, mu_s) =
            get_alphas_betas_zetas_mus::<_, _, P>(log_m, transcript);
        let challenges_end = challenges_start.elapsed();

        // Step 2: Compute g polynomial and sumcheck on it
        let g_start = Instant::now();
        // Setup f_hat_mle for later evaluation of thetas
        let f_hat_mles = w_s
            .iter()
            .map(|w| DenseMultilinearExtension::from_evaluations_slice(log_m, &w.f_hat))
            .collect::<Vec<DenseMultilinearExtension<NTT>>>();

        let zis = cm_i_s
            .iter()
            .zip(w_s.iter())
            .map(|(cm_i, w_i)| cm_i.get_z_vector(&w_i.w_ccs))
            .collect::<Vec<_>>();
        let ris = cm_i_s.iter().map(|cm_i| cm_i.r.clone()).collect::<Vec<_>>();

        // Setup matrix_mles for later evaluation of etas
        // Review creation of this Mi*z mles
        let Mz_mles_vec: Vec<Vec<DenseMultilinearExtension<NTT>>> = zis
            .iter()
            .map(|zi| {
                let Mz_mle = ccs
                    .M
                    .iter()
                    .map(|M| Ok(dense_vec_to_dense_mle(log_m, &mat_vec_mul(M, zi)?)))
                    .collect::<Result<_, FoldingError<_>>>()?;
                Ok(Mz_mle)
            })
            .collect::<Result<_, FoldingError<_>>>()?;

        let g = create_sumcheck_polynomial::<_, P>(
            log_m,
            &f_hat_mles,
            &alpha_s,
            &Mz_mles_vec,
            &zeta_s,
            &ris,
            &beta_s,
            &mu_s,
        )?;

        let (sum_check_proof, prover_state) = MLSumcheck::prove_as_subprotocol(transcript, &g);
        let g_end = g_start.elapsed();
        let r_0 = prover_state
            .randomness
            .into_iter()
            .map(|x| x.into())
            .collect::<Vec<NTT>>();

        // Step 3: Evaluate thetas and etas
        let theta_s_start = Instant::now();
        let theta_s: Vec<NTT> = f_hat_mles
            .iter()
            .map(|f_hat_mle| {
                f_hat_mle
                    .evaluate(&r_0)
                    .ok_or(FoldingError::<NTT>::EvaluationError("f_hat".to_string()))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let theta_s_end = theta_s_start.elapsed();

        let eta_s_start = Instant::now();
        let eta_s: Vec<Vec<NTT>> = Mz_mles_vec
            .iter()
            .map(|Mz_mles| {
                Mz_mles
                    .iter()
                    .map(|mle| {
                        mle.evaluate(r_0.as_slice())
                            .ok_or(FoldingError::<NTT>::EvaluationError("Mz".to_string()))
                    })
                    .collect::<Result<Vec<_>, _>>()
            })
            .collect::<Result<Vec<_>, _>>()?;
        transcript.absorb_slice(&theta_s);
        eta_s.iter().for_each(|etas| transcript.absorb_slice(etas));
        let eta_s_end = eta_s_start.elapsed();

        // Step 5 ge    t rho challenges
        let rho_s_start = Instant::now();
        let rho_s = get_rhos::<_, _, P>(transcript);
        let rho_s_end = rho_s_start.elapsed();

        // Step 6 compute v0, u0, y0, x_w0
        let v0_u0_x0_cm0_start = Instant::now();    
        let (v_0, cm_0, u_0, x_0) = compute_v0_u0_x0_cm_0(&rho_s, &theta_s, cm_i_s, &eta_s, ccs);
        let v0_u0_x0_cm0_end = v0_u0_x0_cm0_start.elapsed();
        // Step 7: Compute f0 and Witness_0
        let f0_wit0_start = Instant::now();
        let h = x_0.last().copied().ok_or(FoldingError::IncorrectLength)?;
        let lcccs = LCCCS {
            r: r_0,
            v: v_0,
            cm: cm_0,
            u: u_0,
            x_w: x_0[0..x_0.len() - 1].to_vec(),
            h,
        };

        let f_0: Vec<NTT> =
            rho_s
                .iter()
                .zip(w_s)
                .fold(vec![NTT::ZERO; w_s[0].f.len()], |acc, (&rho_i, w_i)| {
                    let rho_i: NTT = rho_i.into();

                    acc.into_iter()
                        .zip(w_i.f.iter())
                        .map(|(acc_j, w_ij)| acc_j + rho_i * w_ij)
                        .collect()
                });

        let w_0 = Witness::from_f::<P>(f_0);
        let f0_wit0_end = f0_wit0_start.elapsed();

        let folding_end = folding_start.elapsed();
        println!("Folding Prover:\n\
        folding: {:?} ({:.2}%)\n\
        challenges: {:?} ({:.2}%)\n\
        g: {:?} ({:.2}%)\n\
        theta_s: {:?} ({:.2}%)\n\
        eta_s: {:?} ({:.2}%)\n\
        f0_wit0: {:?} ({:.2}%)",
            folding_end,
            (folding_end.as_micros() as f64 / folding_end.as_micros() as f64) * 100.0,
            challenges_end,
            (challenges_end.as_micros() as f64 / folding_end.as_micros() as f64) * 100.0,
            g_end,
            (g_end.as_micros() as f64 / folding_end.as_micros() as f64) * 100.0,
            theta_s_end,
            (theta_s_end.as_micros() as f64 / folding_end.as_micros() as f64) * 100.0,
            eta_s_end,
            (eta_s_end.as_micros() as f64 / folding_end.as_micros() as f64) * 100.0,
            f0_wit0_end,
            (f0_wit0_end.as_micros() as f64 / folding_end.as_micros() as f64) * 100.0
        );

        let folding_proof = FoldingProof {
            pointshift_sumcheck_proof: sum_check_proof,
            theta_s,
            eta_s,
        };

        Ok((lcccs, w_0, folding_proof))
    }
}

impl<NTT: SuitableRing, T: TranscriptWithSmallChallenges<NTT>> FoldingVerifier<NTT, T>
    for LFFoldingVerifier<NTT, T>
{
    fn verify<const C: usize, P: DecompositionParams>(
        cm_i_s: &[LCCCS<C, NTT>],
        proof: &FoldingProof<NTT>,
        transcript: &mut impl TranscriptWithSmallChallenges<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<LCCCS<C, NTT>, FoldingError<NTT>> {
        if cm_i_s.len() != 2 * P::K {
            return Err(FoldingError::IncorrectLength);
        }

        let log_m = ccs.s;

        // Step 1: Generate alpha, zeta, mu, beta challenges
        let (alpha_s, beta_s, zeta_s, mu_s) =
            get_alphas_betas_zetas_mus::<_, _, P>(log_m, transcript);

        let poly_info = VPAuxInfo::new(log_m, 2 * P::B_SMALL);

        let ris = cm_i_s.iter().map(|cm_i| cm_i.r.clone()).collect::<Vec<_>>();
        let vs = cm_i_s.iter().map(|cm_i| cm_i.v).collect::<Vec<NTT>>();
        let us = cm_i_s.iter().map(|cm_i| cm_i.u.clone()).collect::<Vec<_>>();

        let claim_g1: NTT = alpha_s
            .iter()
            .zip(vs.iter())
            .map(|(&alpha_i, &v_i)| alpha_i * v_i)
            .sum();
        let claim_g3: NTT = zeta_s
            .iter()
            .zip(us.iter())
            .map(|(&zeta_i, ui)| {
                successors(Some(zeta_i), |&zeta| Some(zeta * zeta_i))
                    .zip(ui.iter())
                    .map(|(pow_of_zeta, u_i_j)| pow_of_zeta * u_i_j)
                    .sum::<NTT>()
            })
            .sum();

        //Step 2: The sumcheck.

        // Verify the sumcheck proof.
        let sub_claim = MLSumcheck::verify_as_subprotocol(
            transcript,
            &poly_info,
            claim_g1 + claim_g3,
            &proof.pointshift_sumcheck_proof,
        )?;

        let r_0 = sub_claim
            .point
            .into_iter()
            .map(|x| x.into())
            .collect::<Vec<NTT>>();

        let e_asterisk = eq_eval(&beta_s, &r_0)?;
        let e_s: Vec<NTT> = ris
            .iter()
            .map(|r_i: &Vec<NTT>| eq_eval(r_i, &r_0))
            .collect::<Result<Vec<_>, _>>()?;

        let should_equal_s: NTT = compute_sumcheck_claim_expected_value::<NTT, P>(
            &alpha_s,
            &mu_s,
            &proof.theta_s,
            e_asterisk,
            &e_s,
            &zeta_s,
            &proof.eta_s,
        );

        if should_equal_s != sub_claim.expected_evaluation {
            return Err(FoldingError::SumCheckError(SumCheckFailed(
                should_equal_s,
                sub_claim.expected_evaluation,
            )));
        }

        // Step 5
        transcript.absorb_slice(&proof.theta_s);
        proof
            .eta_s
            .iter()
            .for_each(|etas| transcript.absorb_slice(etas));
        let rho_s = get_rhos::<_, _, P>(transcript);

        // Step 6
        let (v_0, cm_0, u_0, x_0) =
            compute_v0_u0_x0_cm_0(&rho_s, &proof.theta_s, cm_i_s, &proof.eta_s, ccs);

        // Step 7: Compute f0 and Witness_0

        let h = x_0.last().copied().ok_or(FoldingError::IncorrectLength)?;

        Ok(LCCCS {
            r: r_0,
            v: v_0,
            cm: cm_0,
            u: u_0,
            x_w: x_0[0..x_0.len() - 1].to_vec(),
            h,
        })
    }
}

fn check_ring_modulus_128_bits_security(
    ring_modulus: &BigUint,
    kappa: usize,
    degree: usize,
    num_cols: usize,
    b: u128,
    l: usize,
) -> bool {
    // Calculate the logarithm of stark_modulus
    let ring_modulus_log2 = ring_modulus.bits() as f64;
    let ring_modulus_half = ring_modulus / 2u32;

    // Calculate the left side of the inequality
    let bound_l2 = 2f64.powf(
        2.0 * (1.0045f64.ln() / 2f64.ln()).sqrt()
            * (degree as f64 * kappa as f64 * ring_modulus_log2).sqrt(),
    );
    let bound_l2_ceil = bound_l2.ceil() as u64; // Ceil and convert to u64
    let bound_l2_bigint = BigUint::from(bound_l2_ceil); // Convert to BigUint
    let bound_l2_check = bound_l2_bigint < ring_modulus_half;
    // Calculate bound_inf
    let bound_inf = bound_l2 / ((degree as f64 * num_cols as f64).sqrt());

    let b_check = b.to_f64().unwrap() < bound_inf;
    // Calculate the right side of the inequality
    // Check if b^l > stark_modulus/2
    let b_bigint = BigUint::from(b);
    let b_pow_l = b_bigint.pow(l as u32);
    let b_pow_l_check = b_pow_l > ring_modulus_half;

    // Return the result of the condition
    bound_l2_check && b_check && b_pow_l_check
}

#[cfg(test)]
mod tests {
    use rand::thread_rng;

    use crate::{
        arith::{r1cs::tests::get_test_z_split, tests::get_test_ccs, Witness, CCCS},
        commitment::AjtaiCommitmentScheme,
        nifs::{
            decomposition::{
                DecompositionProver, DecompositionVerifier, LFDecompositionProver,
                LFDecompositionVerifier,
            },
            folding::{FoldingProver, FoldingVerifier, LFFoldingProver, LFFoldingVerifier},
            linearization::{
                LFLinearizationProver, LFLinearizationVerifier, LinearizationProver,
                LinearizationVerifier,
            },
        },
        parameters::DecompositionParams,
        transcript::poseidon::PoseidonTranscript,
    };
    use cyclotomic_rings::{StarkChallengeSet, StarkRingNTT};

    // Boilerplate code to generate values needed for testing
    type R = StarkRingNTT;
    type CS = StarkChallengeSet;
    type T = PoseidonTranscript<StarkRingNTT, CS>;

    #[derive(Clone)]
    struct PP;

    impl DecompositionParams for PP {
        const B: u128 = 1_024;
        const L: usize = 1;
        const B_SMALL: usize = 2;
        const K: usize = 10;
    }

    #[test]
    fn test_folding() {
        const WIT_LEN: usize = 4; // 4 is the length of witness in this (Vitalik's) example
        const W: usize = WIT_LEN * PP::L; // the number of columns of the Ajtai matrix

        let ccs = get_test_ccs::<R>(W);
        let (_, x_ccs, w_ccs) = get_test_z_split::<R>(3);
        let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());
        let wit: Witness<R> = Witness::from_w_ccs::<PP>(&w_ccs);
        let cm_i: CCCS<4, R> = CCCS {
            cm: wit.commit::<4, 4, PP>(&scheme).unwrap(),
            x_ccs,
        };

        let mut prover_transcript = PoseidonTranscript::<R, CS>::default();
        let mut verifier_transcript = PoseidonTranscript::<R, CS>::default();

        let (_, linearization_proof) =
            LFLinearizationProver::<_, T>::prove(&cm_i, &wit, &mut prover_transcript, &ccs)
                .unwrap();

        let lcccs = LFLinearizationVerifier::<_, PoseidonTranscript<R, CS>>::verify(
            &cm_i,
            &linearization_proof,
            &mut verifier_transcript,
            &ccs,
        )
        .unwrap();

        let (_, vec_wit, decomposition_proof) = LFDecompositionProver::<_, T>::prove::<4, 4, PP>(
            &lcccs,
            &wit,
            &mut prover_transcript,
            &ccs,
            &scheme,
        )
        .unwrap();

        let vec_lcccs = LFDecompositionVerifier::<_, T>::verify::<4, PP>(
            &lcccs,
            &decomposition_proof,
            &mut verifier_transcript,
            &ccs,
        )
        .unwrap();
        let (lcccs, wit_s) = {
            let mut lcccs = vec_lcccs.clone();
            let mut lcccs_r = vec_lcccs;
            lcccs.append(&mut lcccs_r);

            let mut wit_s = vec_wit.clone();
            let mut wit_s_r = vec_wit;
            wit_s.append(&mut wit_s_r);

            (lcccs, wit_s)
        };
        let (lcccs_prover, _, folding_proof) =
            LFFoldingProver::<_, T>::prove::<4, PP>(&lcccs, &wit_s, &mut prover_transcript, &ccs)
                .unwrap();

        let lcccs_verifier = LFFoldingVerifier::<_, T>::verify::<4, PP>(
            &lcccs,
            &folding_proof,
            &mut verifier_transcript,
            &ccs,
        )
        .unwrap();

        assert_eq!(lcccs_prover, lcccs_verifier);
    }

    #[test]
    fn test_failing_folding_prover() {
        const WIT_LEN: usize = 4; // 4 is the length of witness in this (Vitalik's) example
        const W: usize = WIT_LEN * PP::L; // the number of columns of the Ajtai matrix

        let ccs = get_test_ccs::<R>(W);
        let (_, x_ccs, w_ccs) = get_test_z_split::<R>(3);
        let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());
        let wit: Witness<R> = Witness::from_w_ccs::<PP>(&w_ccs);
        let cm_i: CCCS<4, R> = CCCS {
            cm: wit.commit::<4, 4, PP>(&scheme).unwrap(),
            x_ccs,
        };

        let mut prover_transcript = PoseidonTranscript::<R, CS>::default();
        let mut verifier_transcript = PoseidonTranscript::<R, CS>::default();

        let (_, linearization_proof) =
            LFLinearizationProver::<_, T>::prove(&cm_i, &wit, &mut prover_transcript, &ccs)
                .unwrap();

        let lcccs = LFLinearizationVerifier::<_, PoseidonTranscript<R, CS>>::verify(
            &cm_i,
            &linearization_proof,
            &mut verifier_transcript,
            &ccs,
        )
        .unwrap();

        let (_, mut vec_wit, decomposition_proof) =
            LFDecompositionProver::<_, T>::prove::<4, 4, PP>(
                &lcccs,
                &wit,
                &mut prover_transcript,
                &ccs,
                &scheme,
            )
            .unwrap();

        let vec_lcccs = LFDecompositionVerifier::<_, T>::verify::<4, PP>(
            &lcccs,
            &decomposition_proof,
            &mut verifier_transcript,
            &ccs,
        )
        .unwrap();

        vec_wit[0] = Witness::<R>::from_w_ccs::<PP>(&w_ccs);

        let res = LFFoldingProver::<_, T>::prove::<4, PP>(
            &vec_lcccs,
            &vec_wit,
            &mut prover_transcript,
            &ccs,
        );

        assert!(res.is_err())
    }
}

#[cfg(test)]
mod tests_stark {
    use ark_ff::UniformRand;
    use lattirust_poly::mle::DenseMultilinearExtension;
    use lattirust_ring::{
        cyclotomic_ring::models::stark_prime::{Fq, RqNTT},
        PolyRing,
    };
    use num_bigint::BigUint;
    use rand::thread_rng;

    use crate::arith::tests::{get_test_ccs, get_test_dummy_ccs};
    use crate::{
        arith::r1cs::tests::{get_test_dummy_z_split, get_test_z_split},
        nifs::{
            decomposition::{
                DecompositionProver, DecompositionVerifier, LFDecompositionProver,
                LFDecompositionVerifier,
            },
            folding::{
                check_ring_modulus_128_bits_security, FoldingProver, FoldingVerifier,
                LFFoldingProver, LFFoldingVerifier,
            },
            linearization::LinearizationProver,
        },
    };
    use crate::{
        arith::{Witness, CCCS},
        commitment::AjtaiCommitmentScheme,
        nifs::linearization::{
            LFLinearizationProver, LFLinearizationVerifier, LinearizationVerifier,
        },
        parameters::DecompositionParams,
        transcript::poseidon::PoseidonTranscript,
    };
    use cyclotomic_rings::StarkChallengeSet;

    #[test]
    fn test_dummy_folding() {
        #[cfg(feature = "dhat-heap")]
        #[global_allocator]
        static ALLOC: dhat::Alloc = dhat::Alloc;

        type R = RqNTT;
        type CS = StarkChallengeSet;
        type T = PoseidonTranscript<R, CS>;

        let stark_modulus = BigUint::parse_bytes(
            b"3618502788666131000275863779947924135206266826270938552493006944358698582017",
            10,
        )
        .expect("Failed to parse stark_modulus");

        if check_ring_modulus_128_bits_security(&stark_modulus, C, 16, W, PP::B, PP::L) {
            println!(" Bound condition satisfied");
        } else {
            println!("Bound condition not satisfied");
        }

        #[derive(Clone)]
        struct PP;
        impl DecompositionParams for PP {
            const B: u128 = 3010936384;
            const L: usize = 8;
            const B_SMALL: usize = 38;
            const K: usize = 6;
        }

        const C: usize = 15;
        const IO: usize = 1;
        const WIT_LEN: usize = 512;
        const W: usize = WIT_LEN * PP::L; // the number of columns of the Ajtai matrix
        let r1cs_rows_size = IO + WIT_LEN + 1; // Let's have a square matrix

        #[cfg(feature = "dhat-heap")]
        let _profiler = dhat::Profiler::new_heap(); // Move a round to measure specific parts

        let ccs = get_test_dummy_ccs::<R, IO, WIT_LEN, W>(r1cs_rows_size);
        let (_, x_ccs, w_ccs) = get_test_dummy_z_split::<R, IO, WIT_LEN>();
        let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());

        let wit = Witness::from_w_ccs::<PP>(&w_ccs);
        let cm_i = CCCS {
            cm: wit.commit::<C, W, PP>(&scheme).unwrap(),
            x_ccs,
        };

        let mut prover_transcript = PoseidonTranscript::<R, CS>::default();

        let linearization_proof =
            LFLinearizationProver::<_, T>::prove(&cm_i, &wit, &mut prover_transcript, &ccs);

        let mut verifier_transcript = PoseidonTranscript::<R, CS>::default();

        let linearization_verification = match linearization_proof {
            Ok(res) => {
                println!("Linearization proof generated with success");
                LFLinearizationVerifier::<_, PoseidonTranscript<R, CS>>::verify(
                    &cm_i,
                    &res.1,
                    &mut verifier_transcript,
                    &ccs,
                )
            }
            Err(e) => panic!("Linearization proof generation error: {:?}", e),
        };

        match &linearization_verification {
            Ok(_) => println!("Linearization verified with success"),
            Err(ref e) => println!("Linearization Verification error: {:?}", e),
        };

        let lcccs = linearization_verification.unwrap();

        let decomposition_prover = LFDecompositionProver::<_, T>::prove::<W, C, PP>(
            &lcccs,
            &wit,
            &mut prover_transcript,
            &ccs,
            &scheme,
        );

        let decomposition_proof = match decomposition_prover {
            Ok(res) => {
                println!("Decomposition proof generated with success");
                res
            }
            Err(e) => panic!("Decomposition proof generation error: {:?}", e),
        };

        let decomposition_verification = LFDecompositionVerifier::<_, T>::verify::<C, PP>(
            &lcccs,
            &decomposition_proof.2,
            &mut verifier_transcript,
            &ccs,
        );

        match decomposition_verification {
            Ok(_) => println!("Decomposition verified with success"),
            Err(ref e) => println!("Decomposition Verification error: {:?}", e),
        };

        let lcccs = decomposition_verification.unwrap();

        #[cfg(feature = "dhat-heap")]
        let _profiler = dhat::Profiler::new_heap();
        let (lcccs, wit_s) = {
            let mut lcccs = lcccs.clone();
            let mut lcccs_r = lcccs.clone();
            lcccs.append(&mut lcccs_r);

            let mut wit_s = decomposition_proof.1.clone();
            let mut wit_s_r = decomposition_proof.1;
            wit_s.append(&mut wit_s_r);

            (lcccs, wit_s)
        };
        let folding_prover =
            LFFoldingProver::<_, T>::prove::<C, PP>(&lcccs, &wit_s, &mut prover_transcript, &ccs);

        let folding_proof = match folding_prover {
            Ok(res) => {
                println!("Folding proof generated with success");
                res
            }
            Err(e) => panic!("Folding proof generation error: {:?}", e),
        };

        let folding_verification = LFFoldingVerifier::<_, T>::verify::<C, PP>(
            &lcccs,
            &folding_proof.2,
            &mut verifier_transcript,
            &ccs,
        );

        match folding_verification {
            Ok(_) => println!("Folding verified with success"),
            Err(ref e) => println!("Folding Verification error: {:?}", e),
        };
    }
}
