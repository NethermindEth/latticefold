#![allow(non_snake_case)]
use ark_std::iterable::Iterable;
use ark_std::log2;
use ark_std::marker::PhantomData;
use ark_std::sync::Arc;
use cyclotomic_rings::{rot_sum, SuitableRing};
use lattirust_ring::OverField;
use utils::get_alphas_betas_zetas_mus;

use super::error::FoldingError;
use crate::commitment::Commitment;
use crate::transcript::TranscriptWithSmallChallenges;
use crate::utils::sumcheck::{MLSumcheck, SumCheckError::SumCheckFailed};
use crate::{
    arith::{utils::mat_vec_mul, Instance, Witness, CCS, LCCCS},
    parameters::DecompositionParams,
    transcript::Transcript,
    utils::{mle::dense_vec_to_dense_mle, sumcheck},
};

use lattirust_poly::{
    mle::DenseMultilinearExtension,
    polynomials::{build_eq_x_r, eq_eval, VPAuxInfo, VirtualPolynomial},
};
use lattirust_ring::PolyRing;
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
        assert_eq!(cm_i_s.len(), 2 * P::K);
        let (m, log_m) = (ccs.m, ccs.s);

        // Step 1: Generate alpha, zeta, mu, beta challenges
        // TODO: Get challenges from big set but as NTT
        let (alpha_s, beta_s, zeta_s, mu_s) =
            get_alphas_betas_zetas_mus::<_, _, P>(log_m, transcript);

        // Step 2: Compute g polynomial and sumcheck on it
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
        let vs = cm_i_s.iter().map(|cm_i| cm_i.v).collect::<Vec<NTT>>();
        let us = cm_i_s.iter().map(|cm_i| cm_i.u.clone()).collect::<Vec<_>>();

        // Setup matrix_mles for later evaluation of etas
        // Review creation of this Mi*z mles
        let Mz_mles_vec: Vec<Vec<DenseMultilinearExtension<NTT>>> = zis
            .iter()
            .map(|zi| {
                let Mz_mle = ccs
                    .M
                    .iter()
                    .map(|M| Ok(dense_vec_to_dense_mle(log_m, &mat_vec_mul(&M, &zi)?)))
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

        // let claim_g1: NTT = alpha_s
        //     .iter()
        //     .zip(vs.iter())
        //     .map(|(&alpha_i, &v_i)| alpha_i * v_i)
        //     .sum();
        // let claim_g2 = zeta_s
        //     .iter()
        //     .zip(us.iter())
        //     .fold(NTT::zero(), |acc, (zeta, ui)| {
        //         let mut zeta_i = NTT::one();
        //         let ui_sum = ui.iter().fold(NTT::zero(), |acc, &u_i_t| {
        //             zeta_i = zeta_i * zeta;
        //             acc + (u_i_t * zeta_i)
        //         });
        //         acc + ui_sum
        //     });

        let (sum_check_proof, prover_state) = MLSumcheck::prove_as_subprotocol(transcript, &g);
        let r_0 = prover_state
            .randomness
            .into_iter()
            .map(|x| x.into())
            .collect::<Vec<NTT>>();

        // Step 3: Evaluate thetas and etas
        let thetas: Vec<NTT> = f_hat_mles
            .iter()
            .map(|f_hat_mle| {
                f_hat_mle
                    .evaluate(&r_0)
                    .ok_or(FoldingError::<NTT>::EvaluationError("f_hat".to_string()))
            })
            .collect::<Result<Vec<_>, _>>()?;

        let etas: Vec<Vec<NTT>> = Mz_mles_vec
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

        transcript.absorb_slice(&thetas);
        etas.iter().for_each(|etas| transcript.absorb_slice(&etas));

        // Step 5 get rho challenges
        let rho_s = get_rhos::<_, _, P>(transcript);

        // Step 6 compute v0, u0, y0, x_w0
        let v_0: NTT = rho_s
            .iter()
            .zip(thetas.iter())
            .map(|(&rho_i, theta_i)| NTT::from(rot_sum::<NTT>(rho_i, theta_i.coeffs())))
            .sum();

        let cm_0: Commitment<C, NTT> = rho_s
            .iter()
            .zip(cm_i_s.iter())
            .map(|(&rho_i, cm_i)| cm_i.cm.clone() * NTT::from(rho_i))
            .sum();

        let u_0: Vec<NTT> = rho_s
            .iter()
            .zip(etas.iter())
            .map(|(&rho_i, etas_i)| {
                etas_i
                    .iter()
                    .map(|etas_i_j| NTT::from(rho_i) * etas_i_j)
                    .collect::<Vec<NTT>>()
            })
            .fold(vec![NTT::zero(); ccs.l], |mut acc, rho_i_times_etas_i| {
                acc.iter_mut()
                    .zip(rho_i_times_etas_i)
                    .for_each(|(acc_j, rho_i_times_etas_i_j)| {
                        *acc_j += rho_i_times_etas_i_j;
                    });

                acc
            });

        let x_0: Vec<NTT> = rho_s
            .iter()
            .zip(cm_i_s.iter())
            .map(|(&rho_i, cm_i)| {
                cm_i.x_w
                    .iter()
                    .map(|x_w_i| NTT::from(rho_i) * x_w_i)
                    .collect::<Vec<NTT>>()
            })
            .fold(vec![NTT::zero(); ccs.n], |mut acc, rho_i_times_x_w_i| {
                acc.iter_mut()
                    .zip(rho_i_times_x_w_i)
                    .for_each(|(acc_j, rho_i_times_x_w_i)| {
                        *acc_j += rho_i_times_x_w_i;
                    });

                acc
            });

        // Step 7: Compute f0 and Witness_0

        let h = x_0.last().copied().ok_or(FoldingError::IncorrectLength)?;
        let lcccs = LCCCS {
            r: r_0,
            v: v_0,
            cm: cm_0,
            u: u_0,
            x_w: x_0,
            h,
        };

        let f_0: Vec<NTT> = rho_s.iter().zip(w_s).fold(
            vec![NTT::ZERO; w_s[0].f.len()],
            |mut acc, (&rho_i, w_i)| {
                let rho_i: NTT = rho_i.into();

                acc.into_iter()
                    .zip(w_i.f.iter())
                    .map(|(acc_j, w_ij)| acc_j + rho_i * w_ij)
                    .collect()
            },
        );

        let w_0 = Witness::from_f::<P>(f_0);

        let folding_proof = FoldingProof {
            pointshift_sumcheck_proof: sum_check_proof,
            theta_s: thetas,
            eta_s: etas,
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
        _transcript: &mut impl TranscriptWithSmallChallenges<NTT>,
        _ccs: &CCS<NTT>,
    ) -> Result<LCCCS<C, NTT>, FoldingError<NTT>> {
        let m = _ccs.m;
        let log_m = log2(m) as usize;

        todo!()

        // // Step 1: Generate alpha, zeta, mu, beta challenges
        // cm_i_s.iter().for_each(|lcccs| {
        //     _transcript.absorb_ring_vec(&lcccs.r);
        //     _transcript.absorb_ring(&lcccs.v);
        //     // _transcript.absorb_ring_vec(&lcccs.cm); Not absorbed by transcript?
        //     _transcript.absorb_ring_vec(&lcccs.u);
        //     _transcript.absorb_ring_vec(&lcccs.x_w);
        // });
        // // TODO: Get challenges from big set but as NTT
        // let alpha_s = _transcript.get_small_challenges(2 * P::K);
        // let zeta_s = _transcript.get_small_challenges(2 * P::K);
        // let mu_s = _transcript.get_small_challenges((2 * P::K) - 1); // Note is one challenge less
        // let beta_s = _transcript.get_small_challenges(log_m);

        // let poly_info = VPAuxInfo {
        //     max_degree: _ccs.d + 1,
        //     num_variables: log_m,
        //     phantom: std::marker::PhantomData,
        // };
        // let ris = cm_i_s.iter().map(|cm_i| cm_i.r.clone()).collect::<Vec<_>>();
        // let vs = cm_i_s.iter().map(|cm_i| cm_i.v).collect::<Vec<NTT>>();
        // let us = cm_i_s.iter().map(|cm_i| cm_i.u.clone()).collect::<Vec<_>>();

        // let claim_g1 = alpha_s
        //     .iter()
        //     .zip(vs.iter())
        //     .fold(NTT::zero(), |acc, (&alpha, &vi)| acc + (alpha * vi));
        // let claim_g2 = zeta_s
        //     .iter()
        //     .zip(us.iter())
        //     .fold(NTT::zero(), |acc, (zeta, ui)| {
        //         let mut zeta_i = NTT::one();
        //         let ui_sum = ui.iter().fold(NTT::zero(), |acc, &u_i_t| {
        //             zeta_i = zeta_i * zeta;
        //             acc + (u_i_t * zeta_i)
        //         });
        //         acc + ui_sum
        //     });

        // //Step 2: The sumcheck.

        // // Verify the sumcheck proof.
        // let sub_claim = MLSumcheck::verify_as_subprotocol(
        //     _transcript,
        //     &poly_info,
        //     claim_g1 + claim_g2,
        //     &proof.pointshift_sumcheck_proof,
        // )?;

        // let point_r = sub_claim
        //     .point
        //     .into_iter()
        //     .map(|x| NTT::field_to_base_ring(&x).into())
        //     .collect::<Vec<NTT>>();

        // let e_asterisk = eq_eval(&beta_s, &point_r).unwrap();
        // let e_i_s: Vec<NTT> = ris
        //     .iter()
        //     .map(|r| eq_eval(r.as_slice(), &point_r).unwrap())
        //     .collect::<Vec<_>>();
        // let s = sub_claim.expected_evaluation.clone();

        // let mut should_equal_s =
        //     mu_s.iter()
        //         .zip(proof.theta_s.iter())
        //         .fold(NTT::zero(), |acc, (mu_i, &theta_i)| {
        //             let mut thetas_mul = theta_i;
        //             for j in 1..P::B_SMALL {
        //                 thetas_mul = thetas_mul * (theta_i - NTT::from(j));
        //                 thetas_mul = thetas_mul * (theta_i + NTT::from(j));
        //             }
        //             acc + (thetas_mul * mu_i)
        //         });
        // let last_theta = proof
        //     .theta_s
        //     .last()
        //     .cloned()
        //     .ok_or(FoldingError::IncorrectLength)?;
        // // Recall last mu = 1
        // let mut theta_mul = last_theta;
        // for j in 1..P::B_SMALL {
        //     theta_mul = theta_mul * (last_theta - NTT::from(j));
        //     theta_mul = theta_mul * (last_theta + NTT::from(j));
        // }
        // should_equal_s = should_equal_s + theta_mul;
        // should_equal_s = should_equal_s * e_asterisk;

        // alpha_s
        //     .iter()
        //     .zip(zeta_s)
        //     .zip(proof.theta_s.iter())
        //     .zip(proof.eta_s.iter())
        //     .zip(e_i_s.iter())
        //     .for_each(|((((&alpha_i, zeta_i), theta_i), eta_i), e_i)| {
        //         let mut zeta_i_t = NTT::one();
        //         let combined_eta_i_t = eta_i.iter().fold(NTT::zero(), |acc, &eta_t| {
        //             zeta_i_t = zeta_i_t * zeta_i;
        //             acc + (eta_t * zeta_i_t)
        //         });
        //         let combined_theta_and_eta = (alpha_i * theta_i) + combined_eta_i_t;
        //         should_equal_s += combined_theta_and_eta * e_i;
        //     });

        // match should_equal_s == s {
        //     true => {}
        //     false => {
        //         return Err(FoldingError::SumCheckError(SumCheckFailed(
        //             should_equal_s,
        //             s,
        //         )));
        //     }
        // }
        // // check claim and output o, u0, x0, y0
        // let rhos = _transcript.get_small_challenges((2 * P::K) - 1); // Note that we are missing the first element

        // // Step 6 compute v0, u0, y0, x_w0
        // let v_0: NTT = rhos.iter().zip(proof.theta_s.iter().skip(1)).fold(
        //     proof.theta_s[0],
        //     |acc, (rho_i, theta_i)| {
        //         // acc + rho_i.rot_sum(theta_i) // Note that theta_i is already in NTT form
        //         todo!() // Add WithRot to OverField in lattirust
        //     },
        // ); // Do INTT here

        // let (y_0, u_0, x_0) = rhos
        //     .iter()
        //     .zip(cm_i_s.iter().zip(proof.eta_s.iter()).skip(1))
        //     .fold(
        //         (
        //             cm_i_s[0].cm.clone(),
        //             proof.eta_s[0].clone(),
        //             cm_i_s[0].x_w.clone(),
        //         ),
        //         |(acc_y, acc_u, acc_x), (rho_i, (cm_i, eta))| {
        //             let y = acc_y + (cm_i.cm.clone() * rho_i);
        //             let u = acc_u
        //                 .iter()
        //                 .zip(eta.iter())
        //                 .map(|(&a, &e)| a + (e * rho_i))
        //                 .collect();
        //             let x = acc_x
        //                 .iter()
        //                 .zip(cm_i.x_w.iter())
        //                 .map(|(&a, &x)| a + (x * rho_i))
        //                 .collect();
        //             (y, u, x)
        //         },
        //     );

        // let h = x_0.last().cloned().ok_or(FoldingError::IncorrectLength)?;
        // Ok(LCCCS {
        //     r: point_r,
        //     v: v_0,
        //     cm: y_0,
        //     u: u_0,
        //     x_w: x_0,
        //     h,
        // })
    }
}
// fn create_matrix_mle<NTT: OverField>(
//     log_m: usize,
//     Mi: &Vec<Vec<NTT>>,
//     zi: &Vec<NTT>,
// ) -> DenseMultilinearExtension<NTT> {
//     let zero_vector = usize_to_binary_vector::<NTT>(0, log2(Mi.len()) as usize);
//     let mle_z_ccs_b = mle_val_from_vector(&zi, &zero_vector);
//     let evaluations: Vec<NTT> = mle_matrix_to_val_eval_second(&Mi, &zero_vector)
//         .iter()
//         .map(|val| *val * mle_z_ccs_b)
//         .collect();
//     let mle = DenseMultilinearExtension::from_evaluations_vec(log_m, evaluations);

//     let matrix_mle = (1..Mi.len())
//         .into_iter()
//         .map(|i| usize_to_binary_vector::<NTT>(i, log2(Mi.len()) as usize))
//         .fold(mle, |acc, b| {
//             let mle_z_ccs_b = mle_val_from_vector(&zi, &b);
//             let evaluations: Vec<NTT> = utils::mle_matrix_to_val_eval_second(&Mi, &b)
//                 .iter()
//                 .map(|val| *val * mle_z_ccs_b)
//                 .collect();
//             let mle = DenseMultilinearExtension::from_evaluations_vec(log_m, evaluations);
//             acc + mle
//         });
//     matrix_mle
// }

#[cfg(test)]
mod tests {
    use lattirust_arithmetic::{
        challenge_set::latticefold_challenge_set::BinarySmallSet,
        ring::{Pow2CyclotomicPolyRing, Pow2CyclotomicPolyRingNTT, Zq},
    };
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

    // Boilerplate code to generate values needed for testing
    const Q: u64 = 17; // Replace with an appropriate modulus
    const N: usize = 8;
    type CR = Pow2CyclotomicPolyRing<Zq<Q>, N>;
    type NTT = Pow2CyclotomicPolyRingNTT<Q, N>;
    type CS = BinarySmallSet<Q, N>;
    type T = PoseidonTranscript<Pow2CyclotomicPolyRingNTT<Q, N>, CS>;

    #[derive(Clone)]
    struct PP;

    impl DecompositionParams for PP {
        const B: u128 = 1_024;
        const L: usize = 1;
        const B_SMALL: u128 = 2;
        const K: usize = 10;
    }

    #[test]
    fn test_folding() {
        let ccs = get_test_ccs::<NTT>();
        let (_, x_ccs, w_ccs) = get_test_z_split::<NTT>(3);
        let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());
        let wit: Witness<NTT> = Witness::from_w_ccs::<CR, PP>(w_ccs);
        let cm_i: CCCS<4, NTT> = CCCS {
            cm: wit.commit::<4, 4, CR, PP>(&scheme).unwrap(),
            x_ccs,
        };

        let mut prover_transcript = PoseidonTranscript::<NTT, CS>::default();
        let mut verifier_transcript = PoseidonTranscript::<NTT, CS>::default();

        let (_, linearization_proof) =
            LFLinearizationProver::<_, T>::prove(&cm_i, &wit, &mut prover_transcript, &ccs)
                .unwrap();

        let lcccs = LFLinearizationVerifier::<_, PoseidonTranscript<NTT, CS>>::verify(
            &cm_i,
            &linearization_proof,
            &mut verifier_transcript,
            &ccs,
        )
        .unwrap();

        let (_, vec_wit, decomposition_proof) =
            LFDecompositionProver::<_, T>::prove::<4, 4, CR, PP>(
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
        let (_, _, folding_proof) = LFFoldingProver::<_, T>::prove::<4, CR, PP>(
            &lcccs,
            &wit_s,
            &mut prover_transcript,
            &ccs,
        )
        .unwrap();

        let res = LFFoldingVerifier::<_, T>::verify::<4, PP>(
            &lcccs,
            &folding_proof,
            &mut verifier_transcript,
            &ccs,
        );

        assert!(res.is_ok())
    }

    #[test]
    fn test_failing_folding() {
        let ccs = get_test_ccs::<NTT>();
        let (_, x_ccs, w_ccs) = get_test_z_split::<NTT>(3);
        let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());
        let wit: Witness<NTT> = Witness::from_w_ccs::<CR, PP>(w_ccs.clone());
        let cm_i: CCCS<4, NTT> = CCCS {
            cm: wit.commit::<4, 4, CR, PP>(&scheme).unwrap(),
            x_ccs,
        };

        let mut prover_transcript = PoseidonTranscript::<NTT, CS>::default();
        let mut verifier_transcript = PoseidonTranscript::<NTT, CS>::default();

        let (_, linearization_proof) =
            LFLinearizationProver::<_, T>::prove(&cm_i, &wit, &mut prover_transcript, &ccs)
                .unwrap();

        let lcccs = LFLinearizationVerifier::<_, PoseidonTranscript<NTT, CS>>::verify(
            &cm_i,
            &linearization_proof,
            &mut verifier_transcript,
            &ccs,
        )
        .unwrap();

        let (_, mut vec_wit, decomposition_proof) =
            LFDecompositionProver::<_, T>::prove::<4, 4, CR, PP>(
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

        vec_wit[0] = Witness::<NTT>::from_w_ccs::<CR, PP>(w_ccs);

        let (_, _, folding_proof) = LFFoldingProver::<_, T>::prove::<4, CR, PP>(
            &vec_lcccs,
            &vec_wit,
            &mut prover_transcript,
            &ccs,
        )
        .unwrap();

        let res = LFFoldingVerifier::<_, T>::verify::<4, PP>(
            &vec_lcccs,
            &folding_proof,
            &mut verifier_transcript,
            &ccs,
        );

        assert!(res.is_err())
    }
}
