#![allow(non_snake_case)]

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::cfg_iter;
use ark_std::iter::successors;
use ark_std::iterable::Iterable;
use ark_std::marker::PhantomData;
use cyclotomic_rings::rings::SuitableRing;
use lattirust_ring::{cyclotomic_ring::CRT, OverField};
use utils::get_alphas_betas_zetas_mus;

use super::common_helpers::{to_mles, to_mles_err};
use super::error::FoldingError;
use crate::transcript::TranscriptWithShortChallenges;
use crate::utils::sumcheck::{MLSumcheck, SumCheckError::SumCheckFailed};
use crate::{
    arith::{utils::mat_vec_mul, Instance, Witness, CCS, LCCCS},
    decomposition_parameters::DecompositionParams,
    utils::sumcheck,
};

use lattirust_poly::{
    mle::DenseMultilinearExtension,
    polynomials::{eq_eval, VPAuxInfo},
};
use utils::*;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

mod utils;

#[derive(Clone, Debug, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct FoldingProof<NTT: OverField> {
    // Step 2.
    pub pointshift_sumcheck_proof: sumcheck::Proof<NTT>,
    // Step 3
    pub theta_s: Vec<Vec<NTT>>,
    pub eta_s: Vec<Vec<NTT>>,
}

pub trait FoldingProver<NTT: SuitableRing, T: TranscriptWithShortChallenges<NTT>> {
    fn prove<const C: usize, P: DecompositionParams>(
        cm_i_s: &[LCCCS<C, NTT>],
        w_s: &[Witness<NTT>],
        transcript: &mut impl TranscriptWithShortChallenges<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(LCCCS<C, NTT>, Witness<NTT>, FoldingProof<NTT>), FoldingError<NTT>>;
}

pub trait FoldingVerifier<NTT: SuitableRing, T: TranscriptWithShortChallenges<NTT>> {
    fn verify<const C: usize, P: DecompositionParams>(
        cm_i_s: &[LCCCS<C, NTT>],
        proof: &FoldingProof<NTT>,
        transcript: &mut impl TranscriptWithShortChallenges<NTT>,
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

impl<NTT: SuitableRing, T: TranscriptWithShortChallenges<NTT>> FoldingProver<NTT, T>
    for LFFoldingProver<NTT, T>
{
    fn prove<const C: usize, P: DecompositionParams>(
        cm_i_s: &[LCCCS<C, NTT>],
        w_s: &[Witness<NTT>],
        transcript: &mut impl TranscriptWithShortChallenges<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(LCCCS<C, NTT>, Witness<NTT>, FoldingProof<NTT>), FoldingError<NTT>> {
        if cm_i_s.len() != 2 * P::K {
            return Err(FoldingError::IncorrectLength);
        }

        let log_m = ccs.s;

        // Step 1: Generate alpha, zeta, mu, beta challenges
        let (alpha_s, beta_s, zeta_s, mu_s) =
            get_alphas_betas_zetas_mus::<_, _, P>(log_m, transcript);

        // Step 2: Compute g polynomial and sumcheck on it
        // Setup f_hat_mle for later evaluation of thetas
        let f_hat_mles = w_s
            .iter()
            .map(|w| to_mles::<_, _, FoldingError<_>>(log_m, w.f_hat.iter())) // propagate errors using `?`
            .collect::<Result<Vec<Vec<DenseMultilinearExtension<NTT>>>, _>>()?;

        let zis = cm_i_s
            .iter()
            .zip(w_s.iter())
            .map(|(cm_i, w_i)| cm_i.get_z_vector(&w_i.w_ccs))
            .collect::<Vec<_>>();
        let ris = cm_i_s.iter().map(|cm_i| cm_i.r.clone()).collect::<Vec<_>>();

        let Mz_mles_vec: Vec<Vec<DenseMultilinearExtension<NTT>>> = zis
            .iter()
            .map(|zi| {
                let Mz_mle = to_mles_err::<_, _, FoldingError<NTT>, _>(
                    log_m,
                    ccs.M.iter().map(|M| mat_vec_mul(&M, &zi)),
                )?;
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

        // Step 5: Run sum check prover
        let (sum_check_proof, prover_state) = MLSumcheck::prove_as_subprotocol(transcript, &g);

        let r_0 = prover_state
            .randomness
            .into_iter()
            .map(|x| x.into())
            .collect::<Vec<NTT>>();

        // Step 3: Evaluate thetas and etas
        let theta_s: Vec<Vec<NTT>> = cfg_iter!(f_hat_mles)
            .map(|f_hat_row| {
                f_hat_row
                    .iter()
                    .map(|f_hat_mle| {
                        f_hat_mle
                            .evaluate(&r_0)
                            .ok_or(FoldingError::<NTT>::EvaluationError(
                                super::error::MleEvaluationError::IncorrectLength(
                                    r_0.len(),
                                    f_hat_mle.evaluations.len(),
                                ),
                            ))
                    })
                    .collect::<Result<Vec<_>, _>>()
            })
            .collect::<Result<Vec<_>, _>>()?;

        let eta_s: Vec<Vec<NTT>> = cfg_iter!(Mz_mles_vec)
            .map(|Mz_mles| {
                Mz_mles
                    .iter()
                    .map(|mle| {
                        mle.evaluate(r_0.as_slice())
                            .ok_or(FoldingError::<NTT>::EvaluationError(
                                super::error::MleEvaluationError::IncorrectLength(
                                    r_0.len(),
                                    mle.evaluations.len(),
                                ),
                            ))
                    })
                    .collect::<Result<Vec<_>, _>>()
            })
            .collect::<Result<Vec<_>, _>>()?;

        theta_s
            .iter()
            .for_each(|thetas| transcript.absorb_slice(thetas));
        eta_s.iter().for_each(|etas| transcript.absorb_slice(etas));

        // Step 5 get rho challenges
        let rho_s = get_rhos::<_, _, P>(transcript);

        // Step 6 compute v0, u0, y0, x_w0
        let (v_0, cm_0, u_0, x_0) = compute_v0_u0_x0_cm_0(&rho_s, &theta_s, cm_i_s, &eta_s, ccs);

        // Step 7: Compute f0 and Witness_0
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
                    let rho_i: NTT = rho_i.crt();

                    acc.into_iter()
                        .zip(w_i.f.iter())
                        .map(|(acc_j, w_ij)| acc_j + rho_i * w_ij)
                        .collect()
                });

        let w_0 = Witness::from_f::<P>(f_0);

        let folding_proof = FoldingProof {
            pointshift_sumcheck_proof: sum_check_proof,
            theta_s,
            eta_s,
        };

        Ok((lcccs, w_0, folding_proof))
    }
}

impl<NTT: SuitableRing, T: TranscriptWithShortChallenges<NTT>> FoldingVerifier<NTT, T>
    for LFFoldingVerifier<NTT, T>
{
    fn verify<const C: usize, P: DecompositionParams>(
        cm_i_s: &[LCCCS<C, NTT>],
        proof: &FoldingProof<NTT>,
        transcript: &mut impl TranscriptWithShortChallenges<NTT>,
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
        let vs = cm_i_s
            .iter()
            .map(|cm_i| cm_i.v.clone())
            .collect::<Vec<Vec<NTT>>>();
        let us = cm_i_s.iter().map(|cm_i| cm_i.u.clone()).collect::<Vec<_>>();

        let claim_g1: NTT = alpha_s
            .iter()
            .zip(vs.iter())
            .map(|(&alpha_i, v_i)| {
                successors(Some(alpha_i), |&alpha| Some(alpha * alpha_i))
                    .zip(v_i.iter())
                    .map(|(pow_of_alpha, v_i_j)| pow_of_alpha * v_i_j)
                    .sum::<NTT>()
            })
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
        proof
            .theta_s
            .iter()
            .for_each(|thetas| transcript.absorb_slice(thetas));
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

#[cfg(test)]
mod tests {
    use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Compress, Validate};
    use ark_std::io::Cursor;
    use rand::thread_rng;

    use crate::nifs::folding::FoldingProof;
    use crate::{
        arith::{r1cs::get_test_z_split, tests::get_test_ccs, Witness, CCCS},
        commitment::AjtaiCommitmentScheme,
        decomposition_parameters::{test_params::PPL1, DecompositionParams},
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
        transcript::poseidon::PoseidonTranscript,
    };
    use cyclotomic_rings::rings::{StarkChallengeSet, StarkRingNTT};

    // Boilerplate code to generate values needed for testing
    type R = StarkRingNTT;
    type CS = StarkChallengeSet;
    type T = PoseidonTranscript<StarkRingNTT, CS>;

    #[test]
    fn test_folding() {
        const WIT_LEN: usize = 4; // 4 is the length of witness in this (Vitalik's) example
        const W: usize = WIT_LEN * PPL1::L; // the number of columns of the Ajtai matrix

        let ccs = get_test_ccs::<R>(W);
        let (_, x_ccs, w_ccs) = get_test_z_split::<R>(3);
        let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());
        let wit: Witness<R> = Witness::from_w_ccs::<PPL1>(w_ccs);
        let cm_i: CCCS<4, R> = CCCS {
            cm: wit.commit::<4, 4, PPL1>(&scheme).unwrap(),
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

        let (_, vec_wit, decomposition_proof) = LFDecompositionProver::<_, T>::prove::<4, 4, PPL1>(
            &lcccs,
            &wit,
            &mut prover_transcript,
            &ccs,
            &scheme,
        )
        .unwrap();

        let vec_lcccs = LFDecompositionVerifier::<_, T>::verify::<4, PPL1>(
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
            LFFoldingProver::<_, T>::prove::<4, PPL1>(&lcccs, &wit_s, &mut prover_transcript, &ccs)
                .unwrap();

        let lcccs_verifier = LFFoldingVerifier::<_, T>::verify::<4, PPL1>(
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
        const W: usize = WIT_LEN * PPL1::L; // the number of columns of the Ajtai matrix

        let ccs = get_test_ccs::<R>(W);
        let (_, x_ccs, w_ccs) = get_test_z_split::<R>(3);
        let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());
        let wit: Witness<R> = Witness::from_w_ccs::<PPL1>(w_ccs.clone());
        let cm_i: CCCS<4, R> = CCCS {
            cm: wit.commit::<4, 4, PPL1>(&scheme).unwrap(),
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
            LFDecompositionProver::<_, T>::prove::<4, 4, PPL1>(
                &lcccs,
                &wit,
                &mut prover_transcript,
                &ccs,
                &scheme,
            )
            .unwrap();

        let vec_lcccs = LFDecompositionVerifier::<_, T>::verify::<4, PPL1>(
            &lcccs,
            &decomposition_proof,
            &mut verifier_transcript,
            &ccs,
        )
        .unwrap();

        vec_wit[0] = Witness::<R>::from_w_ccs::<PPL1>(w_ccs);

        let res = LFFoldingProver::<_, T>::prove::<4, PPL1>(
            &vec_lcccs,
            &vec_wit,
            &mut prover_transcript,
            &ccs,
        );

        assert!(res.is_err())
    }

    #[test]
    fn test_folding_proof_serialization() {
        const WIT_LEN: usize = 4; // 4 is the length of witness in this (Vitalik's) example
        const W: usize = WIT_LEN * PPL1::L; // the number of columns of the Ajtai matrix

        let ccs = get_test_ccs::<R>(W);
        let (_, x_ccs, w_ccs) = get_test_z_split::<R>(3);
        let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());
        let wit: Witness<R> = Witness::from_w_ccs::<PPL1>(w_ccs);
        let cm_i: CCCS<4, R> = CCCS {
            cm: wit.commit::<4, 4, PPL1>(&scheme).unwrap(),
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

        let (_, vec_wit, decomposition_proof) = LFDecompositionProver::<_, T>::prove::<4, 4, PPL1>(
            &lcccs,
            &wit,
            &mut prover_transcript,
            &ccs,
            &scheme,
        )
        .unwrap();

        let vec_lcccs = LFDecompositionVerifier::<_, T>::verify::<4, PPL1>(
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
        let (_, _, folding_proof) =
            LFFoldingProver::<_, T>::prove::<4, PPL1>(&lcccs, &wit_s, &mut prover_transcript, &ccs)
                .unwrap();

        let mut serialized = Vec::new();
        folding_proof
            .serialize_with_mode(&mut serialized, Compress::Yes)
            .expect("Failed to serialize proof");

        let mut cursor = Cursor::new(&serialized);
        assert_eq!(
            folding_proof,
            FoldingProof::deserialize_with_mode(&mut cursor, Compress::Yes, Validate::Yes)
                .expect("Failed to deserialize proof")
        );
    }
}

#[cfg(test)]
mod tests_goldilocks {
    use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Compress, Validate};
    use ark_std::io::Cursor;
    use rand::thread_rng;

    use crate::nifs::folding::FoldingProof;
    use crate::{
        arith::{r1cs::get_test_z_split, tests::get_test_ccs, Witness, CCCS},
        commitment::AjtaiCommitmentScheme,
        decomposition_parameters::{test_params::PPL1, DecompositionParams},
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
        transcript::poseidon::PoseidonTranscript,
    };
    use cyclotomic_rings::rings::{GoldilocksChallengeSet, GoldilocksRingNTT};

    // Boilerplate code to generate values needed for testing
    type R = GoldilocksRingNTT;
    type CS = GoldilocksChallengeSet;
    type T = PoseidonTranscript<GoldilocksRingNTT, CS>;

    #[test]
    fn test_folding() {
        const WIT_LEN: usize = 4; // 4 is the length of witness in this (Vitalik's) example
        const W: usize = WIT_LEN * PPL1::L; // the number of columns of the Ajtai matrix

        let ccs = get_test_ccs::<R>(W);
        let (_, x_ccs, w_ccs) = get_test_z_split::<R>(3);
        let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());
        let wit: Witness<R> = Witness::from_w_ccs::<PPL1>(w_ccs);
        let cm_i: CCCS<4, R> = CCCS {
            cm: wit.commit::<4, 4, PPL1>(&scheme).unwrap(),
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

        let (_, vec_wit, decomposition_proof) = LFDecompositionProver::<_, T>::prove::<4, 4, PPL1>(
            &lcccs,
            &wit,
            &mut prover_transcript,
            &ccs,
            &scheme,
        )
        .unwrap();

        let vec_lcccs = LFDecompositionVerifier::<_, T>::verify::<4, PPL1>(
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
            LFFoldingProver::<_, T>::prove::<4, PPL1>(&lcccs, &wit_s, &mut prover_transcript, &ccs)
                .unwrap();

        let lcccs_verifier = LFFoldingVerifier::<_, T>::verify::<4, PPL1>(
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
        const W: usize = WIT_LEN * PPL1::L; // the number of columns of the Ajtai matrix

        let ccs = get_test_ccs::<R>(W);
        let (_, x_ccs, w_ccs) = get_test_z_split::<R>(3);
        let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());
        let wit: Witness<R> = Witness::from_w_ccs::<PPL1>(w_ccs.clone());
        let cm_i: CCCS<4, R> = CCCS {
            cm: wit.commit::<4, 4, PPL1>(&scheme).unwrap(),
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
            LFDecompositionProver::<_, T>::prove::<4, 4, PPL1>(
                &lcccs,
                &wit,
                &mut prover_transcript,
                &ccs,
                &scheme,
            )
            .unwrap();

        let vec_lcccs = LFDecompositionVerifier::<_, T>::verify::<4, PPL1>(
            &lcccs,
            &decomposition_proof,
            &mut verifier_transcript,
            &ccs,
        )
        .unwrap();

        vec_wit[0] = Witness::<R>::from_w_ccs::<PPL1>(w_ccs);

        let res = LFFoldingProver::<_, T>::prove::<4, PPL1>(
            &vec_lcccs,
            &vec_wit,
            &mut prover_transcript,
            &ccs,
        );

        assert!(res.is_err())
    }

    #[test]
    fn test_folding_proof_serialization() {
        const WIT_LEN: usize = 4; // 4 is the length of witness in this (Vitalik's) example
        const W: usize = WIT_LEN * PPL1::L; // the number of columns of the Ajtai matrix

        let ccs = get_test_ccs::<R>(W);
        let (_, x_ccs, w_ccs) = get_test_z_split::<R>(3);
        let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());
        let wit: Witness<R> = Witness::from_w_ccs::<PPL1>(w_ccs);
        let cm_i: CCCS<4, R> = CCCS {
            cm: wit.commit::<4, 4, PPL1>(&scheme).unwrap(),
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

        let (_, vec_wit, decomposition_proof) = LFDecompositionProver::<_, T>::prove::<4, 4, PPL1>(
            &lcccs,
            &wit,
            &mut prover_transcript,
            &ccs,
            &scheme,
        )
        .unwrap();

        let vec_lcccs = LFDecompositionVerifier::<_, T>::verify::<4, PPL1>(
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
        let (_, _, folding_proof) =
            LFFoldingProver::<_, T>::prove::<4, PPL1>(&lcccs, &wit_s, &mut prover_transcript, &ccs)
                .unwrap();

        let mut serialized = Vec::new();
        folding_proof
            .serialize_with_mode(&mut serialized, Compress::Yes)
            .expect("Failed to serialize proof");

        let mut cursor = Cursor::new(&serialized);
        assert_eq!(
            folding_proof,
            FoldingProof::deserialize_with_mode(&mut cursor, Compress::Yes, Validate::Yes)
                .expect("Failed to deserialize proof")
        );
    }
}

#[cfg(test)]
mod tests_stark {
    use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Compress, Validate};
    use ark_std::io::Cursor;
    use lattirust_ring::cyclotomic_ring::models::stark_prime::RqNTT;
    use num_bigint::BigUint;
    use rand::thread_rng;

    use crate::nifs::folding::FoldingProof;
    use crate::{
        arith::r1cs::get_test_dummy_z_split,
        nifs::{
            decomposition::{
                DecompositionProver, DecompositionVerifier, LFDecompositionProver,
                LFDecompositionVerifier,
            },
            folding::{FoldingProver, FoldingVerifier, LFFoldingProver, LFFoldingVerifier},
            linearization::LinearizationProver,
        },
    };
    use crate::{
        arith::tests::get_test_dummy_ccs,
        utils::security_check::check_ring_modulus_128_bits_security,
    };
    use crate::{
        arith::{Witness, CCCS},
        commitment::AjtaiCommitmentScheme,
        decomposition_parameters::{test_params::PP_STARK_FOLDING, DecompositionParams},
        nifs::linearization::{
            LFLinearizationProver, LFLinearizationVerifier, LinearizationVerifier,
        },
        transcript::poseidon::PoseidonTranscript,
    };
    use cyclotomic_rings::rings::StarkChallengeSet;

    #[test]
    fn test_dummy_folding() {
        #[cfg(feature = "dhat-heap")]
        #[global_allocator]
        static ALLOC: dhat::Alloc = dhat::Alloc;

        type R = RqNTT;
        type CS = StarkChallengeSet;
        type T = PoseidonTranscript<R, CS>;

        const C: usize = 15;
        const X_LEN: usize = 1;
        const WIT_LEN: usize = 512;
        const W: usize = WIT_LEN * PP_STARK_FOLDING::L; // the number of columns of the Ajtai matrix
        let r1cs_rows_size = X_LEN + WIT_LEN + 1; // Let's have a square matrix

        #[cfg(feature = "dhat-heap")]
        let _profiler = dhat::Profiler::new_heap(); // Move a round to measure specific parts

        let ccs = get_test_dummy_ccs::<R, X_LEN, WIT_LEN, W>(r1cs_rows_size);
        let (_, x_ccs, w_ccs) = get_test_dummy_z_split::<R, X_LEN, WIT_LEN>();
        let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());

        let wit = Witness::from_w_ccs::<PP_STARK_FOLDING>(w_ccs);

        // Make bound and securitty checks
        let witness_within_bound = wit.within_bound(PP_STARK_FOLDING::B);
        let stark_modulus = BigUint::parse_bytes(
            b"3618502788666131000275863779947924135206266826270938552493006944358698582017",
            10,
        )
        .expect("Failed to parse stark_modulus");

        if check_ring_modulus_128_bits_security(
            &stark_modulus,
            C,
            16,
            W,
            PP_STARK_FOLDING::B,
            PP_STARK_FOLDING::L,
            witness_within_bound,
        ) {
            println!(" Bound condition satisfied for 128 bits security");
        } else {
            println!("Bound condition not satisfied for 128 bits security");
        }

        let cm_i = CCCS {
            cm: wit.commit::<C, W, PP_STARK_FOLDING>(&scheme).unwrap(),
            x_ccs,
        };

        let mut prover_transcript = PoseidonTranscript::<R, CS>::default();

        let linearization_proof =
            LFLinearizationProver::<_, T>::prove(&cm_i, &wit, &mut prover_transcript, &ccs);

        let mut verifier_transcript = PoseidonTranscript::<R, CS>::default();

        let linearization_verification = LFLinearizationVerifier::<_, T>::verify(
            &cm_i,
            &linearization_proof
                .expect("Linearization proof generation error")
                .1,
            &mut verifier_transcript,
            &ccs,
        )
        .expect("Linearization Verification error");

        let lcccs = linearization_verification;

        let decomposition_prover = LFDecompositionProver::<_, T>::prove::<W, C, PP_STARK_FOLDING>(
            &lcccs,
            &wit,
            &mut prover_transcript,
            &ccs,
            &scheme,
        );

        let decomposition_proof =
            decomposition_prover.expect("Decomposition proof generation error");

        let decomposition_verification =
            LFDecompositionVerifier::<_, T>::verify::<C, PP_STARK_FOLDING>(
                &lcccs,
                &decomposition_proof.2,
                &mut verifier_transcript,
                &ccs,
            );

        let lcccs = decomposition_verification.expect("Decomposition Verification error");

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
        let folding_prover = LFFoldingProver::<_, T>::prove::<C, PP_STARK_FOLDING>(
            &lcccs,
            &wit_s,
            &mut prover_transcript,
            &ccs,
        );

        let folding_proof = folding_prover.expect("Folding proof generation error");

        let folding_verification = LFFoldingVerifier::<_, T>::verify::<C, PP_STARK_FOLDING>(
            &lcccs,
            &folding_proof.2,
            &mut verifier_transcript,
            &ccs,
        );

        folding_verification.expect("Folding Verification error");
    }

    #[test]
    fn test_folding_proof_serialization() {
        type R = RqNTT;
        type CS = StarkChallengeSet;
        type T = PoseidonTranscript<R, CS>;

        const C: usize = 15;
        const X_LEN: usize = 1;
        const WIT_LEN: usize = 512;
        const W: usize = WIT_LEN * PP_STARK_FOLDING::L; // the number of columns of the Ajtai matrix
        let r1cs_rows_size = X_LEN + WIT_LEN + 1; // Let's have a square matrix

        let ccs = get_test_dummy_ccs::<R, X_LEN, WIT_LEN, W>(r1cs_rows_size);
        let (_, x_ccs, w_ccs) = get_test_dummy_z_split::<R, X_LEN, WIT_LEN>();
        let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());

        let wit = Witness::from_w_ccs::<PP_STARK_FOLDING>(w_ccs);

        // Make bound and security checks
        let witness_within_bound = wit.within_bound(PP_STARK_FOLDING::B);
        let stark_modulus = BigUint::parse_bytes(
            b"3618502788666131000275863779947924135206266826270938552493006944358698582017",
            10,
        )
        .expect("Failed to parse stark_modulus");

        if check_ring_modulus_128_bits_security(
            &stark_modulus,
            C,
            16,
            W,
            PP_STARK_FOLDING::B,
            PP_STARK_FOLDING::L,
            witness_within_bound,
        ) {
            println!(" Bound condition satisfied for 128 bits security");
        } else {
            println!("Bound condition not satisfied for 128 bits security");
        }

        let cm_i = CCCS {
            cm: wit.commit::<C, W, PP_STARK_FOLDING>(&scheme).unwrap(),
            x_ccs,
        };

        let mut prover_transcript = PoseidonTranscript::<R, CS>::default();

        let linearization_proof =
            LFLinearizationProver::<_, T>::prove(&cm_i, &wit, &mut prover_transcript, &ccs);

        let mut verifier_transcript = PoseidonTranscript::<R, CS>::default();

        let linearization_verification = LFLinearizationVerifier::<_, T>::verify(
            &cm_i,
            &linearization_proof
                .expect("Linearization proof generation error")
                .1,
            &mut verifier_transcript,
            &ccs,
        )
        .expect("Linearization Verification error");

        let lcccs = linearization_verification;

        let decomposition_prover = LFDecompositionProver::<_, T>::prove::<W, C, PP_STARK_FOLDING>(
            &lcccs,
            &wit,
            &mut prover_transcript,
            &ccs,
            &scheme,
        );

        let decomposition_proof =
            decomposition_prover.expect("Decomposition proof generation error");

        let decomposition_verification =
            LFDecompositionVerifier::<_, T>::verify::<C, PP_STARK_FOLDING>(
                &lcccs,
                &decomposition_proof.2,
                &mut verifier_transcript,
                &ccs,
            );

        let lcccs = decomposition_verification.expect("Decomposition Verification error");

        let (lcccs, wit_s) = {
            let mut lcccs = lcccs.clone();
            let mut lcccs_r = lcccs.clone();
            lcccs.append(&mut lcccs_r);

            let mut wit_s = decomposition_proof.1.clone();
            let mut wit_s_r = decomposition_proof.1;
            wit_s.append(&mut wit_s_r);

            (lcccs, wit_s)
        };
        let folding_prover = LFFoldingProver::<_, T>::prove::<C, PP_STARK_FOLDING>(
            &lcccs,
            &wit_s,
            &mut prover_transcript,
            &ccs,
        );

        let folding_proof = folding_prover.expect("Folding proof generation error").2;

        let mut serialized = Vec::new();
        folding_proof
            .serialize_with_mode(&mut serialized, Compress::Yes)
            .expect("Failed to serialize proof");

        let mut cursor = Cursor::new(&serialized);
        assert_eq!(
            folding_proof,
            FoldingProof::deserialize_with_mode(&mut cursor, Compress::Yes, Validate::Yes)
                .expect("Failed to deserialize proof")
        );
    }
}
