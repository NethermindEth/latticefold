use std::marker::PhantomData;
use std::sync::Arc;

use crate::utils::sumcheck::MLSumcheck;
use crate::utils::sumcheck::SumCheckError::SumCheckFailed;
use crate::{
    arith::{utils::mat_vec_mul, Instance, Witness, CCS, LCCCS},
    parameters::DecompositionParams,
    transcript::Transcript,
    utils::{mle::dense_vec_to_dense_mle, sumcheck},
};

use lattirust_arithmetic::{
    challenge_set::latticefold_challenge_set::OverField,
    mle::DenseMultilinearExtension,
    polynomials::{build_eq_x_r, eq_eval, VPAuxInfo, VirtualPolynomial},
    ring::PolyRing,
};

use super::error::FoldingError;

use ark_std::iterable::Iterable;
use ark_std::log2;

#[derive(Clone)]
pub struct FoldingProof<NTT: OverField> {
    // Step 2.
    pub pointshift_sumcheck_proof: sumcheck::Proof<NTT>,
    // Step 3
    pub theta_s: Vec<NTT>,
    pub eta_s: Vec<Vec<NTT>>,
}

pub trait FoldingProver<NTT: OverField, T: Transcript<NTT>> {
    fn prove<const C: usize, CR: PolyRing, P: DecompositionParams>(
        cm_i_s: &[LCCCS<C, NTT>],
        w_s: &[Witness<NTT>],
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(LCCCS<C, NTT>, Witness<NTT>, FoldingProof<NTT>), FoldingError<NTT>>;
}

pub trait FoldingVerifier<NTT: OverField, T: Transcript<NTT>> {
    fn verify<const C: usize, P: DecompositionParams>(
        cm_i_s: &[LCCCS<C, NTT>],
        proof: &FoldingProof<NTT>,
        transcript: &mut impl Transcript<NTT>,
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

impl<NTT: OverField, T: Transcript<NTT>> FoldingProver<NTT, T> for LFFoldingProver<NTT, T> {
    fn prove<const C: usize, CR: PolyRing, P: DecompositionParams>(
        cm_i_s: &[LCCCS<C, NTT>],
        w_s: &[Witness<NTT>],
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(LCCCS<C, NTT>, Witness<NTT>, FoldingProof<NTT>), FoldingError<NTT>> {
        assert_eq!(cm_i_s.len(), 2 * P::K);
        let m = ccs.m;
        let log_m = log2(m) as usize;

        // Step 1: Generate alpha, zeta, mu, beta challenges
        cm_i_s.iter().for_each(|lcccs| {
            transcript.absorb_ring_vec(&lcccs.r);
            transcript.absorb_ring(&lcccs.v);
            // _transcript.absorb_ring_vec(&lcccs.cm); Not absorbed by transcript?
            transcript.absorb_ring_vec(&lcccs.u);
            transcript.absorb_ring_vec(&lcccs.x_w);
        });
        // TODO: Get challenges from big set but as NTT
        let alpha_s = transcript.get_small_challenges(2 * P::K);
        let zeta_s = transcript.get_small_challenges(2 * P::K);
        let mu_s = transcript.get_small_challenges((2 * P::K) - 1); // Note is one challenge less
        let beta_s = transcript.get_small_challenges(log_m);

        // Step 2: Compute g polynomial and sumcheck on it
        // Setup f_hat_mle for later evaluation of thetas
        let f_hat_mles = w_s
            .iter()
            .map(|w| {
                let f_i = w.f.clone();
                DenseMultilinearExtension::from_evaluations_vec(log_m, f_i)
            })
            .collect::<Vec<_>>();

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

        let g = create_sumcheck_polynomial::<NTT, P>(
            log_m,
            &f_hat_mles,
            &alpha_s,
            &Mz_mles_vec,
            &zeta_s,
            ris,
            &beta_s,
            &mu_s,
        );

        let claim_g1 = alpha_s
            .iter()
            .zip(vs.iter())
            .fold(NTT::zero(), |acc, (&alpha, &vi)| acc + (alpha * vi));
        let claim_g2 = zeta_s
            .iter()
            .zip(us.iter())
            .fold(NTT::zero(), |acc, (zeta, ui)| {
                let mut zeta_i = NTT::one();
                let ui_sum = ui.iter().fold(NTT::zero(), |acc, &u_i_t| {
                    zeta_i = zeta_i * zeta;
                    acc + (u_i_t * zeta_i)
                });
                acc + ui_sum
            });
        let (sum_check_proof, prover_state) = MLSumcheck::prove_as_subprotocol(transcript, &g);
        let r_0 = prover_state
            .randomness
            .into_iter()
            .map(|x| NTT::field_to_base_ring(&x).into())
            .collect::<Vec<NTT>>();

        // Step 3: Evaluate thetas and etas
        let thetas = f_hat_mles
            .iter()
            .map(|f_hat_mle| f_hat_mle.evaluate(r_0.as_slice()).unwrap())
            .collect::<Vec<_>>();
        drop(f_hat_mles);
        let etas: Vec<Vec<NTT>> = Mz_mles_vec
            .iter()
            .map(|Mz_mles| {
                Mz_mles
                    .iter()
                    .map(|mle| mle.evaluate(r_0.as_slice()).unwrap())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        drop(Mz_mles_vec);

        transcript.absorb_ring_vec(&thetas);
        etas.iter()
            .for_each(|etas| transcript.absorb_ring_vec(&etas));

        // Step 5 get rho challenges
        let rhos = transcript.get_small_challenges((2 * P::K) - 1); // Note that we are missing the first element

        // Step 6 compute v0, u0, y0, x_w0
        let v_0: NTT =
            rhos.iter()
                .zip(thetas.iter().skip(1))
                .fold(thetas[0], |acc, (rho_i, theta_i)| {
                    // acc + rho_i.rot_sum(theta_i) // Note that theta_i is already in NTT form
                    todo!() // Add WithRot to OverField in lattirust
                }); // Do INTT here

        let (y_0, u_0, x_0) = rhos
            .iter()
            .zip(cm_i_s.iter().zip(etas.iter()).skip(1))
            .fold(
                (cm_i_s[0].cm.clone(), etas[0].clone(), cm_i_s[0].x_w.clone()),
                |(acc_y, acc_u, acc_x), (rho_i, (cm_i, eta))| {
                    let y = acc_y + (cm_i.cm.clone() * rho_i);
                    let u = acc_u
                        .iter()
                        .zip(eta.iter())
                        .map(|(&a, &e)| a + (e * rho_i))
                        .collect();
                    let x = acc_x
                        .iter()
                        .zip(cm_i.x_w.iter())
                        .map(|(&a, &x)| a + (x * rho_i))
                        .collect();
                    (y, u, x)
                },
            );

        // Step 7: Compute f0 and Witness_0

        let h = x_0.last().cloned().ok_or(FoldingError::IncorrectLength)?;
        let lcccs = LCCCS {
            r: r_0,
            v: v_0,
            cm: y_0,
            u: u_0,
            x_w: x_0,
            h,
        };

        let f_0 =
            rhos.iter()
                .zip(w_s.iter().skip(1))
                .fold(w_s[0].f.clone(), |acc, (rho, w_i_s)| {
                    let mut f_i = w_i_s.f.clone();
                    f_i.iter_mut().for_each(|c| *c = *c * rho);
                    acc.iter().zip(f_i.iter()).map(|(a, f)| *a + f).collect()
                });

        let folding_proof = FoldingProof {
            pointshift_sumcheck_proof: sum_check_proof,
            theta_s: thetas,
            eta_s: etas,
        };

        let wit_0 = Witness::<NTT>::from_f::<NTT, P>(f_0);
        Ok((lcccs, wit_0, folding_proof))
    }
}

impl<NTT: OverField, T: Transcript<NTT>> FoldingVerifier<NTT, T> for LFFoldingVerifier<NTT, T> {
    fn verify<const C: usize, P: DecompositionParams>(
        cm_i_s: &[LCCCS<C, NTT>],
        proof: &FoldingProof<NTT>,
        _transcript: &mut impl Transcript<NTT>,
        _ccs: &CCS<NTT>,
    ) -> Result<LCCCS<C, NTT>, FoldingError<NTT>> {
        let m = _ccs.m;
        let log_m = log2(m) as usize;

        // Step 1: Generate alpha, zeta, mu, beta challenges
        cm_i_s.iter().for_each(|lcccs| {
            _transcript.absorb_ring_vec(&lcccs.r);
            _transcript.absorb_ring(&lcccs.v);
            // _transcript.absorb_ring_vec(&lcccs.cm); Not absorbed by transcript?
            _transcript.absorb_ring_vec(&lcccs.u);
            _transcript.absorb_ring_vec(&lcccs.x_w);
        });
        // TODO: Get challenges from big set but as NTT
        let alpha_s = _transcript.get_small_challenges(2 * P::K);
        let zeta_s = _transcript.get_small_challenges(2 * P::K);
        let mu_s = _transcript.get_small_challenges((2 * P::K) - 1); // Note is one challenge less
        let beta_s = _transcript.get_small_challenges(log_m);

        let poly_info = VPAuxInfo {
            max_degree: _ccs.d + 1,
            num_variables: log_m,
            phantom: std::marker::PhantomData,
        };
        let ris = cm_i_s.iter().map(|cm_i| cm_i.r.clone()).collect::<Vec<_>>();
        let vs = cm_i_s.iter().map(|cm_i| cm_i.v).collect::<Vec<NTT>>();
        let us = cm_i_s.iter().map(|cm_i| cm_i.u.clone()).collect::<Vec<_>>();

        let claim_g1 = alpha_s
            .iter()
            .zip(vs.iter())
            .fold(NTT::zero(), |acc, (&alpha, &vi)| acc + (alpha * vi));
        let claim_g2 = zeta_s
            .iter()
            .zip(us.iter())
            .fold(NTT::zero(), |acc, (zeta, ui)| {
                let mut zeta_i = NTT::one();
                let ui_sum = ui.iter().fold(NTT::zero(), |acc, &u_i_t| {
                    zeta_i = zeta_i * zeta;
                    acc + (u_i_t * zeta_i)
                });
                acc + ui_sum
            });

        //Step 2: The sumcheck.

        // Verify the sumcheck proof.
        let sub_claim = MLSumcheck::verify_as_subprotocol(
            _transcript,
            &poly_info,
            claim_g1 + claim_g2,
            &proof.pointshift_sumcheck_proof,
        )?;

        let point_r = sub_claim
            .point
            .into_iter()
            .map(|x| NTT::field_to_base_ring(&x).into())
            .collect::<Vec<NTT>>();

        let e_asterisk = eq_eval(&beta_s, &point_r).unwrap();
        let e_i_s: Vec<NTT> = ris
            .iter()
            .map(|r| eq_eval(r.as_slice(), &point_r).unwrap())
            .collect::<Vec<_>>();
        let s = sub_claim.expected_evaluation.clone();

        let mut should_equal_s =
            mu_s.iter()
                .zip(proof.theta_s.iter())
                .fold(NTT::zero(), |acc, (mu_i, &theta_i)| {
                    let mut thetas_mul = theta_i;
                    for j in 1..P::B_SMALL {
                        thetas_mul = thetas_mul * (theta_i - NTT::from(j));
                        thetas_mul = thetas_mul * (theta_i + NTT::from(j));
                    }
                    acc + (thetas_mul * mu_i)
                });
        let last_theta = proof
            .theta_s
            .last()
            .cloned()
            .ok_or(FoldingError::IncorrectLength)?;
        // Recall last mu = 1
        let mut theta_mul = last_theta;
        for j in 1..P::B_SMALL {
            theta_mul = theta_mul * (last_theta - NTT::from(j));
            theta_mul = theta_mul * (last_theta + NTT::from(j));
        }
        should_equal_s = should_equal_s + theta_mul;
        should_equal_s = should_equal_s * e_asterisk;

        alpha_s
            .iter()
            .zip(zeta_s)
            .zip(proof.theta_s.iter())
            .zip(proof.eta_s.iter())
            .zip(e_i_s.iter())
            .for_each(|((((&alpha_i, zeta_i), theta_i), eta_i), e_i)| {
                let mut zeta_i_t = NTT::one();
                let combined_eta_i_t = eta_i.iter().fold(NTT::zero(), |acc, &eta_t| {
                    zeta_i_t = zeta_i_t * zeta_i;
                    acc + (eta_t * zeta_i_t)
                });
                let combined_theta_and_eta = (alpha_i * theta_i) + combined_eta_i_t;
                should_equal_s += combined_theta_and_eta * e_i;
            });

        match should_equal_s == s {
            true => {}
            false => {
                return Err(FoldingError::SumCheckError(SumCheckFailed(
                    should_equal_s,
                    s,
                )));
            }
        }
        // check claim and output o, u0, x0, y0
        let rhos = _transcript.get_small_challenges((2 * P::K) - 1); // Note that we are missing the first element

        // Step 6 compute v0, u0, y0, x_w0
        let v_0: NTT = rhos.iter().zip(proof.theta_s.iter().skip(1)).fold(
            proof.theta_s[0],
            |acc, (rho_i, theta_i)| {
                // acc + rho_i.rot_sum(theta_i) // Note that theta_i is already in NTT form
                todo!() // Add WithRot to OverField in lattirust
            },
        ); // Do INTT here

        let (y_0, u_0, x_0) = rhos
            .iter()
            .zip(cm_i_s.iter().zip(proof.eta_s.iter()).skip(1))
            .fold(
                (
                    cm_i_s[0].cm.clone(),
                    proof.eta_s[0].clone(),
                    cm_i_s[0].x_w.clone(),
                ),
                |(acc_y, acc_u, acc_x), (rho_i, (cm_i, eta))| {
                    let y = acc_y + (cm_i.cm.clone() * rho_i);
                    let u = acc_u
                        .iter()
                        .zip(eta.iter())
                        .map(|(&a, &e)| a + (e * rho_i))
                        .collect();
                    let x = acc_x
                        .iter()
                        .zip(cm_i.x_w.iter())
                        .map(|(&a, &x)| a + (x * rho_i))
                        .collect();
                    (y, u, x)
                },
            );

        let h = x_0.last().cloned().ok_or(FoldingError::IncorrectLength)?;
        Ok(LCCCS {
            r: point_r,
            v: v_0,
            cm: y_0,
            u: u_0,
            x_w: x_0,
            h,
        })
    }
}
fn create_matrix_mle<NTT: OverField>(
    log_m: usize,
    Mi: &Vec<Vec<NTT>>,
    zi: &Vec<NTT>,
) -> DenseMultilinearExtension<NTT> {
    let zero_vector = usize_to_binary_vector::<NTT>(0, log2(Mi.len()) as usize);
    let mle_z_ccs_b = mle_val_from_vector(&zi, &zero_vector);
    let evaluations: Vec<NTT> = mle_matrix_to_val_eval_second(&Mi, &zero_vector)
        .iter()
        .map(|val| *val * mle_z_ccs_b)
        .collect();
    let mle = DenseMultilinearExtension::from_evaluations_vec(log_m, evaluations);

    let matrix_mle = (1..Mi.len())
        .into_iter()
        .map(|i| usize_to_binary_vector::<NTT>(i, log2(Mi.len()) as usize))
        .fold(mle, |acc, b| {
            let mle_z_ccs_b = mle_val_from_vector(&zi, &b);
            let evaluations: Vec<NTT> = mle_matrix_to_val_eval_second(&Mi, &b)
                .iter()
                .map(|val| *val * mle_z_ccs_b)
                .collect();
            let mle = DenseMultilinearExtension::from_evaluations_vec(log_m, evaluations);
            acc + mle
        });
    matrix_mle
}

fn create_sumcheck_polynomial<NTT: OverField, DP: DecompositionParams>(
    log_m: usize,
    f_hat_mles: &Vec<DenseMultilinearExtension<NTT>>,
    alpha_is: &Vec<NTT>,
    matrix_mles: &Vec<Vec<DenseMultilinearExtension<NTT>>>,
    zeta_is: &Vec<NTT>,
    ris: Vec<Vec<NTT>>,
    Beta: &Vec<NTT>,
    mus: &Vec<NTT>,
) -> VirtualPolynomial<NTT> {
    let k = DP::K;
    let mut g = VirtualPolynomial::new(2 * k);
    let mut g1_plus_g3 = VirtualPolynomial::new(2 * k);
    for i in 0..2 * k {
        let gi_1 = prepare_g1_i_mle(&f_hat_mles[i], alpha_is[i]);
        let gi_3 = prepare_g3_i_mle(&matrix_mles[i], zeta_is[i]);
        let gi_1_plus_gi_3 = Arc::from(gi_1 + gi_3);
        let mut g1_and_g3_virtual = VirtualPolynomial::new_from_mle(&gi_1_plus_gi_3, NTT::ONE);
        let eq_r_i = build_eq_x_r(ris[i].as_slice()).unwrap();
        g1_and_g3_virtual.mul_by_mle(eq_r_i, NTT::one());
        g1_plus_g3 = &g1_plus_g3 + &g1_and_g3_virtual;
    }

    let b = DP::B_SMALL; // Get this from the decomposition step
    let mut g2 = prepare_g2_i_mle(log_m, &f_hat_mles[0], b, mus[0]);
    for i in 1..2 * k - 1 {
        let gi_2 = prepare_g2_i_mle(log_m, &f_hat_mles[i], b, mus[i]);
        g2 = &g2 + &gi_2;
    }
    // Recall that the last mu is 1
    let gi_2 = prepare_g2_i_mle(log_m, &f_hat_mles[2 * k - 1], b, NTT::one());
    g2 = &g2 + &gi_2;

    let eq_beta = build_eq_x_r::<NTT>(Beta.as_slice()).unwrap();
    g2.mul_by_mle(eq_beta, NTT::one());
    g = &g1_plus_g3 + &g2;

    g
}

fn prepare_g1_i_mle<NTT: OverField>(
    fi_mle: &DenseMultilinearExtension<NTT>,
    alpha_i: NTT,
) -> DenseMultilinearExtension<NTT> {
    let mut mle = fi_mle.clone(); // remove clone
    mle.evaluations.iter_mut().for_each(|e| *e = *e * alpha_i);
    mle
}

fn prepare_g2_i_mle<NTT: OverField>(
    log_m: usize,
    fi_mle: &DenseMultilinearExtension<NTT>,
    b: u128,
    mu_i: NTT,
) -> VirtualPolynomial<NTT> {
    let mut mle_list: Vec<Arc<DenseMultilinearExtension<NTT>>> = Vec::new();
    let mle_zero = fi_mle.clone();
    mle_list.push(Arc::from(mle_zero));
    for i in 0..b {
        let mut mle_j = fi_mle.clone();
        mle_j
            .evaluations
            .iter_mut()
            .for_each(|e| *e = *e - NTT::from(i)); // There should be a better way than sub every
                                                   // eval
        mle_list.push(Arc::from(mle_j));
    }
    for i in 0..b {
        let mut mle_j = fi_mle.clone();
        mle_j
            .evaluations
            .iter_mut()
            .for_each(|e| *e = *e + NTT::from(i));
        mle_list.push(Arc::from(mle_j));
    }
    let mut gi_2 = VirtualPolynomial::new(log_m);
    gi_2.add_mle_list(mle_list, mu_i);
    gi_2
}

fn prepare_g3_i_mle<NTT: OverField>(
    matrix_mle: &Vec<DenseMultilinearExtension<NTT>>,
    zeta_i: NTT,
) -> DenseMultilinearExtension<NTT> {
    let (first_mle, mles) = matrix_mle.split_first().unwrap();
    let first_mle = first_mle.clone();
    let mut mle = mles
        .into_iter()
        .fold(first_mle, |acc, mle_i_t| acc + mle_i_t.clone())
        .to_owned();
    let mut zeta = NTT::one();
    mle.evaluations.iter_mut().for_each(|e| {
        zeta = zeta * zeta_i;
        *e = *e * zeta_i
    });
    mle
}

fn mle_val_from_vector<NTT: OverField>(vector: &Vec<NTT>, values: &Vec<NTT>) -> NTT {
    assert_eq!(values.len(), log2(vector.len()) as usize);
    let mle = DenseMultilinearExtension::from_evaluations_vec(values.len(), vector.clone());
    mle.evaluate(values.as_slice()).unwrap()
}

// Convert a bivariate MLE to a univariate MLE by evaluating the second vector
fn mle_matrix_to_val_eval_second<NTT: OverField>(
    matrix: &Vec<Vec<NTT>>,
    values_y: &Vec<NTT>,
) -> Vec<NTT> {
    assert_eq!(values_y.len(), log2(matrix.len()) as usize);
    (0..matrix[0].len())
        .into_iter()
        .map(|i| mle_val_from_vector(&matrix.iter().map(|col| col[i]).collect(), values_y))
        .collect()
}

fn usize_to_binary_vector<NTT: OverField>(n: usize, length: usize) -> Vec<NTT> {
    let mut bits = Vec::new();
    let mut current = n;

    // Extract bits from the number
    while current > 0 {
        bits.push((current & 1) as u8);
        current >>= 1;
    }

    // Reverse to get the bits in correct order
    bits.reverse();

    // Pad with leading zeros if necessary
    if bits.len() < length {
        let padding = length - bits.len();
        bits.splice(0..0, std::iter::repeat(0).take(padding));
    }

    // Convert to the target field elements
    bits.into_iter()
        .map(|bit| if bit == 1 { NTT::one() } else { NTT::zero() })
        .collect()
}

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
