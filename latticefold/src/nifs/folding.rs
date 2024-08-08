use std::sync::Arc;

use crate::arith::utils::mat_vec_mul;
use crate::arith::Instance;
use crate::commitment::AjtaiParams;
use crate::utils::mle::dense_vec_to_dense_mle;
use ark_std::iterable::Iterable;
use lattirust_arithmetic::ring::PolyRing;
use lattirust_arithmetic::{
    challenge_set::latticefold_challenge_set::OverField,
    mle::DenseMultilinearExtension,
    polynomials::{build_eq_x_r, eq_eval, VPAuxInfo, VirtualPolynomial},
};
use num_traits::zero;

use crate::{
    arith::{Witness, CCS, LCCCS},
    transcript::Transcript,
    utils::sumcheck::{
        prover::SumCheckProver, verifier::SumCheckVerifier, SumCheckIP, SumCheckProof,
    },
};

use super::{decomposition::DecompositionParams, error::FoldingError, NIFSProver, NIFSVerifier};
use libm::log2;

#[derive(Clone)]
pub struct FoldingProof<NTT: OverField> {
    // Step 2.
    pub pointshift_sumcheck_proof: SumCheckProof<NTT>,
    // Step 3
    pub theta_s: Vec<NTT>,
    pub eta_s: Vec<Vec<NTT>>,
}

pub trait FoldingProver<CR: PolyRing, NTT: OverField, P: AjtaiParams, T: Transcript<NTT>> {
    type Proof: Clone;
    type Error: std::error::Error;

    fn prove(
        cm_i_s: &[LCCCS<NTT, P>],
        w_s: &[Witness<NTT>],
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(LCCCS<NTT, P>, Witness<NTT>, Self::Proof), Self::Error>;
}

pub trait FoldingVerifier<CR: PolyRing, NTT: OverField, P: AjtaiParams, T: Transcript<NTT>> {
    type Prover: FoldingProver<CR, NTT, P, T>;
    type Error = <Self::Prover as FoldingProver<CR, NTT, P, T>>::Error;

    fn verify(
        cm_i_s: &[LCCCS<NTT, P>],
        proof: &<Self::Prover as FoldingProver<CR, NTT, P, T>>::Proof,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<LCCCS<NTT, P>, Self::Error>;
}

impl<
        CR: PolyRing + From<NTT> + Into<NTT>,
        NTT: OverField,
        P: AjtaiParams,
        DP: DecompositionParams,
        T: Transcript<NTT>,
    > FoldingProver<CR, NTT, P, T> for NIFSProver<CR, NTT, P, DP, T>
{
    type Proof = FoldingProof<NTT>;
    type Error = FoldingError<NTT>;

    fn prove(
        _cm_i_s: &[LCCCS<NTT, P>],
        _w_s: &[Witness<NTT>],
        _transcript: &mut impl Transcript<NTT>,
        _ccs: &CCS<NTT>,
    ) -> Result<(LCCCS<NTT, P>, Witness<NTT>, FoldingProof<NTT>), FoldingError<NTT>> {
        assert_eq!(_cm_i_s.len(), 2 * DP::K);
        let m = _ccs.m;
        let log_m = log2(m as f64) as usize;

        // Step 1: Generate alpha, zeta, mu, beta challenges
        _cm_i_s.iter().for_each(|lcccs| {
            _transcript.absorb_ring_vec(&lcccs.r);
            _transcript.absorb_ring(&lcccs.v);
            // _transcript.absorb_ring_vec(&lcccs.cm); Not absorbed by transcript?
            _transcript.absorb_ring_vec(&lcccs.u);
            _transcript.absorb_ring_vec(&lcccs.x_w);
        });
        // TODO: Get challenges from big set but as NTT
        let alpha_s = _transcript.get_small_challenges(2 * DP::K);
        let zeta_s = _transcript.get_small_challenges(2 * DP::K);
        let mu_s = _transcript.get_small_challenges((2 * DP::K) - 1); // Note is one challenge less
        let beta_s = _transcript.get_small_challenges(log_m);

        // Step 2: Compute g polynomial and sumcheck on it
        // Setup f_hat_mle for later evaluation of thetas
        let f_hat_mles = _w_s
            .iter()
            .map(|w| {
                let f_i = w.f.clone();
                DenseMultilinearExtension::from_evaluations_vec(log_m, f_i)
            })
            .collect::<Vec<_>>();

        let zis = _cm_i_s
            .iter()
            .zip(_w_s.iter())
            .map(|(cm_i, w_i)| cm_i.get_z_vector(&w_i.w_ccs))
            .collect::<Vec<_>>();
        let ris = _cm_i_s
            .iter()
            .map(|cm_i| cm_i.r.clone())
            .collect::<Vec<_>>();
        let vs = _cm_i_s.iter().map(|cm_i| cm_i.v).collect::<Vec<NTT>>();
        let us = _cm_i_s
            .iter()
            .map(|cm_i| cm_i.u.clone())
            .collect::<Vec<_>>();

        // Setup matrix_mles for later evaluation of etas
        // Review creation of this Mi*z mles
        let Mz_mles_vec: Vec<Vec<DenseMultilinearExtension<NTT>>> = zis
            .iter()
            .map(|zi| {
                let Mz_mle = _ccs
                    .M
                    .iter()
                    .map(|M| {
                        let M_times_z = mat_vec_mul(&M, &zi)?;
                        let Mz_mle = dense_vec_to_dense_mle(log_m, &M_times_z);
                        Ok(Mz_mle)
                    })
                    .collect()?;
                Ok(Mz_mle)
            })
            .collect::<Result<_, _>>()?;
        let Mz_mle = Mz_mles_vec.iter().map(|Mz_mles| {});

        let g = create_sumcheck_polynomial::<NTT, DP>(
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
                ui.iter()
                    .fold(NTT::zero(), |acc, &u_i_t| acc + (u_i_t * zeta))
            });

        let prover = SumCheckProver {
            polynomial: g,
            claimed_sum: claim_g1 + claim_g2,
            _marker: std::marker::PhantomData::default(),
        };

        // Run sumcheck
        let (sum_check_proof, subclaim) = prover.prove(_transcript).unwrap();
        let r0 = subclaim.point;

        // Step 3: Evaluate thetas and etas
        let thetas = f_hat_mles
            .iter()
            .map(|f_hat_mle| f_hat_mle.evaluate(r0.as_slice()).unwrap())
            .collect::<Vec<_>>();
        drop(f_hat_mles);
        let etas: Vec<Vec<NTT>> = Mz_mles_vec
            .iter()
            .map(|Mz_mles| {
                Mz_mles
                    .iter()
                    .map(|mle| mle.evaluate(r0.as_slice()).unwrap())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        drop(Mz_mles_vec);

        _transcript.absorb_ring_vec(&thetas);
        etas.iter()
            .for_each(|etas| _transcript.absorb_ring_vec(&etas));

        // Step 5 get rho challenges
        let rhos = _transcript.get_small_challenges((2 * DP::K) - 1); // Note that we are missing the first element

        // Step 6 compute v0, u0, y0, x_w0
        let v0: NTT = rhos
            .iter()
            .zip(thetas.iter().skip(1))
            .fold(thetas[0], |acc, (rho_i, theta_i)| {
                // acc + rho_i.rot_sum(theta_i) // Note that theta_i is already in NTT form
                todo!() // Add WithRot to OverField in lattirust
            }); // Do INTT here

        let y_o = rhos.iter().zip(_cm_i_s.iter().skip(1))
            .fold(_cm_i_s[0].cm, |acc, (rho, cm_i)| {
                acc + (cm_i.cm * rho)
            });

        let u_0 = rhos
            .iter()
            .zip(etas.iter().skip(1)) // Skip the first eta because rho = 1
            .fold(etas[0], |acc, (&rho, &eta)| {
                let new_eta = eta.iter().map(|e| rho * e).collect::<Vec<_>>();
                acc.iter()
                    .zip(new_eta.iter())
                    .map(|(&a, &e)| a + e)
                    .collect()
            });

        let x_0 = rhos.iter().zip(_cm_i_s.iter().skip(1)) // Skip the first x_w because rho = 1
            .fold(_cm_i_s[0].x_w, |acc, (rho, cm_i)| {
                let mut x_w = cm_i.x_w.clone();
                x_w.iter_mut().map(|&mut x| x * rho);
                acc.iter().zip(x_w.iter())
                .map(|(a, x)| *a + *x).collect()
            });

        // Step 7: Compute f0 and Witness_0
        let f_0 = rhos
            .iter()
            .zip(_w_s.iter().skip(1))
            .fold(_w_s[0].f, |acc, (rho, w_i_s)| {
                let mut f_i = w_i_s.f.clone();
                f_i.iter_mut().for_each(|c| *c = *c * rho);
                acc.iter().zip(f_i.iter()).map(|(a, f)| *a + f).collect()
            });
        let w0 = Witness::from_f(f_0);

        let lcccs = LCCCS {
            r: r0,
            v: v0[0],
            cm: y_o,
            u: u_0,
            x_w: x_0,
            h: todo!()
        };

        let folding_proof = FoldingProof {
            pointshift_sumcheck_proof: sum_check_proof,
            theta_s: thetas,
            eta_s: etas,
        };
        Ok((lcccs, w0, folding_proof))
    }
}

impl<
        CR: PolyRing + From<NTT> + Into<NTT>,
        NTT: OverField,
        P: AjtaiParams,
        DP: DecompositionParams,
        T: Transcript<NTT>,
    > FoldingVerifier<CR, NTT, P, T> for NIFSVerifier<CR, NTT, P, DP, T>
{
    type Prover = NIFSProver<CR, NTT, P, DP, T>;

    fn verify(
        _cm_i_s: &[LCCCS<NTT, P>],
        _proof: &<Self::Prover as FoldingProver<CR, NTT, P, T>>::Proof,
        _transcript: &mut impl Transcript<NTT>,
        _ccs: &CCS<NTT>,
    ) -> Result<LCCCS<NTT, P>, FoldingError<NTT>> {
        let m = _ccs.m;
        let k_times_2 = _cm_i_s.len();
        let log_m = log2(m as f64) as usize;
        // Generate challenges
        // Note: absorb commits
        let alphas: Vec<NTT> = (0..k_times_2)
            .map(|_| _transcript.get_big_challenge().into())
            .collect::<Vec<_>>();
        let zetas: Vec<NTT> = (0..k_times_2)
            .map(|_| _transcript.get_big_challenge().into())
            .collect::<Vec<_>>();
        let mut mus: Vec<NTT> = (0..k_times_2 - 1)
            .map(|_| _transcript.get_big_challenge().into())
            .collect::<Vec<_>>();
        let beta_s: Vec<NTT> = (0..log_m)
            .map(|_| _transcript.get_big_challenge().into())
            .collect::<Vec<_>>();

        let poly_info = VPAuxInfo {
            max_degree: _ccs.d + 1,
            num_variables: log_m,
            phantom: std::marker::PhantomData,
        };
        let zis: Vec<Vec<NTT>> = Vec::with_capacity(_ccs.M.len()); // Grab zis from decomposition step
        let ris: Vec<Vec<NTT>> = Vec::new(); // Grab ris from decomposition step
        let vs = _cm_i_s.iter().map(|cm_i| cm_i.v).collect::<Vec<NTT>>();
        let us: Vec<NTT> = Vec::new(); // Grab us from the decomposition step
        let claim_g1 = alphas
            .iter()
            .zip(vs.iter())
            .fold(NTT::zero(), |acc, (&alpha, &vi)| acc + (alpha * vi));
        let claim_g2 = zetas
            .iter()
            .zip(us.iter())
            .fold(NTT::zero(), |acc, (&zeta, &ui)| acc + (zeta * ui));
        let protocol = SumCheckIP {
            claimed_sum: claim_g1 + claim_g2,
            poly_info,
        };
        let verifier = SumCheckVerifier::new(protocol);
        let sub_claim = verifier
            .verify(&_proof.pointshift_sumcheck_proof, _transcript)
            .unwrap();
        let e_asterisk = eq_eval(&beta_s, &sub_claim.point).unwrap();
        let e_i_s: Vec<NTT> = ris
            .iter()
            .map(|r| eq_eval(r.as_slice(), &sub_claim.point).unwrap())
            .collect::<Vec<_>>();
        let s = sub_claim.expected_evaluation.clone();

        let b = 2 as u64; // Get this from decomposition step and also remove from create_sumcheck_poly
        let mut should_equal_s = NTT::one();
        for i in 0..mus.len() {
            let res = _proof.theta_s[i].clone();
            should_equal_s = (0..b).fold(res, |acc, j| {
                let j_ring = NTT::from(j);
                acc * (_proof.theta_s[i] - j_ring)
            });
            should_equal_s = (0..b).fold(should_equal_s, |acc, j| {
                let j_ring = NTT::from(j);
                acc * (_proof.theta_s[i] + j_ring)
            });
            should_equal_s = should_equal_s * mus[i];
        }
        should_equal_s = should_equal_s * e_asterisk;
        for i in 0..e_i_s.len() {
            should_equal_s = should_equal_s + (alphas[i] * e_i_s[i] * _proof.theta_s[i]);
        }
        for i in 0..e_i_s.len() {
            should_equal_s = should_equal_s + (zetas[i] * e_i_s[i] * _proof.eta_s[i]);
        }
        match should_equal_s == s {
            true => {}
            false => {
                return Err(FoldingError::SumCheckError(
                    crate::utils::sumcheck::SumCheckError::SumCheckFailed(should_equal_s, s),
                ));
            }
        }

        let mut rhos = Vec::with_capacity(k_times_2); // need to absorb here as well
        rhos.push(NTT::one());
        for _ in 1..k_times_2 {
            rhos.push(_transcript.get_small_challenge());
        }

        // get y0, u0, v0 and x_w0

        todo!()
    }
}

fn create_matrix_mle<NTT: OverField>(
    log_m: usize,
    Mi: &Vec<Vec<NTT>>,
    zi: &Vec<NTT>,
) -> DenseMultilinearExtension<NTT> {
    let zero_vector = usize_to_binary_vector::<NTT>(0, log2(Mi.len() as f64) as usize);
    let mle_z_ccs_b = mle_val_from_vector(&zi, &zero_vector);
    let evaluations: Vec<NTT> = mle_matrix_to_val_eval_second(&Mi, &zero_vector)
        .iter()
        .map(|val| *val * mle_z_ccs_b)
        .collect();
    let mle = DenseMultilinearExtension::from_evaluations_vec(log_m, evaluations);

    let matrix_mle = (1..Mi.len())
        .into_iter()
        .map(|i| usize_to_binary_vector::<NTT>(i, log2(Mi.len() as f64) as usize))
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
        let gi_1 = create_g1_i_polynomial(&f_hat_mles[i], alpha_is[i]);
        let gi_3 = create_g3_i_polynomial(&matrix_mles[i], zeta_is[i]);
        let gi_1_plus_gi_3 = Arc::from(gi_1 + gi_3);
        let mut g1_and_g3_virtual = VirtualPolynomial::new_from_mle(&gi_1_plus_gi_3, NTT::ONE);
        let eq_r_i = build_eq_x_r(ris[i].as_slice()).unwrap();
        g1_and_g3_virtual.mul_by_mle(eq_r_i, NTT::one());
        g1_plus_g3 = &g1_plus_g3 + &g1_and_g3_virtual;
    }

    let b = DP::SMALL_B; // Get this from the decomposition step
    let mut g2 = create_g2_i_polynomial(log_m, &f_hat_mles[0], b, mus[0]);
    for i in 1..2 * k - 1 {
        let gi_2 = create_g2_i_polynomial(log_m, &f_hat_mles[i], b, mus[i]);
        g2 = &g2 + &gi_2;
    }
    // Recall that the last mu is 1
    let gi_2 = create_g2_i_polynomial(log_m, &f_hat_mles[2 * k - 1], b, NTT::one());
    g2 = &g2 + &gi_2;

    let eq_beta = build_eq_x_r::<NTT>(Beta.as_slice()).unwrap();
    g2.mul_by_mle(eq_beta, NTT::one());
    g = &g1_plus_g3 + &g2;

    g
}

fn create_g1_i_polynomial<NTT: OverField>(
    fi_mle: &DenseMultilinearExtension<NTT>,
    alpha_i: NTT,
) -> DenseMultilinearExtension<NTT> {
    let mut mle = fi_mle.clone(); // remove clone
    mle.evaluations.iter_mut().for_each(|e| *e = *e * alpha_i);
    mle
}

fn create_g2_i_polynomial<NTT: OverField>(
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

fn create_g3_i_polynomial<NTT: OverField>(
    matrix_mle: &Vec<DenseMultilinearExtension<NTT>>,
    zeta_i: NTT,
) -> DenseMultilinearExtension<NTT> {
    let (first_mle, mles) = matrix_mle.split_first().unwrap();
    let first_mle = first_mle.clone();
    let mut mle = mles
        .into_iter()
        .fold(first_mle, |acc, mle_i_t| acc + mle_i_t.clone())
        .to_owned();
    mle.evaluations.iter_mut().for_each(|e| *e = *e * zeta_i);
    mle
}

fn mle_val_from_vector<NTT: OverField>(vector: &Vec<NTT>, values: &Vec<NTT>) -> NTT {
    assert_eq!(values.len(), log2(vector.len() as f64) as usize);
    let mle = DenseMultilinearExtension::from_evaluations_vec(values.len(), vector.clone());
    mle.evaluate(values.as_slice()).unwrap()
}

// Convert a bivariate MLE to a univariate MLE by evaluating the second vector
fn mle_matrix_to_val_eval_second<NTT: OverField>(
    matrix: &Vec<Vec<NTT>>,
    values_y: &Vec<NTT>,
) -> Vec<NTT> {
    assert_eq!(values_y.len(), log2(matrix.len() as f64) as usize);
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
