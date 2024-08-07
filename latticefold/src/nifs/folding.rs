use std::sync::Arc;

use ark_std::iterable::Iterable;
use lattirust_arithmetic::{
    challenge_set::latticefold_challenge_set::OverField,
    mle::{self, DenseMultilinearExtension},
    polynomials::{build_eq_x_r, eq_eval, VPAuxInfo, VirtualPolynomial},
};
use crate::commitment::AjtaiParams;
use lattirust_arithmetic::challenge_set::latticefold_challenge_set::OverField;
use lattirust_arithmetic::ring::PolyRing;

use crate::{
    arith::{Witness, CCCS, CCS, LCCCS},
    transcript::Transcript,
    utils::sumcheck::{prover::SumCheckProver, verifier::SumCheckVerifier, SumCheckIP, SumCheckProof},
};

use libm::log2;
use super::{decomposition::DecompositionParams, error::FoldingError, NIFSProver, NIFSVerifier};

#[derive(Clone)]
pub struct FoldingProof<NTT: OverField> {
    // Step 2.
    pub pointshift_sumcheck_proof: SumCheckProof<NTT>,
    // Step 3
    pub theta_s: Vec<NTT>,
    pub eta_s: Vec<NTT>,
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
        _cm_i_s: &[LCCCS<R>],
        _w_s: &[Witness<R>],
        _transcript: &mut impl Transcript<R>,
        _ccs: &CCS<R>,
    ) -> Result<(LCCCS<R>, Witness<R>, FoldingProof<R>), FoldingError<R>> {
        let m = _ccs.m;
        let log_m = log2(m as f64) as usize;
        let k = 10000000; // this should come from the decomposition step

        // Generate challenges
        // Note: absorb commits
        let alphas: Vec<R> = (0..2 * k)
            .map(|_| _transcript.get_big_challenge().into())
            .collect::<Vec<_>>();
        let zetas: Vec<R> = (0..2 * k)
            .map(|_| _transcript.get_big_challenge().into())
            .collect::<Vec<_>>();
        let mut mus: Vec<R> = (0..2 * k - 1)
            .map(|_| _transcript.get_big_challenge().into())
            .collect::<Vec<_>>();
        mus.push(R::one());
        let Beta: Vec<R> = (0..log_m)
            .map(|_| _transcript.get_big_challenge().into())
            .collect::<Vec<_>>();

        let f_hat_mles = _w_s
            .iter()
            .map(|w| {
                let f_i = w.f_arr.clone();
                DenseMultilinearExtension::from_evaluations_vec(log_m, f_i)
            })
            .collect::<Vec<_>>();

        let zis: Vec<Vec<R>> = Vec::with_capacity(_ccs.M.len()); // Grab zis from decomposition step
        let ris: Vec<Vec<R>> = Vec::new(); // Grab ris from decomposition step
        let vs = _cm_i_s.iter().map(|cm_i| cm_i.v).collect::<Vec<R>>();
        let us: Vec<R> = Vec::new(); // Grab us from the decomposition step

        let matrix_mles = _ccs
            .M
            .iter()
            .zip(zis.iter())
            .map(|(Mi, zi)| create_matrix_mle(log_m, &Mi, &zi))
            .collect::<Vec<_>>();

        let g = create_sumcheck_polynomial(
            k,
            log_m,
            &f_hat_mles,
            &alphas,
            &matrix_mles,
            &zetas,
            ris,
            &Beta,
            &mus,
        );


        let claim_g1 = alphas
            .iter()
            .zip(vs.iter())
            .fold(R::zero(), |acc, (&alpha, &vi)| acc + (alpha * vi));
        let claim_g2 = zetas
            .iter()
            .zip(us.iter())
            .fold(R::zero(), |acc, (&zeta, &ui)| acc + (zeta * ui));

        let prover = SumCheckProver {
            polynomial: g,
            claimed_sum: claim_g1 + claim_g2,
            _marker: std::marker::PhantomData::default(),
        };

        // Run sum check prover
        let (_, sum_check_proof, subclaim) = prover.prove(_transcript).unwrap();
        let r0: Vec<R> = Vec::new(); //  take r0 from transcript

        let thetas = f_hat_mles
            .iter()
            .map(|f_hat_mle| f_hat_mle.evaluate(r0.as_slice()).unwrap())
            .collect::<Vec<_>>();
        drop(f_hat_mles);
        let etas = matrix_mles
            .iter()
            .map(|matrix_mle| matrix_mle.evaluate(r0.as_slice()).unwrap())
            .collect::<Vec<_>>();
        drop(matrix_mles);

        let mut rhos = Vec::with_capacity(2*k); // need to absorb here as well
        rhos.push(R::one());
        for _ in 1..2*k {
            rhos.push(_transcript.get_small_challenge());
        }
        
        // let yi_s = _cm_i_s.iter().map(|cm_i| cm_i.y);
        let u_0 = rhos.iter().zip(etas.iter()).fold(R::zero(), |acc, (&rho, &eta)| {
            acc + (rho * eta)
        });
        // let x_w_len = _cm_i_s[0].x_w.len();
        // let x_0 = rhos.iter().zip(_cm_i_s.iter())
        //     .fold(vec![R::zero(); x_w_len], |acc, (rho, cm_i)| {
        //         let mut x_w = cm_i.x_w.clone();
        //         x_w.iter_mut().map(|&mut x| x * rho);
        //         acc.iter().zip(x_w.iter())
        //         .map(|(a, x)| *a + *x).collect()
        //     });
        // let y_o = rhos.iter().zip(_cm_i_s.iter())
        //     .fold(vec![R::zero(); x_w_len], |acc, (rho, cm_i)| {
        //         let mut y_i = cm_i.y.clone();
        //         y_i.iter_mut().map(|&mut x| x * rho);
        //         acc.iter().zip(y_i.iter())
        //         .map(|(a, x)| *a + *x).collect()
        //     });

        let f_0 = rhos.iter().zip(_w_s.iter())
            .fold(vec![R::zero();4], |acc, (rho, w_i_s)| {
                let mut f_i = w_i_s.f_arr.clone();
                f_i.iter_mut().for_each(|c| *c = *c * rho);
                acc.iter().zip(f_i.iter()).map(|(a, f)| *a + f).collect()
            }
        );

        let folding_proof = FoldingProof {
            pointshift_sumcheck_proof: sum_check_proof,
            theta_s: thetas,
            eta_s: etas,
        };


        todo!()
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
        _cm_i_s: &[LCCCS<R>],
        _proof: &<Self::Prover as FoldingProver<R, T>>::Proof,
        _transcript: &mut impl Transcript<R>,
        _ccs: &CCS<R>,
    ) -> Result<LCCCS<R>, FoldingError<R>> {
        let m = _ccs.m;
        let k_times_2 = _cm_i_s.len();
        let log_m = log2(m as f64) as usize;
        // Generate challenges
        // Note: absorb commits
        let alphas: Vec<R> = (0..k_times_2)
            .map(|_| _transcript.get_big_challenge().into())
            .collect::<Vec<_>>();
        let zetas: Vec<R> = (0..k_times_2)
            .map(|_| _transcript.get_big_challenge().into())
            .collect::<Vec<_>>();
        let mut mus: Vec<R> = (0..k_times_2 - 1)
            .map(|_| _transcript.get_big_challenge().into())
            .collect::<Vec<_>>();
        mus.push(R::one());
        let Beta: Vec<R> = (0..log_m)
            .map(|_| _transcript.get_big_challenge().into())
            .collect::<Vec<_>>();

        let poly_info = VPAuxInfo {
            max_degree: _ccs.d + 1,
            num_variables: log_m,
            phantom: std::marker::PhantomData,
        };
        let zis: Vec<Vec<R>> = Vec::with_capacity(_ccs.M.len()); // Grab zis from decomposition step
        let ris: Vec<Vec<R>> = Vec::new(); // Grab ris from decomposition step
        let vs = _cm_i_s.iter().map(|cm_i| cm_i.v).collect::<Vec<R>>();
        let us: Vec<R> = Vec::new(); // Grab us from the decomposition step
        let claim_g1 = alphas
            .iter()
            .zip(vs.iter())
            .fold(R::zero(), |acc, (&alpha, &vi)| acc + (alpha * vi));
        let claim_g2 = zetas
            .iter()
            .zip(us.iter())
            .fold(R::zero(), |acc, (&zeta, &ui)| acc + (zeta * ui));
        let protocol = SumCheckIP {
            claimed_sum: claim_g1 + claim_g2,
            poly_info,
        };
        let verifier = SumCheckVerifier::new(protocol);
        let sub_claim = verifier
            .verify(&_proof.pointshift_sumcheck_proof, _transcript).unwrap();
        let e_asterisk = eq_eval(&Beta, &sub_claim.point).unwrap();
        let e_i_s: Vec<R> = ris.iter().map(|r| eq_eval(r, &sub_claim.point).unwrap()).collect::<Vec<_>>();
        let s = sub_claim.expected_evaluation.clone();

        let b = 2 as u64; // Get this from decomposition step and also remove from create_sumcheck_poly
        let mut should_equal_s = R::one();
        for i in 0..mus.len() {
            let res = _proof.theta_s[i].clone();
            should_equal_s = (0..b).fold(res, |acc, j| {
                let j_ring = R::from(j);
                acc*(_proof.theta_s[i] - j_ring)
            });
            should_equal_s = (0..b).fold(should_equal_s, |acc, j| {
                let j_ring = R::from(j);
                acc*(_proof.theta_s[i] + j_ring)
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
            true => {},
            false => {
                return Err(FoldingError::SumCheckError(crate::utils::sumcheck::SumCheckError::SumCheckFailed(should_equal_s, s)));
            }
        }

        let mut rhos = Vec::with_capacity(k_times_2); // need to absorb here as well
        rhos.push(R::one());
        for _ in 1..k_times_2 {
            rhos.push(_transcript.get_small_challenge());
        }

        // get y0, u0, v0 and x_w0

        todo!()
    }
}

fn create_matrix_mle<R: OverField>(
    log_m: usize,
    Mi: &Vec<Vec<R>>,
    zi: &Vec<R>,
) -> DenseMultilinearExtension<R> {
    let zero_vector = usize_to_binary_vector::<R>(0, log2(Mi.len() as f64) as usize);
    let mle_z_ccs_b = mle_val_from_vector(&zi, &zero_vector);
    let evaluations: Vec<R> = mle_matrix_to_val_eval_second(&Mi, &zero_vector)
        .iter()
        .map(|val| *val * mle_z_ccs_b)
        .collect();
    let mle = DenseMultilinearExtension::from_evaluations_vec(log_m, evaluations);

    let matrix_mle = (1..Mi.len())
        .into_iter()
        .map(|i| usize_to_binary_vector::<R>(i, log2(Mi.len() as f64) as usize))
        .fold(mle, |acc, b| {
            let mle_z_ccs_b = mle_val_from_vector(&zi, &b);
            let evaluations: Vec<R> = mle_matrix_to_val_eval_second(&Mi, &b)
                .iter()
                .map(|val| *val * mle_z_ccs_b)
                .collect();
            let mle = DenseMultilinearExtension::from_evaluations_vec(log_m, evaluations);
            acc + mle
        });
    matrix_mle
}

fn create_sumcheck_polynomial<R: OverField>(
    k: usize,
    log_m: usize,
    f_hat_mles: &Vec<DenseMultilinearExtension<R>>,
    alpha_is: &Vec<R>,
    matrix_mles: &Vec<DenseMultilinearExtension<R>>,
    zeta_is: &Vec<R>,
    ris: Vec<Vec<R>>,
    Beta: &Vec<R>,
    mus: &Vec<R>,
) -> VirtualPolynomial<R> {
    let mut g = VirtualPolynomial::new(2 * k);
    let mut g1_plus_g3 = VirtualPolynomial::new(2 * k);
    for i in 0..2 * k {
        let gi_1 = create_g1_i_polynomial(&f_hat_mles[i], alpha_is[i]);
        let gi_3 = create_g3_i_polynomial(&matrix_mles[i], zeta_is[i]);
        let gi_1_plus_gi_3 = Arc::from(gi_1 + gi_3);
        let mut g1_and_g3_virtual = VirtualPolynomial::new_from_mle(&gi_1_plus_gi_3, R::ONE);
        let eq_r_i = build_eq_x_r(ris[i].as_slice()).unwrap();
        g1_and_g3_virtual.mul_by_mle(eq_r_i, R::one());
        g1_plus_g3 = &g1_plus_g3 + &g1_and_g3_virtual;
    }

    let b = 2; // Get this from the decomposition step
    let mut g2 = create_g2_i_polynomial(log_m, &f_hat_mles[0], b, mus[0]);
    for i in 1..2 * k {
        let gi_2 = create_g2_i_polynomial(log_m, &f_hat_mles[i], b, mus[i]);
        g2 = &g2 + &gi_2;
    }
    let eq_beta = build_eq_x_r::<R>(Beta.as_slice()).unwrap();
    g2.mul_by_mle(eq_beta, R::one());
    g = &g1_plus_g3 + &g2;

    g
}

fn create_g1_i_polynomial<R: OverField>(
    fi_mle: &DenseMultilinearExtension<R>,
    alpha_i: R,
) -> DenseMultilinearExtension<R> {
    let mut mle = fi_mle.clone(); // remove clone
    mle.evaluations.iter_mut().for_each(|e| *e = *e * alpha_i);
    mle
}

fn create_g2_i_polynomial<R: OverField>(
    log_m: usize,
    fi_mle: &DenseMultilinearExtension<R>,
    b: u64,
    mu_i: R,
) -> VirtualPolynomial<R> {
    let mut mle_list: Vec<Arc<DenseMultilinearExtension<R>>> = Vec::new();
    let mle_zero = fi_mle.clone();
    mle_list.push(Arc::from(mle_zero));
    for i in 0..b {
        let mut mle_j = fi_mle.clone();
        mle_j
            .evaluations
            .iter_mut()
            .for_each(|e| *e = *e - R::from(i)); // There should be a better way than sub every
                                                 // eval
        mle_list.push(Arc::from(mle_j));
    }
    for i in 0..b {
        let mut mle_j = fi_mle.clone();
        mle_j
            .evaluations
            .iter_mut()
            .for_each(|e| *e = *e + R::from(i));
        mle_list.push(Arc::from(mle_j));
    }
    let mut gi_2 = VirtualPolynomial::new(log_m);
    gi_2.add_mle_list(mle_list, mu_i);
    gi_2
}

fn create_g3_i_polynomial<R: OverField>(
    matrix_mle: &DenseMultilinearExtension<R>,
    zeta_i: R,
) -> DenseMultilinearExtension<R> {
    let mut mle = matrix_mle.clone();
    mle.evaluations.iter_mut().for_each(|e| *e = *e * zeta_i);
    mle
}

fn mle_val_from_vector<R: OverField>(vector: &Vec<R>, values: &Vec<R>) -> R {
    assert_eq!(values.len(), log2(vector.len() as f64) as usize);
    let mle = DenseMultilinearExtension::from_evaluations_vec(values.len(), vector.clone());
    mle.evaluate(values.as_slice()).unwrap()
}

// Convert a bivariate MLE to a univariate MLE by evaluating the second vector
fn mle_matrix_to_val_eval_second<R: OverField>(matrix: &Vec<Vec<R>>, values_y: &Vec<R>) -> Vec<R> {
    assert_eq!(values_y.len(), log2(matrix.len() as f64) as usize);
    (0..matrix[0].len())
        .into_iter()
        .map(|i| mle_val_from_vector(&matrix.iter().map(|col| col[i]).collect(), values_y))
        .collect()
}

fn usize_to_binary_vector<R: OverField>(n: usize, length: usize) -> Vec<R> {
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
        .map(|bit| if bit == 1 { R::one() } else { R::zero() })
        .collect()
}
