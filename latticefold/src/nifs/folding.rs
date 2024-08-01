use std::sync::Arc;

use lattirust_arithmetic::{
    challenge_set::latticefold_challenge_set::OverField,
    mle::{self, DenseMultilinearExtension},
    polynomials::{build_eq_x_r, VirtualPolynomial},
};

use crate::{
    arith::{Witness, CCCS, CCS, LCCCS},
    transcript::Transcript,
    utils::sumcheck::{prover::SumCheckProver, SumCheckProof},
};

use super::{error::FoldingError, NIFSProver, NIFSVerifier};
use libm::log2;

#[derive(Clone)]
pub struct FoldingProof<R: OverField> {
    // Step 2.
    pub pointshift_sumcheck_proof: SumCheckProof<R>,
    // Step 3
    pub theta_s: Vec<R>,
    pub eta_s: Vec<R>,
}

pub trait FoldingProver<R: OverField, T: Transcript<R>> {
    type Proof: Clone;
    type Error: std::error::Error;

    fn prove(
        cm_i_s: &[LCCCS<R>],
        w_s: &[Witness<R>],
        transcript: &mut impl Transcript<R>,
        ccs: &CCS<R>,
    ) -> Result<(LCCCS<R>, Witness<R>, Self::Proof), Self::Error>;
}

pub trait FoldingVerifier<R: OverField, T: Transcript<R>> {
    type Prover: FoldingProver<R, T>;
    type Error = <Self::Prover as FoldingProver<R, T>>::Error;

    fn verify(
        cm_i_s: &[LCCCS<R>],
        proof: &<Self::Prover as FoldingProver<R, T>>::Proof,
        transcript: &mut impl Transcript<R>,
        ccs: &CCS<R>,
    ) -> Result<LCCCS<R>, Self::Error>;
}

impl<R: OverField, T: Transcript<R>> FoldingProver<R, T> for NIFSProver<R, T> {
    type Proof = FoldingProof<R>;
    type Error = FoldingError<R>;

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

        let zis: Vec<Vec<R>> = Vec::new(); // Grab zis from decomposition step
        let ris: Vec<Vec<R>> = Vec::new(); // Grab ris from decomposition step
        let vs: Vec<R> = Vec::new(); // Grab ris from decomposition step
        let us: Vec<R> = Vec::new(); // Grab us from the decomposition step

        let g = create_sumcheck_polynomial(
            k,
            log_m,
            &f_hat_mles,
            &alphas,
            &_ccs.M,
            zis.as_slice(),
            &zetas,
            ris,
            &Beta,
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
            .map(|f_hat_mle| f_hat_mle.evaluate(r0.as_slice()))
            .collect::<Vec<_>>();
        let etas: Vec<R> = Vec::new();

        todo!()
    }
}

impl<R: OverField, T: Transcript<R>> FoldingVerifier<R, T> for NIFSVerifier<R, T> {
    type Prover = NIFSProver<R, T>;

    fn verify(
        _cm_i_s: &[LCCCS<R>],
        _proof: &<Self::Prover as FoldingProver<R, T>>::Proof,
        _transcript: &mut impl Transcript<R>,
        _ccs: &CCS<R>,
    ) -> Result<LCCCS<R>, FoldingError<R>> {
        todo!()
    }
}

fn create_sumcheck_polynomial<R: OverField>(
    k: usize,
    log_m: usize,
    f_hat_mles: &Vec<DenseMultilinearExtension<R>>,
    alpha_is: &Vec<R>,
    Mis: &[Vec<Vec<R>>],
    zis: &[Vec<R>],
    zeta_is: &Vec<R>,
    ris: Vec<Vec<R>>,
    Beta: &Vec<R>,
) -> VirtualPolynomial<R> {
    let mut g = VirtualPolynomial::new(2 * k);
    let mut g1_plus_g3 = VirtualPolynomial::new(2 * k);
    for i in 0..2 * k {
        let gi_1 = create_g1_i_polynomial(&f_hat_mles[i], alpha_is[i]);
        let gi_3 = create_g3_i_polynomial(log_m, Mis[i].clone(), zis[i].clone(), zeta_is[i]);
        let gi_1_plus_gi_3 = Arc::from(gi_1 + gi_3);
        let mut g1_and_g3_virtual = VirtualPolynomial::new_from_mle(&gi_1_plus_gi_3, R::ONE);
        let eq_r_i = build_eq_x_r(ris[i].as_slice()).unwrap();
        g1_and_g3_virtual.mul_by_mle(eq_r_i, R::one());
        g1_plus_g3 = &g1_plus_g3 + &g1_and_g3_virtual;
    }

    let mut g2_dense = create_g2_i_polynomial();
    for i in 1..2 * k {
        let gi_2 = create_g2_i_polynomial();
        g2_dense = g2_dense + gi_2;
    }
    let g2_dense_arc = Arc::from(g2_dense);
    let mut g2 = VirtualPolynomial::new_from_mle(&g2_dense_arc, R::one());
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

fn create_g2_i_polynomial<R: OverField>() -> DenseMultilinearExtension<R> {
    todo!()
}

fn create_g3_i_polynomial<R: OverField>(
    log_m: usize,
    Mi: Vec<Vec<R>>,
    zi: Vec<R>,
    zeta_i: R,
) -> DenseMultilinearExtension<R> {
    let zero_vector = usize_to_binary_vector::<R>(0, log2(Mi.len() as f64) as usize);
    let mle_z_ccs_b = mle_val_from_vector(&zi, &zero_vector);
    let evaluations: Vec<R> = mle_matrix_to_val_eval_second(&Mi, &zero_vector)
        .iter()
        .map(|val| *val * mle_z_ccs_b)
        .collect();
    let mut mle = DenseMultilinearExtension::from_evaluations_vec(log_m, evaluations);
    mle.evaluations.iter_mut().for_each(|e| *e = *e * zeta_i);

    let matrix_mle = (1..Mi.len())
        .into_iter()
        .map(|i| usize_to_binary_vector::<R>(i, log2(Mi.len() as f64) as usize))
        .fold(mle, |acc, b| {
            let mle_z_ccs_b = mle_val_from_vector(&zi, &b);
            let evaluations: Vec<R> = mle_matrix_to_val_eval_second(&Mi, &b)
                .iter()
                .map(|val| *val * mle_z_ccs_b)
                .collect();
            let mut mle = DenseMultilinearExtension::from_evaluations_vec(log_m, evaluations);
            mle.evaluations.iter_mut().for_each(|e| *e = *e * zeta_i);
            acc + mle
        });
    matrix_mle
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
