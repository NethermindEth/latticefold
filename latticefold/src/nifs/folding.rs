use std::sync::Arc;

use lattirust_arithmetic::{
    challenge_set::latticefold_challenge_set::OverField,
    mle::{self, DenseMultilinearExtension},
    polynomials::{build_eq_x_r, VirtualPolynomial},
};

use crate::{
    arith::{Witness, CCCS, CCS, LCCCS},
    transcript::Transcript,
    utils::sumcheck::SumCheckProof,
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
        let k = 0; // this should come from the decomposition step

        // Generate challenges
        // Note: absorb commits
        let alphas = (0..2 * k)
            .map(|_| _transcript.get_big_challenge())
            .collect::<Vec<_>>();
        let zetas = (0..2 * k)
            .map(|_| _transcript.get_big_challenge())
            .collect::<Vec<_>>();
        let mut mus = (0..2 * k - 1)
            .map(|_| _transcript.get_big_challenge())
            .collect::<Vec<_>>();
        let Beta = (0..log_m)
            .map(|_| _transcript.get_big_challenge())
            .collect::<Vec<_>>();
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
    fis: &[Vec<R>],
    alpha_is: Vec<R>,
    Mis: &[Vec<Vec<R>>],
    zis: &[Vec<R>],
    zeta_is: Vec<R>,
    ris: Vec<Vec<R>>,
    Beta: &Vec<R>,
) -> VirtualPolynomial<R> {
    let mut g = VirtualPolynomial::new(2 * k);
    let mut g1_plus_g3 = VirtualPolynomial::new(2 * k);
    for i in 0..2 * k {
        let gi_1 = create_g1_i_polynomial(log_m, fis[i].as_slice(), alpha_is[i]);
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

fn create_gi_polynomial<R: OverField>(
    log_m: usize,
    Mi: Vec<Vec<R>>,
    fi: Vec<R>,
    alpha_i: R,
    zeta_i: R,
    zi: Vec<R>,
    ri: &[R],
) -> VirtualPolynomial<R> {
    let g1 = create_g1_i_polynomial(log_m, &fi, alpha_i);
    let g2 = create_g2_i_polynomial::<R>();
    let g3 = create_g3_i_polynomial(log_m, Mi, zi, zeta_i);
    let g1_g3 = Arc::from(g1 + g3);
    let mut g1_and_g3_virtual = VirtualPolynomial::new_from_mle(&g1_g3, R::ONE);
    let eq_r_i = build_eq_x_r(ri).unwrap();
    g1_and_g3_virtual.mul_by_mle(eq_r_i, R::one());
    g1_and_g3_virtual =
        &g1_and_g3_virtual + &VirtualPolynomial::new_from_mle(&Arc::from(g2), R::one());
    todo!()
}

fn create_g1_i_polynomial<R: OverField>(
    log_m: usize,
    fi: &[R],
    alpha_i: R,
) -> DenseMultilinearExtension<R> {
    let evals = fi.iter().map(|&e| e * alpha_i).collect();
    DenseMultilinearExtension::from_evaluations_vec(log_m, evals)
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
        .map(|val| zeta_i * *val * mle_z_ccs_b)
        .collect();
    let mle = DenseMultilinearExtension::from_evaluations_vec(log_m, evaluations);

    let matrix_mle = (1..Mi.len())
        .into_iter()
        .map(|i| usize_to_binary_vector::<R>(i, log2(Mi.len() as f64) as usize))
        .fold(mle, |acc, b| {
            let mle_z_ccs_b = mle_val_from_vector(&zi, &b);
            let evaluations: Vec<R> = mle_matrix_to_val_eval_second(&Mi, &b)
                .iter()
                .map(|val| zeta_i * *val * mle_z_ccs_b)
                .collect();
            let mle = DenseMultilinearExtension::from_evaluations_vec(log_m, evaluations);
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
