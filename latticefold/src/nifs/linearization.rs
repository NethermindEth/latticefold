use crate::{
    arith::{Witness, CCCS, CCS, LCCCS},
    transcript::Transcript,
    utils::sumcheck::SumCheckProof,
};
use ark_ff::Field;
use ark_std::iterable::Iterable;
use lattirust_arithmetic::{
    challenge_set::latticefold_challenge_set::OverField, mle::DenseMultilinearExtension,
    ring::PolyRing,
};
use libm::log2;

use super::{error::LinearizationError, NIFSProver, NIFSVerifier};

#[derive(Clone)]
pub struct LinearizationProof<R: OverField> {
    // Sent in the step 2. of the linearization subprotocol
    pub linearization_sumcheck: SumCheckProof<R>,
    // Sent in the step 3.
    pub v: R,
    pub u: Vec<R>,
}

pub trait LinearizationProver<R: OverField, T: Transcript<R>> {
    type Proof: Clone;
    type Error: std::error::Error;

    fn prove(
        cm_i: &CCCS<R>,
        wit: &Witness<R>,
        transcript: &mut impl Transcript<R>,
        ccs: &CCS<R>,
    ) -> Result<(LCCCS<R>, Self::Proof), Self::Error>;
}

pub trait LinearizationVerifier<R: OverField, T: Transcript<R>> {
    type Prover: LinearizationProver<R, T>;
    type Error = <Self::Prover as LinearizationProver<R, T>>::Error;

    fn verify(
        cm_i: &CCCS<R>,
        proof: &<Self::Prover as LinearizationProver<R, T>>::Proof,
        transcript: &mut impl Transcript<R>,
        ccs: &CCS<R>,
    ) -> Result<LCCCS<R>, Self::Error>;
}

impl<R: OverField, T: Transcript<R>> LinearizationProver<R, T> for NIFSProver<R, T> {
    type Proof = LinearizationProof<R>;
    type Error = LinearizationError<R>;

    fn prove(
        _cm_i: &CCCS<R>,
        _wit: &Witness<R>,
        _transcript: &mut impl Transcript<R>,
        _ccs: &CCS<R>,
    ) -> Result<(LCCCS<R>, LinearizationProof<R>), LinearizationError<R>> {
        todo!()
    }
}

impl<R: OverField, T: Transcript<R>> LinearizationVerifier<R, T> for NIFSVerifier<R, T> {
    type Prover = NIFSProver<R, T>;

    fn verify(
        _cm_i: &CCCS<R>,
        _proof: &<Self::Prover as LinearizationProver<R, T>>::Proof,
        _transcript: &mut impl Transcript<R>,
        _ccs: &CCS<R>,
    ) -> Result<LCCCS<R>, LinearizationError<R>> {
        todo!()
    }
}

fn get_challenge_vector<R: OverField, T: Transcript<R>>(
    _cm_i: &CCCS<R>,
    _transcript: &mut T,
    len: usize,
) -> Vec<<R as PolyRing>::BaseRing> {
    (0..len)
        .map(|_| {
            let challenge = _transcript.get_big_challenge();
            _transcript.absorb_ring_vec(&_cm_i.cm);
            challenge
        })
        .collect()
}

// Outputs 1 if vectors are the same, zero otherwise
fn eq<R: OverField>(b_arr: &Vec<R>, x_arr: &Vec<R>) -> R {
    assert_eq!(
        b_arr.len(),
        x_arr.len(),
        "Eq function takes two vectors of the same length!"
    );
    b_arr.iter().zip(x_arr).fold(R::one(), |ret_value, (b, x)| {
        ret_value * ((R::one() - b) * (R::one() - x) + *b * x)
    })
}

fn usize_to_binary_vector<R: OverField>(n: usize) -> Vec<R> {
    let mut bits = Vec::new();
    let mut current = n;

    while current > 0 {
        bits.push((current & 1) as u8);
        current >>= 1;
    }

    bits.reverse();
    bits.iter()
        .map(|bit| if *bit == 1 { R::one() } else { R::zero() })
        .collect()
}

fn mle_val_from_vector<R: OverField>(vector: &Vec<R>, values: &Vec<R>) -> R {
    assert_eq!(values.len(), log2(vector.len() as f64) as usize);
    let mle = DenseMultilinearExtension::from_evaluations_vec(values.len(), vector.clone());
    mle.evaluate(values.as_slice()).unwrap()
}

fn mle_val_from_matrix<R: OverField>(
    matrix: &Vec<Vec<R>>,
    values_x: &Vec<R>,
    values_y: &Vec<R>,
) -> R {
    assert_eq!(values_y.len(), log2(matrix.len() as f64) as usize);
    assert_eq!(values_x.len(), log2(matrix[0].len() as f64) as usize);
    let univariate_mle_evaluations = matrix
        .into_iter()
        .map(|col| mle_val_from_vector(col, values_x))
        .collect();
    mle_val_from_vector(&univariate_mle_evaluations, values_y)
}
fn mle_matrix_to_val_eval_first<R: OverField>(matrix: &Vec<Vec<R>>, values_x: &Vec<R>) -> Vec<R> {
    assert_eq!(values_x.len(), log2(matrix[0].len() as f64) as usize);
    matrix
        .into_iter()
        .map(|col| mle_val_from_vector(col, values_x))
        .collect()
}
// Convert a bivariate MLE to a univariate MLE by evaluating the second vector
fn mle_matrix_to_val_eval_second<R: OverField>(matrix: &Vec<Vec<R>>, values_y: &Vec<R>) -> Vec<R> {
    assert_eq!(values_y.len(), log2(matrix.len() as f64) as usize);
    (0..matrix[0].len())
        .into_iter()
        .map(|i| mle_val_from_vector(&matrix.iter().map(|col| col[i]).collect(), values_y))
        .collect()
}

#[cfg(test)]
mod tests {
    use ark_ff::{One, Zero};
    use lattirust_arithmetic::ring::{Pow2CyclotomicPolyRingNTT, Zq};

    use crate::nifs::linearization::{
        eq, mle_matrix_to_val_eval_first, mle_matrix_to_val_eval_second, mle_val_from_matrix,
        mle_val_from_vector, usize_to_binary_vector,
    };

    // Boilerplate code to generate values needed for testing
    const Q: u64 = 17; // Replace with an appropriate modulus
    const N: usize = 8;

    fn generate_coefficient_i(index: usize) -> Zq<Q> {
        Zq::<Q>::from(index as u64) // Simple example: use the index as the coefficient value
    }
    fn poly_ntt() -> Pow2CyclotomicPolyRingNTT<Q, N> {
        Pow2CyclotomicPolyRingNTT::<Q, N>::from_fn(generate_coefficient_i)
    }
    fn one() -> Pow2CyclotomicPolyRingNTT<Q, N> {
        Pow2CyclotomicPolyRingNTT::<Q, N>::one()
    }
    fn zero() -> Pow2CyclotomicPolyRingNTT<17, 8> {
        Pow2CyclotomicPolyRingNTT::<Q, N>::zero()
    }

    // Actual Tests
    #[test]
    fn test_utils() {
        // Test evaluation of mle from a vector
        let evaluation_vector = vec![poly_ntt(), zero()];
        assert_eq!(
            mle_val_from_vector(&evaluation_vector, &vec![one()]),
            zero()
        );
        assert_ne!(
            mle_val_from_vector(&evaluation_vector, &vec![one()]),
            poly_ntt()
        );
        assert_eq!(
            mle_val_from_vector(&evaluation_vector, &vec![zero()]),
            poly_ntt()
        );
        assert_ne!(
            mle_val_from_vector(&evaluation_vector, &vec![zero()]),
            zero()
        );

        // Test evaluation of bivariate MLE from a matrix
        let evaluation_matrix = vec![vec![poly_ntt(), zero()], vec![one(), poly_ntt()]];
        assert_eq!(
            mle_val_from_matrix(&evaluation_matrix, &vec![zero()], &vec![one()]),
            one()
        );
        assert_ne!(
            mle_val_from_matrix(&evaluation_matrix, &vec![zero()], &vec![one()]),
            poly_ntt()
        );
        assert_eq!(
            mle_val_from_matrix(&evaluation_matrix, &vec![zero()], &vec![zero()]),
            poly_ntt()
        );
        assert_ne!(
            mle_val_from_matrix(&evaluation_matrix, &vec![zero()], &vec![zero()]),
            zero()
        );

        // Test the eq function
        let vector_one = vec![zero(), one(), one(), zero()];
        let vector_two = vec![zero(), one(), one(), zero()];
        let vector_three = vec![zero(), one(), one(), one()];

        assert_eq!(eq(&vector_one, &vector_two), one());
        assert_ne!(eq(&vector_one, &vector_two), zero());
        assert_eq!(eq(&vector_one, &vector_three), zero());
        assert_ne!(eq(&vector_one, &vector_three), one());

        assert_eq!(
            usize_to_binary_vector::<Pow2CyclotomicPolyRingNTT<Q, N>>(4),
            vec![one(), zero(), zero()]
        );
        assert_eq!(
            usize_to_binary_vector::<Pow2CyclotomicPolyRingNTT<Q, N>>(5),
            vec![one(), zero(), one()]
        );
        // Test the conversion of Bivariate MLE to univariate MLE by evaluating first values
        let bivariate_mle = vec![
            vec![poly_ntt(), poly_ntt(), one(), zero()],
            vec![zero(), poly_ntt(), zero(), one()],
        ];
        assert_eq!(
            mle_matrix_to_val_eval_first(&bivariate_mle, &vec![zero(), zero()]),
            vec![poly_ntt(), zero()]
        );
        assert_eq!(
            mle_matrix_to_val_eval_first(&bivariate_mle, &vec![one(), zero()]),
            vec![poly_ntt(), poly_ntt()]
        );

        // Test the conversion of Bivariate MLE to univariate MLE by evaluating second values
        assert_eq!(
            mle_matrix_to_val_eval_second(&bivariate_mle, &vec![one()]),
            vec![zero(), poly_ntt(), zero(), one()]
        );
        assert_eq!(
            mle_matrix_to_val_eval_second(&bivariate_mle, &vec![zero()]),
            vec![poly_ntt(), poly_ntt(), one(), zero()]
        );
    }
}
