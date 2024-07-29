use crate::{
    arith::{ Witness, CCCS, LCCCS },
    transcript::Transcript,
    utils::sumcheck::SumCheckProof,
};
use ark_crypto_primitives::sponge::Absorb;
use ark_ff::Field;
use ark_std::iterable::Iterable;
use lattirust_arithmetic::{
    challenge_set::latticefold_challenge_set::OverField,
    mle::DenseMultilinearExtension,
    ring::PolyRing,
};
use libm::log2;

use super::{ error::LinearizationError, NIFSProver, NIFSVerifier };

#[derive(Clone)]
pub struct LinearizationProof<F: Field, R: OverField<F>> where F: Absorb {
    // Sent in the step 2. of the linearization subprotocol
    pub linearization_sumcheck: SumCheckProof<F, R>,
    // Sent in the step 3.
    pub v: R,
    pub u: Vec<R>,
}

pub trait LinearizationProver<F: Field, R: OverField<F>, T: Transcript<F, R>> {
    type Proof: Clone;
    type Error: std::error::Error;

    fn prove(
        cm_i: &CCCS<R>,
        wit: &Witness<R>,
        transcript: &mut impl Transcript<F, R>
    ) -> Result<(LCCCS<R>, Self::Proof), Self::Error>;
}

pub trait LinearizationVerifier<F: Field, R: OverField<F>, T: Transcript<F, R>> {
    type Prover: LinearizationProver<F, R, T>;
    type Error = <Self::Prover as LinearizationProver<F, R, T>>::Error;

    fn verify(
        cm_i: &CCCS<R>,
        proof: &<Self::Prover as LinearizationProver<F, R, T>>::Proof,
        transcript: &mut impl Transcript<F, R>
    ) -> Result<LCCCS<R>, Self::Error>;
}

impl<F: Field, R: OverField<F>, T: Transcript<F, R>> LinearizationProver<F, R, T>
    for NIFSProver<F, R, T>
    where F: Absorb
{
    type Proof = LinearizationProof<F, R>;
    type Error = LinearizationError<R>;

    fn prove(
        _cm_i: &CCCS<R>,
        _wit: &Witness<R>,
        _transcript: &mut impl Transcript<F, R>
    ) -> Result<(LCCCS<R>, LinearizationProof<F, R>), LinearizationError<R>> {
        // Step 1: Get the public coin randomness
        // Step 2:
        todo!()
    }
}

impl<F: Field, R: OverField<F>, T: Transcript<F, R>> LinearizationVerifier<F, R, T>
    for NIFSVerifier<F, R, T>
    where F: Absorb
{
    type Prover = NIFSProver<F, R, T>;

    fn verify(
        _cm_i: &CCCS<R>,
        _proof: &<Self::Prover as LinearizationProver<F, R, T>>::Proof,
        _transcript: &mut impl Transcript<F, R>
    ) -> Result<LCCCS<R>, LinearizationError<R>> {
        todo!()
    }
}

fn get_challenge_vector<F: Field, R: OverField<F>, T: Transcript<F, R>>(
    _cm_i: &CCCS<R>,
    _transcript: &mut T,
    len: usize
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
fn eq<F: Field, R: OverField<F>>(b_arr: &Vec<R>, x_arr: &Vec<R>) -> R {
    assert_eq!(b_arr.len(), x_arr.len(), "Eq function takes two vectors of the same length!");
    b_arr
        .iter()
        .zip(x_arr)
        .fold(R::one(), |ret_value, (b, x)| {
            ret_value * ((R::one() - b) * (R::one() - x) + *b * x)
        })
}

fn usize_to_binary_vector<F: Field, R: OverField<F>>(n: usize) -> Vec<R> {
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

fn mle_val_from_vector<F: Field, R: OverField<F>>(vector: &Vec<R>, values: &Vec<R>) -> R {
    assert_eq!(values.len(), log2(vector.len() as f64) as usize);
    let mle = DenseMultilinearExtension::from_evaluations_vec(values.len(), vector.clone());
    mle.evaluate(values.as_slice()).unwrap()
}
fn mle_val_from_matrix<F: Field, R: OverField<F>>(
    matrix: &Vec<Vec<R>>,
    values_x: &Vec<R>,
    values_y: &Vec<R>
) -> R {
    assert_eq!(values_y.len(), log2(matrix.len() as f64) as usize);
    assert_eq!(values_x.len(), log2(matrix[0].len() as f64) as usize);
    let univariate_mle_evaluations = matrix
        .into_iter()
        .map(|row| mle_val_from_vector(row, values_x))
        .collect();
    mle_val_from_vector(&univariate_mle_evaluations, values_y)
}

#[cfg(test)]
mod tests {
    use ark_ff::{ One, Zero };
    use lattirust_arithmetic::ring::{ Pow2CyclotomicPolyRingNTT, Zq };

    use crate::nifs::linearization::{
        eq,
        mle_val_from_matrix,
        mle_val_from_vector,
        usize_to_binary_vector,
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
        assert_eq!(mle_val_from_vector(&evaluation_vector, &vec![one()]), zero());
        assert_ne!(mle_val_from_vector(&evaluation_vector, &vec![one()]), poly_ntt());
        assert_eq!(mle_val_from_vector(&evaluation_vector, &vec![zero()]), poly_ntt());
        assert_ne!(mle_val_from_vector(&evaluation_vector, &vec![zero()]), zero());

        // Test evaluation of bivariate MLE from a matrix
        let evaluation_matrix = vec![vec![poly_ntt(), zero()], vec![one(), poly_ntt()]];
        assert_eq!(mle_val_from_matrix(&evaluation_matrix, &vec![zero()], &vec![one()]), one());
        assert_ne!(
            mle_val_from_matrix(&evaluation_matrix, &vec![zero()], &vec![one()]),
            poly_ntt()
        );
        assert_eq!(
            mle_val_from_matrix(&evaluation_matrix, &vec![zero()], &vec![zero()]),
            poly_ntt()
        );
        assert_ne!(mle_val_from_matrix(&evaluation_matrix, &vec![zero()], &vec![zero()]), zero());

        // Test the eq function
        let vector_one = vec![zero(), one(), one(), zero()];
        let vector_two = vec![zero(), one(), one(), zero()];
        let vector_three = vec![zero(), one(), one(), one()];

        assert_eq!(eq(&vector_one, &vector_two), one());
        assert_ne!(eq(&vector_one, &vector_two), zero());
        assert_eq!(eq(&vector_one, &vector_three), zero());
        assert_ne!(eq(&vector_one, &vector_three), one());

        assert_eq!(
            usize_to_binary_vector::<Zq<Q>, Pow2CyclotomicPolyRingNTT<Q, N>>(4),
            vec![one(), zero(), zero()]
        );
        assert_eq!(
            usize_to_binary_vector::<Zq<Q>, Pow2CyclotomicPolyRingNTT<Q, N>>(5),
            vec![one(), zero(), one()]
        )
    }
}
