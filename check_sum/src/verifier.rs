use ark_crypto_primitives::sponge::Absorb;
use ark_ff::PrimeField;
use lattirust_arithmetic::challenge_set::latticefold_challenge_set::{
    LatticefoldChallengeSet,
    OverField,
};

use crate::sum_check_transcript::SumCheckIP;
pub struct SumCheckVerifier<F: PrimeField, R: OverField<F>, CS: LatticefoldChallengeSet<F, R>>
    where F: Absorb {
    _marker: std::marker::PhantomData<(F, CS)>,
    transcript: SumCheckIP<F, R, CS>,
}

impl<F: PrimeField, R: OverField<F>, CS: LatticefoldChallengeSet<F, R>> SumCheckVerifier<F, R, CS>
    where F: Absorb
{
    pub fn new(transcript: SumCheckIP<F, R, CS>) -> SumCheckVerifier<F, R, CS> {
        SumCheckVerifier {
            _marker: std::marker::PhantomData,
            transcript,
        }
    }
    // Verifies the transcript
    // Returns true if the transcript represents a true claim
    // False otherwise
    pub fn verify(&self) -> bool {
        let mut check_sum = self.transcript.claimed_sum;
        let poly = &self.transcript.polynomial;
        let rounds_len = self.transcript.rounds.len();
        for (i, round) in self.transcript.rounds.iter().enumerate() {
            let j = rounds_len - 1 - i;

            let eval1 = round.unipoly.evaluate(&[R::one()]);
            let eval0 = round.unipoly.evaluate(&[R::zero()]);
            match (&eval1, &eval0) {
                (Ok(_), Ok(_)) => {}
                _ => {
                    return false;
                }
            }
            let sum = eval1.unwrap() + eval0.unwrap();
            if sum != check_sum {
                return false;
            } else {
                // if round.unipoly.degree() > poly.multi_degree()[j] {
                //     return false;
                // } else {
                check_sum = round.unipoly.evaluate(&[round.challenge]).unwrap();
                // }
            }
        }
        let challenge_vector: Vec<R> = self.transcript.rounds
            .iter()
            .map(|x| x.challenge)
            .collect();
        let challenges: &[R] = &challenge_vector;
        let oracle_evaluation = self.transcript.polynomial.evaluate(&challenges).unwrap();
        return oracle_evaluation == check_sum;
    }
}
#[cfg(test)]
mod test {
    use std::sync::Arc;

    use crate::{ prover::SumCheckProver, verifier::SumCheckVerifier };

    use ark_ff::{ One, Zero };
    use lattirust_arithmetic::{
        challenge_set::latticefold_challenge_set::BinarySmallSet,
        mle::DenseMultilinearExtension,
        polynomials::VirtualPolynomial,
        ring::{ Pow2CyclotomicPolyRingNTT, Zq },
    };

    #[test]
    fn test_sumcheck_protocol() {
        // Define the modulus Q and the dimension N
        const Q: u64 = 17; // Replace with an appropriate modulus
        const N: usize = 8; // Replace with an appropriate dimension

        // Example function to generate coefficients
        fn generate_coefficient(index: usize) -> Zq<Q> {
            Zq::<Q>::from(index as u64) // Simple example: use the index as the coefficient value
        }

        // Create an instance of Pow2CyclotomicPolyRingNTT using from_fn
        let poly_ntt = Pow2CyclotomicPolyRingNTT::<Q, N>::from_fn(generate_coefficient);

        let mle = DenseMultilinearExtension::from_evaluations_slice(
            1,
            &[Pow2CyclotomicPolyRingNTT::<Q, N>::zero(), poly_ntt]
        );
        let polynomial = VirtualPolynomial::new_from_mle(
            &Arc::from(mle),
            Pow2CyclotomicPolyRingNTT::<Q, N>::one()
        );
        // Define the claimed sum for testing
        let claimed_sum = poly_ntt; // Example sum

        // Create an instance of the prover
        let prover = SumCheckProver::<
            Zq<Q>,
            Pow2CyclotomicPolyRingNTT<Q, N>,
            BinarySmallSet<Q, N>
        > {
            _marker: std::marker::PhantomData,
            polynomial: polynomial.clone(),
            claimed_sum: Into::into(claimed_sum),
        };

        // Prove the statement
        let transcript = prover.prove();

        // Create an instance of the verifier
        let verifier = SumCheckVerifier::<
            Zq<Q>,
            Pow2CyclotomicPolyRingNTT<Q, N>,
            BinarySmallSet<Q, N>
        >::new(transcript);

        // Verify the transcript
        let result = verifier.verify();
        assert_eq!(result, true)
    }
}
