use ark_crypto_primitives::sponge::Absorb;
use ark_ff::PrimeField;
use lattirust_arithmetic::challenge_set::latticefold_challenge_set::{
    LatticefoldChallengeSet,
    OverField,
};

use crate::sum_check_transcript::SumCheckTranscript;
pub struct SumCheckVerifier<F: PrimeField, R: OverField<F>, CS: LatticefoldChallengeSet<F, R>>
    where F: Absorb {
    _marker: std::marker::PhantomData<(F, CS)>,
    transcript: SumCheckTranscript<F, R, CS>,
}

impl<F: PrimeField, R: OverField<F>, CS: LatticefoldChallengeSet<F, R>> SumCheckVerifier<F, R, CS>
    where F: Absorb
{
    pub fn new(transcript: SumCheckTranscript<F, R, CS>) -> SumCheckVerifier<F, R, CS> {
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

            let sum = round.unipoly.at_one() + round.unipoly.at_zero();
            if sum != check_sum {
                return false;
            } else {
                if round.unipoly.degree() > poly.multi_degree()[j] {
                    return false;
                } else {
                    check_sum = round.unipoly.eval(&round.challenge);
                }
            }
        }
        let challenges = self.transcript.rounds
            .iter()
            .map(|x| x.challenge)
            .collect();
        let oracle_evaluation = self.transcript.polynomial.eval_poly(&challenges);
        return oracle_evaluation == check_sum;
    }
}
#[cfg(test)]
mod test {
    use crate::{ poly_utils::MultiPoly, prover::SumCheckProver, verifier::SumCheckVerifier };

    use lattirust_arithmetic::{
        challenge_set::latticefold_challenge_set::BinarySmallSet,
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

        // Create a MultiPoly instance
        let poly = MultiPoly {
            terms: vec![(poly_ntt, vec![(0, 1)])],
        };

        // Define the claimed sum for testing
        let claimed_sum = poly_ntt; // Example sum

        // Create an instance of the prover
        let prover = SumCheckProver::<
            Zq<Q>,
            Pow2CyclotomicPolyRingNTT<Q, N>,
            BinarySmallSet<Q, N>
        > {
            _marker: std::marker::PhantomData,
            polynomial: poly.clone(),
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
