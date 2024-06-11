use lattirust_arithmetic::ring::Ring;

use crate::transcript::{ self, SumCheckTranscript };
pub struct SumCheckVerifier<R: Ring> {
    transcript: SumCheckTranscript<R>,
}

impl<R: Ring> SumCheckVerifier<R> {
    pub fn new(transcript: SumCheckTranscript<R>) -> SumCheckVerifier<R> {
        SumCheckVerifier {
            transcript,
        }
    }
    // Verifies the transcript
    // Returns true if the transcript represents a true claim
    // False otherwise
    pub fn verify(&self) -> bool {
        let zero = R::zero();
        let one = R::one();
        let mut check_sum = self.transcript.claimed_sum;
        let mut poly = &self.transcript.polynomial;
        let rounds_len = self.transcript.rounds.len();
        for (i, round) in self.transcript.rounds.iter().enumerate() {
            let j = rounds_len - 1 - i;

            let sum = round.unipoly.eval(&one) + round.unipoly.eval(&zero);
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
    use crate::poly_utils::MultiPoly;
    use crate::prover::SumCheckProver;

    use lattirust_arithmetic::ring::Z2_64;
    use crate::transcript;

    use super::SumCheckVerifier;
    #[test]
    fn test_prover() {
        // Define a polynomial 3 + 2x1^2 + x1*x2 + 5x2 in Z2_64
        let poly = MultiPoly {
            terms: vec![
                (Z2_64::from(3 as u8), vec![]),
                (Z2_64::from(2 as u8), vec![(0, 2)]),
                (Z2_64::from(1 as u8), vec![(0, 1), (1, 1)]),
                (Z2_64::from(5 as u8), vec![(1, 1)])
            ],
        };
        let claimed_sum = Z2_64::from(27u8);

        let prover = SumCheckProver::<Z2_64> {
            polynomial: poly.clone().into(),
            cyclotomic_ring_modulus: 2,
            claimed_sum,
        };

        let transcript = prover.prove();
        let verifier = SumCheckVerifier::new(transcript);
        assert_eq!(verifier.verify(), true);
    }
}
