use lattirust_arithmetic::challenge_set::latticefold_challenge_set::{
    LatticefoldChallengeSet, OverField,
};

use crate::transcript::Transcript;

use super::SumCheckError;
use super::SumCheckIP;
use super::SumCheckProof;
use super::SumCheckSubClaim;

pub struct SumCheckVerifier<R: OverField, CS: LatticefoldChallengeSet<R>> {
    _marker: std::marker::PhantomData<CS>,
    protocol: SumCheckIP<R>,
}

impl<R: OverField, CS: LatticefoldChallengeSet<R>> SumCheckVerifier<R, CS> {
    pub fn new(protocol: SumCheckIP<R>) -> SumCheckVerifier<R, CS> {
        SumCheckVerifier {
            _marker: std::marker::PhantomData,
            protocol,
        }
    }
    // Verifies the transcript
    // Returns evaluation subclaim,
    // returns an error otherwise.
    pub fn verify(
        &self,
        proof: &SumCheckProof<R>,
        transcript: &mut impl Transcript<R, ChallengeSet = CS>,
    ) -> Result<SumCheckSubClaim<R>, SumCheckError<R>> {
        let mut check_sum = self.protocol.claimed_sum;

        let mut challenge_vector: Vec<R> = Vec::with_capacity(proof.rounds.len());

        for round in proof.rounds.iter() {
            let eval1 = round
                .unipoly
                .coeffs
                .iter()
                .fold(R::zero(), |sum, x| sum + x);
            let eval0 = round.unipoly.coeffs[0];
            let sum = eval1 + eval0;

            if sum != check_sum {
                return Err(SumCheckError::SumCheckFailed(sum, check_sum));
            } else if round.unipoly.degree() > self.protocol.poly_info.max_degree {
                return Err(SumCheckError::MaxDegreeExceeded);
            }

            transcript.absorb_ring_vec(&round.unipoly.coeffs);

            let challenge = transcript.get_big_challenge().into();

            check_sum = round.unipoly.evaluate(challenge);
            challenge_vector.push(challenge);
        }

        Ok(SumCheckSubClaim {
            expected_evaluation: check_sum,
            point: challenge_vector,
        })
    }
}
#[cfg(test)]
mod test {
    use std::sync::Arc;

    use super::SumCheckVerifier;
    use crate::{
        transcript::poseidon::PoseidonTranscript,
        utils::sumcheck::{prover::SumCheckProver, SumCheckIP},
    };

    use ark_ff::{One, Zero};
    use lattirust_arithmetic::{
        challenge_set::latticefold_challenge_set::BinarySmallSet,
        mle::DenseMultilinearExtension,
        polynomials::{VPAuxInfo, VirtualPolynomial},
        ring::{Pow2CyclotomicPolyRingNTT, Zq},
    };

    #[test]
    fn test_sumcheck_protocol() {
        // Define the modulus Q and the dimension N
        const Q: u64 = 17; // Replace with an appropriate modulus
        const N: usize = 8; // Replace with an appropriate dimension

        // Example function to generate coefficients
        fn generate_coefficient_i(index: usize) -> Zq<Q> {
            Zq::<Q>::from(index as u64) // Simple example: use the index as the coefficient value
        }

        // Create an instance of Pow2CyclotomicPolyRingNTT using from_fn
        let poly_ntt = Pow2CyclotomicPolyRingNTT::<Q, N>::from_fn(generate_coefficient_i);

        let one = Pow2CyclotomicPolyRingNTT::<Q, N>::one();
        let zero = Pow2CyclotomicPolyRingNTT::<Q, N>::zero();
        let mle1 = DenseMultilinearExtension::from_evaluations_slice(
            3,
            &[poly_ntt, poly_ntt, zero, zero, zero, zero, zero, zero],
        );
        let mle2 = DenseMultilinearExtension::from_evaluations_slice(
            3,
            &[zero, zero, poly_ntt, poly_ntt, zero, zero, zero, zero],
        );
        let mle3 = DenseMultilinearExtension::from_evaluations_slice(
            3,
            &[zero, zero, zero, zero, poly_ntt, poly_ntt, zero, zero],
        );
        let mle4 = DenseMultilinearExtension::from_evaluations_slice(
            3,
            &[zero, zero, zero, zero, zero, zero, poly_ntt, poly_ntt],
        );
        let mut polynomial = VirtualPolynomial::new(3);
        polynomial
            .add_mle_list(
                vec![
                    Arc::from(mle1.clone()),
                    Arc::from(mle2),
                    Arc::from(mle3),
                    Arc::from(mle4),
                ],
                one,
            )
            .unwrap();
        polynomial.add_mle_list(vec![Arc::from(mle1)], one).unwrap();
        // Define the claimed sum for testing
        let claimed_sum = poly_ntt + poly_ntt; // Example sum

        let protocol = SumCheckIP::new(claimed_sum, VPAuxInfo::new(3, 4));
        // Create an instance of the prover
        let prover = SumCheckProver::<Pow2CyclotomicPolyRingNTT<Q, N>, BinarySmallSet<Q, N>>::new(
            polynomial.clone(),
            claimed_sum,
        );

        let mut transcript = PoseidonTranscript::default();

        // Prove the statement
        let (proof, _subclaim) = prover.prove(&mut transcript).unwrap();

        let mut transcript = PoseidonTranscript::default();
        // Create an instance of the verifier
        let verifier =
            SumCheckVerifier::<Pow2CyclotomicPolyRingNTT<Q, N>, BinarySmallSet<Q, N>>::new(
                protocol,
            );

        // Verify the transcript
        let subclaim = verifier.verify(&proof, &mut transcript).unwrap();
        assert!(
            polynomial.evaluate(&subclaim.point).unwrap() == subclaim.expected_evaluation,
            "wrong subclaim"
        );
    }
    #[test]
    fn test_failing_sumcheck_protocol() {
        // Define the modulus Q and the dimension N
        const Q: u64 = 17; // Replace with an appropriate modulus
        const N: usize = 8; // Replace with an appropriate dimension

        // Example function to generate coefficients
        fn generate_coefficient_i(index: usize) -> Zq<Q> {
            Zq::<Q>::from(index as u64) // Simple example: use the index as the coefficient value
        }

        // Create an instance of Pow2CyclotomicPolyRingNTT using from_fn
        let poly_ntt = Pow2CyclotomicPolyRingNTT::<Q, N>::from_fn(generate_coefficient_i);

        let one = Pow2CyclotomicPolyRingNTT::<Q, N>::one();
        let zero = Pow2CyclotomicPolyRingNTT::<Q, N>::zero();
        let mle1 = DenseMultilinearExtension::from_evaluations_slice(
            3,
            &[poly_ntt, poly_ntt, zero, zero, zero, zero, zero, zero],
        );
        let mle2 = DenseMultilinearExtension::from_evaluations_slice(
            3,
            &[zero, zero, poly_ntt, poly_ntt, zero, zero, zero, zero],
        );
        let mle3 = DenseMultilinearExtension::from_evaluations_slice(
            3,
            &[zero, zero, zero, zero, poly_ntt, poly_ntt, zero, zero],
        );
        let mle4 = DenseMultilinearExtension::from_evaluations_slice(
            3,
            &[zero, zero, zero, zero, zero, zero, poly_ntt, poly_ntt],
        );
        let mut polynomial = VirtualPolynomial::new(3);
        polynomial
            .add_mle_list(
                vec![
                    Arc::from(mle1.clone()),
                    Arc::from(mle2),
                    Arc::from(mle3),
                    Arc::from(mle4),
                ],
                one,
            )
            .unwrap();
        polynomial.add_mle_list(vec![Arc::from(mle1)], one).unwrap();
        // Define the claimed sum for testing
        let claimed_sum = poly_ntt; // Example sum

        let mut transcript = PoseidonTranscript::default();
        // Create an instance of the prover

        let protocol = SumCheckIP::new(claimed_sum, VPAuxInfo::new(3, 4));

        let prover = SumCheckProver::<Pow2CyclotomicPolyRingNTT<Q, N>, BinarySmallSet<Q, N>>::new(
            polynomial,
            claimed_sum,
        );
        // Prove the statement
        let (proof, _) = prover.prove(&mut transcript).unwrap();

        // Create an instance of the verifier
        let verifier =
            SumCheckVerifier::<Pow2CyclotomicPolyRingNTT<Q, N>, BinarySmallSet<Q, N>>::new(
                protocol,
            );

        let mut transcript = PoseidonTranscript::default();
        // Verify the transcript
        let subclaim = verifier.verify(&proof, &mut transcript);
        assert!(
            subclaim.is_err(),
            "Wrong sumcheck claim should return error"
        );
    }
}
