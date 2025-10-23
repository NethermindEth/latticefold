use ark_crypto_primitives::sponge::{
    poseidon::{PoseidonConfig, PoseidonSponge},
    CryptographicSponge,
};
use ark_ff::Field;
use latticefold_rings::{
    challenge_set::ChallengeSet,
    rings::{GetPoseidonParams, SuitableRing},
};
use stark_rings::OverField;

use super::{Transcript, TranscriptWithSmallChallenges};
use crate::ark_base::*;

/// PoseidonTranscript implements the Transcript trait using the Poseidon hash
#[derive(Clone)]
pub struct PoseidonTranscript<R: OverField> {
    sponge: PoseidonSponge<<R::BaseRing as Field>::BasePrimeField>,
}

impl<R: SuitableRing> Default for PoseidonTranscript<R> {
    fn default() -> Self {
        Self::new(&R::PoseidonParams::get_poseidon_config())
    }
}

impl<R: OverField> Transcript<R> for PoseidonTranscript<R> {
    type TranscriptConfig = PoseidonConfig<<R::BaseRing as Field>::BasePrimeField>;

    fn new(config: &Self::TranscriptConfig) -> Self {
        let sponge = PoseidonSponge::<<R::BaseRing as Field>::BasePrimeField>::new(config);
        Self { sponge }
    }

    fn absorb(&mut self, v: &R) {
        self.sponge.absorb(
            &v.coeffs()
                .iter()
                .flat_map(|x| x.to_base_prime_field_elements())
                .collect::<Vec<_>>(),
        );
    }

    fn get_challenge(&mut self) -> R::BaseRing {
        let extension_degree = R::BaseRing::extension_degree();
        let c = self
            .sponge
            .squeeze_field_elements(extension_degree as usize);
        self.sponge.absorb(&c);
        <R::BaseRing as Field>::from_base_prime_field_elems(&c)
            .expect("something went wrong: c does not contain extension_degree elements")
    }

    fn squeeze_bytes(&mut self, n: usize) -> Vec<u8> {
        self.sponge.squeeze_bytes(n)
    }
}

impl<R: SuitableRing> TranscriptWithSmallChallenges<R> for PoseidonTranscript<R> {
    type ChallengeSet = R::ChallengeSet;

    fn get_small_challenge(&mut self) -> R::CoefficientRepresentation {
        let random_bytes = self.sponge.squeeze_bytes(Self::ChallengeSet::BYTES_NEEDED);

        Self::ChallengeSet::small_challenge_from_random_bytes(&random_bytes)
            .expect("not enough bytes to get a small challenge")
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::BigInt;
    use latticefold_rings::rings::{GoldilocksRingNTT, GoldilocksRingPoly};
    use stark_rings::cyclotomic_ring::models::goldilocks::{Fq, Fq3};

    use super::*;

    #[test]
    fn test_get_big_challenge() {
        let mut transcript = PoseidonTranscript::<GoldilocksRingNTT>::default();

        transcript
            .sponge
            .absorb(&Fq::from(BigInt::<1>::from(0xFFu32)));

        let expected: Fq3 = Fq3::new(
            Fq::new(BigInt([10462816198028961279])),
            Fq::new(BigInt([17217694161994925895])),
            Fq::new(BigInt([6163269596856181508])),
        );

        assert_eq!(expected, transcript.get_challenge())
    }

    #[test]
    fn test_get_small_challenge() {
        let mut transcript = PoseidonTranscript::<GoldilocksRingNTT>::default();

        transcript
            .sponge
            .absorb(&Fq::from(BigInt::<1>::from(0xFFu32)));

        let expected_coeffs: Vec<Fq> = vec![
            Fq::new(BigInt([31])),
            Fq::new(BigInt([18446744069414584312])),
            Fq::new(BigInt([18446744069414584291])),
            Fq::new(BigInt([14])),
            Fq::new(BigInt([18446744069414584306])),
            Fq::new(BigInt([18446744069414584312])),
            Fq::new(BigInt([30])),
            Fq::new(BigInt([18446744069414584313])),
            Fq::new(BigInt([19])),
            Fq::new(BigInt([18446744069414584317])),
            Fq::new(BigInt([20])),
            Fq::new(BigInt([18446744069414584306])),
            Fq::new(BigInt([18446744069414584295])),
            Fq::new(BigInt([4])),
            Fq::new(BigInt([18446744069414584320])),
            Fq::new(BigInt([7])),
            Fq::new(BigInt([18446744069414584298])),
            Fq::new(BigInt([18446744069414584295])),
            Fq::new(BigInt([18446744069414584304])),
            Fq::new(BigInt([18446744069414584290])),
            Fq::new(BigInt([3])),
            Fq::new(BigInt([18446744069414584304])),
            Fq::new(BigInt([25])),
            Fq::new(BigInt([18446744069414584304])),
        ];

        let expected = GoldilocksRingPoly::from(expected_coeffs);

        assert_eq!(expected, transcript.get_small_challenge())
    }
}
