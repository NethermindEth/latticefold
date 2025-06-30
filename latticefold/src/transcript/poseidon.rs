use ark_crypto_primitives::sponge::{
    poseidon::{PoseidonConfig, PoseidonSponge},
    CryptographicSponge,
};
use ark_ff::Field;
use ark_std::marker::PhantomData;
use cyclotomic_rings::{
    challenge_set::LatticefoldChallengeSet,
    rings::{GetPoseidonParams, SuitableRing},
};
use stark_rings::OverField;

use super::{Transcript, TranscriptWithShortChallenges};
use crate::ark_base::*;

/// PoseidonTranscript implements the Transcript trait using the Poseidon hash
#[derive(Clone)]
pub struct PoseidonTranscript<R: OverField, CS> {
    _marker: PhantomData<CS>,
    sponge: PoseidonSponge<<R::BaseRing as Field>::BasePrimeField>,
}

impl<R: SuitableRing, CS: LatticefoldChallengeSet<R>> Default for PoseidonTranscript<R, CS> {
    fn default() -> Self {
        Self::new(&R::PoseidonParams::get_poseidon_config())
    }
}

#[derive(Clone)]
pub struct PoseidonTS<R: OverField> {
    sponge: PoseidonSponge<<R::BaseRing as Field>::BasePrimeField>,
}
impl<R: OverField> PoseidonTS<R> {
    pub fn default<P: GetPoseidonParams<<<R>::BaseRing as Field>::BasePrimeField>>() -> Self {
        Self::new(&P::get_poseidon_config())
    }
}
impl<R: OverField> Transcript<R> for PoseidonTS<R> {
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
}

impl<R: OverField, CS> Transcript<R> for PoseidonTranscript<R, CS> {
    type TranscriptConfig = PoseidonConfig<<R::BaseRing as Field>::BasePrimeField>;

    fn new(config: &Self::TranscriptConfig) -> Self {
        let sponge = PoseidonSponge::<<R::BaseRing as Field>::BasePrimeField>::new(config);
        Self {
            sponge,
            _marker: PhantomData,
        }
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
}

impl<R: SuitableRing, CS: LatticefoldChallengeSet<R>> TranscriptWithShortChallenges<R>
    for PoseidonTranscript<R, CS>
{
    type ChallengeSet = CS;

    fn get_short_challenge(&mut self) -> R::CoefficientRepresentation {
        let random_bytes = self.sponge.squeeze_bytes(Self::ChallengeSet::BYTES_NEEDED);

        Self::ChallengeSet::short_challenge_from_random_bytes(&random_bytes)
            .expect("not enough bytes to get a small challenge")
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::BigInt;
    use cyclotomic_rings::rings::{GoldilocksChallengeSet, GoldilocksRingNTT, GoldilocksRingPoly};
    use stark_rings::cyclotomic_ring::models::goldilocks::{Fq, Fq3};

    use super::*;

    #[test]
    fn test_get_big_challenge() {
        let mut transcript =
            PoseidonTranscript::<GoldilocksRingNTT, GoldilocksChallengeSet>::default();

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
        let mut transcript =
            PoseidonTranscript::<GoldilocksRingNTT, GoldilocksChallengeSet>::default();

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

        assert_eq!(expected, transcript.get_short_challenge())
    }
}
