use ark_crypto_primitives::sponge::{
    poseidon::{PoseidonConfig, PoseidonSponge},
    CryptographicSponge,
};
use ark_ff::{BigInteger, Field, PrimeField};
use ark_std::marker::PhantomData;
use lattirust_ring::{OverField, PolyRing};

use super::{Transcript, TranscriptWithSmallChallenges};
use cyclotomic_rings::{challenge_set::LatticefoldChallengeSet, GetPoseidonParams, SuitableRing};

/// PoseidonTranscript implements the Transcript trait using the Poseidon hash
pub struct PoseidonTranscript<R: OverField, CS> {
    _marker: PhantomData<CS>,
    sponge: PoseidonSponge<<R::BaseRing as Field>::BasePrimeField>,
}

impl<R: SuitableRing, CS: LatticefoldChallengeSet<R>> Default for PoseidonTranscript<R, CS> {
    fn default() -> Self {
        Self::new(&R::PoseidonParams::get_poseidon_config())
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
                .into_iter()
                .map(|x| x.to_base_prime_field_elements()).flatten().collect::<Vec<_>>(),
        );
    }

    fn absorb_slice(&mut self, v: &[R]) {
        for ring in v {
            self.absorb(ring);
        }
    }

    fn get_big_challenge(&mut self) -> R::BaseRing {
        let extension_degree = R::BaseRing::extension_degree();
        let c = self
            .sponge
            .squeeze_field_elements(extension_degree as usize);
        self.sponge.absorb(&c);
        <R::BaseRing as Field>::from_base_prime_field_elems(&c).expect("something went wrong: c does not contain extension_degree elements")
    }
}

impl<R: SuitableRing, CS: LatticefoldChallengeSet<R>> TranscriptWithSmallChallenges<R>
    for PoseidonTranscript<R, CS>
{
    type ChallengeSet = CS;

    fn get_small_challenge(&mut self) -> R::CoefficientRepresentation {
        todo!()
        // let c: Vec<R::BaseRing> = self.sponge.squeeze_field_elements(1);
        // self.sponge.absorb(&c);
        // Self::ChallengeSet::small_challenge_from_random_bytes(&c[0].into_bigint().to_bytes_be())
    }
}
