use ark_crypto_primitives::sponge::{
    poseidon::{PoseidonConfig, PoseidonSponge},
    CryptographicSponge,
};
use ark_ff::{BigInteger, PrimeField};
use ark_std::marker::PhantomData;
use lattirust_ring::{OverField, PolyRing};

use super::Transcript;
use cyclotomic_rings::{challenge_set::LatticefoldChallengeSet, GetPoseidonParams, SuitableRing};

/// PoseidonTranscript implements the Transcript trait using the Poseidon hash
pub struct PoseidonTranscript<R: OverField, CS: LatticefoldChallengeSet<R>>
where
    <R as PolyRing>::BaseRing: PrimeField,
{
    _marker: PhantomData<CS>,
    sponge: PoseidonSponge<R::BaseRing>,
}

impl<R: SuitableRing, CS: LatticefoldChallengeSet<R>> Default for PoseidonTranscript<R, CS>
where
    <R as PolyRing>::BaseRing: PrimeField,
{
    fn default() -> Self {
        Self::new(&R::PoseidonParams::get_poseidon_config())
    }
}

impl<R: OverField, CS: LatticefoldChallengeSet<R>> Transcript<R> for PoseidonTranscript<R, CS>
where
    <R as PolyRing>::BaseRing: PrimeField,
{
    type TranscriptConfig = PoseidonConfig<R::BaseRing>;

    type ChallengeSet = CS;

    fn new(config: &Self::TranscriptConfig) -> Self {
        let sponge = PoseidonSponge::<R::BaseRing>::new(config);
        Self {
            sponge,
            _marker: PhantomData,
        }
    }

    fn absorb(&mut self, v: &R) {
        self.sponge.absorb(&v.coeffs());
    }

    fn absorb_slice(&mut self, v: &[R]) {
        for ring in v {
            self.absorb(ring);
        }
    }

    fn get_big_challenge(&mut self) -> R::BaseRing {
        let c: Vec<R::BaseRing> = self.sponge.squeeze_field_elements(1);
        self.sponge.absorb(&c);
        c[0]
    }

    fn get_small_challenge(&mut self) -> R {
        let c: Vec<R::BaseRing> = self.sponge.squeeze_field_elements(1);
        self.sponge.absorb(&c);
        Self::ChallengeSet::small_challenge_from_random_bytes(&c[0].into_bigint().to_bytes_be())
    }
}
