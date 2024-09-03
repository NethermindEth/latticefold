use ark_crypto_primitives::sponge::{
    poseidon::{PoseidonConfig, PoseidonSponge},
    CryptographicSponge,
};
use ark_ff::{BigInteger, PrimeField, Zero};
use ark_std::marker::PhantomData;
use lattirust_arithmetic::{
    challenge_set::latticefold_challenge_set::LatticefoldChallengeSet, ring::OverField,
};

use super::Transcript;

/// PoseidonTranscript implements the Transcript trait using the Poseidon hash
pub struct PoseidonTranscript<R: OverField, CS: LatticefoldChallengeSet<R>> {
    _marker: PhantomData<CS>,
    sponge: PoseidonSponge<R::BaseRing>,
}

impl<R: OverField, CS: LatticefoldChallengeSet<R>> Default for PoseidonTranscript<R, CS> {
    fn default() -> Self {
        let config = PoseidonConfig {
            full_rounds: 8, // Example values, adjust according to your needs
            partial_rounds: 57,
            alpha: 5,
            ark: vec![vec![R::BaseRing::zero(); 3]; 8 + 57], // Adjust to actual ark parameters
            mds: vec![vec![R::BaseRing::zero(); 3]; 3], // Adjust to actual MDS matrix parameters
            rate: 2,
            capacity: 1,
        };

        Self::new(&config)
    }
}

impl<R: OverField, CS: LatticefoldChallengeSet<R>> Transcript<R> for PoseidonTranscript<R, CS> {
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
