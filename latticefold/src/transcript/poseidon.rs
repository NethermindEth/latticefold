use ark_crypto_primitives::sponge::{
    poseidon::{PoseidonConfig, PoseidonSponge},
    Absorb, CryptographicSponge,
};
use ark_ff::{BigInteger, PrimeField};
use lattirust_arithmetic::{
    challenge_set::latticefold_challenge_set::{LatticefoldChallengeSet, OverField},
    ring::UnsignedRepresentative,
};

use super::Transcript;

/// PoseidonTranscript implements the Transcript trait using the Poseidon hash
pub struct PoseidonTranscript<F: PrimeField, R: OverField<F>, CS: LatticefoldChallengeSet<F, R>>
where
    F: Absorb,
{
    _marker: std::marker::PhantomData<(R, CS)>,
    sponge: PoseidonSponge<F>,
}

impl<F: PrimeField, R: OverField<F>, CS: LatticefoldChallengeSet<F, R>> Transcript<F, R>
    for PoseidonTranscript<F, R, CS>
where
    F: Absorb,
{
    type TranscriptConfig = PoseidonConfig<F>;

    type ChallengeSet = CS;

    fn new(config: &Self::TranscriptConfig) -> Self {
        let sponge = PoseidonSponge::<F>::new(config);
        Self {
            sponge,
            _marker: std::marker::PhantomData,
        }
    }

    fn absorb(&mut self, v: &F) {
        self.sponge.absorb(&v);
    }

    fn absorb_vec(&mut self, v: &[F]) {
        self.sponge.absorb(&v);
    }

    fn absorb_ring(&mut self, v: &R) {
        self.sponge.absorb(&ring_to_field(v));
    }

    fn absorb_ring_vec(&mut self, v: &[R]) {
        for ring in v {
            self.absorb_ring(ring);
        }
    }

    fn get_big_challenge(&mut self) -> <R>::BaseRing {
        let c: Vec<F> = self.sponge.squeeze_field_elements(1);
        self.sponge.absorb(&c);
        Self::ChallengeSet::big_challenge_from_field(&c[0])
    }

    fn get_small_challenge(&mut self) -> R {
        let c: Vec<F> = self.sponge.squeeze_field_elements(1);
        self.sponge.absorb(&c);
        Self::ChallengeSet::small_challenge(&c[0].into_bigint().to_bytes_be())
    }
}

fn ring_to_field<F: PrimeField, R: OverField<F>>(x: &R) -> Vec<F> {
    x.coeffs()
        .into_iter()
        .map(|coeff| F::from(<<R>::BaseRing as Into<UnsignedRepresentative>>::into(coeff).0))
        .collect()
}
