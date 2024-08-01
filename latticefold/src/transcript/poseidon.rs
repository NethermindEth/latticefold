use super::Transcript;
use ark_crypto_primitives::sponge::{
    poseidon::{PoseidonConfig, PoseidonSponge},
    CryptographicSponge,
};
use ark_ff::BigInteger;
use ark_ff::PrimeField;
use ark_ff::Zero;
use lattirust_arithmetic::{
    challenge_set::latticefold_challenge_set::{LatticefoldChallengeSet, OverField},
    ring::UnsignedRepresentative,
};

/// PoseidonTranscript implements the Transcript trait using the Poseidon hash
pub struct PoseidonTranscript<R: OverField, CS: LatticefoldChallengeSet<R>> {
    _marker: std::marker::PhantomData<(R, CS)>,
    sponge: PoseidonSponge<R::F>,
}

impl<R: OverField, CS: LatticefoldChallengeSet<R>> Default for PoseidonTranscript<R, CS> {
    fn default() -> Self {
        let config = PoseidonConfig {
            full_rounds: 8, // Example values, adjust according to your needs
            partial_rounds: 57,
            alpha: 5,
            ark: vec![vec![R::F::zero(); 3]; 8 + 57], // Adjust to actual ark parameters
            mds: vec![vec![R::F::zero(); 3]; 3],      // Adjust to actual MDS matrix parameters
            rate: 2,
            capacity: 1,
        };

        Self::new(&config)
    }
}

impl<R: OverField, CS: LatticefoldChallengeSet<R>> Transcript<R> for PoseidonTranscript<R, CS> {
    type TranscriptConfig = PoseidonConfig<R::F>;

    type ChallengeSet = CS;

    fn new(config: &Self::TranscriptConfig) -> Self {
        let sponge = PoseidonSponge::<R::F>::new(config);
        Self {
            sponge,
            _marker: std::marker::PhantomData,
        }
    }

    fn absorb(&mut self, v: &R::F) {
        self.sponge.absorb(&v);
    }

    fn absorb_vec(&mut self, v: &[R::F]) {
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
        let c: Vec<R::F> = self.sponge.squeeze_field_elements(1);
        self.sponge.absorb(&c);
        Self::ChallengeSet::big_challenge_from_field(&c[0])
    }

    fn get_small_challenge(&mut self) -> R {
        let c: Vec<R::F> = self.sponge.squeeze_field_elements(1);
        self.sponge.absorb(&c);
        Self::ChallengeSet::small_challenge(&c[0].into_bigint().to_bytes_be())
    }

    fn get_small_challenges(&mut self, n: usize) -> Vec<R> {
        let mut challenges: Vec<R> = Vec::with_capacity(n);
        (0..n).for_each(|_| {
            challenges.push(self.get_small_challenge());
        });
        challenges
    }
    fn get_big_challenges(&mut self, n: usize) -> Vec<R::BaseRing> {
        let mut challenges: Vec<R::BaseRing> = Vec::with_capacity(n);
        (0..n).for_each(|_| {
            challenges.push(self.get_big_challenge());
        });
        challenges
    }
}

fn ring_to_field<R: OverField>(x: &R) -> Vec<R::F> {
    x.coeffs()
        .into_iter()
        .map(|coeff| R::F::from(<<R>::BaseRing as Into<UnsignedRepresentative>>::into(coeff).0))
        .collect()
}
