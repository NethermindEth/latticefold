use ark_crypto_primitives::sponge::{
    poseidon::{PoseidonConfig, PoseidonSponge},
    CryptographicSponge,
};
use ark_ff::{Field, PrimeField};
use cyclotomic_rings::rings::GetPoseidonParams;
use latticefold::transcript::Transcript;
use stark_rings::OverField;

/// Transcript using Poseidon sponge
///
/// Same implementation as LatticeFold's `PoseidonTranscript` though without the challenge set
/// generic / requirement on `SuitableRing`.
#[derive(Clone)]
pub struct PoseidonTranscript<R: OverField> {
    sponge: PoseidonSponge<<R::BaseRing as Field>::BasePrimeField>,
}

impl<R: OverField> PoseidonTranscript<R> {
    pub fn empty<P: GetPoseidonParams<<<R>::BaseRing as Field>::BasePrimeField>>() -> Self {
        Self::new(&P::get_poseidon_config())
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

pub fn squeeze_challenges<R: OverField>(
    transcript: &mut impl Transcript<R>,
    name: &str,
    n: usize,
) -> Vec<R::BaseRing> {
    transcript.absorb_field_element(&<R::BaseRing as Field>::from_base_prime_field(
        <R::BaseRing as Field>::BasePrimeField::from_be_bytes_mod_order(name.as_bytes()),
    ));

    transcript.get_challenges(n)
}

pub fn squeeze_rchallenges<R: OverField>(
    transcript: &mut impl Transcript<R>,
    name: &str,
    n: usize,
) -> Vec<R> {
    squeeze_challenges(transcript, name, n)
        .into_iter()
        .map(|z| R::from(z))
        .collect::<Vec<R>>()
}
