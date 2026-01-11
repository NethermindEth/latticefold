use ark_crypto_primitives::sponge::{
    poseidon::{PoseidonConfig, PoseidonSponge},
    CryptographicSponge,
};
use ark_ff::{Field, PrimeField};
use cyclotomic_rings::rings::GetPoseidonParams;
use latticefold::transcript::Transcript;
use stark_rings::OverField;

/// Lightweight transcript metrics to estimate Poseidon sponge work.
///
/// Counts are in units of the Poseidon sponge's **base prime field** elements.
#[derive(Clone, Copy, Debug, Default)]
pub struct PoseidonTranscriptMetrics {
    /// Number of base-prime-field elements absorbed into the sponge.
    pub absorbed_elems: u64,
    /// Number of base-prime-field elements squeezed as challenges (`get_challenge`).
    pub squeezed_field_elems: u64,
    /// Number of bytes squeezed via `squeeze_bytes`.
    pub squeezed_bytes: u64,
}

impl PoseidonTranscriptMetrics {
    /// Estimate number of Poseidon permutations based on sponge rate.
    /// 
    /// For a rate-r sponge:
    /// - Each permutation absorbs r field elements
    /// - Each permutation produces r field elements (squeeze)
    /// 
    /// This is a lower bound; actual count depends on interleaving.
    /// 
    /// `bytes_per_field_elem`: For Goldilocks (~64-bit), use 8. For larger fields, adjust.
    pub fn estimated_permutations(&self, rate: usize, bytes_per_field_elem: usize) -> u64 {
        let absorb_perms = (self.absorbed_elems + rate as u64 - 1) / rate as u64;
        let squeeze_perms = (self.squeezed_field_elems + rate as u64 - 1) / rate as u64;
        // Convert squeezed bytes to field elements
        let byte_field_elems = (self.squeezed_bytes + bytes_per_field_elem as u64 - 1) 
            / bytes_per_field_elem as u64;
        let byte_perms = (byte_field_elems + rate as u64 - 1) / rate as u64;
        absorb_perms + squeeze_perms + byte_perms
    }
}

/// Transcript using Poseidon sponge with metrics tracking.
///
/// Same implementation as LatticeFold's `PoseidonTranscript` though without the challenge set
/// generic / requirement on `SuitableRing`.
#[derive(Clone)]
pub struct PoseidonTranscript<R: OverField> {
    sponge: PoseidonSponge<<R::BaseRing as Field>::BasePrimeField>,
    metrics: PoseidonTranscriptMetrics,
}

impl<R: OverField> PoseidonTranscript<R> {
    pub fn empty<P: GetPoseidonParams<<<R>::BaseRing as Field>::BasePrimeField>>() -> Self {
        Self::new(&P::get_poseidon_config())
    }

    /// Get current transcript metrics.
    pub fn metrics(&self) -> PoseidonTranscriptMetrics {
        self.metrics
    }

    /// Print a summary of transcript work.
    /// 
    /// Uses Goldilocks-like parameters (64-bit field, rate=8 for Poseidon2).
    pub fn print_metrics(&self) {
        let m = &self.metrics;
        println!("=== LF+ Transcript Metrics ===");
        println!("  Absorbed base-field elems: {}", m.absorbed_elems);
        println!("  Squeezed field elems:      {}", m.squeezed_field_elems);
        println!("  Squeezed bytes:            {}", m.squeezed_bytes);
        // Goldilocks-style: rate=8 for Poseidon2, 8 bytes per field elem
        let perms_p2 = m.estimated_permutations(8, 8);
        // Conservative: rate=2 for standard Poseidon
        let perms_p1 = m.estimated_permutations(2, 8);
        println!("  Est. Poseidon2 permutations (rate=8): {}", perms_p2);
        println!("  Est. Poseidon permutations (rate=2):  {}", perms_p1);
        println!("==============================");
    }
}

impl<R: OverField> Transcript<R> for PoseidonTranscript<R> {
    type TranscriptConfig = PoseidonConfig<<R::BaseRing as Field>::BasePrimeField>;

    fn new(config: &Self::TranscriptConfig) -> Self {
        let sponge = PoseidonSponge::<<R::BaseRing as Field>::BasePrimeField>::new(config);
        Self { 
            sponge,
            metrics: PoseidonTranscriptMetrics::default(),
        }
    }

    fn absorb(&mut self, v: &R) {
        let elems: Vec<_> = v.coeffs()
            .iter()
            .flat_map(|x| x.to_base_prime_field_elements())
            .collect();
        self.metrics.absorbed_elems += elems.len() as u64;
        self.sponge.absorb(&elems);
    }

    fn absorb_field_element(&mut self, v: &R::BaseRing) {
        // IMPORTANT (encoding):
        // Absorb base-ring field elements directly into the sponge (as base-prime-field elems),
        // instead of converting to a constant-coeff ring element and absorbing `d` coefficients.
        //
        // This reduces transcript IO (and thus WE-gate Poseidon constraints) by ~`R::dimension()`
        // for scalar absorbs, while keeping ring-element absorbs unchanged.
        let elems: Vec<_> = v.to_base_prime_field_elements().collect();
        self.metrics.absorbed_elems += elems.len() as u64;
        self.sponge.absorb(&elems);
    }

    fn get_challenge(&mut self) -> R::BaseRing {
        let extension_degree = R::BaseRing::extension_degree();
        let c = self
            .sponge
            .squeeze_field_elements(extension_degree as usize);
        self.metrics.squeezed_field_elems += extension_degree as u64;
        // Re-absorb squeezed elements (Fiat-Shamir)
        self.metrics.absorbed_elems += c.len() as u64;
        self.sponge.absorb(&c);
        <R::BaseRing as Field>::from_base_prime_field_elems(&c)
            .expect("something went wrong: c does not contain extension_degree elements")
    }

    fn squeeze_bytes(&mut self, n: usize) -> Vec<u8> {
        self.metrics.squeezed_bytes += n as u64;
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
