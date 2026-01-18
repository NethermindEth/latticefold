use ark_crypto_primitives::sponge::{
    poseidon::{PoseidonConfig, PoseidonSponge},
    CryptographicSponge,
};
use ark_ff::{Field, PrimeField};
use cyclotomic_rings::rings::GetPoseidonParams;
use latticefold::transcript::Transcript;
use stark_rings::OverField;

use crate::transcript::PoseidonTranscriptMetrics;

/// Poseidon sponge transcript operation trace (in the sponge's base prime field).
#[derive(Clone, Debug)]
pub enum PoseidonTraceOp<BF: PrimeField> {
    Absorb(Vec<BF>),
    SqueezeField(Vec<BF>),
    SqueezeBytes { n: usize, out: Vec<u8> },
}

/// Full transcript trace in terms of the Poseidon sponge's base prime field.
#[derive(Clone, Debug, Default)]
pub struct PoseidonTranscriptTrace<BF: PrimeField> {
    pub ops: Vec<PoseidonTraceOp<BF>>,
    pub absorbed: Vec<BF>,
    pub squeezed_field: Vec<BF>,
    pub squeezed_bytes: Vec<u8>,
}

/// Poseidon transcript that records a full operation trace.
///
/// This is intended for WE/DPP arithmetization frontends: the prover can record a trace and later
/// provide it as part of a witness, with constraints enforcing that it matches the Poseidon
/// permutation schedule.
#[derive(Clone)]
pub struct TracePoseidonTranscript<R: OverField> {
    sponge: PoseidonSponge<<R::BaseRing as Field>::BasePrimeField>,
    metrics: PoseidonTranscriptMetrics,
    scratch: Vec<<R::BaseRing as Field>::BasePrimeField>,
    trace: PoseidonTranscriptTrace<<R::BaseRing as Field>::BasePrimeField>,
}

impl<R: OverField> TracePoseidonTranscript<R> {
    pub fn empty<P: GetPoseidonParams<<<R>::BaseRing as Field>::BasePrimeField>>() -> Self {
        Self::new(&P::get_poseidon_config())
    }

    pub fn metrics(&self) -> PoseidonTranscriptMetrics {
        self.metrics
    }

    pub fn trace(&self) -> &PoseidonTranscriptTrace<<R::BaseRing as Field>::BasePrimeField> {
        &self.trace
    }

    fn absorb_base_prime_field_elems_vec(
        &mut self,
        elems: Vec<<R::BaseRing as Field>::BasePrimeField>,
    ) {
        self.metrics.absorbed_elems += elems.len() as u64;
        self.sponge.absorb(&elems);
        self.trace.absorbed.extend_from_slice(&elems);
        self.trace.ops.push(PoseidonTraceOp::Absorb(elems));
    }
}

impl<R: OverField> Transcript<R> for TracePoseidonTranscript<R> {
    type TranscriptConfig = PoseidonConfig<<R::BaseRing as Field>::BasePrimeField>;

    fn new(config: &Self::TranscriptConfig) -> Self {
        let sponge = PoseidonSponge::<<R::BaseRing as Field>::BasePrimeField>::new(config);
        Self {
            sponge,
            metrics: PoseidonTranscriptMetrics::default(),
            scratch: Vec::with_capacity(64),
            trace: PoseidonTranscriptTrace::default(),
        }
    }

    fn absorb(&mut self, v: &R) {
        self.scratch.clear();
        for c in v.coeffs() {
            self.scratch.extend(c.to_base_prime_field_elements());
        }
        let elems = self.scratch.clone();
        self.absorb_base_prime_field_elems_vec(elems);
    }

    fn absorb_field_element(&mut self, v: &R::BaseRing) {
        // IMPORTANT (encoding):
        // Absorb base-ring field elements directly (as base-prime-field elems), instead of
        // converting to a constant-coeff ring element and absorbing `d` coefficients.
        self.absorb_base_prime_field_elems_vec(v.to_base_prime_field_elements().collect());
    }

    fn get_challenge(&mut self) -> R::BaseRing {
        let extension_degree = R::BaseRing::extension_degree();
        let c = self
            .sponge
            .squeeze_field_elements(extension_degree as usize);
        self.metrics.squeezed_field_elems += c.len() as u64;
        self.trace.squeezed_field.extend_from_slice(&c);
        self.trace.ops.push(PoseidonTraceOp::SqueezeField(c.clone()));

        // `get_challenge` re-absorbs the squeezed elements to evolve the sponge state.
        self.absorb_base_prime_field_elems_vec(c.clone());

        <R::BaseRing as Field>::from_base_prime_field_elems(&c)
            .expect("TracePoseidonTranscript: wrong extension_degree")
    }

    fn squeeze_bytes(&mut self, n: usize) -> Vec<u8> {
        self.metrics.squeezed_bytes += n as u64;
        let out = self.sponge.squeeze_bytes(n);
        self.trace.squeezed_bytes.extend_from_slice(&out);
        self.trace
            .ops
            .push(PoseidonTraceOp::SqueezeBytes { n, out: out.clone() });
        out
    }
}

