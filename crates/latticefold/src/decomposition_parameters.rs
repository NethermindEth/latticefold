//! Provides utility for decomposition parameters.
//!
//! Decomposition parameters dictate how higher-bound witness are
//! decomposed into lower-bound witnesses.

use ark_std::fmt::Display;

/// Decomposition parameters.
/// Contains both gadget matrix data and Latticefold decomposition step data.
#[derive(Clone, Debug, PartialEq)]
pub struct DecompositionParams {
    pub B: u128,
    pub l: usize,
    pub b: usize,
    pub k: usize,
}

impl DecompositionParams {
    pub fn new(B: u128, l: usize, b: usize, k: usize) -> Self {
        Self { B, l, b, k }
    }
}

impl From<DecompositionParams> for DecompositionParamData {
    fn from(dparams: DecompositionParams) -> Self {
        {
            Self {
                b: dparams.B,
                l: dparams.l,
            }
        }
    }
}

// Nice representation of parameters for printing out in benchmarks.
#[derive(Clone, Copy)]
pub struct DecompositionParamData {
    // The MSIS bound.
    b: u128,
    // The ring modulus should be < B^L.
    l: usize,
}

impl Display for DecompositionParamData {
    fn fmt(&self, f: &mut ark_std::fmt::Formatter<'_>) -> ark_std::fmt::Result {
        write!(f, "B={}, l={}", self.b, self.l,)
    }
}

#[cfg(test)]
pub mod test_params {
    use super::DecompositionParams;

    pub fn dp() -> DecompositionParams {
        DecompositionParams::new(1024, 2, 2, 10)
    }

    pub fn dp_l1() -> DecompositionParams {
        DecompositionParams::new(1024, 1, 2, 10)
    }

    pub fn dp_stark() -> DecompositionParams {
        DecompositionParams::new(10485760000, 8, 320, 4)
    }

    pub fn dp_goldilocks() -> DecompositionParams {
        DecompositionParams::new(1 << 15, 5, 2, 15)
    }

    pub fn dp_babybear() -> DecompositionParams {
        DecompositionParams::new(1 << 8, 4, 2, 8)
    }

    pub fn dp_stark_folding() -> DecompositionParams {
        DecompositionParams::new(3010936384, 8, 38, 6)
    }

    pub fn dp_frog() -> DecompositionParams {
        DecompositionParams::new(1 << 8, 8, 2, 10)
    }
}
