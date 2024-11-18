use ark_std::fmt::Display;

/// Decomposition parameters.
/// Convenient to enforce them compile-time.
/// Contains both gadget matrix data and Latticefold decomposition step data.
pub trait DecompositionParams: Clone {
    /// The MSIS bound.
    const B: u128;
    /// The ring modulus should be < B^L.
    const L: usize;
    /// The small b from the decomposition step of LF.
    const B_SMALL: usize;
    /// K = log_b B.
    const K: usize;
}

impl<P: DecompositionParams> From<P> for DecompositionParamData {
    fn from(_: P) -> Self {
        {
            Self { b: P::B, l: P::L }
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

#[allow(non_camel_case_types)]
#[cfg(test)]
pub mod test_params {
    use super::DecompositionParams;

    #[derive(Clone)]
    pub struct PP;

    #[cfg(test)]
    impl DecompositionParams for PP {
        const B: u128 = 1024;
        const L: usize = 2;
        const B_SMALL: usize = 2;
        const K: usize = 10;
    }

    #[derive(Clone)]
    pub struct PPL1;

    #[cfg(test)]
    impl DecompositionParams for PPL1 {
        const B: u128 = 1024;
        const L: usize = 1;
        const B_SMALL: usize = 2;
        const K: usize = 10;
    }
    #[derive(Clone)]
    pub struct PP_STARK;
    impl DecompositionParams for PP_STARK {
        const B: u128 = 10485760000;
        const L: usize = 8;
        const B_SMALL: usize = 320;
        const K: usize = 4;
    }
    #[derive(Clone)]
    pub struct PP_STARK_FOLDING;
    impl DecompositionParams for PP_STARK_FOLDING {
        const B: u128 = 3010936384;
        const L: usize = 8;
        const B_SMALL: usize = 38;
        const K: usize = 6;
    }
}
