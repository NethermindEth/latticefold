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
