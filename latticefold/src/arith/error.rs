//! Provides error functionality for constraint systems.
use crate::ark_base::*;
use thiserror::Error;

/// Errors that can arise in constraint system calculations
#[derive(Debug, Error)]
pub enum CSError {
    /// The constraint system is not satisfied by the provided witness
    #[error("constraint system is not satisfied")]
    NotSatisfied,

    /// The constraint system is not of length 2.
    ///
    /// More to the point, the witness length will not be a power of 2,
    /// so we cannot use it as a MLE.
    #[error("constraint system matrices rows length (m) not a power of 2: {0}")]
    MatricesRowsLengthNotPowerOf2(usize),

    /// This error occurs if the CCS instance is not correctly padded
    ///
    /// See [definition 4.3](https://eprint.iacr.org/2024/257.pdf#page=40) of LatticeFold paper.
    #[error("constraint system has invalid size bounds: m = {0}, n = {1}, L = {2}")]
    InvalidSizeBounds(usize, usize, usize),

    /// This error occurs when performing operations on vectors of differing lengths.
    #[error("vectors {0} and {1} have different lengths: {0} and {1}")]
    LengthsNotEqual(String, String, usize, usize),
}
