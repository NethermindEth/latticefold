use crate::ark_base::*;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum CSError {
    #[error("constraint system is not satisfied")]
    NotSatisfied,
    #[error("constraint system matrices rows length (m) not a power of 2: {0}")]
    MatricesRowsLengthNotPowerOf2(usize),
    #[error("constraint system has invalid size bounds: m = {0}, n = {1}, L = {2}")]
    InvalidSizeBounds(usize, usize, usize),
    #[error("vectors {0} and {1} have different lengths: {0} and {1}")]
    LengthsNotEqual(String, String, usize, usize),
}
