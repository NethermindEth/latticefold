use crate::ark_base::*;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum CSError {
    #[error("constraint system is not satisfied")]
    NotSatisfied,
    #[error("vectors {0} and {1} have different lengths: {0} and {1}")]
    LengthsNotEqual(String, String, usize, usize),
}
