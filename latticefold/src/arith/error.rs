use crate::ark_base::*;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum CSError {
    #[error("constraint system is not satisfied")]
    NotSatisfied,
    #[error("vectors {0} and {1} have different lengths: {0} and {1}")]
    LengthsNotEqual(String, String, usize, usize),
    #[error("original length {0} is greater than padded length {1}")]
    IncorrectPadding(usize, usize),
}
