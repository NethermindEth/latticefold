use thiserror::Error;

use crate::ark_base::*;

/// Short challenge generation error.
#[derive(Debug, Error)]
pub enum ChallengeSetError {
    /// An error meaning there is not enough bytes to generate
    /// a short challenge.
    #[error("too few bytes: got {0}, expected {1}")]
    TooFewBytes(usize, usize),
}
