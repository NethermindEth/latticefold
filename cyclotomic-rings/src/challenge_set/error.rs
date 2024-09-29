use lattirust_poly::polynomials::ArithErrors;
use lattirust_ring::Ring;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum ChallengeSetError {
    #[error("too few bytes: got {0}, expected {1}")]
    TooFewBytes(usize, usize),
}
