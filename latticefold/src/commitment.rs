use thiserror::Error;

mod commitment_scheme;
mod homomorphic_commitment;
#[macro_use]
mod operations;
pub use commitment_scheme::*;
pub use homomorphic_commitment::*;

#[derive(Debug, Error)]
pub enum CommitmentError {
    #[error("Wrong length of the witness: {0}, expected: {1}")]
    WrongWitnessLength(usize, usize),
    #[error("Wrong length of the commitment: {0}, expected: {1}")]
    WrongCommitmentLength(usize, usize),
    #[error("Ajtai matrix has dimensions: {0}x{1}, expected: {2}x{3}")]
    WrongAjtaiMatrixDimensions(usize, usize, usize, usize),
}
