use lattirust_arithmetic::ring::Ring;
use thiserror::Error;

use crate::{arith::error::NotSatisfiedError, utils::sumcheck::SumCheckError};

#[derive(Debug, Error)]
pub enum LatticefoldError<R: Ring> {
    #[error("linearization failed: {0}")]
    LinearizationError(#[from] LinearizationError<R>),
    #[error("decomposition failed: {0}")]
    DecompositionError(#[from] DecompositionError<R>),
    #[error("folding failed: {0}")]
    FoldingError(#[from] FoldingError<R>),
}

#[derive(Debug, Error)]
pub enum LinearizationError<R: Ring> {
    #[error("sum check failed at linearization step: {0}")]
    SumCheckError(#[from] SumCheckError<R>),
}

#[derive(Debug, Error)]
pub enum DecompositionError<R: Ring> {
    #[error("phantom decomposition error constructor")]
    PhantomRRemoveThisLater(R),
    #[error("input vectors have incorrect length")]
    IncorrectLength,
    #[error("verification failed: {0}")]
    VerificationFailed(#[from] NotSatisfiedError),
}

#[derive(Debug, Error)]
pub enum FoldingError<R: Ring> {
    #[error("phantom folding error")]
    PhantomRRemoveThisLater(R),
}
