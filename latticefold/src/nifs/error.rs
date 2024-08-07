use lattirust_arithmetic::{polynomials::ArithErrors, ring::Ring};
use thiserror::Error;

use crate::{arith::error::CSError, commitment::CommitmentError, utils::sumcheck::SumCheckError};

#[derive(Debug, Error)]
pub enum LatticefoldError<R: Ring> {
    #[error("linearization failed: {0}")]
    LinearizationError(#[from] LinearizationError<R>),
    #[error("decomposition failed: {0}")]
    DecompositionError(#[from] DecompositionError),
    #[error("folding failed: {0}")]
    FoldingError(#[from] FoldingError<R>),
}

#[derive(Debug, Error)]
pub enum LinearizationError<R: Ring> {
    #[error("sum check failed at linearization step: {0}")]
    SumCheckError(#[from] SumCheckError<R>),
    #[error("parameters error: {0}")]
    ParametersError(String),
    #[error("constraint system related error: {0}")]
    ConstraintSystemError(#[from] CSError),
    #[error("Arithmetic error: {0}")]
    ArithmeticError(#[from] ArithErrors),
}

#[derive(Debug, Error)]
pub enum DecompositionError {
    #[error("input vectors have incorrect length")]
    IncorrectLength,
    #[error("ajtai commitment error: {0}")]
    CommitmentError(#[from] CommitmentError),
    #[error("failed to evaluate witness MLE")]
    WitnessMleEvalFail,
    #[error("constraint system related error: {0}")]
    ConstraintSystemError(#[from] CSError),
    #[error("recomposing proof checked failed")]
    RecomposedError,
}

#[derive(Debug, Error)]
pub enum FoldingError<R: Ring> {
    #[error("sum check failed at folding step: {0}")]
    SumCheckError(#[from] SumCheckError<R>),
}
