use lattirust_arithmetic::{polynomials::ArithErrors, ring::Ring};
use thiserror::Error;

use crate::{arith::error::CSError, utils::sumcheck::SumCheckError};

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
    #[error("parameters error: {0}")]
    ParametersError(String),
    #[error("constraint system related error: {0}")]
    ConstraintSystemError(#[from] CSError),
    #[error("Arithmetic error: {0}")]
    ArithmeticError(#[from] ArithErrors),
}

#[derive(Debug, Error)]
pub enum DecompositionError<R: Ring> {
    #[error("phantom decomposition error constructor")]
    PhantomRRemoveThisLater(R),
    #[error("input vectors have incorrect length")]
    IncorrectLength,
}

#[derive(Debug, Error)]
pub enum FoldingError<R: Ring> {
    #[error("phantom folding error")]
    PhantomRRemoveThisLater(R),
}
