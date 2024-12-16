use stark_rings::Ring;
use stark_rings_poly::polynomials::ArithErrors;
use thiserror::Error;

use crate::{
    arith::error::CSError, ark_base::*, commitment::CommitmentError,
    utils::mle_helpers::MleEvaluationError, utils::sumcheck::SumCheckError,
};

#[derive(Debug, Error)]
pub enum LatticefoldError<R: Ring> {
    #[error("linearization failed: {0}")]
    LinearizationError(#[from] LinearizationError<R>),
    #[error("decomposition failed: {0}")]
    DecompositionError(#[from] DecompositionError),
    #[error("folding failed: {0}")]
    FoldingError(#[from] FoldingError<R>),
    #[error("constraint system related error: {0}")]
    ConstraintSystemError(#[from] CSError),
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
    #[error("mle evaluation failed: {0}")]
    EvaluationError(#[from] MleEvaluationError),
}

#[derive(Debug, Error)]
pub enum DecompositionError {
    #[error("input vectors have incorrect length")]
    IncorrectLength,
    #[error("ajtai commitment error: {0}")]
    CommitmentError(#[from] CommitmentError),
    #[error("constraint system related error: {0}")]
    ConstraintSystemError(#[from] CSError),
    #[error("recomposing proof checked failed")]
    RecomposedError,
    #[error("mle evaluation failed: {0}")]
    EvaluationError(#[from] MleEvaluationError),
}

#[derive(Debug, Error)]
pub enum FoldingError<R: Ring> {
    #[error("input vectors have incorrect length")]
    IncorrectLength,
    #[error("sum check failed at folding step: {0}")]
    SumCheckError(#[from] SumCheckError<R>),
    #[error("constraint system related error: {0}")]
    ConstraintSystemError(#[from] CSError),
    #[error("virtual polynomial error: {0}")]
    ArithError(#[from] ArithErrors),
    #[error("mle evaluation failed: {0}")]
    EvaluationError(#[from] MleEvaluationError),
    #[error("sumcheck challenge point were not generate correctly")]
    SumcheckChallengeError,
}
