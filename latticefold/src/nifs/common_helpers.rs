//!
//! Helper function used by all three subprotocols.
//!  

use ark_std::cfg_into_iter;

use lattirust_poly::mle::DenseMultilinearExtension;
use lattirust_ring::Ring;

use super::error::MleEvaluationError;

#[cfg(feature = "parallel")]
use rayon::iter::{IntoParallelIterator, ParallelIterator};

pub trait Evaluate<R: Ring> {
    fn evaluate(self, point: &[R]) -> Result<R, MleEvaluationError>;
}

impl<R: Ring> Evaluate<R> for Vec<R> {
    fn evaluate(self, point: &[R]) -> Result<R, MleEvaluationError> {
        let evals_len = self.len();

        DenseMultilinearExtension::from_evaluations_vec(point.len(), self)
            .evaluate(point)
            .ok_or(MleEvaluationError::IncorrectLength(point.len(), evals_len))
    }
}

impl<'a, R: Ring> Evaluate<R> for &'a [R] {
    fn evaluate(self, point: &[R]) -> Result<R, MleEvaluationError> {
        let evals_len = self.len();

        DenseMultilinearExtension::from_evaluations_slice(point.len(), self)
            .evaluate(point)
            .ok_or(MleEvaluationError::IncorrectLength(point.len(), evals_len))
    }
}

impl<'a, R: Ring> Evaluate<R> for &'a DenseMultilinearExtension<R> {
    fn evaluate(self, point: &[R]) -> Result<R, MleEvaluationError> {
        DenseMultilinearExtension::<R>::evaluate(self, point).ok_or(
            MleEvaluationError::IncorrectLength(point.len(), self.evaluations.len()),
        )
    }
}

#[cfg(not(feature = "parallel"))]
pub fn evaluate_mles<R, V, I, E>(mle_s: I, point: &[R]) -> Result<Vec<R>, E>
where
    R: Ring,
    V: Evaluate<R>,
    I: IntoIterator<Item = V>,
    E: From<MleEvaluationError>,
{
    cfg_into_iter!(mle_s)
        .map(|evals| evals.evaluate(point).map_err(From::from))
        .collect()
}

#[cfg(feature = "parallel")]
pub fn evaluate_mles<R, V, I, E>(mle_s: I, point: &[R]) -> Result<Vec<R>, E>
where
    R: Ring,
    V: Evaluate<R>,
    I: IntoParallelIterator<Item = V>,
    E: From<MleEvaluationError> + Send + Sync,
{
    cfg_into_iter!(mle_s)
        .map(|evals| evals.evaluate(point).map_err(From::from))
        .collect()
}

#[cfg(not(feature = "parallel"))]
pub fn to_mles<'a, I, R, E>(n_vars: usize, mle_s: I) -> Result<Vec<DenseMultilinearExtension<R>>, E>
where
    I: Iterator<Item = &'a Vec<R>>,
    R: Ring,
    E: From<MleEvaluationError> + Sync + Send,
{
    mle_s
        .into_iter()
        .map(|M| Ok(DenseMultilinearExtension::from_slice(n_vars, &M)))
        .collect::<Result<_, E>>()
}

#[cfg(not(feature = "parallel"))]
pub fn to_mles_owned<I, R, E>(
    n_vars: usize,
    mle_s: I,
) -> Result<Vec<DenseMultilinearExtension<R>>, E>
where
    I: IntoIterator<Item = Vec<R>>,
    R: Ring,
    E: From<MleEvaluationError> + Sync + Send,
{
    mle_s
        .into_iter()
        .map(|M| Ok(DenseMultilinearExtension::from_slice(n_vars, &M)))
        .collect::<Result<_, E>>()
}

#[cfg(not(feature = "parallel"))]
pub fn to_mles_err<I, R, E, E1>(
    n_vars: usize,
    mle_s: I,
) -> Result<Vec<DenseMultilinearExtension<R>>, E>
where
    I: IntoIterator<Item = Result<Vec<R>, E1>>,
    R: Ring,
    E: From<MleEvaluationError>,
    E1: Into<E>,
{
    mle_s
        .into_iter()
        .map(|M| match M {
            Ok(M) => Ok(DenseMultilinearExtension::from_slice(n_vars, &M)),
            Err(M) => Err(M.into()),
        })
        .collect::<Result<_, E>>()
}

#[cfg(feature = "parallel")]
pub fn to_mles<'a, I, R, E>(n_vars: usize, mle_s: I) -> Result<Vec<DenseMultilinearExtension<R>>, E>
where
    I: Iterator<Item = &'a Vec<R>>,
    R: Ring,
    E: From<MleEvaluationError> + Sync + Send,
{
    mle_s
        .into_par_iter()
        .map(|M| Ok(DenseMultilinearExtension::from_slice(n_vars, &M)))
        .collect::<Result<_, E>>()
}
#[cfg(feature = "parallel")]
pub fn to_mles_owned<I, R, E>(
    num_vars: usize,
    mle_s: I,
) -> Result<Vec<DenseMultilinearExtension<R>>, E>
where
    I: IntoParallelIterator<Item = Vec<R>>,
    R: Ring,
    E: From<MleEvaluationError> + Sync + Send,
{
    mle_s
        .into_par_iter()
        .map(|M| Ok(DenseMultilinearExtension::from_slice(n_vars, &M)))
        .collect::<Result<_, E>>()
}

#[cfg(feature = "parallel")]
pub fn to_mles_err<I, R, E, E1>(
    n_vars: usize,
    mle_s: I,
) -> Result<Vec<DenseMultilinearExtension<R>>, E>
where
    I: IntoParallelIterator<Item = Result<Vec<R>, E1>>,
    R: Ring,
    E: From<MleEvaluationError>,
    E1: Into<E>,
{
    mle_s
        .into_par_iter()
        .map(|M| match M {
            Ok(M) => Ok(DenseMultilinearExtension::from_slice(n_vars, &M)),
            Err(M) => Err(M.into()),
        })
        .collect::<Result<_, E>>()
}
