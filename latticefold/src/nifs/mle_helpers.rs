//!
//! Helper function used by all three subprotocols.
//!  

use ark_std::{cfg_into_iter, vec::Vec};

use lattirust_poly::mle::DenseMultilinearExtension;
use lattirust_ring::Ring;

use super::error::MleEvaluationError;
use crate::arith::{error::CSError, utils::mat_vec_mul_with_padding, CCS};

#[cfg(feature = "parallel")]
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};

pub trait Evaluate<R: Ring> {
    fn evaluate(self, point: &[R]) -> Result<R, MleEvaluationError>;
}

impl<R: Ring> Evaluate<R> for Vec<R> {
    fn evaluate(self, point: &[R]) -> Result<R, MleEvaluationError> {
        if self.len() != 1 << point.len() {
            return Err(MleEvaluationError::IncorrectLength(point.len(), self.len()));
        }

        DenseMultilinearExtension::from_evaluations_vec(point.len(), self).evaluate(point)
    }
}

impl<R: Ring> Evaluate<R> for &[R] {
    fn evaluate(self, point: &[R]) -> Result<R, MleEvaluationError> {
        if self.len() != 1 << point.len() {
            return Err(MleEvaluationError::IncorrectLength(point.len(), self.len()));
        }

        DenseMultilinearExtension::from_evaluations_slice(point.len(), self).evaluate(point)
    }
}

impl<R: Ring> Evaluate<R> for &Vec<R> {
    fn evaluate(self, point: &[R]) -> Result<R, MleEvaluationError> {
        if self.len() != 1 << point.len() {
            return Err(MleEvaluationError::IncorrectLength(point.len(), self.len()));
        }

        DenseMultilinearExtension::from_evaluations_slice(point.len(), self).evaluate(point)
    }
}

impl<R: Ring> Evaluate<R> for &DenseMultilinearExtension<R> {
    fn evaluate(self, point: &[R]) -> Result<R, MleEvaluationError> {
        DenseMultilinearExtension::<R>::evaluate(self, point).ok_or(
            MleEvaluationError::IncorrectLength(point.len(), self.evaluations.len()),
        )
    }
}

impl<R: Ring> Evaluate<R> for DenseMultilinearExtension<R> {
    fn evaluate(self, point: &[R]) -> Result<R, MleEvaluationError> {
        DenseMultilinearExtension::<R>::evaluate(&self, point).ok_or(
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
    I: IntoIterator<Item = &'a Vec<R>>,
    R: Ring,
    E: From<MleEvaluationError> + Sync + Send,
{
    mle_s
        .into_iter()
        .map(|m| if 1 << n_vars != m.len() {
            Err(MleEvaluationError::IncorrectLength(1 << n_vars, m.len()).into())
        } else {
            Ok(DenseMultilinearExtension::from_slice(n_vars, &m))
        })
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
    E: From<MleEvaluationError> + From<E1>,
{
    mle_s
        .into_iter()
        .map(|m| {
            let m = m?;
            if 1 << n_vars != m.len() {
                Err(MleEvaluationError::IncorrectLength(1 << n_vars, m.len()).into())
            } else {
                Ok(DenseMultilinearExtension::from_slice(n_vars, &m))
            }
        })
        .collect::<Result<_, E>>()
}

#[cfg(feature = "parallel")]
pub fn to_mles<'a, I, R, E>(n_vars: usize, mle_s: I) -> Result<Vec<DenseMultilinearExtension<R>>, E>
where
    I: IntoParallelIterator<Item = &'a Vec<R>>,
    R: Ring,
    E: From<MleEvaluationError> + Sync + Send,
{
    mle_s
        .into_par_iter()
        .map(|m| {
            if 1 << n_vars != m.len() {
                Err(MleEvaluationError::IncorrectLength(1 << n_vars, m.len()).into())
            } else {
                Ok(DenseMultilinearExtension::from_slice(n_vars, &m))
            }
        })
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
    E: From<MleEvaluationError> + Sync + Send + From<E1>,
{
    mle_s
        .into_par_iter()
        .map(|m| {
            let m = m?;
            if 1 << n_vars != m.len() {
                Err(MleEvaluationError::IncorrectLength(1 << n_vars, m.len()).into())
            } else {
                Ok(DenseMultilinearExtension::from_slice(n_vars, &m))
            }
        })
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
        .map(|m| if 1 << n_vars != m.len() {
            Err(MleEvaluationError::IncorrectLength(1 << n_vars, m.len()).into())
        } else {
            Ok(DenseMultilinearExtension::from_evaluations_vec(n_vars, m))
        })
        .collect::<Result<_, E>>()
}

#[cfg(feature = "parallel")]
pub fn to_mles_owned<I, R, E>(
    n_vars: usize,
    mle_s: I,
) -> Result<Vec<DenseMultilinearExtension<R>>, E>
where
    I: IntoParallelIterator<Item = Vec<R>>,
    R: Ring,
    E: From<MleEvaluationError> + Sync + Send,
{
    mle_s
        .into_par_iter()
        .map(|M| if 1 << n_vars != m.len() {
            Err(MleEvaluationError::IncorrectLength(1 << n_vars, m.len()).into())
        } else {
            Ok(DenseMultilinearExtension::from_evaluations_vec(n_vars, m))
        })
        .collect::<Result<_, E>>()
}


pub fn calculate_Mz_mles<R: Ring, E: From<CSError> + From<MleEvaluationError> + Send + Sync>(
    ccs: &CCS<R>,
    z_ccs: &[R],
) -> Result<Vec<DenseMultilinearExtension<R>>, E> {
    to_mles_err::<_, _, E, _>(
        ccs.s,
        cfg_iter!(ccs.M)
            .map(|M| mat_vec_mul_with_padding(M, z_ccs, ccs.padded_length_of_decomposed_witness)),
    )
}

pub fn prepare_f_hat_mle<R: Ring, E: From<MleEvaluationError> + Send + Sync>(
    ccs: &CCS<R>,
    f_hat: &[Vec<R>],
) -> Result<Vec<DenseMultilinearExtension<R>>, E> {
    to_mles_owned(
        ccs.padded_length_of_decomposed_witness_log,
        cfg_iter!(f_hat).map(|f_hat_row| {
            let mut f_hat_row_cloned: Vec<R> =
                Vec::with_capacity(ccs.padded_length_of_decomposed_witness);

            f_hat_row_cloned.extend_from_slice(f_hat_row);
            f_hat_row_cloned.resize(ccs.padded_length_of_decomposed_witness, R::zero());

            f_hat_row_cloned
        }),
    )
}

pub fn calculate_f_hat_mle_evaluations<R: Ring, E: From<MleEvaluationError> + Send + Sync>(
    ccs: &CCS<R>,
    f_hat: &[Vec<R>],
    point_r: &[R],
) -> Result<Vec<R>, E> {
    evaluate_mles(
        cfg_iter!(f_hat).map(|f_hat_row| {
            let mut f_hat_row_cloned: Vec<R> =
                Vec::with_capacity(ccs.padded_length_of_decomposed_witness);

            f_hat_row_cloned.extend_from_slice(f_hat_row);
            f_hat_row_cloned.resize(ccs.padded_length_of_decomposed_witness, R::zero());

            f_hat_row_cloned
        }),
        point_r,
    )
}
