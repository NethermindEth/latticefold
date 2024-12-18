//!
//! Helper function used by all three subprotocols.
//!  

use ark_std::{cfg_into_iter, cfg_iter, vec::Vec};
use thiserror::Error;

use stark_rings::Ring;
use stark_rings_poly::mle::DenseMultilinearExtension;

use crate::arith::{error::CSError, utils::mat_vec_mul, CCS};
use cyclotomic_rings::rings::SuitableRing;

#[cfg(feature = "parallel")]
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};

#[derive(Debug, Error)]
pub enum MleEvaluationError {
    #[error("lengths of evaluation point and evaluations are not consistent: 1 << {0} != {1}")]
    IncorrectLength(usize, usize),
}

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

impl<R: Ring> Evaluate<R> for &[R] {
    fn evaluate(self, point: &[R]) -> Result<R, MleEvaluationError> {
        let evals_len = self.len();

        DenseMultilinearExtension::from_evaluations_slice(point.len(), self)
            .evaluate(point)
            .ok_or(MleEvaluationError::IncorrectLength(point.len(), evals_len))
    }
}

impl<R: Ring> Evaluate<R> for &DenseMultilinearExtension<R> {
    fn evaluate(self, point: &[R]) -> Result<R, MleEvaluationError> {
        DenseMultilinearExtension::<R>::evaluate(self, point)
            .ok_or(MleEvaluationError::IncorrectLength(point.len(), self.elen))
    }
}

impl<R: Ring> Evaluate<R> for &Vec<R> {
    fn evaluate(self, point: &[R]) -> Result<R, MleEvaluationError> {
        if self.len() != 1 << point.len() {
            return Err(MleEvaluationError::IncorrectLength(point.len(), self.len()));
        }

        DenseMultilinearExtension::from_evaluations_slice(point.len(), self)
            .evaluate(point)
            .ok_or(MleEvaluationError::IncorrectLength(point.len(), self.len()))
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
            if 1 << n_vars < m.len() {
                Err(MleEvaluationError::IncorrectLength(1 << n_vars, m.len()).into())
            } else {
                Ok(DenseMultilinearExtension::from_evaluations_vec(n_vars, m))
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
            if 1 << n_vars < m.len() {
                Err(MleEvaluationError::IncorrectLength(1 << n_vars, m.len()).into())
            } else {
                Ok(DenseMultilinearExtension::from_evaluations_vec(n_vars, m))
            }
        })
        .collect::<Result<_, E>>()
}

// Prepare MLE's of the form mle[M_i \cdot z_ccs](x), a.k.a. \sum mle[M_i](x, b) * mle[z_ccs](b).
pub fn calculate_Mz_mles<NTT, E>(
    ccs: &CCS<NTT>,
    z_ccs: &[NTT],
) -> Result<Vec<DenseMultilinearExtension<NTT>>, E>
where
    NTT: SuitableRing,
    E: From<MleEvaluationError> + From<CSError> + Sync + Send,
{
    to_mles_err::<_, _, E, CSError>(ccs.s, cfg_iter!(ccs.M).map(|M| mat_vec_mul(M, z_ccs)))
}
