#![allow(non_snake_case, clippy::upper_case_acronyms)]
use crate::{
    arith::{error::CSError, utils::mat_vec_mul, Witness, CCS, LCCCS},
    ark_base::*,
    commitment::{AjtaiCommitmentScheme, Commitment, CommitmentError},
    decomposition_parameters::DecompositionParams,
    nifs::{
        error::DecompositionError,
        mle_helpers::{evaluate_mles, to_mles_err},
    },
    transcript::Transcript,
};
use cyclotomic_rings::rings::SuitableRing;
use lattirust_linear_algebra::SparseMatrix;
use lattirust_poly::polynomials::DenseMultilinearExtension;

use lattirust_ring::OverField;
use utils::{decompose_B_vec_into_k_vec, decompose_big_vec_into_k_vec_and_compose_back};

use ark_std::{cfg_into_iter, cfg_iter};
#[cfg(feature = "parallel")]
use rayon::prelude::*;

pub use structs::*;

use super::mle_helpers::to_mles;
mod structs;
#[cfg(test)]
mod tests;
mod utils;

impl<NTT: SuitableRing, T: Transcript<NTT>> LFDecompositionProver<NTT, T> {
    fn decompose_witness<P: DecompositionParams>(wit: &Witness<NTT>) -> Vec<Witness<NTT>> {
        let f_s = decompose_B_vec_into_k_vec::<NTT, P>(&wit.f_coeff);
        cfg_into_iter!(f_s)
            .map(|f| Witness::from_f_coeff::<P>(f))
            .collect()
    }

    fn compute_x_s<P: DecompositionParams>(mut x_w: Vec<NTT>, h: NTT) -> Vec<Vec<NTT>> {
        x_w.push(h);
        decompose_big_vec_into_k_vec_and_compose_back::<NTT, P>(x_w)
    }

    fn commit_from_witnesses<const C: usize, const W: usize, P: DecompositionParams>(
        wit_s: &Vec<Witness<NTT>>,
        scheme: &AjtaiCommitmentScheme<C, W, NTT>,
    ) -> Result<Vec<Commitment<C, NTT>>, CommitmentError> {
        cfg_iter!(wit_s)
            .map(|wit| wit.commit::<C, W, P>(scheme))
            .collect::<Result<Vec<_>, _>>()
    }

    fn compute_v_s(
        wit_s: &Vec<Witness<NTT>>,
        mle_length: usize,
        point_r: &[NTT],
    ) -> Result<Vec<Vec<NTT>>, DecompositionError> {
        cfg_iter!(wit_s)
            .map(|wit| {
                evaluate_mles::<NTT, _, _, DecompositionError>(
                    &to_mles::<_, _, DecompositionError>(mle_length, &wit.f_hat)?,
                    point_r,
                )
            })
            .collect::<Result<Vec<_>, _>>()
    }

    fn compute_u_s(
        wit_s: &Vec<Witness<NTT>>,
        M: &[SparseMatrix<NTT>],
        decomposed_statements: &[Vec<NTT>],
        point_r: &[NTT],
        num_mle_vars: usize,
    ) -> Result<Vec<Vec<NTT>>, DecompositionError> {
        cfg_iter!(wit_s)
            .enumerate()
            .map(|(i, wit)| {
                let z: Vec<_> = {
                    let mut z =
                        Vec::with_capacity(decomposed_statements[i].len() + wit.w_ccs.len());

                    z.extend_from_slice(&decomposed_statements[i]);
                    z.extend_from_slice(&wit.w_ccs);

                    z
                };

                let mles = to_mles_err::<_, _, DecompositionError, _>(
                    num_mle_vars,
                    cfg_iter!(M).map(|M| mat_vec_mul(M, &z)),
                )?;

                let u_s_for_i =
                    evaluate_mles::<NTT, &DenseMultilinearExtension<_>, _, DecompositionError>(
                        &mles, point_r,
                    )?;

                Ok(u_s_for_i)
            })
            .collect::<Result<Vec<Vec<NTT>>, DecompositionError>>()
    }
}
impl<NTT: SuitableRing, T: Transcript<NTT>> DecompositionProver<NTT, T>
    for LFDecompositionProver<NTT, T>
{
    fn prove<const W: usize, const C: usize, P: DecompositionParams>(
        cm_i: &LCCCS<C, NTT>,
        wit: &Witness<NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
        scheme: &AjtaiCommitmentScheme<C, W, NTT>,
    ) -> Result<
        (
            Vec<LCCCS<C, NTT>>,
            Vec<Witness<NTT>>,
            DecompositionProof<C, NTT>,
        ),
        DecompositionError,
    > {
        sanity_check::<NTT, P>(ccs)?;
        let log_m = ccs.s;

        let wit_s: Vec<Witness<NTT>> = Self::decompose_witness::<P>(wit);

        let x_s = Self::compute_x_s::<P>(cm_i.x_w.clone(), cm_i.h);

        let y_s: Vec<Commitment<C, NTT>> = Self::commit_from_witnesses::<C, W, P>(&wit_s, scheme)?;

        let v_s: Vec<Vec<NTT>> = Self::compute_v_s(&wit_s, log_m, &cm_i.r)?;

        let u_s = Self::compute_u_s(&wit_s, &ccs.M, &x_s, &cm_i.r, log_m)?;

        let mut lcccs_s = Vec::with_capacity(P::K);

        for (((x, y), u), v) in x_s.iter().zip(&y_s).zip(&u_s).zip(&v_s) {
            transcript.absorb_slice(x);
            transcript.absorb_slice(y.as_ref());
            transcript.absorb_slice(u);
            transcript.absorb_slice(v);

            let h = x
                .last()
                .cloned()
                .ok_or(DecompositionError::IncorrectLength)?;
            lcccs_s.push(LCCCS {
                r: cm_i.r.clone(),
                v: v.clone(),
                cm: y.clone(),
                u: u.clone(),
                x_w: x[0..x.len() - 1].to_vec(),
                h,
            })
        }

        let proof = DecompositionProof { u_s, v_s, x_s, y_s };

        Ok((lcccs_s, wit_s, proof))
    }
}

impl<NTT: OverField, T: Transcript<NTT>> LFDecompositionVerifier<NTT, T> {
    pub fn recompose_commitment<const C: usize>(
        y_s: &[Commitment<C, NTT>],
        coeffs: &[NTT],
    ) -> Result<Commitment<C, NTT>, DecompositionError> {
        y_s.iter()
            .zip(coeffs)
            .map(|(y_i, b_i)| y_i * b_i)
            .reduce(|acc, bi_part| acc + bi_part)
            .ok_or(DecompositionError::RecomposedError)
    }

    pub fn recompose_u(u_s: &[Vec<NTT>], coeffs: &[NTT]) -> Result<Vec<NTT>, DecompositionError> {
        u_s.iter()
            .zip(coeffs)
            .map(|(u_i, b_i)| u_i.iter().map(|&u| u * b_i).collect())
            .reduce(|acc, u_i_times_b_i: Vec<NTT>| {
                acc.into_iter()
                    .zip(&u_i_times_b_i)
                    .map(|(u0, ui)| u0 + ui)
                    .collect()
            })
            .ok_or(DecompositionError::RecomposedError)
    }

    pub fn recompose_v(v_s: &[Vec<NTT>], coeffs: &[NTT], row: usize) -> NTT {
        v_s.iter()
            .zip(coeffs)
            .map(|(v_i, b_i)| v_i[row] * b_i)
            .sum()
    }

    pub fn recompose_xw_and_h(
        x_s: &[Vec<NTT>],
        coeffs: &[NTT],
    ) -> Result<(Vec<NTT>, NTT), DecompositionError> {
        let mut should_equal_xw = x_s
            .iter()
            .zip(coeffs)
            .map(|(x_i, b_i)| x_i.iter().map(|&x| x * b_i).collect())
            .reduce(|acc, x_i_times_b_i: Vec<NTT>| {
                acc.into_iter()
                    .zip(&x_i_times_b_i)
                    .map(|(x0, xi)| x0 + xi)
                    .collect()
            })
            .ok_or(DecompositionError::RecomposedError)?;

        let should_equal_h = should_equal_xw
            .pop()
            .ok_or(DecompositionError::RecomposedError)?;

        Ok((should_equal_xw, should_equal_h))
    }
}

impl<NTT: OverField, T: Transcript<NTT>> DecompositionVerifier<NTT, T>
    for LFDecompositionVerifier<NTT, T>
{
    fn verify<const C: usize, P: DecompositionParams>(
        cm_i: &LCCCS<C, NTT>,
        proof: &DecompositionProof<C, NTT>,
        transcript: &mut impl Transcript<NTT>,
        _ccs: &CCS<NTT>,
    ) -> Result<Vec<LCCCS<C, NTT>>, DecompositionError> {
        let mut lcccs_s = Vec::<LCCCS<C, NTT>>::with_capacity(P::K);

        for (((x, y), u), v) in proof
            .x_s
            .iter()
            .zip(&proof.y_s)
            .zip(&proof.u_s)
            .zip(&proof.v_s)
        {
            transcript.absorb_slice(x);
            transcript.absorb_slice(y.as_ref());
            transcript.absorb_slice(u);
            transcript.absorb_slice(v);

            let h = x
                .last()
                .cloned()
                .ok_or(DecompositionError::IncorrectLength)?;
            lcccs_s.push(LCCCS {
                r: cm_i.r.clone(),
                v: v.clone(),
                cm: y.clone(),
                u: u.clone(),
                x_w: x[0..x.len() - 1].to_vec(),
                h,
            });
        }

        let b_s: Vec<_> = (0..P::K)
            .map(|i| NTT::from((P::B_SMALL as u128).pow(i as u32)))
            .collect();

        let should_equal_y0 = Self::recompose_commitment::<C>(&proof.y_s, &b_s)?;

        if should_equal_y0 != cm_i.cm {
            return Err(DecompositionError::RecomposedError);
        }

        let should_equal_u0: Vec<NTT> = Self::recompose_u(&proof.u_s, &b_s)?;

        if should_equal_u0 != cm_i.u {
            return Err(DecompositionError::RecomposedError);
        }

        for (row, &cm_i_value) in cm_i.v.iter().enumerate() {
            let should_equal_v0: NTT = Self::recompose_v(&proof.v_s, &b_s, row);

            if should_equal_v0 != cm_i_value {
                return Err(DecompositionError::RecomposedError);
            }
        }

        let (should_equal_xw, should_equal_h) = Self::recompose_xw_and_h(&proof.x_s, &b_s)?;

        if should_equal_h != cm_i.h {
            return Err(DecompositionError::RecomposedError);
        }

        if should_equal_xw != cm_i.x_w {
            return Err(DecompositionError::RecomposedError);
        }

        Ok(lcccs_s)
    }
}

fn sanity_check<NTT: SuitableRing, DP: DecompositionParams>(
    ccs: &CCS<NTT>,
) -> Result<(), DecompositionError> {
    if ccs.m != usize::max((ccs.n - ccs.l - 1) * DP::L, ccs.m).next_power_of_two() {
        return Err(CSError::InvalidSizeBounds(ccs.m, ccs.n, DP::L).into());
    }

    Ok(())
}
