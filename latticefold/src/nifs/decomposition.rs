#![allow(non_snake_case, clippy::upper_case_acronyms)]

use ark_std::{cfg_into_iter, cfg_iter};
use cyclotomic_rings::rings::SuitableRing;
use num_traits::Zero;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use stark_rings::OverField;
use stark_rings_linalg::SparseMatrix;
use stark_rings_poly::polynomials::DenseMultilinearExtension;

pub use self::structs::*;
use self::utils::{decompose_B_vec_into_k_vec, decompose_big_vec_into_k_vec_and_compose_back};
use crate::{
    arith::{error::CSError, utils::mat_vec_mul, Witness, CCS, LCCCS},
    ark_base::*,
    commitment::{AjtaiCommitmentScheme, Commitment, CommitmentError},
    decomposition_parameters::DecompositionParams,
    nifs::error::DecompositionError,
    transcript::Transcript,
    utils::mle_helpers::{evaluate_mles, to_mles_err},
};

mod structs;

#[cfg(test)]
mod tests;

mod utils;

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
            Vec<Vec<DenseMultilinearExtension<NTT>>>,
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

        let y_s: Vec<Commitment<C, NTT>> = Self::commit_witnesses::<C, W, P>(&wit_s, scheme, cm_i)?;

        let v_s: Vec<Vec<NTT>> = Self::compute_v_s(&wit_s, &cm_i.r)?;

        let mz_mles = Self::compute_mz_mles(&wit_s, &ccs.M, &x_s, log_m)?;

        let u_s = Self::compute_u_s(&mz_mles, &cm_i.r)?;

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

        Ok((mz_mles, lcccs_s, wit_s, proof))
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

        let b_s: Vec<_> = Self::calculate_b_s::<P>();

        let y = Self::recompose_commitment::<C>(&proof.y_s, &b_s)?;
        if y != cm_i.cm {
            return Err(DecompositionError::RecomposedError);
        }

        let v = Self::recompose(&proof.v_s, &b_s)?;
        if v != cm_i.v {
            return Err(DecompositionError::RecomposedError);
        }

        let u = Self::recompose(&proof.u_s, &b_s)?;
        if u != cm_i.u {
            return Err(DecompositionError::RecomposedError);
        }

        let mut x_w = Self::recompose(&proof.x_s, &b_s)?;
        let h = x_w.pop().ok_or(DecompositionError::IncorrectLength)?;
        if x_w != cm_i.x_w {
            return Err(DecompositionError::RecomposedError);
        }
        if h != cm_i.h {
            return Err(DecompositionError::RecomposedError);
        }

        Ok(lcccs_s)
    }
}

impl<NTT: SuitableRing, T: Transcript<NTT>> LFDecompositionProver<NTT, T> {
    /// Decomposes a witness `wit` into `P::K` vectors norm `< P::B_SMALL` such that
    /// $$ \text{wit} = \sum\limits_{i=0}^{\text{P::K} - 1} \text{P::B\\_SMALL}^i \cdot \text{wit}_i.$$
    ///
    fn decompose_witness<P: DecompositionParams>(wit: &Witness<NTT>) -> Vec<Witness<NTT>> {
        let f_s = decompose_B_vec_into_k_vec::<NTT, P>(&wit.f_coeff);
        cfg_into_iter!(f_s)
            .map(|f| Witness::from_f_coeff::<P>(f))
            .collect()
    }

    /// Takes the concatenation `x_w || h`, performs gadget decomposition of it,
    /// decomposes the resulting `P::B`-short vector into `P::K` `P::B_SMALL`-vectors
    /// and gadget-composes each of the vectors back to obtain `P::K` vectors in their NTT form.
    fn compute_x_s<P: DecompositionParams>(mut x_w: Vec<NTT>, h: NTT) -> Vec<Vec<NTT>> {
        x_w.push(h);
        decompose_big_vec_into_k_vec_and_compose_back::<NTT, P>(x_w)
    }

    /// Ajtai commits to witnesses `wit_s` using Ajtai commitment scheme `scheme`.
    fn commit_witnesses<const C: usize, const W: usize, P: DecompositionParams>(
        wit_s: &[Witness<NTT>],
        scheme: &AjtaiCommitmentScheme<C, W, NTT>,
        cm_i: &LCCCS<C, NTT>,
    ) -> Result<Vec<Commitment<C, NTT>>, CommitmentError> {
        let b = NTT::from(P::B_SMALL as u128);

        let commitments_k1: Vec<_> = cfg_iter!(wit_s[1..])
            .map(|wit| wit.commit::<C, W, P>(scheme))
            .collect::<Result<_, _>>()?;

        let b_sum = commitments_k1
            .iter()
            .rev()
            .fold(Commitment::zero(), |acc, y_i| (acc + y_i) * b);

        let mut result = Vec::with_capacity(wit_s.len());
        result.push(&cm_i.cm - b_sum);
        result.extend(commitments_k1);

        Ok(result)
    }

    /// Compute f-hat evaluation claims.
    fn compute_v_s(
        wit_s: &[Witness<NTT>],
        point_r: &[NTT],
    ) -> Result<Vec<Vec<NTT>>, DecompositionError> {
        cfg_iter!(wit_s)
            .map(|wit| evaluate_mles::<NTT, _, _, DecompositionError>(&wit.f_hat, point_r))
            .collect::<Result<Vec<_>, _>>()
    }

    /// Compute CCS-linearization evaluation claims.
    fn compute_u_s(
        mz_mles: &[Vec<DenseMultilinearExtension<NTT>>],
        point_r: &[NTT],
    ) -> Result<Vec<Vec<NTT>>, DecompositionError> {
        cfg_iter!(mz_mles)
            .map(|mles| {
                let u_s_for_i =
                    evaluate_mles::<NTT, &DenseMultilinearExtension<_>, _, DecompositionError>(
                        mles, point_r,
                    )?;

                Ok(u_s_for_i)
            })
            .collect::<Result<Vec<Vec<NTT>>, DecompositionError>>()
    }
    fn compute_mz_mles(
        wit_s: &Vec<Witness<NTT>>,
        M: &[SparseMatrix<NTT>],
        decomposed_statements: &[Vec<NTT>],
        num_mle_vars: usize,
    ) -> Result<Vec<Vec<DenseMultilinearExtension<NTT>>>, DecompositionError> {
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

                Ok(mles)
            })
            .collect::<Result<Vec<Vec<_>>, DecompositionError>>()
    }
}

impl<NTT: OverField, T: Transcript<NTT>> LFDecompositionVerifier<NTT, T> {
    /// Recomposes `s`, calculating the linear combination `b[0] * s[0][j] + b[1] * s[1][j] + ... + b[s.len() - 1] * s[s.len() - 1][j]`
    /// for each element indexed at `j`.
    pub fn recompose(s: &[Vec<NTT>], b: &[NTT]) -> Result<Vec<NTT>, DecompositionError> {
        if s.is_empty() {
            return Err(DecompositionError::RecomposedError);
        }
        let len = s[0].len();
        Ok((0..len)
            .map(|j| {
                s.iter()
                    .zip(b)
                    .fold(NTT::zero(), |acc, (s_i, b_i)| acc + s_i[j] * b_i)
            })
            .collect())
    }

    /// Computes the linear combination `coeffs[0] * y_s[0] + coeffs[1] * y_s[1] + ... + coeffs[y_s.len() - 1] * y_s[y_s.len() - 1]`.
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

    fn calculate_b_s<P: DecompositionParams>() -> Vec<NTT> {
        (0..P::K)
            .map(|i| NTT::from((P::B_SMALL as u128).pow(i as u32)))
            .collect()
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
