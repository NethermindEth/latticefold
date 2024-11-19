#![allow(non_snake_case, clippy::upper_case_acronyms)]
use crate::{
    arith::{utils::mat_vec_mul, Witness, CCS, LCCCS},
    commitment::AjtaiCommitmentScheme,
    commitment::Commitment,
    decomposition_parameters::DecompositionParams,
    nifs::error::DecompositionError,
    transcript::Transcript,
};
use cyclotomic_rings::rings::SuitableRing;
use lattirust_poly::polynomials::DenseMultilinearExtension;
use lattirust_ring::OverField;
use utils::{decompose_B_vec_into_k_vec, decompose_big_vec_into_k_vec_and_compose_back};

use ark_std::{cfg_into_iter, cfg_iter};
#[cfg(feature = "parallel")]
use rayon::prelude::*;

pub use structs::*;
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
            Vec<LCCCS<C, NTT>>,
            Vec<Witness<NTT>>,
            DecompositionProof<C, NTT>,
        ),
        DecompositionError,
    > {
        let log_m = ccs.s;

        let wit_s: Vec<Witness<NTT>> = {
            let f_s = decompose_B_vec_into_k_vec::<NTT, P>(&wit.f_coeff);
            cfg_into_iter!(f_s)
                .map(|f| Witness::from_f_coeff::<P>(f))
                .collect()
        };

        let mut cm_i_x_w = cm_i.x_w.clone();
        cm_i_x_w.push(cm_i.h);
        let x_s = decompose_big_vec_into_k_vec_and_compose_back::<NTT, P>(cm_i_x_w);

        let y_s: Vec<Commitment<C, NTT>> = cfg_iter!(wit_s)
            .map(|wit| wit.commit::<C, W, P>(scheme))
            .collect::<Result<Vec<_>, _>>()?;

        let v_s: Vec<Vec<NTT>> = cfg_iter!(wit_s)
            .map(|wit| {
                wit.f_hat
                    .iter()
                    .map(|f_hat_row| {
                        DenseMultilinearExtension::from_slice(log_m, f_hat_row)
                            .evaluate(&cm_i.r)
                            .ok_or(DecompositionError::WitnessMleEvalFail)
                    })
                    .collect::<Result<Vec<_>, _>>()
            })
            .collect::<Result<Vec<_>, _>>()?;

        let u_s = cfg_iter!(wit_s)
            .enumerate()
            .map(|(i, wit)| {
                let mut u_s_for_i = Vec::with_capacity(ccs.t);

                let z: Vec<NTT> = {
                    let mut z = Vec::with_capacity(x_s[i].len() + wit.w_ccs.len());

                    z.extend_from_slice(&x_s[i]);
                    z.extend_from_slice(&wit.w_ccs);

                    z
                };

                for M in &ccs.M {
                    u_s_for_i.push(
                        DenseMultilinearExtension::from_slice(ccs.s, &mat_vec_mul(M, &z)?)
                            .evaluate(&cm_i.r)
                            .ok_or(DecompositionError::WitnessMleEvalFail)?,
                    );
                }

                Ok(u_s_for_i)
            })
            .collect::<Result<Vec<Vec<NTT>>, DecompositionError>>()?;

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

        let should_equal_y0 = proof
            .y_s
            .iter()
            .zip(&b_s)
            .map(|(y_i, b_i)| y_i * b_i)
            .reduce(|acc, bi_part| acc + bi_part)
            .ok_or(DecompositionError::RecomposedError)?;

        if should_equal_y0 != cm_i.cm {
            return Err(DecompositionError::RecomposedError);
        }

        let should_equal_u0: Vec<NTT> = proof
            .u_s
            .iter()
            .zip(&b_s)
            .map(|(u_i, b_i)| u_i.iter().map(|&u| u * b_i).collect())
            .reduce(|acc, u_i_times_b_i: Vec<NTT>| {
                acc.into_iter()
                    .zip(&u_i_times_b_i)
                    .map(|(u0, ui)| u0 + ui)
                    .collect()
            })
            .ok_or(DecompositionError::RecomposedError)?;

        if should_equal_u0 != cm_i.u {
            return Err(DecompositionError::RecomposedError);
        }

        for (i, &cm_i_value) in cm_i.v.iter().enumerate() {
            let should_equal_v0: NTT = proof
                .v_s
                .iter()
                .zip(&b_s)
                .map(|(v_i, b_i)| v_i[i] * b_i)
                .sum();

            if should_equal_v0 != cm_i_value {
                return Err(DecompositionError::RecomposedError);
            }
        }

        let mut should_equal_xw: Vec<NTT> = proof
            .x_s
            .iter()
            .zip(&b_s)
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

        if should_equal_h != cm_i.h {
            return Err(DecompositionError::RecomposedError);
        }

        if should_equal_xw != cm_i.x_w {
            return Err(DecompositionError::RecomposedError);
        }

        Ok(lcccs_s)
    }
}
