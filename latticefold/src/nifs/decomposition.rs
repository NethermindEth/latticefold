#![allow(non_snake_case, clippy::upper_case_acronyms)]

use ark_std::marker::PhantomData;
use lattirust_ring::{
    balanced_decomposition::{decompose_balanced_vec, pad_and_transpose, recompose},
    OverField, Ring,
};

use crate::{
    arith::{utils::mat_vec_mul, Witness, CCS, LCCCS},
    commitment::AjtaiCommitmentScheme,
    commitment::Commitment,
    nifs::error::DecompositionError,
    parameters::DecompositionParams,
    transcript::Transcript,
    utils::mle::dense_vec_to_dense_mle,
};
use cyclotomic_rings::SuitableRing;

#[derive(Clone)]
pub struct DecompositionProof<const C: usize, NTT: Ring> {
    pub u_s: Vec<Vec<NTT>>,
    pub v_s: Vec<NTT>,
    pub x_s: Vec<Vec<NTT>>,
    pub y_s: Vec<Commitment<C, NTT>>,
}
pub trait DecompositionProver<NTT: SuitableRing, T: Transcript<NTT>> {
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
    >;
}

pub trait DecompositionVerifier<NTT: OverField, T: Transcript<NTT>> {
    fn verify<const C: usize, P: DecompositionParams>(
        cm_i: &LCCCS<C, NTT>,
        proof: &DecompositionProof<C, NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<Vec<LCCCS<C, NTT>>, DecompositionError>;
}

pub struct LFDecompositionProver<NTT, T> {
    _ntt: PhantomData<NTT>,
    _t: PhantomData<T>,
}

pub struct LFDecompositionVerifier<NTT, T> {
    _ntt: PhantomData<NTT>,
    _t: PhantomData<T>,
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
        let wit_s: Vec<Witness<NTT>> = {
            let f_s = decompose_B_vec_into_k_vec::<NTT, P>(&wit.f);
            f_s.into_iter().map(|f| Witness::from_f::<P>(f)).collect()
        };

        let mut cm_i_x_w = cm_i.x_w.clone();
        cm_i_x_w.push(cm_i.h);
        let x_s = decompose_big_vec_into_k_vec_and_compose_back::<NTT, P>(&cm_i_x_w);

        let y_s: Vec<Commitment<C, NTT>> = wit_s
            .iter()
            .map(|wit| wit.commit::<C, W, P>(scheme))
            .collect::<Result<Vec<_>, _>>()?;

        let v_s: Vec<NTT> = wit_s
            .iter()
            .map(|wit| {
                dense_vec_to_dense_mle(ccs.s, &wit.f_hat)
                    .evaluate(&cm_i.r)
                    .ok_or(DecompositionError::WitnessMleEvalFail)
            })
            .collect::<Result<Vec<_>, _>>()?;

        let mut u_s: Vec<Vec<NTT>> = Vec::with_capacity(ccs.M.len());

        for (i, wit) in wit_s.iter().enumerate() {
            let mut u_s_for_i = Vec::with_capacity(P::K);
            let z: Vec<NTT> = {
                let mut z = Vec::with_capacity(x_s[i].len() + wit.w_ccs.len());

                z.extend_from_slice(&x_s[i]);
                z.extend_from_slice(&wit.w_ccs);

                z
            };

            for M in &ccs.M {
                u_s_for_i.push(
                    dense_vec_to_dense_mle(ccs.s, &mat_vec_mul(M, &z)?)
                        .evaluate(&cm_i.r)
                        .ok_or(DecompositionError::WitnessMleEvalFail)?,
                );
            }

            u_s.push(u_s_for_i)
        }

        let mut lcccs_s = Vec::with_capacity(P::K);

        for (((x, y), u), v) in x_s.iter().zip(&y_s).zip(&u_s).zip(&v_s) {
            transcript.absorb_slice(x);
            transcript.absorb_slice(y.as_ref());
            transcript.absorb_slice(u);
            transcript.absorb(v);

            let h = x
                .last()
                .cloned()
                .ok_or(DecompositionError::IncorrectLength)?;
            lcccs_s.push(LCCCS {
                r: cm_i.r.clone(),
                v: *v,
                cm: y.clone(),
                u: u.clone(),
                x_w: x.clone(),
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
            transcript.absorb(v);

            let h = x
                .last()
                .cloned()
                .ok_or(DecompositionError::IncorrectLength)?;
            lcccs_s.push(LCCCS {
                r: cm_i.r.clone(),
                v: *v,
                cm: y.clone(),
                u: u.clone(),
                x_w: x.clone(),
                h,
            });
        }

        let b = P::B_SMALL;
        let b_s: Vec<_> = (0..P::K).map(|i| NTT::from(b.pow(i as u32))).collect();

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

        let should_equal_v0: NTT = proof
            .v_s
            .iter()
            .zip(&b_s)
            .map(|(&v_i, b_i)| v_i * b_i)
            .sum();

        if should_equal_v0 != cm_i.v {
            return Err(DecompositionError::RecomposedError);
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

/// Decompose a vector of arbitrary norm in its NTT form into DP::K vectors
/// and applies the gadget-B matrix again.
fn decompose_big_vec_into_k_vec_and_compose_back<NTT: SuitableRing, DP: DecompositionParams>(
    x: &[NTT],
) -> Vec<Vec<NTT>> {
    let coeff_repr: Vec<NTT::CoefficientRepresentation> = x.iter().map(|&x| x.into()).collect();

    // radix-B
    let decomposed_in_B: Vec<NTT::CoefficientRepresentation> =
        pad_and_transpose(decompose_balanced_vec(&coeff_repr, DP::B, Some(DP::L)))
            .into_iter()
            .flatten()
            .collect();

    decompose_balanced_vec(&decomposed_in_B, DP::B_SMALL, Some(DP::K))
        .into_iter()
        .map(|vec| {
            vec.chunks(DP::L)
                .map(|chunk| recompose(chunk, NTT::BaseRing::from(DP::B)).into())
                .collect()
        })
        .collect()
}

/// Decompose a vector of norm B in its NTT form into DP::K small vectors.
fn decompose_B_vec_into_k_vec<NTT: SuitableRing, DP: DecompositionParams>(
    x: &[NTT],
) -> Vec<Vec<NTT>> {
    let coeff_repr: Vec<NTT::CoefficientRepresentation> = x.iter().map(|&x| x.into()).collect();

    decompose_balanced_vec(&coeff_repr, DP::B_SMALL, Some(DP::K))
        .into_iter()
        .map(|vec| vec.into_iter().map(|x| x.into()).collect())
        .collect()
}

#[cfg(test)]
mod tests {
    use lattirust_ring::Pow2CyclotomicPolyRingNTT;
    use rand::thread_rng;

    use crate::{
        arith::{
            r1cs::tests::{get_test_dummy_z_split, get_test_vitalik_z_split},
            tests::{get_test_dummy_ccs, get_test_vitalik_ccs},
            Witness, CCCS,
        },
        commitment::AjtaiCommitmentScheme,
        nifs::{
            decomposition::{
                DecompositionParams, DecompositionProver, DecompositionVerifier,
                LFDecompositionProver, LFDecompositionVerifier,
            },
            linearization::{
                LFLinearizationProver, LFLinearizationVerifier, LinearizationProver,
                LinearizationVerifier,
            },
        },
        transcript::poseidon::PoseidonTranscript,
    };
    use cyclotomic_rings::challenge_set::BinarySmallSet;

    // Boilerplate code to generate values needed for testing
    const Q: u64 = 17; // Replace with an appropriate modulus
    const N: usize = 8;
    type NTT = Pow2CyclotomicPolyRingNTT<Q, N>;
    type CS = BinarySmallSet<Q, N>;
    type T = PoseidonTranscript<Pow2CyclotomicPolyRingNTT<Q, N>, CS>;

    #[derive(Clone)]
    struct PP;

    impl DecompositionParams for PP {
        const B: u128 = 1_024;
        const L: usize = 2;
        const B_SMALL: u128 = 2;
        const K: usize = 10;
    }
    // Actual Tests
    #[test]
    fn test_decomposition() {
        const WIT_LEN: usize = 4; // 4 is the length of witness in this (Vitalik's) example
        const W: usize = WIT_LEN * PP::L; // the number of columns of the Ajtai matrix

        let ccs = get_test_vitalik_ccs::<NTT>(W);
        let (_, x_ccs, w_ccs) = get_test_vitalik_z_split::<NTT>(3);
        let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());
        let wit: Witness<NTT> = Witness::from_w_ccs::<PP>(&w_ccs);
        let cm_i: CCCS<4, NTT> = CCCS {
            cm: wit.commit::<4, W, PP>(&scheme).unwrap(),
            x_ccs,
        };

        let mut prover_transcript = PoseidonTranscript::<NTT, CS>::default();
        let mut verifier_transcript = PoseidonTranscript::<NTT, CS>::default();

        let (_, linearization_proof) =
            LFLinearizationProver::<_, T>::prove(&cm_i, &wit, &mut prover_transcript, &ccs)
                .unwrap();

        let lcccs = LFLinearizationVerifier::<_, PoseidonTranscript<NTT, CS>>::verify(
            &cm_i,
            &linearization_proof,
            &mut verifier_transcript,
            &ccs,
        )
        .unwrap();

        let (_, _, decomposition_proof) = LFDecompositionProver::<_, T>::prove::<W, 4, PP>(
            &lcccs,
            &wit,
            &mut prover_transcript,
            &ccs,
            &scheme,
        )
        .unwrap();

        let res = LFDecompositionVerifier::<_, T>::verify::<4, PP>(
            &lcccs,
            &decomposition_proof,
            &mut verifier_transcript,
            &ccs,
        );

        assert!(res.is_ok());
    }

    #[test]
    fn test_failing_decomposition() {
        const WIT_LEN: usize = 4; // 4 is the length of witness in this (Vitalik's) example
        const W: usize = WIT_LEN * PP::L; // the number of columns of the Ajtai matrix

        let ccs = get_test_vitalik_ccs::<NTT>(W);
        let (_, x_ccs, w_ccs) = get_test_vitalik_z_split::<NTT>(3);
        let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());
        let wit: Witness<NTT> = Witness::from_w_ccs::<PP>(&w_ccs);
        let cm_i: CCCS<4, NTT> = CCCS {
            cm: wit.commit::<4, W, PP>(&scheme).unwrap(),
            x_ccs,
        };

        let mut prover_transcript = PoseidonTranscript::<NTT, CS>::default();
        let mut verifier_transcript = PoseidonTranscript::<NTT, CS>::default();

        let (_, linearization_proof) =
            LFLinearizationProver::<_, T>::prove(&cm_i, &wit, &mut prover_transcript, &ccs)
                .unwrap();

        let lcccs = LFLinearizationVerifier::<_, PoseidonTranscript<NTT, CS>>::verify(
            &cm_i,
            &linearization_proof,
            &mut verifier_transcript,
            &ccs,
        )
        .unwrap();

        let (_, _, w_ccs) = get_test_vitalik_z_split::<NTT>(100);
        let fake_witness = Witness::<NTT>::from_w_ccs::<PP>(&w_ccs);

        let (_, _, decomposition_proof) = LFDecompositionProver::<_, T>::prove::<W, 4, PP>(
            &lcccs,
            &fake_witness,
            &mut prover_transcript,
            &ccs,
            &scheme,
        )
        .unwrap();

        let res = LFDecompositionVerifier::<_, T>::verify::<4, PP>(
            &lcccs,
            &decomposition_proof,
            &mut verifier_transcript,
            &ccs,
        );

        assert!(res.is_err());
    }

    #[test]
    fn test_dummy_decomposition() {
        const C: usize = 10;
        const IO: usize = 1;
        const WIT_LEN: usize = 1 << 12;
        const W: usize = WIT_LEN * PP::L;
        let r1cs_rows = WIT_LEN;

        let ccs = get_test_dummy_ccs::<NTT, IO, WIT_LEN>(r1cs_rows);
        let (_, x_ccs, w_ccs) = get_test_dummy_z_split::<NTT, IO, WIT_LEN>();
        let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());
        let wit: Witness<NTT> = Witness::from_w_ccs::<PP>(&w_ccs);
        let cm_i: CCCS<C, NTT> = CCCS {
            cm: wit.commit::<C, W, PP>(&scheme).unwrap(),
            x_ccs,
        };

        let mut prover_transcript = PoseidonTranscript::<NTT, CS>::default();
        let mut verifier_transcript = PoseidonTranscript::<NTT, CS>::default();

        let (_, linearization_proof) =
            LFLinearizationProver::<_, T>::prove(&cm_i, &wit, &mut prover_transcript, &ccs)
                .unwrap();

        let lcccs = LFLinearizationVerifier::<_, PoseidonTranscript<NTT, CS>>::verify(
            &cm_i,
            &linearization_proof,
            &mut verifier_transcript,
            &ccs,
        )
        .unwrap();

        let (_, _, decomposition_proof) = LFDecompositionProver::<_, T>::prove::<W, C, PP>(
            &lcccs,
            &wit,
            &mut prover_transcript,
            &ccs,
            &scheme,
        )
        .unwrap();

        let res = LFDecompositionVerifier::<_, T>::verify::<C, PP>(
            &lcccs,
            &decomposition_proof,
            &mut verifier_transcript,
            &ccs,
        );

        assert!(res.is_ok());
    }
}
