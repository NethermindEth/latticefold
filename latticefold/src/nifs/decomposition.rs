#![allow(non_snake_case, clippy::upper_case_acronyms)]

use lattirust_arithmetic::{
    balanced_decomposition::{decompose_balanced_slice_polyring, pad_and_transpose, recompose},
    challenge_set::latticefold_challenge_set::OverField,
    ring::PolyRing,
};
use std::marker::PhantomData;

use crate::{
    arith::{utils::mat_vec_mul, Witness, CCS, LCCCS},
    commitment::AjtaiCommitmentScheme,
    commitment::Commitment,
    nifs::error::DecompositionError,
    parameters::DecompositionParams,
    transcript::Transcript,
    utils::mle::dense_vec_to_dense_mle,
};

#[derive(Clone)]
pub struct DecompositionProof<const C: usize, NTT: OverField> {
    pub u_s: Vec<Vec<NTT>>,
    pub v_s: Vec<NTT>,
    pub x_s: Vec<Vec<NTT>>,
    pub y_s: Vec<Commitment<C, NTT>>,
}
pub trait DecompositionProver<NTT: OverField, T: Transcript<NTT>> {
    fn prove<
        const W: usize,
        const C: usize,
        CR: PolyRing<BaseRing = NTT::BaseRing> + Into<NTT> + From<NTT>,
        P: DecompositionParams,
    >(
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

impl<NTT: OverField, T: Transcript<NTT>> DecompositionProver<NTT, T>
    for LFDecompositionProver<NTT, T>
{
    fn prove<
        const W: usize,
        const C: usize,
        CR: PolyRing<BaseRing = NTT::BaseRing> + Into<NTT> + From<NTT>,
        P: DecompositionParams,
    >(
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
            let f_s = decompose_B_vec_into_k_vec::<CR, NTT, P>(&wit.f);
            f_s.into_iter()
                .map(|f| Witness::<NTT>::from_f::<CR, P>(f))
                .collect()
        };

        let mut cm_i_x_w = cm_i.x_w.clone();
        cm_i_x_w.push(cm_i.h);
        let x_s = decompose_big_vec_into_k_vec_and_compose_back::<CR, NTT, P>(&cm_i_x_w);

        let y_s: Vec<Commitment<C, NTT>> = wit_s
            .iter()
            .map(|wit| wit.commit::<C, W, CR, P>(scheme))
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
            transcript.absorb_ring_vec(x);
            transcript.absorb_ring_vec(y.as_ref());
            transcript.absorb_ring_vec(u);
            transcript.absorb_ring(v);

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
            transcript.absorb_ring_vec(x);
            transcript.absorb_ring_vec(y.as_ref());
            transcript.absorb_ring_vec(u);
            transcript.absorb_ring(v);

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

        let mut should_equal_y0 = proof.y_s[0].clone();

        proof.y_s.iter().enumerate().skip(1).for_each(|(i, y)| {
            let bi_part = y * NTT::from(b.pow(i as u32));

            should_equal_y0 = should_equal_y0.clone() + bi_part;
        });

        if should_equal_y0 != cm_i.cm {
            return Err(DecompositionError::RecomposedError);
        }

        let mut should_equal_u0 = proof.u_s[0].clone();
        proof.u_s.iter().enumerate().skip(1).for_each(|(i, u_i)| {
            let bi_part: Vec<NTT> = u_i
                .iter()
                .map(|&u| u * NTT::from(b.pow(i as u32)))
                .collect();
            should_equal_u0 = should_equal_u0
                .iter()
                .zip(&bi_part)
                .map(|(&u0, ui)| u0 + ui)
                .collect();
        });

        if should_equal_u0 != cm_i.u {
            return Err(DecompositionError::RecomposedError);
        }

        let mut should_equal_v0 = proof.v_s[0];
        proof.v_s.iter().enumerate().skip(1).for_each(|(i, &v_i)| {
            let bi_part = v_i * NTT::from(b.pow(i as u32));
            should_equal_v0 += bi_part;
        });

        if should_equal_v0 != cm_i.v {
            return Err(DecompositionError::RecomposedError);
        }

        let mut should_equal_xw = proof.x_s[0].clone();

        proof.x_s.iter().enumerate().skip(1).for_each(|(i, x_i)| {
            let bi_part: Vec<NTT> = x_i
                .iter()
                .map(|&u| u * NTT::from(b.pow(i as u32)))
                .collect();

            should_equal_xw = should_equal_xw
                .iter()
                .zip(&bi_part)
                .map(|(&xw, xwi)| xw + xwi)
                .collect();
        });
        should_equal_xw.pop();
        if should_equal_xw != cm_i.x_w {
            return Err(DecompositionError::RecomposedError);
        }

        Ok(lcccs_s)
    }
}

/// Decompose a vector of arbitrary norm in its NTT form into DP::K vectors
/// and applies the gadget-B matrix again.
fn decompose_big_vec_into_k_vec_and_compose_back<
    CR: PolyRing + From<NTT> + Into<NTT>,
    NTT: OverField,
    DP: DecompositionParams,
>(
    x: &[NTT],
) -> Vec<Vec<NTT>> {
    let coeff_repr: Vec<CR> = x.iter().map(|&x| x.into()).collect();

    // radix-B
    let decomposed_in_B: Vec<CR> = pad_and_transpose(decompose_balanced_slice_polyring(
        &coeff_repr,
        DP::B,
        Some(DP::L),
    ))
    .into_iter()
    .flatten()
    .collect();

    decompose_balanced_slice_polyring(&decomposed_in_B, DP::B_SMALL, Some(DP::K))
        .into_iter()
        .map(|vec| {
            vec.chunks(DP::L)
                .map(|chunk| recompose(chunk, CR::from(DP::B)).into())
                .collect()
        })
        .collect()
}

/// Decompose a vector of norm B in its NTT form into DP::K small vectors.
fn decompose_B_vec_into_k_vec<
    CR: PolyRing + From<NTT> + Into<NTT>,
    NTT: OverField,
    DP: DecompositionParams,
>(
    x: &[NTT],
) -> Vec<Vec<NTT>> {
    let coeff_repr: Vec<CR> = x.iter().map(|&x| x.into()).collect();

    decompose_balanced_slice_polyring(&coeff_repr, DP::B_SMALL, Some(DP::K))
        .into_iter()
        .map(|vec| vec.into_iter().map(|x| x.into()).collect())
        .collect()
}

#[cfg(test)]
mod tests {
    use lattirust_arithmetic::{
        challenge_set::latticefold_challenge_set::BinarySmallSet,
        ring::{Pow2CyclotomicPolyRing, Pow2CyclotomicPolyRingNTT, Zq},
    };
    use rand::thread_rng;

    use crate::{
        arith::{r1cs::tests::get_test_z_split, tests::get_test_ccs, Witness, CCCS},
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

    // Boilerplate code to generate values needed for testing
    const Q: u64 = 17; // Replace with an appropriate modulus
    const N: usize = 8;

    // Actual Tests
    #[test]
    fn test_decomposition() {
        type CR = Pow2CyclotomicPolyRing<Zq<Q>, N>;
        type NTT = Pow2CyclotomicPolyRingNTT<Q, N>;
        type CS = BinarySmallSet<Q, N>;
        type T = PoseidonTranscript<Pow2CyclotomicPolyRingNTT<Q, N>, CS>;
        let ccs = get_test_ccs::<NTT>();
        let (_, x_ccs, w_ccs) = get_test_z_split::<NTT>(3);
        let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());
        #[derive(Clone)]
        struct PP;

        impl DecompositionParams for PP {
            const B: u128 = 1_024;
            const L: usize = 1;
            const B_SMALL: u128 = 2;
            const K: usize = 10;
        }

        let wit: Witness<NTT> = Witness::from_w_ccs::<CR, PP>(&w_ccs);
        let cm_i: CCCS<4, NTT> = CCCS {
            cm: wit.commit::<4, 4, CR, PP>(&scheme).unwrap(),
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

        let (_, _, decomposition_proof) = LFDecompositionProver::<_, T>::prove::<4, 4, CR, PP>(
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
        type CR = Pow2CyclotomicPolyRing<Zq<Q>, N>;
        type NTT = Pow2CyclotomicPolyRingNTT<Q, N>;
        type CS = BinarySmallSet<Q, N>;
        type T = PoseidonTranscript<Pow2CyclotomicPolyRingNTT<Q, N>, CS>;
        let ccs = get_test_ccs::<NTT>();
        let (_, x_ccs, w_ccs) = get_test_z_split::<NTT>(3);
        let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());
        #[derive(Clone)]
        struct PP;

        impl DecompositionParams for PP {
            const B: u128 = 1_024;
            const L: usize = 1;
            const B_SMALL: u128 = 2;
            const K: usize = 10;
        }

        let wit: Witness<NTT> = Witness::from_w_ccs::<CR, PP>(&w_ccs);
        let cm_i: CCCS<4, NTT> = CCCS {
            cm: wit.commit::<4, 4, CR, PP>(&scheme).unwrap(),
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

        let (_, _, w_ccs) = get_test_z_split::<NTT>(100);
        let fake_witness = Witness::<NTT>::from_w_ccs::<CR, PP>(&w_ccs);
        let (_, _, decomposition_proof) = LFDecompositionProver::<_, T>::prove::<4, 4, CR, PP>(
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
}
