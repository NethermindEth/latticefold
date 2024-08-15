use lattirust_arithmetic::balanced_decomposition::{
    decompose_balanced_slice_polyring, pad_and_transpose, recompose,
};
use lattirust_arithmetic::challenge_set::latticefold_challenge_set::OverField;
use lattirust_arithmetic::ring::PolyRing;

use crate::arith::utils::mat_vec_mul;
use crate::commitment::{AjtaiCommitmentScheme, AjtaiParams, Commitment};

use crate::utils::mle::dense_vec_to_dense_mle;
use crate::{
    arith::{Witness, CCS, LCCCS},
    transcript::Transcript,
};

use super::{error::DecompositionError, NIFSProver, NIFSVerifier};

#[derive(Clone)]
pub struct DecompositionProof<NTT: OverField, P: AjtaiParams> {
    pub u_s: Vec<Vec<NTT>>,
    pub v_s: Vec<NTT>,
    pub x_s: Vec<Vec<NTT>>,
    pub y_s: Vec<Commitment<NTT, P>>,
}

pub trait DecompositionParams: Clone {
    type AP: AjtaiParams;
    // the small b such that we decompose a witness into vectors of the norm b.
    const SMALL_B: u128;
    // k such that b^k = B
    const K: usize;
}

pub trait DecompositionProver<
    CR: PolyRing + From<NTT> + Into<NTT>,
    NTT: OverField,
    P: DecompositionParams,
    T: Transcript<NTT>,
>
{
    type Proof: Clone;
    type Error: std::error::Error;

    fn prove(
        cm_i: &LCCCS<NTT, P::AP>,
        wit: &Witness<NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
        ajtai: &AjtaiCommitmentScheme<CR, NTT, P::AP>,
    ) -> Result<(Vec<LCCCS<NTT, P::AP>>, Vec<Witness<NTT>>, Self::Proof), Self::Error>;
}

pub trait DecompositionVerifier<
    CR: PolyRing + From<NTT> + Into<NTT>,
    NTT: OverField,
    P: DecompositionParams,
    T: Transcript<NTT>,
>
{
    type Prover: DecompositionProver<CR, NTT, P, T>;
    type Error = <Self::Prover as DecompositionProver<CR, NTT, P, T>>::Error;

    fn verify(
        cm_i: &LCCCS<NTT, P::AP>,
        proof: &<Self::Prover as DecompositionProver<CR, NTT, P, T>>::Proof,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<Vec<LCCCS<NTT, P::AP>>, Self::Error>;
}

impl<
        CR: PolyRing<BaseRing = NTT::BaseRing> + From<NTT> + Into<NTT>,
        NTT: OverField,
        P: AjtaiParams,
        DP: DecompositionParams<AP = P>,
        T: Transcript<NTT>,
    > DecompositionProver<CR, NTT, DP, T> for NIFSProver<CR, NTT, P, DP, T>
{
    type Proof = DecompositionProof<NTT, P>;
    type Error = DecompositionError;

    fn prove(
        cm_i: &LCCCS<NTT, P>,
        wit: &Witness<NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
        ajtai: &AjtaiCommitmentScheme<CR, NTT, P>,
    ) -> Result<
        (
            Vec<LCCCS<NTT, P>>,
            Vec<Witness<NTT>>,
            DecompositionProof<NTT, P>,
        ),
        DecompositionError,
    > {
        let wit_s: Vec<Witness<NTT>> = {
            let f_s = decompose_B_vec_into_k_vec::<CR, NTT, DP>(&wit.f);
            f_s.into_iter()
                .map(|f| Witness::<NTT>::from_f::<CR, P>(f))
                .collect()
        };

        let mut cm_i_x_w = cm_i.x_w.clone();
        cm_i_x_w.push(cm_i.h);
        let x_s = decompose_big_vec_into_k_vec_and_compose_back::<CR, NTT, DP>(&cm_i_x_w);

        let y_s: Vec<Commitment<NTT, P>> = wit_s
            .iter()
            .map(|wit| wit.commit(ajtai))
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
            let mut u_s_for_i = Vec::with_capacity(DP::K);
            let z: Vec<NTT> = {
                let mut z = Vec::with_capacity(x_s[i].len() + wit.w_ccs.len());

                z.extend_from_slice(&x_s[i]);
                z.extend_from_slice(&wit.w_ccs);

                z
            };

            for M in &ccs.M {
                u_s_for_i.push(
                    dense_vec_to_dense_mle(ccs.m, &mat_vec_mul(M, &z)?)
                        .evaluate(&cm_i.r)
                        .ok_or(DecompositionError::WitnessMleEvalFail)?,
                );
            }

            u_s.push(u_s_for_i)
        }

        let mut lcccs_s = Vec::<LCCCS<NTT, P>>::with_capacity(DP::K);
        for (((x, y), u), v) in x_s.iter().zip(&y_s).zip(&u_s).zip(&v_s) {
            transcript.absorb_ring_vec(x);
            transcript.absorb_ring_vec(y.as_ref());
            transcript.absorb_ring_vec(u);
            transcript.absorb_ring(v);

            let x = x.clone();
            let h = x
                .last()
                .cloned()
                .ok_or(DecompositionError::IncorrectLength)?;
            lcccs_s.push(LCCCS {
                r: cm_i.r.clone(),
                v: *v,
                cm: y.clone(),
                u: u.clone(),
                x_w: x,
                h,
            })
        }

        let proof = DecompositionProof { u_s, v_s, x_s, y_s };

        Ok((lcccs_s, wit_s, proof))
    }
}

impl<
        CR: PolyRing<BaseRing = NTT::BaseRing> + From<NTT> + Into<NTT>,
        NTT: OverField,
        P: AjtaiParams,
        DP: DecompositionParams<AP = P>,
        T: Transcript<NTT>,
    > DecompositionVerifier<CR, NTT, DP, T> for NIFSVerifier<CR, NTT, P, DP, T>
{
    type Prover = NIFSProver<CR, NTT, P, DP, T>;

    fn verify(
        cm_i: &LCCCS<NTT, P>,
        proof: &<Self::Prover as DecompositionProver<CR, NTT, DP, T>>::Proof,
        transcript: &mut impl Transcript<NTT>,
        _ccs: &CCS<NTT>,
    ) -> Result<Vec<LCCCS<NTT, P>>, DecompositionError> {
        let mut lcccs_s = Vec::<LCCCS<NTT, P>>::with_capacity(DP::K);

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

            let x = x.clone();
            let h = x
                .last()
                .cloned()
                .ok_or(DecompositionError::IncorrectLength)?;
            lcccs_s.push(LCCCS {
                r: cm_i.r.clone(),
                v: *v,
                cm: y.clone(),
                u: u.clone(),
                x_w: x,
                h,
            });
        }

        // TODO: Add consistency checks! That small commitments sum up to the big commitment etc.
        let b = NTT::from(DP::SMALL_B);

        let mut should_equal_y0 = proof.y_s[0].clone();
        proof.y_s.iter().enumerate().skip(1).for_each(|(i, y)| {
            let bi_part = y * b.pow([i as u64]);
            should_equal_y0 = should_equal_y0.clone() + bi_part;
        });

        match should_equal_y0 == cm_i.cm {
            true => {}
            false => {
                return Err(DecompositionError::RecomposedError);
            }
        }

        let mut should_equal_u0 = proof.u_s[0].clone();
        proof.u_s.iter().enumerate().skip(1).for_each(|(i, u_i)| {
            let bi_part: Vec<NTT> = u_i.iter().map(|&u| u * b.pow([i as u64])).collect();
            should_equal_u0 = should_equal_u0
                .iter()
                .zip(&bi_part)
                .map(|(&u0, ui)| u0 + ui)
                .collect();
        });

        match should_equal_u0 == cm_i.u {
            true => {}
            false => {
                return Err(DecompositionError::RecomposedError);
            }
        }

        let mut should_equal_v0 = proof.v_s[0];
        proof.v_s.iter().enumerate().skip(1).for_each(|(i, &v_i)| {
            let bi_part = v_i * b.pow([i as u64]);
            should_equal_v0 = should_equal_v0 + bi_part;
        });

        match should_equal_v0 == cm_i.v {
            true => {}
            false => {
                return Err(DecompositionError::RecomposedError);
            }
        }

        let mut should_equal_xw = proof.x_s[0].clone();
        proof.x_s.iter().enumerate().skip(1).for_each(|(i, x_i)| {
            let bi_part: Vec<NTT> = x_i.iter().map(|&u| u * b.pow([i as u64])).collect();
            should_equal_xw = should_equal_xw
                .iter()
                .zip(&bi_part)
                .map(|(&xw, xwi)| xw + xwi)
                .collect();
        });

        match should_equal_xw == cm_i.x_w {
            true => {}
            false => {
                return Err(DecompositionError::RecomposedError);
            }
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
        DP::AP::B,
        Some(DP::AP::L),
    ))
    .into_iter()
    .flatten()
    .collect();

    decompose_balanced_slice_polyring(&decomposed_in_B, DP::SMALL_B, Some(DP::K))
        .into_iter()
        .map(|vec| {
            vec.chunks(DP::AP::L)
                .map(|chunk| recompose(chunk, CR::from(DP::AP::B)).into())
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

    decompose_balanced_slice_polyring(&coeff_repr, DP::SMALL_B, Some(DP::K))
        .into_iter()
        .map(|vec| vec.into_iter().map(|x| x.into()).collect())
        .collect()
}
#[cfg(test)]
mod tests {
    use ark_ff::UniformRand;
    use lattirust_arithmetic::{
        challenge_set::latticefold_challenge_set::BinarySmallSet,
        ring::{Pow2CyclotomicPolyRingNTT, Zq},
    };
    use rand::thread_rng;

    use crate::{
        arith::{r1cs::tests::get_test_z_split, tests::get_test_ccs, Witness, CCCS},
        commitment::{AjtaiCommitmentScheme, AjtaiParams},
        nifs::{
            decomposition::{DecompositionParams, DecompositionProver, DecompositionVerifier},
            linearization::{LinearizationProver, LinearizationVerifier},
            NIFSProver, NIFSVerifier,
        },
        transcript::poseidon::PoseidonTranscript,
    };

    // Boilerplate code to generate values needed for testing
    const Q: u64 = 17; // Replace with an appropriate modulus
    const N: usize = 8;

    fn generate_coefficient_i(_i: usize) -> Zq<Q> {
        let mut rng = thread_rng();
        Zq::<Q>::rand(&mut rng)
    }

    fn generate_a_ring_elem() -> Pow2CyclotomicPolyRingNTT<Q, N> {
        Pow2CyclotomicPolyRingNTT::<Q, N>::from_fn(generate_coefficient_i)
    }

    // Actual Tests
    #[test]
    fn test_decomposition() {
        const Q: u64 = 17;
        const N: usize = 8;
        type NTT = Pow2CyclotomicPolyRingNTT<Q, N>;
        type CR = Pow2CyclotomicPolyRingNTT<Q, N>;
        type CS = BinarySmallSet<Q, N>;
        type T = PoseidonTranscript<Pow2CyclotomicPolyRingNTT<Q, N>, CS>;
        let ccs = get_test_ccs::<NTT>();
        let (_, x_ccs, w_ccs) = get_test_z_split::<NTT>(3);
        let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());
        #[derive(Clone, Eq, PartialEq)]
        struct P;

        #[derive(Clone, Eq, PartialEq)]
        struct DP;

        impl AjtaiParams for P {
            const B: u128 = 1000;
            const L: usize = 1;
            const WITNESS_SIZE: usize = 4;
            const OUTPUT_SIZE: usize = 4;
        }

        impl DecompositionParams for DP {
            type AP = P;
            const SMALL_B: u128 = 10;
            const K: usize = 3;
        }

        let wit: Witness<NTT> = Witness::<NTT>::from_w_ccs::<CR, P>(w_ccs);
        let cm_i: CCCS<NTT, P> = CCCS {
            cm: wit.commit::<NTT, P>(&scheme).unwrap(),
            x_ccs,
        };
        let mut transcript = PoseidonTranscript::<NTT, CS>::default();

        let (_, linearization_proof) = <NIFSProver<CR, NTT, P, DP, T> as LinearizationProver<
            NTT,
            P,
            T,
        >>::prove(&cm_i, &wit, &mut transcript, &ccs)
        .unwrap();

        let mut transcript = PoseidonTranscript::<NTT, CS>::default();

        let lcccs = <NIFSVerifier<CR, NTT, P, DP, T> as LinearizationVerifier<NTT, P, T>>::verify(
            &cm_i,
            &linearization_proof,
            &mut transcript,
            &ccs,
        )
        .unwrap();

        let mut transcript = PoseidonTranscript::<NTT, CS>::default();

        let (_, _, decomposition_proof) = <NIFSProver<CR, NTT, P, DP, T> as DecompositionProver<
            CR,
            NTT,
            DP,
            T,
        >>::prove(
            &lcccs, &wit, &mut transcript, &ccs, &scheme
        )
        .unwrap();

        let mut transcript = PoseidonTranscript::<NTT, CS>::default();

        <NIFSVerifier<CR, NTT, P, DP, T> as DecompositionVerifier<CR, NTT, DP, T>>::verify(
            &lcccs,
            &decomposition_proof,
            &mut transcript,
            &ccs,
        )
        .unwrap();
    }
}
