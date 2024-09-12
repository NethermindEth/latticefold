#![allow(non_snake_case)]
use ark_ff::PrimeField;
use ark_std::{marker::PhantomData, sync::Arc};
use lattirust_poly::{
    mle::DenseMultilinearExtension,
    polynomials::{build_eq_x_r, eq_eval, VPAuxInfo, VirtualPolynomial},
};
use lattirust_ring::OverField;

use super::error::LinearizationError;
use crate::{
    arith::{utils::mat_vec_mul, Instance, Witness, CCCS, CCS, LCCCS},
    transcript::Transcript,
    utils::{
        mle::dense_vec_to_dense_mle,
        sumcheck,
        sumcheck::{MLSumcheck, SumCheckError::SumCheckFailed},
    },
};

#[derive(Clone)]
pub struct LinearizationProof<NTT: OverField> {
    // Sent in the step 2. of the linearization subprotocol
    pub linearization_sumcheck: sumcheck::Proof<NTT>,
    // Sent in the step 3.
    pub v: NTT,
    pub u: Vec<NTT>,
}

pub trait LinearizationProver<NTT: OverField, T: Transcript<NTT>> {
    fn prove<const C: usize>(
        cm_i: &CCCS<C, NTT>,
        wit: &Witness<NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(LCCCS<C, NTT>, LinearizationProof<NTT>), LinearizationError<NTT>>;
}

pub trait LinearizationVerifier<NTT: OverField, T: Transcript<NTT>> {
    fn verify<const C: usize>(
        cm_i: &CCCS<C, NTT>,
        proof: &LinearizationProof<NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<LCCCS<C, NTT>, LinearizationError<NTT>>;
}

pub struct LFLinearizationProver<NTT, T> {
    _ntt: PhantomData<NTT>,
    _t: PhantomData<T>,
}

pub struct LFLinearizationVerifier<NTT, T> {
    _ntt: PhantomData<NTT>,
    _t: PhantomData<T>,
}

impl<NTT: OverField, T: Transcript<NTT>> LinearizationProver<NTT, T>
    for LFLinearizationProver<NTT, T>
{
    fn prove<const C: usize>(
        cm_i: &CCCS<C, NTT>,
        wit: &Witness<NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<(LCCCS<C, NTT>, LinearizationProof<NTT>), LinearizationError<NTT>> {
        let log_m = ccs.s;
        // Step 1: Generate the beta challenges.
        transcript.absorb_field_element(&NTT::BaseRing::from_be_bytes_mod_order(b"beta_s"));
        let beta_s: Vec<NTT> = transcript
            .get_big_challenges(log_m)
            .into_iter()
            .map(|x| x.into())
            .collect();
        // Step 2: Sum check protocol

        // z_ccs vector, i.e. concatenation x || 1 || w.
        let z_ccs: Vec<NTT> = cm_i.get_z_vector(&wit.w_ccs);

        // Prepare MLE's of the form mle[M_i \cdot z_ccs](x), a.k.a. \sum mle[M_i](x, b) * mle[z_ccs](b).
        let Mz_mles: Vec<DenseMultilinearExtension<NTT>> = ccs
            .M
            .iter()
            .map(|M| Ok(dense_vec_to_dense_mle(log_m, &mat_vec_mul(M, &z_ccs)?)))
            .collect::<Result<_, LinearizationError<_>>>()?;

        // The sumcheck polynomial
        let g = prepare_lin_sumcheck_polynomial(log_m, &ccs.c, &Mz_mles, &ccs.S, &beta_s)?;

        // Run sum check prover
        let (sum_check_proof, prover_state) = MLSumcheck::prove_as_subprotocol(transcript, &g);
        // Extract the evaluation point
        let r = prover_state
            .randomness
            .into_iter()
            .map(|x| x.into())
            .collect::<Vec<NTT>>();

        // Step 3: Compute v, u_vector

        let v = dense_vec_to_dense_mle(log_m, &wit.f_hat)
            .evaluate(&r)
            .expect("cannot end up here, because the sumcheck subroutine must yield a point of the length log m");
        let u = compute_u(&Mz_mles, &r)?;

        // Absorbing the prover's messages to the verifier.
        transcript.absorb(&v);
        transcript.absorb_slice(&u);

        // Step 5: Output linearization_proof and lcccs
        let linearization_proof = LinearizationProof {
            linearization_sumcheck: sum_check_proof,
            v,
            u: u.clone(),
        };

        let lcccs = LCCCS {
            r,
            v,
            cm: cm_i.cm.clone(),
            u,
            x_w: cm_i.x_ccs.clone(),
            h: NTT::one(),
        };

        Ok((lcccs, linearization_proof))
    }
}

impl<NTT: OverField, T: Transcript<NTT>> LinearizationVerifier<NTT, T>
    for LFLinearizationVerifier<NTT, T>
{
    fn verify<const C: usize>(
        cm_i: &CCCS<C, NTT>,
        proof: &LinearizationProof<NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<LCCCS<C, NTT>, LinearizationError<NTT>> {
        let log_m = ccs.s;
        // Step 1: Generate the beta challenges.
        transcript.absorb_field_element(&NTT::BaseRing::from_be_bytes_mod_order(b"beta_s"));
        let beta_s: Vec<NTT> = transcript
            .get_big_challenges(log_m)
            .into_iter()
            .map(|x| x.into())
            .collect();

        //Step 2: The sumcheck.
        // The polynomial has degree <= ccs.d + 1 and log_m vars.
        let poly_info = VPAuxInfo::new(log_m, ccs.d + 1);

        // Verify the sumcheck proof.
        let subclaim = MLSumcheck::verify_as_subprotocol(
            transcript,
            &poly_info,
            NTT::zero(),
            &proof.linearization_sumcheck,
        )?;

        // Absorbing the prover's messages to the verifier.
        transcript.absorb(&proof.v);
        transcript.absorb_slice(&proof.u);

        // The final evaluation claim from the sumcheck.
        let s = subclaim.expected_evaluation;

        let point_r = subclaim
            .point
            .into_iter()
            .map(|x| x.into())
            .collect::<Vec<NTT>>();

        // Step 4: reshaping the evaluation claim.
        // eq(beta, r)
        let e = eq_eval(&point_r, &beta_s)?;
        let should_equal_s = e * ccs // e * (\sum c_i * \Pi_{j \in S_i} u_j)
            .c
            .iter()
            .enumerate()
            .map(|(i, &c)| c * ccs.S[i].iter().map(|&j| proof.u[j]).product::<NTT>()) // c_i * \Pi_{j \in S_i} u_j
            .sum::<NTT>(); // \sum c_i * \Pi_{j \in S_i} u_j

        if should_equal_s != s {
            return Err(LinearizationError::SumCheckError(SumCheckFailed(
                should_equal_s,
                s,
            )));
        }

        Ok(LCCCS::<C, NTT> {
            r: point_r,
            v: proof.v,
            cm: cm_i.cm.clone(),
            u: proof.u.clone(),
            x_w: cm_i.x_ccs.clone(),
            h: NTT::one(),
        })
    }
}

/// Batch compute the values of mles at the point r.
fn compute_u<NTT: OverField>(
    Mz_mles: &[DenseMultilinearExtension<NTT>],
    r: &[NTT],
) -> Result<Vec<NTT>, LinearizationError<NTT>> {
    Mz_mles
        .iter()
        .map(|M_i_mle| {
            M_i_mle
                .evaluate(r)
                .ok_or(LinearizationError::ParametersError(format!(
                    "one of the CCS matrices has an incorrect length {}, expected {}",
                    M_i_mle.evaluations.len(),
                    1 << r.len(),
                )))
        })
        .collect()
}

/// Prepare the main linearization polynomial.
fn prepare_lin_sumcheck_polynomial<NTT: OverField>(
    log_m: usize,
    c: &[NTT],
    M_mles: &[DenseMultilinearExtension<NTT>],
    S: &[Vec<usize>],
    beta_s: &[NTT],
) -> Result<VirtualPolynomial<NTT>, LinearizationError<NTT>> {
    let mut g = VirtualPolynomial::new(log_m);

    for (i, coefficient) in c.iter().enumerate().filter(|(_, c)| !c.is_zero()) {
        let mut mle_list: Vec<Arc<DenseMultilinearExtension<NTT>>> = Vec::with_capacity(S[i].len());

        for &j in &S[i] {
            mle_list.push(Arc::new(M_mles[j].clone()));
        }

        g.add_mle_list(mle_list, *coefficient)?;
    }

    g.mul_by_mle(build_eq_x_r(beta_s)?, NTT::one())?;

    Ok(g)
}

#[cfg(test)]
mod tests {
    use ark_ff::UniformRand;
    use lattirust_poly::mle::DenseMultilinearExtension;
    use lattirust_ring::{Pow2CyclotomicPolyRingNTT, Zq};
    use rand::thread_rng;

    use crate::{
        arith::{
            r1cs::tests::{get_test_dummy_z_split, get_test_vitalik_z_split},
            tests::{get_test_dummy_ccs, get_test_vitalik_ccs},
            Witness, CCCS,
        },
        commitment::AjtaiCommitmentScheme,
        nifs::linearization::{
            LFLinearizationProver, LFLinearizationVerifier, LinearizationVerifier,
        },
        parameters::DecompositionParams,
        transcript::poseidon::PoseidonTranscript,
    };
    use cyclotomic_rings::challenge_set::BinarySmallSet;

    use super::{compute_u, LinearizationProver};

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

    #[test]
    fn test_compute_u() {
        let mut mles = Vec::with_capacity(10);

        // generate evals
        for _i in 0..10 {
            let evals: Vec<Pow2CyclotomicPolyRingNTT<Q, N>> =
                (0..8).map(|_| generate_a_ring_elem()).collect();

            mles.push(DenseMultilinearExtension::from_evaluations_slice(3, &evals))
        }

        for b in 0..8_u8 {
            let us: Vec<Pow2CyclotomicPolyRingNTT<Q, N>> = compute_u(
                &mles,
                &[
                    (b & 0x01).into(),
                    ((b & 0x2) >> 1).into(),
                    ((b & 0x4) >> 2).into(),
                ],
            )
            .unwrap();

            for (i, &u) in us.iter().enumerate() {
                assert_eq!(u, mles[i].evaluations[b.to_le() as usize]);
            }
        }
    }

    // Actual Tests
    #[test]
    fn test_linearization() {
        const Q: u64 = 17;
        const N: usize = 8;
        type R = Pow2CyclotomicPolyRingNTT<Q, N>;
        type CS = BinarySmallSet<Q, N>;
        type T = PoseidonTranscript<Pow2CyclotomicPolyRingNTT<Q, N>, CS>;

        impl DecompositionParams for PP {
            const B: u128 = 1_024;
            const L: usize = 2;
            const B_SMALL: u128 = 2;
            const K: usize = 10;
        }

        const WIT_LEN: usize = 4; // 4 is the length of witness in this (Vitalik's) example
        const W: usize = WIT_LEN * PP::L; // the number of columns of the Ajtai matrix

        let ccs = get_test_ccs::<R>(W);
        let (_, x_ccs, w_ccs) = get_test_z_split::<R>(3);
        let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());
        #[derive(Clone)]
        struct PP;

        let wit: Witness<R> = Witness::from_w_ccs::<PP>(&w_ccs);
        let cm_i: CCCS<4, R> = CCCS {
            cm: wit.commit::<4, W, PP>(&scheme).unwrap(),
            x_ccs,
        };
        let mut transcript = PoseidonTranscript::<R, CS>::default();

        let res = LFLinearizationProver::<_, T>::prove(&cm_i, &wit, &mut transcript, &ccs);

        let mut transcript = PoseidonTranscript::<R, CS>::default();

        let res = LFLinearizationVerifier::<_, PoseidonTranscript<R, CS>>::verify(
            &cm_i,
            &res.unwrap().1,
            &mut transcript,
            &ccs,
        );

        res.unwrap();
    }

    #[test]
    fn test_dummy_linearization() {
        const Q: u64 = (1 << 16) + 1;
        const N: usize = 8;
        type R = Pow2CyclotomicPolyRingNTT<Q, N>;
        type CS = BinarySmallSet<Q, N>;
        type T = PoseidonTranscript<R, CS>;
        const C: usize = 10;
        const IO: usize = 1; // io length
        const W: usize = 1 << 10; // witness length
        let ccs = get_test_dummy_ccs::<R, IO, W>();
        let (_, x_ccs, w_ccs) = get_test_dummy_z_split::<R, IO, W>();
        let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());
        #[derive(Clone)]
        struct PP;

        impl DecompositionParams for PP {
            const B: u128 = 1_024;
            const L: usize = 1;
            const B_SMALL: u128 = 2;
            const K: usize = 10;
        }

        let wit: Witness<R> = Witness::from_w_ccs::<PP>(&w_ccs);
        let cm = wit.commit::<C, W, PP>(&scheme).unwrap();
        let cm_i: CCCS<C, R> = CCCS {
            cm,
            x_ccs,
        };

        let mut transcript = PoseidonTranscript::<R, CS>::default();

        let proof = LFLinearizationProver::<_, T>::prove(&cm_i, &wit, &mut transcript, &ccs);

        let mut transcript = PoseidonTranscript::<R, CS>::default();

        let res = LFLinearizationVerifier::<_, PoseidonTranscript<R, CS>>::verify(
            &cm_i,
            &proof.unwrap().1,
            &mut transcript,
            &ccs,
        );

        res.unwrap();
    }
}
