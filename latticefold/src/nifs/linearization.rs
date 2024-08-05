use std::sync::Arc;

use super::{
    error::LinearizationError::{self},
    NIFSProver, NIFSVerifier,
};
use crate::arith::utils::mat_vec_mul;
use crate::{
    arith::Instance,
    utils::{mle::dense_vec_to_dense_mle, sumcheck::SumCheckError::SumCheckFailed},
};
use crate::{
    arith::{Witness, CCCS, CCS, LCCCS},
    transcript::Transcript,
    utils::sumcheck::{
        prover::SumCheckProver, verifier::SumCheckVerifier, SumCheckIP, SumCheckProof,
    },
};
use ark_ff::PrimeField;
use lattirust_arithmetic::polynomials::eq_eval;
use lattirust_arithmetic::{
    challenge_set::latticefold_challenge_set::OverField,
    mle::DenseMultilinearExtension,
    polynomials::{build_eq_x_r, VPAuxInfo, VirtualPolynomial},
};

#[derive(Clone)]
pub struct LinearizationProof<R: OverField> {
    // Sent in the step 2. of the linearization subprotocol
    pub linearization_sumcheck: SumCheckProof<R>,
    // Sent in the step 3.
    pub v: R,
    pub u: Vec<R>,
}

pub trait LinearizationProver<R: OverField, T: Transcript<R>> {
    type Proof: Clone;
    type Error: std::error::Error;

    fn prove(
        cm_i: &CCCS<R>,
        wit: &Witness<R>,
        transcript: &mut impl Transcript<R>,
        ccs: &CCS<R>,
    ) -> Result<(LCCCS<R>, Self::Proof), Self::Error>;
}

pub trait LinearizationVerifier<R: OverField, T: Transcript<R>> {
    type Prover: LinearizationProver<R, T>;
    type Error = <Self::Prover as LinearizationProver<R, T>>::Error;

    fn verify(
        cm_i: &CCCS<R>,
        proof: &<Self::Prover as LinearizationProver<R, T>>::Proof,
        transcript: &mut impl Transcript<R>,
        ccs: &CCS<R>,
    ) -> Result<LCCCS<R>, Self::Error>;
}

impl<R: OverField, T: Transcript<R>> LinearizationProver<R, T> for NIFSProver<R, T> {
    type Proof = LinearizationProof<R>;
    type Error = LinearizationError<R>;

    fn prove(
        cm_i: &CCCS<R>,
        wit: &Witness<R>,
        transcript: &mut impl Transcript<R>,
        ccs: &CCS<R>,
    ) -> Result<(LCCCS<R>, LinearizationProof<R>), LinearizationError<R>> {
        let log_m = ccs.s;
        // Step 1: Generate the public coin randomness
        transcript.absorb(&R::F::from_be_bytes_mod_order(b"beta_s"));
        let beta_s = transcript.get_small_challenges(log_m);

        // Step 2: Sum check protocol
        // z_ccs vector
        let z_ccs: Vec<R> = cm_i.get_z_vector(&wit.w_ccs);

        let Mz_mles: Vec<DenseMultilinearExtension<R>> = ccs
            .M
            .iter()
            .map(|M| Ok(dense_vec_to_dense_mle(log_m, &mat_vec_mul(M, &z_ccs)?)))
            .collect::<Result<_, LinearizationError<_>>>()?;

        // Create polynomial
        let g = prepare_lin_sumcheck_polynomial(log_m, &ccs.c, &Mz_mles, &ccs.S, &beta_s)?;
        // Run sum check prover
        let (sum_check_proof, subclaim) = SumCheckProver::new(g, R::zero()).prove(transcript)?;

        // Step 3: Compute v, u_vector
        let r = subclaim.point;

        let v = dense_vec_to_dense_mle(log_m, &wit.f_hat)
            .evaluate(&r)
            .expect("cannot end up here, because the sumcheck subroutine must yield a point of the length log m");
        let u = compute_u(&Mz_mles, &r)?;

        transcript.absorb_ring(&v);
        transcript.absorb_ring_vec(&u);

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
            h: R::one(),
        };

        Ok((lcccs, linearization_proof))
    }
}

impl<R: OverField, T: Transcript<R>> LinearizationVerifier<R, T> for NIFSVerifier<R, T> {
    type Prover = NIFSProver<R, T>;

    fn verify(
        cm_i: &CCCS<R>,
        proof: &<Self::Prover as LinearizationProver<R, T>>::Proof,
        transcript: &mut impl Transcript<R>,
        ccs: &CCS<R>,
    ) -> Result<LCCCS<R>, LinearizationError<R>> {
        let log_m = ccs.s;

        transcript.absorb(&R::F::from_be_bytes_mod_order(b"beta_s"));
        let beta_s = transcript.get_small_challenges(log_m);

        let poly_info = VPAuxInfo::new(ccs.d + 1, log_m);
        let verifier = SumCheckVerifier::new(SumCheckIP::new(R::zero(), poly_info));

        // Verify the transcript
        let subclaim = verifier.verify(&proof.linearization_sumcheck, transcript)?;

        transcript.absorb_ring(&proof.v);
        transcript.absorb_ring_vec(&proof.u);

        let e = eq_eval(&subclaim.point, &beta_s)?;
        let s = subclaim.expected_evaluation;

        let should_equal_s = e * ccs
            .c
            .iter()
            .enumerate()
            .map(|(i, &c)| c * ccs.S[i].iter().map(|&j| proof.u[j]).product::<R>())
            .sum::<R>();

        if should_equal_s != s {
            return Err(LinearizationError::SumCheckError(SumCheckFailed(
                should_equal_s,
                s,
            )));
        }

        Ok(LCCCS {
            r: subclaim.point,
            v: proof.v,
            cm: cm_i.cm.clone(),
            u: proof.u.clone(),
            x_w: cm_i.x_ccs.clone(),
            h: R::one(),
        })
    }
}

fn compute_u<R: OverField>(
    Mz_mles: &[DenseMultilinearExtension<R>],
    r: &[R],
) -> Result<Vec<R>, LinearizationError<R>> {
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

fn prepare_lin_sumcheck_polynomial<R: OverField>(
    log_m: usize,
    c: &[R],
    M_mles: &[DenseMultilinearExtension<R>],
    S: &[Vec<usize>],
    beta_s: &[R],
) -> Result<VirtualPolynomial<R>, LinearizationError<R>> {
    let mut g = VirtualPolynomial::new(log_m);

    for (i, coefficient) in c.iter().enumerate().filter(|(_, c)| !c.is_zero()) {
        let mut mle_list: Vec<Arc<DenseMultilinearExtension<R>>> = Vec::with_capacity(S[i].len());

        for &j in &S[i] {
            mle_list.push(Arc::new(M_mles[j].clone()));
        }

        g.add_mle_list(mle_list, *coefficient)?;
    }

    g.mul_by_mle(build_eq_x_r(beta_s)?, R::one())?;

    Ok(g)
}

#[cfg(test)]
mod tests {
    use ark_ff::UniformRand;
    use lattirust_arithmetic::{
        mle::DenseMultilinearExtension,
        ring::{Pow2CyclotomicPolyRingNTT, Zq},
    };
    use rand::thread_rng;

    use super::compute_u;

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
    // #[test]
    // fn test_utils() {
    //     // Test evaluation of mle from a vector
    //     let evaluation_vector = vec![generate_a_ring_elem(), zero()];
    //     assert_eq!(
    //         mle_val_from_vector(&evaluation_vector, &vec![one()]).unwrap(),
    //         zero()
    //     );
    //     assert_ne!(
    //         mle_val_from_vector(&evaluation_vector, &vec![one()]).unwrap(),
    //         generate_a_ring_elem()
    //     );
    //     assert_eq!(
    //         mle_val_from_vector(&evaluation_vector, &vec![zero()]).unwrap(),
    //         generate_a_ring_elem()
    //     );
    //     assert_ne!(
    //         mle_val_from_vector(&evaluation_vector, &vec![zero()]).unwrap(),
    //         zero()
    //     );

    //     let evaluation_matrix = vec![vec![generate_a_ring_elem(), zero()], vec![one(), generate_a_ring_elem()]];
    //     assert_eq!(
    //         mle_val_from_matrix(&evaluation_matrix, &vec![zero()], &vec![one()]).unwrap(),
    //         one()
    //     );
    //     assert_ne!(
    //         mle_val_from_matrix(&evaluation_matrix, &vec![zero()], &vec![one()]).unwrap(),
    //         generate_a_ring_elem()
    //     );
    //     assert_eq!(
    //         mle_val_from_matrix(&evaluation_matrix, &vec![zero()], &vec![zero()]).unwrap(),
    //         generate_a_ring_elem()
    //     );
    //     assert_ne!(
    //         mle_val_from_matrix(&evaluation_matrix, &vec![zero()], &vec![zero()]).unwrap(),
    //         zero()
    //     );

    //     assert_eq!(
    //         usize_to_binary_vector::<Pow2CyclotomicPolyRingNTT<Q, N>>(4, 8).unwrap(),
    //         vec![
    //             zero(),
    //             zero(),
    //             zero(),
    //             zero(),
    //             zero(),
    //             one(),
    //             zero(),
    //             zero()
    //         ]
    //     );
    //     assert_eq!(
    //         usize_to_binary_vector::<Pow2CyclotomicPolyRingNTT<Q, N>>(5, 5).unwrap(),
    //         vec![zero(), zero(), one(), zero(), one()]
    //     );
    //     // Test the conversion of Bivariate MLE to univariate MLE by evaluating first values
    //     let bivariate_mle = vec![
    //         vec![generate_a_ring_elem(), generate_a_ring_elem(), one(), zero()],
    //         vec![zero(), generate_a_ring_elem(), zero(), one()],
    //     ];
    //     assert_eq!(
    //         mle_matrix_to_val_eval_first(&bivariate_mle, &vec![zero(), zero()]).unwrap(),
    //         vec![generate_a_ring_elem(), zero()]
    //     );
    //     assert_eq!(
    //         mle_matrix_to_val_eval_first(&bivariate_mle, &vec![one(), zero()]).unwrap(),
    //         vec![generate_a_ring_elem(), generate_a_ring_elem()]
    //     );

    //     // Test the conversion of Bivariate MLE to univariate MLE by evaluating second values
    //     assert_eq!(
    //         mle_matrix_to_val_eval_second(&bivariate_mle, &vec![one()]).unwrap(),
    //         vec![zero(), generate_a_ring_elem(), zero(), one()]
    //     );
    //     assert_eq!(
    //         mle_matrix_to_val_eval_second(&bivariate_mle, &vec![zero()]).unwrap(),
    //         vec![generate_a_ring_elem(), generate_a_ring_elem(), one(), zero()]
    //     );
    // }
}
