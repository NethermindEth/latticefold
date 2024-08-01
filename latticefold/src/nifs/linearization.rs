use std::sync::Arc;

use crate::utils::sumcheck::SumCheckError::SumCheckFailed;
use crate::{
    arith::{Witness, CCCS, CCS, LCCCS},
    transcript::Transcript,
    utils::sumcheck::{
        prover::SumCheckProver, verifier::SumCheckVerifier, SumCheckIP, SumCheckProof,
    },
};
use ark_std::iterable::Iterable;
use ark_std::log2;
use lattirust_arithmetic::{
    challenge_set::latticefold_challenge_set::OverField,
    mle::DenseMultilinearExtension,
    polynomials::{build_eq_x_r, VPAuxInfo, VirtualPolynomial},
};

use super::{
    error::LinearizationError::{self, ParametersError},
    NIFSProver, NIFSVerifier,
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
        let beta_s = transcript.get_small_challenges(log_m);

        // Step 2: Sum check protocol
        // z_ccs vector
        let mut z_ccs: Vec<R> = Vec::new();
        let x_ccs = cm_i.x_ccs.clone();
        let one = vec![R::one()];
        let w_ccs = wit.w_ccs.clone();
        z_ccs.extend(x_ccs);
        z_ccs.extend(one);
        z_ccs.extend(w_ccs);

        // Create polynomial
        let g = create_sumcheck_polynomial(log_m, &ccs.c, &ccs.M, &z_ccs, &beta_s)?;
        let prover = SumCheckProver {
            polynomial: g,
            claimed_sum: R::zero(),
            _marker: std::marker::PhantomData,
        };
        // Run sum check prover
        let (_, sum_check_proof, subclaim) = prover.prove(transcript).unwrap();

        // Step 3: Compute v, u_vector
        let r_arr = subclaim.point;

        let v = mle_val_from_vector(&wit.f_hat, &r_arr)?;
        let u = create_u(ccs.t, &ccs.M, &r_arr, &z_ccs)?;

        // Step 5: Output linearization_proof and lcccs
        let linearization_proof = LinearizationProof {
            linearization_sumcheck: sum_check_proof,
            v,
            u: u.clone(),
        };

        let lcccs = LCCCS {
            r_arr,
            v,
            y: cm_i.cm.clone(),
            u: u.clone(),
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
        let beta_s = transcript.get_small_challenges(log_m);

        let poly_info = VPAuxInfo {
            max_degree: ccs.d + 1,
            num_variables: log_m,
            phantom: std::marker::PhantomData,
        };
        let protocol = SumCheckIP::new(R::zero(), poly_info);
        let verifier = SumCheckVerifier::new(protocol);

        // Verify the transcript
        let subclaim = verifier
            .verify(&proof.linearization_sumcheck, transcript)
            .unwrap();
        let e = eq(&beta_s, &subclaim.point)?;
        let s = subclaim.expected_evaluation;

        let should_equal_s = e * ccs.c.iter().fold(R::zero(), |sum, c| {
            sum + *c * proof.u.iter().fold(R::one(), |product, u_j| product * u_j)
        });

        if should_equal_s != s {
            return Err(LinearizationError::SumCheckError(SumCheckFailed(
                should_equal_s,
                s,
            )));
        }

        Ok(LCCCS {
            r_arr: subclaim.point,
            v: proof.v,
            y: cm_i.cm.clone(),
            u: proof.u.clone(),
            x_w: cm_i.x_ccs.clone(),
            h: R::one(),
        })
    }
}

fn create_u<R: OverField>(
    length: usize,
    M: &Vec<Vec<Vec<R>>>,
    r_arr: &Vec<R>,
    z_ccs: &[R],
) -> Result<Vec<R>, LinearizationError<R>> {
    let mut u: Vec<R> = Vec::with_capacity(length);

    for M_i in M.iter() {
        for i in 1..M_i.len() {
            let b = usize_to_binary_vector::<R>(i, log2(M_i.len()) as usize)?;
            let mle_matrix_val = mle_val_from_matrix(M_i, r_arr, &b)?;
            let mle_vector_val = mle_val_from_vector(z_ccs, &b)?;
            u.push(mle_matrix_val * mle_vector_val);
        }
    }
    Ok(u)
}

fn create_sumcheck_polynomial<R: OverField>(
    log_m: usize,
    c: &Vec<R>,
    M: &Vec<Vec<Vec<R>>>,
    z_ccs: &[R],
    beta_s: &[R],
) -> Result<VirtualPolynomial<R>, LinearizationError<R>> {
    let mut g = VirtualPolynomial::new(log_m);

    for coefficient in c.iter() {
        let mut mle_list: Vec<Arc<DenseMultilinearExtension<R>>> = Vec::new();

        for matrix in M.iter() {
            // Initialize MLE for the zero vector
            let zero_vector = usize_to_binary_vector::<R>(0, log2(matrix.len()) as usize)?;
            let mle_z_ccs_b = mle_val_from_vector(z_ccs, &zero_vector)?;
            let evaluations: Vec<R> = mle_matrix_to_val_eval_second(matrix, &zero_vector)?
                .iter()
                .map(|val| *val * mle_z_ccs_b)
                .collect();
            let mut mle = DenseMultilinearExtension::from_evaluations_vec(log_m, evaluations);

            // Accumulate MLEs for other binary vectors
            for i in 1..matrix.len() {
                let b = usize_to_binary_vector::<R>(i, log2(matrix.len()) as usize)?;
                let mle_z_ccs_b = mle_val_from_vector(z_ccs, &b)?;
                let evaluations: Vec<R> = mle_matrix_to_val_eval_second(matrix, &b)?
                    .iter()
                    .map(|val| *val * mle_z_ccs_b)
                    .collect();
                let mle_next = DenseMultilinearExtension::from_evaluations_vec(log_m, evaluations);
                mle += mle_next;
            }

            mle_list.push(Arc::new(mle));
        }

        let _ = g.add_mle_list(mle_list, *coefficient);
    }

    // Multiply by eq function
    let eq = build_eq_x_r::<R>(beta_s).unwrap();
    let _ = g.mul_by_mle(eq, R::one());

    Ok(g)
}

fn eq<R: OverField>(b_arr: &Vec<R>, x_arr: &Vec<R>) -> Result<R, LinearizationError<R>> {
    if b_arr.len() != x_arr.len() {
        return Err(LinearizationError::ParametersError(String::from(
            "Eq function takes two vectors of the same length!",
        )));
    }
    Ok(b_arr.iter().zip(x_arr).fold(R::one(), |ret_value, (b, x)| {
        ret_value * ((R::one() - b) * (R::one() - x) + *b * x)
    }))
}

fn usize_to_binary_vector<R: OverField>(
    n: usize,
    length: usize,
) -> Result<Vec<R>, LinearizationError<R>> {
    if 1 << length <= n {
        return Err(LinearizationError::ParametersError(String::from(
            "Cannot put number in binary, number is too big",
        )));
    }
    let mut bits = Vec::new();
    let mut current = n;

    // Extract bits from the number
    while current > 0 {
        bits.push((current & 1) as u8);
        current >>= 1;
    }

    // Reverse to get the bits in correct order
    bits.reverse();

    // Pad with leading zeros if necessary
    if bits.len() < length {
        let padding = length - bits.len();
        bits.splice(0..0, std::iter::repeat(0).take(padding));
    }

    // Convert to the target field elements
    Ok(bits
        .into_iter()
        .map(|bit| if bit == 1 { R::one() } else { R::zero() })
        .collect())
}

fn mle_val_from_vector<R: OverField>(
    vector: &[R],
    values: &Vec<R>,
) -> Result<R, LinearizationError<R>> {
    if (1 << values.len()) != vector.len() {
        return Err(LinearizationError::ParametersError(String::from(
            "values and MLE do not match",
        )));
    }
    let mle = DenseMultilinearExtension::from_evaluations_vec(values.len(), vector.to_owned());
    Ok(mle.evaluate(values.as_slice()).unwrap())
}

fn mle_val_from_matrix<R: OverField>(
    matrix: &Vec<Vec<R>>,
    values_x: &Vec<R>,
    values_y: &Vec<R>,
) -> Result<R, LinearizationError<R>> {
    if (1 << values_y.len()) != matrix.len() || (1 << values_x.len()) != matrix[0].len() {
        return Err(LinearizationError::ParametersError(String::from(
            "values and MLE do not match",
        )));
    }

    let univariate_mle_evaluations: Result<Vec<R>, LinearizationError<R>> = matrix
        .iter()
        .map(|col| mle_val_from_vector(col, values_x))
        .collect();

    mle_val_from_vector(&univariate_mle_evaluations?, values_y)
}
#[allow(dead_code)]
fn mle_matrix_to_val_eval_first<R: OverField>(
    matrix: &Vec<Vec<R>>,
    values_x: &Vec<R>,
) -> Result<Vec<R>, LinearizationError<R>> {
    if (1 << values_x.len()) != matrix[0].len() {
        return Err(ParametersError(String::from("values and MLE do not match")));
    }
    matrix
        .iter()
        .map(|col| mle_val_from_vector(col, values_x))
        .collect()
}
// Convert a bivariate MLE to a univariate MLE by evaluating the second vector
fn mle_matrix_to_val_eval_second<R: OverField>(
    matrix: &Vec<Vec<R>>,
    values_y: &Vec<R>,
) -> Result<Vec<R>, LinearizationError<R>> {
    if (1 << values_y.len()) != matrix.len() {
        return Err(ParametersError(String::from("values and MLE do not match")));
    }

    (0..matrix[0].len())
        .map(|i| {
            mle_val_from_vector(
                matrix
                    .iter()
                    .map(|col| col[i])
                    .collect::<Vec<R>>()
                    .as_slice(),
                values_y,
            )
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use ark_ff::{One, Zero};
    use lattirust_arithmetic::ring::{Pow2CyclotomicPolyRingNTT, Zq};

    use crate::nifs::linearization::{
        eq, mle_matrix_to_val_eval_first, mle_matrix_to_val_eval_second, mle_val_from_matrix,
        mle_val_from_vector, usize_to_binary_vector,
    };

    // Boilerplate code to generate values needed for testing
    const Q: u64 = 17; // Replace with an appropriate modulus
    const N: usize = 8;

    fn generate_coefficient_i(index: usize) -> Zq<Q> {
        Zq::<Q>::from(index as u64) // Simple example: use the index as the coefficient value
    }
    fn poly_ntt() -> Pow2CyclotomicPolyRingNTT<Q, N> {
        Pow2CyclotomicPolyRingNTT::<Q, N>::from_fn(generate_coefficient_i)
    }
    fn one() -> Pow2CyclotomicPolyRingNTT<Q, N> {
        Pow2CyclotomicPolyRingNTT::<Q, N>::one()
    }
    fn zero() -> Pow2CyclotomicPolyRingNTT<17, 8> {
        Pow2CyclotomicPolyRingNTT::<Q, N>::zero()
    }

    // Actual Tests
    #[test]
    fn test_utils() {
        // Test evaluation of mle from a vector
        let evaluation_vector = vec![poly_ntt(), zero()];
        assert_eq!(
            mle_val_from_vector(&evaluation_vector, &vec![one()]).unwrap(),
            zero()
        );
        assert_ne!(
            mle_val_from_vector(&evaluation_vector, &vec![one()]).unwrap(),
            poly_ntt()
        );
        assert_eq!(
            mle_val_from_vector(&evaluation_vector, &vec![zero()]).unwrap(),
            poly_ntt()
        );
        assert_ne!(
            mle_val_from_vector(&evaluation_vector, &vec![zero()]).unwrap(),
            zero()
        );

        let evaluation_matrix = vec![vec![poly_ntt(), zero()], vec![one(), poly_ntt()]];
        assert_eq!(
            mle_val_from_matrix(&evaluation_matrix, &vec![zero()], &vec![one()]).unwrap(),
            one()
        );
        assert_ne!(
            mle_val_from_matrix(&evaluation_matrix, &vec![zero()], &vec![one()]).unwrap(),
            poly_ntt()
        );
        assert_eq!(
            mle_val_from_matrix(&evaluation_matrix, &vec![zero()], &vec![zero()]).unwrap(),
            poly_ntt()
        );
        assert_ne!(
            mle_val_from_matrix(&evaluation_matrix, &vec![zero()], &vec![zero()]).unwrap(),
            zero()
        );

        // Test the eq function
        let vector_one = vec![zero(), one(), one(), zero()];
        let vector_two = vec![zero(), one(), one(), zero()];
        let vector_three = vec![zero(), one(), one(), one()];

        assert_eq!(eq(&vector_one, &vector_two).unwrap(), one());
        assert_ne!(eq(&vector_one, &vector_two).unwrap(), zero());
        assert_eq!(eq(&vector_one, &vector_three).unwrap(), zero());
        assert_ne!(eq(&vector_one, &vector_three).unwrap(), one());

        assert_eq!(
            usize_to_binary_vector::<Pow2CyclotomicPolyRingNTT<Q, N>>(4, 8).unwrap(),
            vec![
                zero(),
                zero(),
                zero(),
                zero(),
                zero(),
                one(),
                zero(),
                zero()
            ]
        );
        assert_eq!(
            usize_to_binary_vector::<Pow2CyclotomicPolyRingNTT<Q, N>>(5, 5).unwrap(),
            vec![zero(), zero(), one(), zero(), one()]
        );
        // Test the conversion of Bivariate MLE to univariate MLE by evaluating first values
        let bivariate_mle = vec![
            vec![poly_ntt(), poly_ntt(), one(), zero()],
            vec![zero(), poly_ntt(), zero(), one()],
        ];
        assert_eq!(
            mle_matrix_to_val_eval_first(&bivariate_mle, &vec![zero(), zero()]).unwrap(),
            vec![poly_ntt(), zero()]
        );
        assert_eq!(
            mle_matrix_to_val_eval_first(&bivariate_mle, &vec![one(), zero()]).unwrap(),
            vec![poly_ntt(), poly_ntt()]
        );

        // Test the conversion of Bivariate MLE to univariate MLE by evaluating second values
        assert_eq!(
            mle_matrix_to_val_eval_second(&bivariate_mle, &vec![one()]).unwrap(),
            vec![zero(), poly_ntt(), zero(), one()]
        );
        assert_eq!(
            mle_matrix_to_val_eval_second(&bivariate_mle, &vec![zero()]).unwrap(),
            vec![poly_ntt(), poly_ntt(), one(), zero()]
        );
    }
}
