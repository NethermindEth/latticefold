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
use lattirust_arithmetic::{
    challenge_set::latticefold_challenge_set::OverField,
    mle::DenseMultilinearExtension,
    polynomials::{build_eq_x_r, VPAuxInfo, VirtualPolynomial},
};
use libm::log2;

use super::{error::LinearizationError, NIFSProver, NIFSVerifier};

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
        _cm_i: &CCCS<R>,
        _wit: &Witness<R>,
        _transcript: &mut impl Transcript<R>,
        _ccs: &CCS<R>,
    ) -> Result<(LCCCS<R>, LinearizationProof<R>), LinearizationError<R>> {
        // Define some constants
        // TODO decide which ones of these are best as explicit
        let m = _ccs.m;
        let log_m = log2(m as f64) as usize;

        // Step 1: Generate the public coin randomness
        let Beta = get_challenge_vector(_cm_i, _transcript, log_m);

        // Step 2: Sum check protocol
        // z_ccs vector
        let mut z_ccs: Vec<R> = Vec::new();
        let x_ccs = _cm_i.x_ccs.clone();
        let one = vec![R::one()];
        let w_ccs = _wit.w_ccs.clone();
        z_ccs.extend(x_ccs);
        z_ccs.extend(one);
        z_ccs.extend(w_ccs);

        // Create polynomial
        let g = create_sumcheck_polynomial(log_m, &_ccs.c, &_ccs.M, &z_ccs, &Beta);
        let prover = SumCheckProver {
            polynomial: g,
            claimed_sum: R::zero(),
            _marker: std::marker::PhantomData::default(),
        };
        // Run sum check prover
        let (_, sum_check_proof, subclaim) = prover.prove(_transcript).unwrap();

        // Step 3: Compute v, u_vector
        let r_arr = subclaim.point;

        let v = mle_val_from_vector(&_wit.f_hat, &r_arr);
        let u = create_u(_ccs.t, &_ccs.M, &r_arr, &z_ccs);

        // Step 5: Output linearization_proof and lcccs
        let linearization_proof = LinearizationProof {
            linearization_sumcheck: sum_check_proof,
            v,
            u: u.clone(),
        };

        let lcccs = LCCCS {
            r_arr,
            v,
            y: _cm_i.cm.clone(),
            u: u.clone(),
            x_w: _cm_i.x_ccs.clone(),
            h: R::one(),
        };

        Ok((lcccs, linearization_proof))
    }
}

impl<R: OverField, T: Transcript<R>> LinearizationVerifier<R, T> for NIFSVerifier<R, T> {
    type Prover = NIFSProver<R, T>;

    fn verify(
        _cm_i: &CCCS<R>,
        _proof: &<Self::Prover as LinearizationProver<R, T>>::Proof,
        _transcript: &mut impl Transcript<R>,
        _ccs: &CCS<R>,
    ) -> Result<LCCCS<R>, LinearizationError<R>> {
        let m = _ccs.m;
        let log_m = log2(m as f64) as usize;
        let Beta = get_challenge_vector(_cm_i, _transcript, log_m);

        let poly_info = VPAuxInfo {
            max_degree: _ccs.d + 1,
            num_variables: log_m,
            phantom: std::marker::PhantomData,
        };
        let protocol = SumCheckIP {
            claimed_sum: R::zero(),
            poly_info,
        };
        let verifier = SumCheckVerifier::new(protocol);

        // Verify the transcript
        let subclaim = verifier
            .verify(&_proof.linearization_sumcheck, _transcript)
            .unwrap();
        let e = eq(&Beta, &subclaim.point);
        let s = subclaim.expected_evaluation.clone();

        let should_equal_s = e * _ccs.c.iter().fold(R::zero(), |sum, c| {
            sum + *c * _proof.u.iter().fold(R::one(), |product, u_j| product * u_j)
        });

        match should_equal_s == s {
            true => {}
            false => {
                return Err(LinearizationError::SumCheckError(SumCheckFailed(
                    should_equal_s,
                    s,
                )));
            }
        }
        Ok(LCCCS {
            r_arr: subclaim.point,
            v: _proof.v,
            y: _cm_i.cm.clone(),
            u: _proof.u.clone(),
            x_w: _cm_i.x_ccs.clone(),
            h: R::one(),
        })
    }
}

fn create_u<R: OverField>(
    length: usize,
    M: &Vec<Vec<Vec<R>>>,
    r_arr: &Vec<R>,
    z_ccs: &Vec<R>,
) -> Vec<R> {
    let mut u: Vec<R> = Vec::with_capacity(length);
    M.iter().for_each(|M_i| {
        (1..M_i.len())
            .into_iter()
            .map(|i| usize_to_binary_vector::<R>(i, log2(M_i.len() as f64) as usize))
            .for_each(|b| {
                let u_i = mle_val_from_matrix(&M_i, &r_arr, &b) * mle_val_from_vector(&z_ccs, &b);
                u.push(u_i);
            })
    });
    u
}

fn create_sumcheck_polynomial<R: OverField>(
    log_m: usize,
    c: &Vec<R>,
    M: &Vec<Vec<Vec<R>>>,
    z_ccs: &Vec<R>,
    Beta: &Vec<R>,
) -> VirtualPolynomial<R> {
    let mut g = VirtualPolynomial::new(log_m);
    c.iter().for_each(|coefficient| {
        let mut mle_list: Vec<Arc<DenseMultilinearExtension<R>>> = Vec::new();
        M.iter().for_each(|matrix| {
            //Initialise MLE
            let zero_vector = usize_to_binary_vector::<R>(0, log2(matrix.len() as f64) as usize);
            let mle_z_ccs_b = mle_val_from_vector(&z_ccs, &zero_vector);
            let evaluations: Vec<R> = mle_matrix_to_val_eval_second(&matrix, &zero_vector)
                .iter()
                .map(|val| *val * mle_z_ccs_b)
                .collect();
            let mle = DenseMultilinearExtension::from_evaluations_vec(log_m, evaluations);

            let matrix_mle = (1..matrix.len())
                .into_iter()
                .map(|i| usize_to_binary_vector::<R>(i, log2(matrix.len() as f64) as usize))
                .fold(mle, |acc, b| {
                    let mle_z_ccs_b = mle_val_from_vector(&z_ccs, &b);
                    let evaluations: Vec<R> = mle_matrix_to_val_eval_second(&matrix, &b)
                        .iter()
                        .map(|val| *val * mle_z_ccs_b)
                        .collect();
                    let mle = DenseMultilinearExtension::from_evaluations_vec(log_m, evaluations);
                    acc + mle
                });
            mle_list.push(Arc::from(matrix_mle));
        });
        g.add_mle_list(mle_list, *coefficient);
    });

    // multiply by eq function
    let eq = build_eq_x_r::<R>(Beta.as_slice()).unwrap();
    g.mul_by_mle(eq, R::one());
    g
}
fn get_challenge_vector<R: OverField, T: Transcript<R>>(
    _cm_i: &CCCS<R>,
    _transcript: &mut T,
    len: usize,
) -> Vec<R> {
    (0..len)
        .map(|_| {
            let challenge = _transcript.get_big_challenge().into();
            _transcript.absorb_ring_vec(&_cm_i.cm);
            challenge
        })
        .collect()
}

fn eq<R: OverField>(b_arr: &Vec<R>, x_arr: &Vec<R>) -> R {
    assert_eq!(
        b_arr.len(),
        x_arr.len(),
        "Eq function takes two vectors of the same length!"
    );
    b_arr.iter().zip(x_arr).fold(R::one(), |ret_value, (b, x)| {
        ret_value * ((R::one() - b) * (R::one() - x) + *b * x)
    })
}

fn usize_to_binary_vector<R: OverField>(n: usize, length: usize) -> Vec<R> {
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
    bits.into_iter()
        .map(|bit| if bit == 1 { R::one() } else { R::zero() })
        .collect()
}

fn mle_val_from_vector<R: OverField>(vector: &Vec<R>, values: &Vec<R>) -> R {
    assert_eq!(values.len(), log2(vector.len() as f64) as usize);
    let mle = DenseMultilinearExtension::from_evaluations_vec(values.len(), vector.clone());
    mle.evaluate(values.as_slice()).unwrap()
}

fn mle_val_from_matrix<R: OverField>(
    matrix: &Vec<Vec<R>>,
    values_x: &Vec<R>,
    values_y: &Vec<R>,
) -> R {
    assert_eq!(values_y.len(), log2(matrix.len() as f64) as usize);
    assert_eq!(values_x.len(), log2(matrix[0].len() as f64) as usize);
    let univariate_mle_evaluations = matrix
        .into_iter()
        .map(|col| mle_val_from_vector(col, values_x))
        .collect();
    mle_val_from_vector(&univariate_mle_evaluations, values_y)
}
fn mle_matrix_to_val_eval_first<R: OverField>(matrix: &Vec<Vec<R>>, values_x: &Vec<R>) -> Vec<R> {
    assert_eq!(values_x.len(), log2(matrix[0].len() as f64) as usize);
    matrix
        .into_iter()
        .map(|col| mle_val_from_vector(col, values_x))
        .collect()
}
// Convert a bivariate MLE to a univariate MLE by evaluating the second vector
fn mle_matrix_to_val_eval_second<R: OverField>(matrix: &Vec<Vec<R>>, values_y: &Vec<R>) -> Vec<R> {
    assert_eq!(values_y.len(), log2(matrix.len() as f64) as usize);
    (0..matrix[0].len())
        .into_iter()
        .map(|i| mle_val_from_vector(&matrix.iter().map(|col| col[i]).collect(), values_y))
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
            mle_val_from_vector(&evaluation_vector, &vec![one()]),
            zero()
        );
        assert_ne!(
            mle_val_from_vector(&evaluation_vector, &vec![one()]),
            poly_ntt()
        );
        assert_eq!(
            mle_val_from_vector(&evaluation_vector, &vec![zero()]),
            poly_ntt()
        );
        assert_ne!(
            mle_val_from_vector(&evaluation_vector, &vec![zero()]),
            zero()
        );

        let evaluation_matrix = vec![vec![poly_ntt(), zero()], vec![one(), poly_ntt()]];
        assert_eq!(
            mle_val_from_matrix(&evaluation_matrix, &vec![zero()], &vec![one()]),
            one()
        );
        assert_ne!(
            mle_val_from_matrix(&evaluation_matrix, &vec![zero()], &vec![one()]),
            poly_ntt()
        );
        assert_eq!(
            mle_val_from_matrix(&evaluation_matrix, &vec![zero()], &vec![zero()]),
            poly_ntt()
        );
        assert_ne!(
            mle_val_from_matrix(&evaluation_matrix, &vec![zero()], &vec![zero()]),
            zero()
        );

        // Test the eq function
        let vector_one = vec![zero(), one(), one(), zero()];
        let vector_two = vec![zero(), one(), one(), zero()];
        let vector_three = vec![zero(), one(), one(), one()];

        assert_eq!(eq(&vector_one, &vector_two), one());
        assert_ne!(eq(&vector_one, &vector_two), zero());
        assert_eq!(eq(&vector_one, &vector_three), zero());
        assert_ne!(eq(&vector_one, &vector_three), one());

        assert_eq!(
            usize_to_binary_vector::<Pow2CyclotomicPolyRingNTT<Q, N>>(4, 8),
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
            usize_to_binary_vector::<Pow2CyclotomicPolyRingNTT<Q, N>>(5, 5),
            vec![zero(), zero(), one(), zero(), one()]
        );
        // Test the conversion of Bivariate MLE to univariate MLE by evaluating first values
        let bivariate_mle = vec![
            vec![poly_ntt(), poly_ntt(), one(), zero()],
            vec![zero(), poly_ntt(), zero(), one()],
        ];
        assert_eq!(
            mle_matrix_to_val_eval_first(&bivariate_mle, &vec![zero(), zero()]),
            vec![poly_ntt(), zero()]
        );
        assert_eq!(
            mle_matrix_to_val_eval_first(&bivariate_mle, &vec![one(), zero()]),
            vec![poly_ntt(), poly_ntt()]
        );

        // Test the conversion of Bivariate MLE to univariate MLE by evaluating second values
        assert_eq!(
            mle_matrix_to_val_eval_second(&bivariate_mle, &vec![one()]),
            vec![zero(), poly_ntt(), zero(), one()]
        );
        assert_eq!(
            mle_matrix_to_val_eval_second(&bivariate_mle, &vec![zero()]),
            vec![poly_ntt(), poly_ntt(), one(), zero()]
        );
    }
}
