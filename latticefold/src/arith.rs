//! The arith module provides utility for operating with arithmetic constraint systems.

#![allow(non_snake_case)]

use core::mem;

use ark_ff::Field;
use ark_std::log2;
use cyclotomic_rings::rings::SuitableRing;
use r1cs::R1CS;
use stark_rings::{
    balanced_decomposition::{GadgetDecompose, GadgetRecompose},
    cyclotomic_ring::{CRT, ICRT},
    PolyRing, Ring,
};
use stark_rings_linalg::SparseMatrix;
use stark_rings_poly::mle::DenseMultilinearExtension;

use self::{
    error::CSError as Error,
    utils::{hadamard, mat_vec_mul, vec_add, vec_scalar_mul},
};
use crate::{
    ark_base::*,
    commitment::{AjtaiCommitmentScheme, Commitment, CommitmentError},
    decomposition_parameters::DecompositionParams,
};

pub mod ccs;
pub mod error;
pub mod r1cs;
pub mod utils;

/// A trait for defining the behaviour of an arithmetic constraint system.
///
/// ## Type Parameters
///
///  * `R: Ring` - the ring algebra over which the constraint system operates
pub trait Arith<R: Ring> {
    /// Checks that the given Arith structure is satisfied by a z vector. Used only for testing.
    fn check_relation(&self, z: &[R]) -> Result<(), Error>;

    /// Returns the bytes that represent the parameters, that is, the matrices sizes, the amount of
    /// public inputs, etc, without the matrices/polynomials values.
    fn params_to_le_bytes(&self) -> Vec<u8>;
}

/// CCS represents the Customizable Constraint Systems structure defined in
/// the [CCS paper](https://eprint.iacr.org/2023/552)
#[derive(Debug, Clone, PartialEq)]
pub struct CCS<R: Ring> {
    /// m: number of rows in M_i (such that M_i \in F^{m, n})
    pub m: usize,
    /// n = |z|, number of cols in M_i
    pub n: usize,
    /// l = |io|, size of public input/output
    pub l: usize,
    /// t = |M|, number of matrices
    pub t: usize,
    /// q = |c| = |S|, number of multisets
    pub q: usize,
    /// d: max degree in each variable
    pub d: usize,
    /// s = log(m), dimension of x
    pub s: usize,
    /// s_prime = log(n), dimension of y
    pub s_prime: usize,
    /// vector of matrices
    pub M: Vec<SparseMatrix<R>>,
    /// vector of multisets
    pub S: Vec<Vec<usize>>,
    /// vector of coefficients
    pub c: Vec<R>,
}

impl<R: Ring> Arith<R> for CCS<R> {
    /// check that a CCS structure is satisfied by a z vector. Only for testing.
    fn check_relation(&self, z: &[R]) -> Result<(), Error> {
        let mut result: Vec<R> = vec![R::zero(); self.m];

        for i in 0..self.q {
            // extract the needed M_j matrices out of S_i
            let vec_M_j: Vec<&SparseMatrix<R>> = self.S[i].iter().map(|j| &self.M[*j]).collect();

            // complete the hadamard chain
            let mut hadamard_result = vec![R::one(); self.m];
            for M_j in vec_M_j.into_iter() {
                let mut res = mat_vec_mul(M_j, z)?;
                res.resize(self.m, R::ZERO);
                hadamard_result = hadamard(&hadamard_result, &res)?;
            }

            // multiply by the coefficient of this step
            let c_M_j_z = vec_scalar_mul(&hadamard_result, &self.c[i]);

            // add it to the final vector
            result = vec_add(&result, &c_M_j_z)?;
        }

        // make sure the final vector is all zeroes
        result
            .iter()
            .all(|item| item.is_zero())
            .then_some(())
            .ok_or(Error::NotSatisfied)
    }

    fn params_to_le_bytes(&self) -> Vec<u8> {
        [
            self.l.to_le_bytes(),
            self.m.to_le_bytes(),
            self.n.to_le_bytes(),
            self.t.to_le_bytes(),
            self.q.to_le_bytes(),
            self.d.to_le_bytes(),
        ]
        .concat()
    }
}

impl<R: Ring> CCS<R> {
    /// Constructs a [`CCS`] instance from a [`R1CS`].
    pub fn from_r1cs(r1cs: R1CS<R>, W: usize) -> Self {
        let m = W;
        let n = r1cs.A.ncols();

        CCS {
            m,
            n,
            l: r1cs.l,
            s: log2(m) as usize,
            s_prime: log2(n) as usize,
            t: 3,
            q: 2,
            d: 2,

            S: vec![vec![0, 1], vec![2]],
            c: vec![R::one(), R::one().neg()],
            M: vec![r1cs.A, r1cs.B, r1cs.C],
        }
    }

    /// Constructs a [`CCS`] instance from a [`R1CS`]. The CCS instance matrices rows are then padded
    /// to max { `(n - l - 1) * L`, `m` } next power of 2.
    pub fn from_r1cs_padded(r1cs: R1CS<R>, W: usize, L: usize) -> Self {
        let mut ccs = Self::from_r1cs(r1cs, W);
        let len = usize::max((ccs.n - ccs.l - 1) * L, ccs.m).next_power_of_two();
        ccs.pad_rows_to(len);
        ccs
    }

    ///Constructs a [`R1CS`] instance from a [`CCS`] instance.
    // TODO check that this is actually valid R1CS and make this return a result
    pub fn to_r1cs(self) -> R1CS<R> {
        R1CS::<R> {
            l: self.l,
            A: self.M[0].clone(),
            B: self.M[1].clone(),
            C: self.M[2].clone(),
        }
    }

    fn pad_rows_to(&mut self, size: usize) {
        let size = size.next_power_of_two();
        if size > self.m {
            self.m = size;
            self.s = log2(size) as usize;

            // Update matrices
            self.M.iter_mut().for_each(|mat| mat.pad_rows(size));
        }
    }
}

/// A representation of a CCS witness commitment and statement.
///
/// # Type Parameters
/// - `C`: The length of the commitment vector.
/// - `R`: The ring in which the CCS is operating.
///
#[derive(Debug, Clone, PartialEq)]
pub struct CCCS<const C: usize, R: Ring> {
    /// A commitment to the B-decomposed CCS witness.
    pub cm: Commitment<C, R>,
    /// A CCS statement
    pub x_ccs: Vec<R>,
}

/// A representation a linearized CCS witness commitment and statement.
///
/// # Type Parameters
/// - `C`: The length of the commitment vector.
/// - `R`: The ring in which the CCS is operating.
///
#[derive(Debug, Clone, PartialEq)]
pub struct LCCCS<const C: usize, R: Ring> {
    /// The linearization sumcheck challenge vector
    pub r: Vec<R>,
    /// The evaluation of the linearized CCS commitment at `r`.
    pub v: Vec<R>,
    /// The commitment to the B-decomposed CCS witness.
    pub cm: Commitment<C, R>,
    /// The evaluation of the MLEs of $\\{ M_j \mathbf{z} \mid j = 1, 2, \dots, t \\}$ at `r`.
    pub u: Vec<R>,
    /// The CCS statement
    pub x_w: Vec<R>,
    /// Constant term of the z-vector
    pub h: R,
}

/// A representation of a CCS witness.
///
/// # Type Parameters
/// - `NTT`: The ring in which the CCS is operating.
///
#[derive(Debug, Clone, PartialEq)]
pub struct Witness<NTT: SuitableRing> {
    /// `w_ccs` is the original CCS witness.
    pub w_ccs: Vec<NTT>,
    /// `f` is B-decomposed CCS witness in NTT form
    pub f: Vec<NTT>,
    /// `f_coeff` is a gadget-decomposed witness slice in the coefficient form.
    pub f_coeff: Vec<NTT::CoefficientRepresentation>,
    /// See full description of f_hat [here](crate::arith::Witness::get_fhat).
    pub f_hat: Vec<DenseMultilinearExtension<NTT>>,
}

impl<NTT: SuitableRing> Witness<NTT> {
    /// Create a [`Witness`] from a ccs witness.
    ///
    /// The main operations that need to be done are decomposing the ccs witness.
    /// We can then construct [`f_hat`](crate::arith::Witness::get_fhat).
    pub fn from_w_ccs<P: DecompositionParams>(w_ccs: Vec<NTT>) -> Self {
        // iNTT
        let w_coeff: Vec<NTT::CoefficientRepresentation> = ICRT::elementwise_icrt(w_ccs.clone());

        // decompose radix-B
        let f_coeff: Vec<NTT::CoefficientRepresentation> = w_coeff.gadget_decompose(P::B, P::L);

        // NTT(coef_repr_decomposed)
        let f: Vec<NTT> = CRT::elementwise_crt(f_coeff.clone());
        // coef_repr_decomposed -> coefs -> NTT = coeffs.
        let f_hat: Vec<DenseMultilinearExtension<NTT>> = Self::get_fhat(&f_coeff);

        Self {
            f,
            f_coeff,
            f_hat,
            w_ccs,
        }
    }

    /// Given a gadget-decomposed witness slice `f` returns the f-hat matrix from the Latticefold paper,
    /// i.e. a matrix of dimension `tau x f.len()`, where `tau = NTT::CoefficientRepresentation::dimension() / NTT::dimension()`,
    /// such that the `j`th row of the matrix is obtained as (in pseudocode)
    /// ```text
    ///   [
    ///      NTT(f[0][NTT::dimension() * j .. NTT::dimension() * (j+1)]),
    ///      NTT(f[1][NTT::dimension() * j .. NTT::dimension() * (j+1)]),
    ///      ...,
    ///      NTT(f[f.len() - 1][NTT::dimension() * j .. NTT::dimension() * (j+1)])
    ///   ].
    /// ```
    /// Note, that our definition of the hat-matrix is equivalent up to transposing to the
    /// one in the Latticefold paper, since it is more often the case that we need to iterate
    /// over the columns (of the original latticefold f-hat) rather than the rows.
    ///
    /// # Arguments
    ///
    /// * `f` - A gadget-decomposed witness slice in the coefficient form.
    ///
    /// # Returns
    ///
    /// The hat matrix `Vec<Vec<NTT>>`.
    ///
    fn get_fhat(f: &[NTT::CoefficientRepresentation]) -> Vec<DenseMultilinearExtension<NTT>> {
        let num_vars = log2(f.len().next_power_of_two()) as usize;

        let mut fhat =
            vec![
                DenseMultilinearExtension::from_evaluations_vec(num_vars, Vec::<NTT>::new());
                NTT::CoefficientRepresentation::dimension() / NTT::dimension()
            ];

        for (i, f_i) in f.iter().enumerate() {
            for (j, coeff_chunk_j) in f_i.coeffs().chunks(NTT::dimension()).enumerate() {
                for (&coeff_from_chunk, fhat_j_i_coeff) in
                    coeff_chunk_j.iter().zip(fhat[j][i].coeffs_mut().iter_mut())
                {
                    *fhat_j_i_coeff =
                        <NTT::BaseRing as Field>::from_base_prime_field(coeff_from_chunk);
                }
            }
        }

        // Make any possible MLEs smaller
        fhat.iter_mut().for_each(|mle| mle.truncate_lnze());

        fhat
    }

    pub(crate) fn from_f<P: DecompositionParams>(f: Vec<NTT>) -> Self {
        let f_coeff: Vec<NTT::CoefficientRepresentation> = ICRT::elementwise_icrt(f.clone());
        let f_hat: Vec<DenseMultilinearExtension<NTT>> = Self::get_fhat(&f_coeff);
        // Reconstruct the original CCS witness from the Ajtai witness
        // Ajtai witness has bound B
        // WE multiply by the base B gadget matrix to reconstruct w_ccs
        let w_ccs = f.gadget_recompose(P::B, P::L);

        Self {
            f,
            f_coeff,
            f_hat,
            w_ccs,
        }
    }

    #[allow(dead_code)]
    fn from_f_slice<P: DecompositionParams>(f: &[NTT]) -> Self {
        Self::from_f::<P>(f.into())
    }

    /// Reconstruct the original CCS witness from the Ajtai witness
    ///
    /// Assume that Ajtai witness has bound B.
    /// We can multiply by the base B gadget matrix to reconstruct w_ccs.
    pub fn from_f_coeff<P: DecompositionParams>(
        f_coeff: Vec<NTT::CoefficientRepresentation>,
    ) -> Self {
        let f: Vec<NTT> = CRT::elementwise_crt(f_coeff.clone());
        let f_hat: Vec<DenseMultilinearExtension<NTT>> = Self::get_fhat(&f_coeff);

        let w_ccs = f.gadget_recompose(P::B, P::L);

        Self {
            f,
            f_coeff,
            f_hat,
            w_ccs,
        }
    }

    /// Generates a random witness by firstly generating a random
    /// vector of arbitrary norm and then computing the rest of the data
    /// needed for a witness.
    ///
    /// # Arguments
    /// * `rng` is a mutable reference to the random number generator.
    /// * `w_ccs_len` is the length of the non-decomposed witness (a.k.a. the CCS witness).
    pub fn rand<Rng: rand::Rng + ?Sized, P: DecompositionParams>(
        rng: &mut Rng,
        w_ccs_len: usize,
    ) -> Self {
        Self::from_w_ccs::<P>((0..w_ccs_len).map(|_| NTT::rand(rng)).collect())
    }

    /// Produces a commitment from a witness
    ///
    /// Ajtai commitments are produced by multiplying an Ajtai matrix by the witness vector
    pub fn commit<const C: usize, const W: usize, P: DecompositionParams>(
        &self,
        ajtai: &AjtaiCommitmentScheme<C, W, NTT>,
    ) -> Result<Commitment<C, NTT>, CommitmentError> {
        ajtai.commit_ntt(&self.f)
    }

    /// Takes the `f_hat` value.
    ///
    /// Leaves the value in the struct as `None`.
    pub fn take_f_hat(&mut self) -> Vec<DenseMultilinearExtension<NTT>> {
        mem::take(&mut self.f_hat)
    }

    #[allow(dead_code)]
    fn within_bound(&self, b: u128) -> bool {
        // TODO consider using signed representatives instead

        let coeffs_repr: Vec<NTT::CoefficientRepresentation> =
            ICRT::elementwise_icrt(self.f.clone());

        // linf_norm should be used in CyclotomicGeneral not in specific ring
        let b = <<NTT as PolyRing>::BaseRing as Field>::BasePrimeField::from(b);
        let all_under_bound = coeffs_repr.iter().all(|ele| {
            let coeffs = ele.coeffs();
            coeffs.iter().all(|x| x < &b)
        });

        all_under_bound
    }
}

/// A trait for defining the behaviour of a satisfying instance of a constraint system
///
/// # Types
///  - `R: Ring` - the ring in which the constraint system is operating.
///
pub trait Instance<R: Ring> {
    /// Given a witness vector, produce a concatonation of the statement and the witness
    fn get_z_vector(&self, w: &[R]) -> Vec<R>;
}

impl<const C: usize, R: Ring> Instance<R> for CCCS<C, R> {
    fn get_z_vector(&self, w: &[R]) -> Vec<R> {
        let mut z: Vec<R> = Vec::with_capacity(self.x_ccs.len() + w.len() + 1);

        z.extend_from_slice(&self.x_ccs);
        z.push(R::one());
        z.extend_from_slice(w);

        z
    }
}

impl<const C: usize, R: Ring> Instance<R> for LCCCS<C, R> {
    fn get_z_vector(&self, w: &[R]) -> Vec<R> {
        let mut z: Vec<R> = Vec::with_capacity(self.x_w.len() + w.len() + 1);

        z.extend_from_slice(&self.x_w);
        z.push(self.h);
        z.extend_from_slice(w);

        z
    }
}

#[cfg(test)]
pub mod tests {
    use ark_ff::{One, Zero};
    use cyclotomic_rings::rings::{
        BabyBearRingNTT, GoldilocksRingNTT, GoldilocksRingPoly, StarkRingNTT,
    };
    use stark_rings::cyclotomic_ring::models::goldilocks::{Fq, Fq3};

    use super::*;
    use crate::{
        arith::r1cs::{get_test_r1cs, get_test_z as r1cs_get_test_z},
        decomposition_parameters::test_params::{BabyBearDP, GoldilocksDP, StarkDP},
    };

    pub(crate) fn get_test_ccs<R: Ring>(W: usize, L: usize) -> CCS<R> {
        let r1cs = get_test_r1cs::<R>();
        CCS::<R>::from_r1cs_padded(r1cs, W, L)
    }

    pub(crate) fn get_test_z<R: Ring>(input: usize) -> Vec<R> {
        r1cs_get_test_z(input)
    }

    /// Test that a basic CCS relation is satisfied by a witness
    #[test]
    fn test_ccs_relation() {
        let ccs = get_test_ccs::<BabyBearRingNTT>(4, 1);
        let z = get_test_z(3);

        ccs.check_relation(&z).unwrap();
    }

    #[test]
    fn test_get_fhat() {
        let mut f_1_coeffs = vec![Fq::from(1), Fq::from(2), Fq::from(3)];
        let mut f_2_coeffs = vec![Fq::from(4), Fq::from(5), Fq::from(6)];

        f_1_coeffs.extend((0..GoldilocksRingPoly::dimension() - 3).map(|_| Fq::zero()));
        f_2_coeffs.extend((0..GoldilocksRingPoly::dimension() - 3).map(|_| Fq::one()));

        let f: Vec<GoldilocksRingPoly> = vec![
            GoldilocksRingPoly::from(f_1_coeffs),
            GoldilocksRingPoly::from(f_2_coeffs),
        ];

        let fhat: Vec<Vec<_>> = Witness::<GoldilocksRingNTT>::get_fhat(&f)
            .into_iter()
            .map(|f_hat_mle| f_hat_mle.evaluations)
            .collect();

        assert_eq!(
            fhat,
            vec![
                vec![
                    GoldilocksRingNTT::from(vec![
                        Fq3::from_base_prime_field(Fq::from(1)),
                        Fq3::from_base_prime_field(Fq::from(2)),
                        Fq3::from_base_prime_field(Fq::from(3)),
                        Fq3::zero(),
                        Fq3::zero(),
                        Fq3::zero(),
                        Fq3::zero(),
                        Fq3::zero(),
                    ]),
                    GoldilocksRingNTT::from(vec![
                        Fq3::from_base_prime_field(Fq::from(4)),
                        Fq3::from_base_prime_field(Fq::from(5)),
                        Fq3::from_base_prime_field(Fq::from(6)),
                        Fq3::one(),
                        Fq3::one(),
                        Fq3::one(),
                        Fq3::one(),
                        Fq3::one(),
                    ]),
                ],
                vec![GoldilocksRingNTT::zero(), GoldilocksRingNTT::one()],
                vec![GoldilocksRingNTT::zero(), GoldilocksRingNTT::one()]
            ]
        );
    }

    impl<NTT: SuitableRing> Witness<NTT> {
        fn check_data<P: DecompositionParams>(&self) -> bool {
            let w_coeff = ICRT::elementwise_icrt(self.w_ccs.clone());

            (self.f_coeff == w_coeff.gadget_decompose(P::B, P::L))
                && (CRT::elementwise_crt(self.f_coeff.clone()) == self.f)
                && (self.f_hat == Self::get_fhat(&self.f_coeff))
        }
    }

    const WIT_LEN: usize = 1024;

    #[test]
    fn test_from_w_ccs() {
        let mut rng = ark_std::test_rng();

        let random_witness =
            Witness::<GoldilocksRingNTT>::rand::<_, GoldilocksDP>(&mut rng, WIT_LEN);
        let recreated_witness = Witness::from_w_ccs::<GoldilocksDP>(random_witness.w_ccs.clone());

        assert!(recreated_witness.check_data::<GoldilocksDP>());
        assert_eq!(recreated_witness, random_witness);
    }

    #[test]
    fn test_from_f() {
        let mut rng = ark_std::test_rng();

        let random_witness = Witness::<BabyBearRingNTT>::rand::<_, BabyBearDP>(&mut rng, WIT_LEN);
        let recreated_witness = Witness::from_f::<BabyBearDP>(random_witness.f.clone());

        assert!(recreated_witness.check_data::<BabyBearDP>());
        assert_eq!(recreated_witness, random_witness);
    }

    #[test]
    fn test_from_f_coeff() {
        let mut rng = ark_std::test_rng();

        let random_witness = Witness::<StarkRingNTT>::rand::<_, StarkDP>(&mut rng, WIT_LEN);
        let recreated_witness = Witness::from_f_coeff::<StarkDP>(random_witness.f_coeff.clone());

        assert!(recreated_witness.check_data::<StarkDP>());
        assert_eq!(recreated_witness, random_witness);
    }
}
