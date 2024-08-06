#![allow(non_snake_case)]

pub mod error;
pub mod r1cs;
pub mod utils;

use crate::arith::utils::hadamard;
use crate::arith::utils::mat_vec_mul;
use crate::arith::utils::vec_add;
use crate::arith::utils::vec_scalar_mul;
use crate::commitment::AjtaiParams;
use crate::commitment::Commitment;
use ark_std::log2;
use error::CSError as Error;
use lattirust_arithmetic::linear_algebra::SparseMatrix;
use lattirust_arithmetic::ring::Ring;
use r1cs::R1CS;

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
                hadamard_result = hadamard(&hadamard_result, &mat_vec_mul(M_j, z)?)?;
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
    pub fn from_r1cs(r1cs: R1CS<R>) -> Self {
        let m = r1cs.A.nrows();
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

    pub fn to_r1cs(self) -> R1CS<R> {
        R1CS::<R> {
            l: self.l,
            A: self.M[0].clone(),
            B: self.M[1].clone(),
            C: self.M[2].clone(),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct CCCS<R: Ring, P: AjtaiParams> {
    pub cm: Commitment<R, P>,
    pub x_ccs: Vec<R>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct LCCCS<R: Ring, P: AjtaiParams> {
    pub r: Vec<R>,
    pub v: R,
    pub cm: Commitment<R, P>,
    pub u: Vec<R>,
    pub x_w: Vec<R>,
    pub h: R,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Witness<NTT: Ring> {
    // F is B-decomposed ccs witness
    pub f: Vec<NTT>,
    // f_hat = vec(CR repr of f)
    pub f_hat: Vec<NTT>,
    pub w_ccs: Vec<NTT>,
}

pub trait Instance<R: Ring> {
    fn get_z_vector(&self, w: &[R]) -> Vec<R>;
}

impl<R: Ring, P: AjtaiParams> Instance<R> for CCCS<R, P> {
    fn get_z_vector(&self, w: &[R]) -> Vec<R> {
        let mut z: Vec<R> = Vec::with_capacity(self.x_ccs.len() + w.len() + 1);

        z.extend_from_slice(&self.x_ccs);
        z.push(R::one());
        z.extend_from_slice(w);

        z
    }
}

impl<R: Ring, P: AjtaiParams> Instance<R> for LCCCS<R, P> {
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
    use super::*;
    use crate::arith::r1cs::tests::{get_test_r1cs, get_test_z as r1cs_get_test_z};
    use lattirust_arithmetic::ring::Pow2CyclotomicPolyRingNTT;

    pub fn get_test_ccs<R: Ring>() -> CCS<R> {
        let r1cs = get_test_r1cs::<R>();
        CCS::<R>::from_r1cs(r1cs)
    }
    pub fn get_test_z<R: Ring>(input: usize) -> Vec<R> {
        r1cs_get_test_z(input)
    }

    /// Test that a basic CCS relation can be satisfied
    #[test]
    fn test_ccs_relation() {
        let ccs = get_test_ccs::<Pow2CyclotomicPolyRingNTT<101u64, 64>>();
        let z = get_test_z(3);

        ccs.check_relation(&z).unwrap();
    }
}
