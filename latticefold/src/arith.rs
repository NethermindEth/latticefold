#![allow(non_snake_case)]

pub mod error;
pub mod r1cs;
pub mod utils;

use ark_std::log2;
use error::NotSatisfiedError;
use lattirust_arithmetic::ring::Ring;
use r1cs::R1CS;
use utils::hadamard_vec;

/// CCS represents the Customizable Constraint Systems structure defined in
/// the [CCS paper](https://eprint.iacr.org/2023/552)
#[derive(Debug, Clone, Eq, PartialEq)]
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
    // TODO Make sure this is right, d in paper is maximum size of multiset
    pub d: usize,
    /// s = log(m), dimension of x
    pub s: usize,
    /// s_prime = log(n), dimension of y
    pub s_prime: usize,

    /// vector of matrices
    pub M: Vec<Vec<Vec<R>>>,
    /// vector of multisets
    pub S: Vec<Vec<usize>>,
    /// vector of coefficients
    pub c: Vec<R>,
}

impl<R: Ring> CCS<R> {
    /// check that a CCS structure is satisfied by a z vector. Only for testing.
    pub fn check_relation(&self, z: &[R]) -> Result<(), NotSatisfiedError> {
        let mut result = vec![R::zero(); self.m];
        //  Calculates \sum_{i=1}^{n_r} c_i \cdot \bigg( \bigcirc_{j \in S_i} (M_j \cdot \vec{z}) \bigg)
        for i in 0..self.q {
            // Extract the needed M_j matrices out of S_i
            let vec_M_j: Vec<&Vec<Vec<R>>> = self.S[i].iter().map(|&j| &self.M[j]).collect();

            // complete the hadamard chain
            let mut hadamard_result = vec![R::one(); self.m];
            for M_j in vec_M_j.into_iter() {
                let M_j_z: Vec<R> = utils::mat_by_vec(M_j, z);

                hadamard_result = hadamard_vec(&hadamard_result, &M_j_z);
            }

            // multiply by the coefficient of this step
            let c_M_j_z = utils::vec_value_mul(&hadamard_result, &self.c[i]);

            // Add it to the final vector
            result = utils::vec_add(&result, &c_M_j_z);
        }

        if result.iter().all(|x| x == &R::zero()) {
            Ok(())
        } else {
            Err(NotSatisfiedError)
        }
    }
}

// CCS is a generalisation of R1CS, so we can convert from R1CS to CSS
impl<R: Ring> CCS<R> {
    pub fn from_r1cs(r1cs: R1CS<R>) -> Self {
        let m = r1cs.A.len();
        let n = r1cs.A[0].len();
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
pub struct CCCS<R: Ring> {
    pub cm: Vec<R>,
    pub x_ccs: Vec<R>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct LCCCS<R: Ring> {
    pub r_arr: Vec<R>,
    pub v: R,
    pub y: Vec<R>,
    pub u: Vec<R>,
    pub x_w: Vec<R>,
    pub h: R,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Witness<R: Ring> {
    pub f_arr: Vec<R>,
    // NTT(f_hat) = f_arr
    pub f_hat: Vec<R>,
    pub w_ccs: Vec<R>,
}
