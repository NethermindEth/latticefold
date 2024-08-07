use super::error::CSError as Error;
use super::utils::mat_vec_mul;
use super::utils::vec_add;
use super::utils::vec_scalar_mul;
use lattirust_arithmetic::linear_algebra::SparseMatrix;
use lattirust_arithmetic::ring::Ring;

use crate::arith::hadamard;

#[derive(Debug, Clone, PartialEq)]
pub struct R1CS<R: Ring> {
    pub l: usize,
    pub A: SparseMatrix<R>,
    pub B: SparseMatrix<R>,
    pub C: SparseMatrix<R>,
}

impl<R: Ring> R1CS<R> {
    // returns a tuple containing (w, x) (witness and public inputs respectively)
    pub fn split_z(&self, z: &[R]) -> (Vec<R>, Vec<R>) {
        (z[self.l + 1..].to_vec(), z[1..self.l + 1].to_vec())
    }
    // check that a R1CS structure is satisfied by a z vector. Only for testing.
    pub fn check_relation(&self, z: &[R]) -> Result<(), Error> {
        let Az = mat_vec_mul(&self.A, z)?;
        let Bz = mat_vec_mul(&self.B, z)?;

        let Cz = mat_vec_mul(&self.C, z)?;
        let AzBz = hadamard(&Az, &Bz)?;

        if AzBz != Cz {
            Err(Error::NotSatisfied)
        } else {
            Ok(())
        }
    }
    // converts the R1CS instance into a RelaxedR1CS as described in
    // [Nova](https://eprint.iacr.org/2021/370.pdf) section 4.1.
    pub fn relax(self) -> RelaxedR1CS<R> {
        RelaxedR1CS::<R> {
            l: self.l,
            E: vec![R::zero(); self.A.nrows()],
            A: self.A,
            B: self.B,
            C: self.C,
            u: R::one(),
        }
    }
}
#[derive(Debug, Clone, PartialEq)]
pub struct RelaxedR1CS<R: Ring> {
    pub l: usize, // io len
    pub A: SparseMatrix<R>,
    pub B: SparseMatrix<R>,
    pub C: SparseMatrix<R>,
    pub u: R,
    pub E: Vec<R>,
}

impl<R: Ring> RelaxedR1CS<R> {
    /// check that a RelaxedR1CS structure is satisfied by a z vector.
    pub fn check_relation(&self, z: &[R]) -> Result<(), Error> {
        let Az = mat_vec_mul(&self.A, z)?;
        let Bz = mat_vec_mul(&self.B, z)?;
        let Cz = mat_vec_mul(&self.C, z)?;

        let uCz = vec_scalar_mul(&Cz, &self.u);
        let uCzE = vec_add(&uCz, &self.E)?;
        let AzBz = hadamard(&Az, &Bz)?;
        if AzBz != uCzE {
            Err(Error::NotSatisfied)
        } else {
            Ok(())
        }
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;
    use lattirust_arithmetic::ring::{Pow2CyclotomicPolyRingNTT, Ring};

    pub fn to_F_matrix<R: Ring>(M: Vec<Vec<usize>>) -> SparseMatrix<R> {
        to_F_dense_matrix::<R>(M).as_slice().into()
    }

    pub fn to_F_dense_matrix<R: Ring>(M: Vec<Vec<usize>>) -> Vec<Vec<R>> {
        M.iter()
            .map(|m| m.iter().map(|r| R::from(*r as u64)).collect())
            .collect()
    }
    pub fn to_F_vec<R: Ring>(z: Vec<usize>) -> Vec<R> {
        z.iter().map(|c| R::from(*c as u64)).collect()
    }

    pub fn get_test_r1cs<R: Ring>() -> R1CS<R> {
        // R1CS for: x^3 + x + 5 = y (example from article
        // https://www.vitalik.ca/general/2016/12/10/qap.html )
        let A = to_F_matrix::<R>(vec![
            vec![1, 0, 0, 0, 0, 0],
            vec![0, 0, 0, 1, 0, 0],
            vec![1, 0, 0, 0, 1, 0],
            vec![0, 5, 0, 0, 0, 1],
        ]);
        let B = to_F_matrix::<R>(vec![
            vec![1, 0, 0, 0, 0, 0],
            vec![1, 0, 0, 0, 0, 0],
            vec![0, 1, 0, 0, 0, 0],
            vec![0, 1, 0, 0, 0, 0],
        ]);
        let C = to_F_matrix::<R>(vec![
            vec![0, 0, 0, 1, 0, 0],
            vec![0, 0, 0, 0, 1, 0],
            vec![0, 0, 0, 0, 0, 1],
            vec![0, 0, 1, 0, 0, 0],
        ]);

        R1CS::<R> { l: 1, A, B, C }
    }

    pub fn get_test_z<R: Ring>(input: usize) -> Vec<R> {
        // z = (1, io, w)
        to_F_vec(vec![
            1,
            input,                             // io
            input * input * input + input + 5, // x^3 + x + 5
            input * input,                     // x^2
            input * input * input,             // x^2 * x
            input * input * input + input,     // x^3 + x
        ])
    }

    pub fn get_test_z_split<R: Ring>(input: usize) -> (R, Vec<R>, Vec<R>) {
        // z = (1, io, w)
        (
            R::one(),
            to_F_vec(vec![
                input, // io
            ]),
            to_F_vec(vec![
                input * input * input + input + 5, // x^3 + x + 5
                input * input,                     // x^2
                input * input * input,             // x^2 * x
                input * input * input + input,     // x^3 + x
            ]),
        )
    }

    #[test]
    fn test_check_relation() {
        let r1cs = get_test_r1cs::<Pow2CyclotomicPolyRingNTT<101, 16>>();
        let z = get_test_z(5);

        r1cs.check_relation(&z).unwrap();
        r1cs.relax().check_relation(&z).unwrap();
    }
}
