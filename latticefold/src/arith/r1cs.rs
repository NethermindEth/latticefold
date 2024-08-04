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
        let Cz = mat_vec_mul(&self.B, z)?;

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

// TODO Find an example to test on
