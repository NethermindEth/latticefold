use lattirust_arithmetic::ring::Ring;
use crate::utils::{ self, mat_by_vec };
use crate::utils::hadamard_vec;
use crate::NotSatisfiedError;

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct R1CS<R: Ring> {
    pub l: usize,
    pub A: Vec<Vec<R>>,
    pub B: Vec<Vec<R>>,
    pub C: Vec<Vec<R>>,
}

impl<R: Ring> R1CS<R> {
    // returns a tuple containing (w, x) (witness and public inputs respectively)
    pub fn split_z(&self, z: &[R]) -> (Vec<R>, Vec<R>) {
        (z[self.l + 1..].to_vec(), z[1..self.l + 1].to_vec())
    }
    // check that a R1CS structure is satisfied by a z vector. Only for testing.
    pub fn check_relation(&self, z: &[R]) -> Result<(), NotSatisfiedError> {
        let Az = utils::mat_by_vec(&self.A, &z.to_vec());
        let Bz = utils::mat_by_vec(&self.B, &z.to_vec());

        let Cz = utils::mat_by_vec(&self.C, &z.to_vec());
        let AzBz = utils::hadamard_vec(&Az, &Bz);

        if AzBz != Cz {
            Err(NotSatisfiedError)
        } else {
            Ok(())
        }
    }
    // converts the R1CS instance into a RelaxedR1CS as described in
    // [Nova](https://eprint.iacr.org/2021/370.pdf) section 4.1.
    pub fn relax(self) -> RelaxedR1CS<R> {
        RelaxedR1CS::<R> {
            l: self.l,
            E: vec![R::zero(); self.A.len()],
            A: self.A,
            B: self.B,
            C: self.C,
            u: R::one(),
        }
    }
}
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct RelaxedR1CS<R: Ring> {
    pub l: usize, // io len
    pub A: Vec<Vec<R>>,
    pub B: Vec<Vec<R>>,
    pub C: Vec<Vec<R>>,
    pub u: R,
    pub E: Vec<R>,
}

impl<R: Ring> RelaxedR1CS<R> {
    /// check that a RelaxedR1CS structure is satisfied by a z vector.
    pub fn check_relation(&self, z: &[R]) -> Result<(), NotSatisfiedError> {
        let Az = mat_by_vec(&self.A, &z.to_vec());
        let Bz = mat_by_vec(&self.B, &z.to_vec());
        let Cz = mat_by_vec(&self.B, &z.to_vec());

        let uCz = utils::vec_value_mul(&Cz, &self.u);
        let uCzE = utils::vec_add(&uCz, &self.E);
        let AzBz = utils::hadamard_vec(&Az, &Bz);
        if AzBz != uCzE {
            Err(NotSatisfiedError)
        } else {
            Ok(())
        }
    }
}

// TODO Find an example to test on
