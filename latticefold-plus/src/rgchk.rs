use ark_std::log2;
use latticefold::{
    transcript::Transcript,
    utils::sumcheck::{
        utils::{build_eq_x_r, eq_eval},
        MLSumcheck, Proof, SumCheckError,
    },
};
use stark_rings::{
    balanced_decomposition::{convertible_ring::ConvertibleRing, DecomposeToVec, GadgetDecompose},
    exp, CoeffRing, OverField, PolyRing, Ring,
};
use stark_rings_linalg::{ops::Transpose, Matrix, SparseMatrix};
use stark_rings_poly::mle::{DenseMultilinearExtension, SparseMultilinearExtension};
use thiserror::Error;

use crate::setchk::{Arity, In, Out};

// D_f: decomposed cf(f), Z n x dk
// M_f: EXP(D_f)

#[derive(Debug)]
pub struct Rg<R> {
    pub nvars: usize,
    pub M_f: Vec<Matrix<R>>, // n x d, k matrices
    pub f: Vec<R>,           // n

    // decomposition
    pub b: u128,
    pub k: u128,
}

#[derive(Debug)]
pub struct Dcom<R> {
    cm_f: Vec<R>,
}

impl<R: CoeffRing> Rg<R>
where
    R::Coeff: ConvertibleRing,
{
    pub fn range_check(&self, transcript: &mut impl Transcript<R>) {
        let outs = self
            .M_f
            .iter()
            .map(|m| {
                In {
                    nvars: self.nvars,
                    M: Arity::Single(SparseMatrix::from_dense(m)),
                }
                .set_check(transcript)
            })
            .collect::<Vec<_>>();

        let cfs = self
            .f
            .iter()
            .map(|r| r.coeffs().to_vec())
            .collect::<Vec<_>>()
            .transpose();
        let v = cfs
            .into_iter()
            .map(|evals| DenseMultilinearExtension::from_evaluations_vec(self.nvars, evals))
            .collect::<Vec<_>>();
        println!("{:?}", v);
    }
}

#[cfg(test)]
mod tests {
    use ark_std::{One, Zero};
    use cyclotomic_rings::rings::GoldilocksPoseidonConfig as PC;
    use latticefold::transcript::poseidon::PoseidonTS;
    use stark_rings::{
        balanced_decomposition::DecomposeToVec, cyclotomic_ring::models::goldilocks::RqPoly as R,
    };
    use stark_rings_linalg::SparseMatrix;

    use super::*;

    #[test]
    fn test_range_check() {
        // f: [
        // 2 + 5X
        // 4 + X^2
        // ]
        let mut f = vec![R::zero(); 2];
        f[0].coeffs_mut()[0] = 2u128.into();
        f[0].coeffs_mut()[1] = 5u128.into();
        f[1].coeffs_mut()[0] = 4u128.into();
        f[1].coeffs_mut()[2] = 1u128.into();

        let n = f.len();
        let b = 2;
        let k = 4;

        let cfs: Matrix<_> = f
            .iter()
            .map(|r| r.coeffs().to_vec())
            .collect::<Vec<Vec<_>>>()
            .into();
        let dec = cfs
            .vals
            .iter()
            .map(|row| row.decompose_to_vec(b, k as usize))
            .collect::<Vec<_>>();

        let mut D_f = vec![Matrix::zero(n, R::dimension()); k as usize];

        // map dec: (Z n x d x k) to D_f: (Z n x d, k matrices)
        dec.iter().enumerate().for_each(|(n_i, drow)| {
            drow.iter().enumerate().for_each(|(d_i, coeffs)| {
                coeffs.iter().enumerate().for_each(|(k_i, coeff)| {
                    D_f[k_i].vals[n_i][d_i] = *coeff;
                });
            });
        });

        let M_f: Vec<Matrix<R>> = D_f
            .iter()
            .map(|m| {
                m.vals
                    .iter()
                    .map(|row| {
                        row.iter()
                            .map(|c| exp::<R>(*c).unwrap())
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>()
                    .into()
            })
            .collect::<Vec<_>>();

        let rg = Rg {
            nvars: log2(n) as usize,
            M_f,
            f,
            b,
            k,
        };

        let mut ts = PoseidonTS::default::<PC>();

        rg.range_check(&mut ts);
    }
}
