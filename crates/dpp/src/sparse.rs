//! Sparse vectors and dot-products for efficient DPP frontends.

use ark_ff::PrimeField;
use rayon::prelude::*;

/// Sparse vector: list of (coefficient, index) terms.
#[derive(Clone, Debug, Default)]
pub struct SparseVec<F: PrimeField> {
    pub terms: Vec<(F, usize)>,
}

impl<F: PrimeField> SparseVec<F> {
    pub fn new(terms: Vec<(F, usize)>) -> Self {
        Self { terms }
    }

    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }

    pub fn dot(&self, v: &[F]) -> F {
        const PAR_SPARSE_DOT_MIN_TERMS: usize = 2048;
        if self.terms.len() >= PAR_SPARSE_DOT_MIN_TERMS {
            self.terms
                .par_iter()
                .map(|(c, idx)| *c * v[*idx])
                .reduce(|| F::ZERO, |acc, t| acc + t)
        } else {
            self.terms
                .iter()
                .fold(F::ZERO, |acc, (c, idx)| acc + (*c * v[*idx]))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::Field;
    use ark_ff::{Fp64, MontBackend, MontConfig};

    #[derive(MontConfig)]
    #[modulus = "18446744069414584321"]
    #[generator = "7"]
    pub struct FBigConfig;
    type FBig = Fp64<MontBackend<FBigConfig, 1>>;

    #[test]
    fn test_sparse_dot_parallel_matches_sequential() {
        // Force the parallel path by using many terms.
        let n = 10_000usize;
        let v = (0..n).map(|i| FBig::from(i as u64)).collect::<Vec<_>>();
        let terms = (0..5000usize).map(|i| (FBig::from((i + 1) as u64), i)).collect::<Vec<_>>();
        let sv = SparseVec::new(terms);

        let seq = sv
            .terms
            .iter()
            .fold(FBig::ZERO, |acc, (c, idx)| acc + (*c * v[*idx]));
        let par = sv.dot(&v);
        assert_eq!(par, seq);
    }
}

