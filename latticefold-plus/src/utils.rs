use ark_std::{iter::once, log2, Zero};
use latticefold::{
    transcript::Transcript,
    utils::sumcheck::{
        utils::{build_eq_x_r, eq_eval},
        MLSumcheck, Proof, SumCheckError,
    },
};
// Import the concrete ring type for testing
#[cfg(test)]
use stark_rings::cyclotomic_ring::models::goldilocks::RqPoly as TestRing;
use stark_rings::{
    balanced_decomposition::{
        convertible_ring::ConvertibleRing, Decompose, DecomposeToVec, GadgetDecompose,
    },
    exp, psi, psi_range_check, CoeffRing, OverField, PolyRing, Ring, Zq,
};
use stark_rings_linalg::{ops::Transpose, Matrix, SparseMatrix};
use stark_rings_poly::mle::{DenseMultilinearExtension, SparseMultilinearExtension};

pub fn split<R: PolyRing>(com: &Matrix<R>, n: usize, b: u128, k: usize) -> Vec<R::BaseRing>
where
    R: Decompose,
{
    let M_prime = com.gadget_decompose(b, k);
    let M_dprime = M_prime.vals.into_iter().fold(vec![], |mut acc, row| {
        // TODO pre-alloc
        acc.extend(row);
        acc
    });
    let mut tau = M_dprime
        .iter()
        .map(|r| r.coeffs().to_vec())
        .into_iter()
        .fold(vec![], |mut acc, row| {
            // TODO pre-alloc
            acc.extend(row);
            acc
        });
    if tau.len() < n {
        // TODO handle when opposite
        tau.resize(n, R::BaseRing::zero());
    } else {
        panic!(
            "small n {} unsupported, must be >= tau unpadded {}",
            n,
            tau.len()
        );
    }
    tau
}

pub fn to_bits(n: usize, num_bits: usize) -> Vec<u8> {
    (0..num_bits).map(|i| ((n >> i) & 1) as u8).collect()
}

pub fn from_bits(bits: &[u8]) -> usize {
    bits.iter()
        .rev()
        .enumerate()
        .fold(0, |acc, (i, &b)| acc | ((b as usize) << i))
}
