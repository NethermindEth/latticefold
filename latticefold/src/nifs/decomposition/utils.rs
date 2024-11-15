use cyclotomic_rings::rings::SuitableRing;
use lattirust_linear_algebra::ops::Transpose;
use lattirust_ring::{
    balanced_decomposition::{decompose_balanced_vec, gadget_decompose, recompose},
    cyclotomic_ring::{CRT, ICRT},
};

use crate::decomposition_parameters::DecompositionParams;

/// Decompose a vector of arbitrary norm in its NTT form into DP::K vectors
/// and applies the gadget-B matrix again.
pub(super) fn decompose_big_vec_into_k_vec_and_compose_back<
    NTT: SuitableRing,
    DP: DecompositionParams,
>(
    x: Vec<NTT>,
) -> Vec<Vec<NTT>> {
    // Allow x to have length m
    let coeff_repr: Vec<NTT::CoefficientRepresentation> = ICRT::elementwise_icrt(x);

    // radix-B
    let decomposed_in_B: Vec<NTT::CoefficientRepresentation> =
        gadget_decompose(&coeff_repr, DP::B, DP::L);

    // We now have a m * l length vector
    // Each element from original vector is mapped to l-length chunk

    decompose_balanced_vec(&decomposed_in_B, DP::B_SMALL as u128, DP::K)
        // We have a k by (m*l) matrix
        .transpose()
        // We have a (m*l) by k matrix
        .into_iter()
        // We recompose to a m * k matrix
        // Where could recompose basis b horizontally to recreate the original vector
        .map(|vec| {
            vec.chunks(DP::L)
                .map(|chunk| recompose(chunk, DP::B).crt())
                .collect()
        })
        .collect()
}

/// Decompose a vector of norm B in its NTT form into DP::K small vectors.
pub(super) fn decompose_B_vec_into_k_vec<NTT: SuitableRing, DP: DecompositionParams>(
    x: &[NTT::CoefficientRepresentation],
) -> Vec<Vec<NTT::CoefficientRepresentation>> {
    decompose_balanced_vec(x, DP::B_SMALL as u128, DP::K).transpose()
}
