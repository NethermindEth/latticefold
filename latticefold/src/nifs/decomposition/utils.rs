use cyclotomic_rings::SuitableRing;
use lattirust_ring::balanced_decomposition::{
    decompose_balanced_vec, pad_and_transpose, recompose,
};

use crate::decomposition_parameters::DecompositionParams;

/// Decompose a vector of arbitrary norm in its NTT form into DP::K vectors
/// and applies the gadget-B matrix again.
pub(super) fn decompose_big_vec_into_k_vec_and_compose_back<
    NTT: SuitableRing,
    DP: DecompositionParams,
>(
    x: &[NTT],
) -> Vec<Vec<NTT>> {
    let coeff_repr: Vec<NTT::CoefficientRepresentation> = x.iter().map(|&x| x.into()).collect();

    // radix-B
    let decomposed_in_B: Vec<NTT::CoefficientRepresentation> =
        pad_and_transpose(decompose_balanced_vec(&coeff_repr, DP::B, Some(DP::L)))
            .into_iter()
            .flatten()
            .collect();

    decompose_balanced_vec(&decomposed_in_B, DP::B_SMALL as u128, Some(DP::K))
        .into_iter()
        .map(|vec| {
            vec.chunks(DP::L)
                .map(|chunk| recompose(chunk, NTT::CoefficientRepresentation::from(DP::B)).into())
                .collect()
        })
        .collect()
}

/// Decompose a vector of norm B in its NTT form into DP::K small vectors.
pub(super) fn decompose_B_vec_into_k_vec<NTT: SuitableRing, DP: DecompositionParams>(
    x: &[NTT],
) -> Vec<Vec<NTT>> {
    let coeff_repr: Vec<NTT::CoefficientRepresentation> = x.iter().map(|&x| x.into()).collect();

    decompose_balanced_vec(&coeff_repr, DP::B_SMALL as u128, Some(DP::K))
        .into_iter()
        .map(|vec| vec.into_iter().map(|x| x.into()).collect())
        .collect()
}
