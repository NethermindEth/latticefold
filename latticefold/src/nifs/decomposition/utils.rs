use cyclotomic_rings::rings::SuitableRing;
use lattirust_linear_algebra::ops::Transpose;
use lattirust_ring::balanced_decomposition::{decompose_balanced_vec, recompose};

use crate::decomposition_parameters::DecompositionParams;

use ark_std::ops::MulAssign;
/// Decompose a vector of arbitrary norm in its NTT form into DP::K vectors
/// and applies the gadget-B matrix again.
pub(super) fn decompose_big_vec_into_k_vec_and_compose_back<
    NTT: SuitableRing,
    DP: DecompositionParams,
>(
    x: &[NTT],
) -> Vec<Vec<NTT>>
where
    for<'a> <NTT as SuitableRing>::CoefficientRepresentation: MulAssign<&'a u128>,
{
    let coeff_repr: Vec<NTT::CoefficientRepresentation> = x.iter().map(|&x| x.into()).collect();

    // radix-B
    let decomposed_in_B: Vec<NTT::CoefficientRepresentation> =
        decompose_balanced_vec(&coeff_repr, DP::B, Some(DP::L))
            .into_iter()
            .flatten()
            .collect();

    decompose_balanced_vec(&decomposed_in_B, DP::B_SMALL as u128, Some(DP::K))
        .transpose()
        .into_iter()
        .map(|vec| {
            vec.chunks(DP::L)
                .map(|chunk| recompose(chunk, DP::B).into())
                .collect()
        })
        .collect()
}

/// Decompose a vector of norm B in its NTT form into DP::K small vectors.
pub(super) fn decompose_B_vec_into_k_vec<NTT: SuitableRing, DP: DecompositionParams>(
    x: &[NTT::CoefficientRepresentation],
) -> Vec<Vec<NTT::CoefficientRepresentation>> {
    decompose_balanced_vec(x, DP::B_SMALL as u128, Some(DP::K)).transpose()
}
