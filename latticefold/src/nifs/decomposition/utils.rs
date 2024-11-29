use cyclotomic_rings::rings::SuitableRing;
use lattirust_linear_algebra::ops::Transpose;
use lattirust_ring::{
    balanced_decomposition::{decompose_balanced_vec, gadget_decompose, recompose},
    cyclotomic_ring::{CRT, ICRT},
};

use crate::ark_base::*;
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

/// Decompose a vector of norm B in its coefficient form into DP::K small vectors.
pub(super) fn decompose_B_vec_into_k_vec<NTT: SuitableRing, DP: DecompositionParams>(
    x: &[NTT::CoefficientRepresentation],
) -> Vec<Vec<NTT::CoefficientRepresentation>> {
    decompose_balanced_vec(x, DP::B_SMALL as u128, DP::K).transpose()
}

#[cfg(test)]
mod tests {
    use cyclotomic_rings::rings::SuitableRing;
    use lattirust_ring::{
        balanced_decomposition::{decompose_balanced_vec, recompose},
        cyclotomic_ring::{
            models::goldilocks::{RqNTT, RqPoly},
            CRT,
        },
        PolyRing,
    };
    use rand::Rng;

    use crate::{
        ark_base::*,
        decomposition_parameters::{test_params::DP, DecompositionParams},
        nifs::decomposition::utils::{
            decompose_B_vec_into_k_vec, decompose_big_vec_into_k_vec_and_compose_back,
        },
    };

    fn draw_ring_below_bound<RqPoly, const B: u128>(rng: &mut impl Rng) -> RqPoly
    where
        RqPoly: PolyRing + CRT,
    {
        let degree = <RqPoly as PolyRing>::dimension();
        let mut coeffs = Vec::with_capacity(degree);
        for _ in 0..degree {
            let random_coeff = rng.gen_range(0..B);
            coeffs.push(<RqPoly as PolyRing>::BaseRing::from(random_coeff));
        }
        RqPoly::from(coeffs)
    }
    fn test_decompose_B_vec_into_k_vec<RqNTT, RqPoly>()
    where
        RqNTT: SuitableRing<CoefficientRepresentation = RqPoly>,
        RqPoly: PolyRing + CRT,
    {
        // Create a test vector
        const N: usize = 32;
        let mut rng = ark_std::test_rng();
        let test_vector: Vec<RqPoly> = (0..N)
            .map(|_| draw_ring_below_bound::<RqPoly, { DP::B }>(&mut rng))
            .collect();

        // Call the function
        let decomposed = decompose_B_vec_into_k_vec::<RqNTT, DP>(&test_vector);

        // Check that we get K vectors back from the decomposition
        assert_eq!(
            decomposed.len(),
            DP::K,
            "Decomposition should output K={} vectors",
            DP::K
        );

        // Check the length of each inner vector
        for vec in &decomposed {
            assert_eq!(vec.len(), N);
        }

        // Check that the decomposition is correct
        for i in 0..N {
            let decomp_i = decomposed.iter().map(|d_j| d_j[i]).collect::<Vec<_>>();
            assert_eq!(
                test_vector[i],
                recompose(&decomp_i, RqPoly::from(DP::B_SMALL as u128))
            );
        }
    }

    fn recompose_from_k_vec_to_big_vec<NTT: SuitableRing>(
        k_vecs: &[Vec<NTT>],
    ) -> Vec<NTT::CoefficientRepresentation> {
        let decomposed_in_b: Vec<Vec<NTT::CoefficientRepresentation>> = k_vecs
            .iter()
            .map(|vec| {
                vec.iter()
                    .flat_map(|&x| decompose_balanced_vec(&[x.icrt()], DP::B, DP::L))
                    .flatten()
                    .collect()
            })
            .collect();

        // Transpose the decomposed vectors
        let mut transposed = vec![vec![]; decomposed_in_b[0].len()];
        for row in &decomposed_in_b {
            for (j, &val) in row.iter().enumerate() {
                transposed[j].push(val);
            }
        }

        // Recompose first with B_SMALL, then with B
        transposed
            .iter()
            .map(|vec| {
                recompose(
                    vec,
                    NTT::CoefficientRepresentation::from(DP::B_SMALL as u128),
                )
            })
            .collect::<Vec<_>>()
            .chunks(DP::L)
            .map(|chunk| recompose(chunk, NTT::CoefficientRepresentation::from(DP::B)))
            .collect()
    }

    fn test_decompose_big_vec_into_k_vec_and_compose_back<RqNTT, RqPoly>()
    where
        RqNTT: SuitableRing<CoefficientRepresentation = RqPoly>,
        RqPoly: PolyRing + CRT,
        Vec<RqNTT>: FromIterator<<RqPoly as CRT>::CRTForm>,
    {
        // Create a test vector
        const N: usize = 32;
        let mut rng = ark_std::test_rng();
        let test_vector: Vec<RqNTT> = (0..N)
            .map(|_| draw_ring_below_bound::<RqPoly, { DP::B }>(&mut rng).crt())
            .collect();
        let decomposed_and_composed_back =
            decompose_big_vec_into_k_vec_and_compose_back::<RqNTT, DP>(test_vector.clone());
        let restore_decomposed =
            recompose_from_k_vec_to_big_vec::<RqNTT>(&decomposed_and_composed_back);

        // Check each entry matches
        for i in 0..N {
            assert_eq!(
                restore_decomposed[i],
                test_vector[i].icrt(),
                "Mismatch at index {}: decomposed_and_composed_back={}, test_vector={}",
                i,
                restore_decomposed[i],
                test_vector[i].icrt()
            );
        }
    }

    #[test]
    fn test_decompose_B_vec_into_k_vec_gold() {
        test_decompose_B_vec_into_k_vec::<RqNTT, RqPoly>();
    }

    #[test]
    fn test_decompose_big_vec_into_k_vec_and_compose_back_gold() {
        test_decompose_big_vec_into_k_vec_and_compose_back::<RqNTT, RqPoly>();
    }
}
