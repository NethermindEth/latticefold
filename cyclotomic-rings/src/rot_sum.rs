use ark_ff::{Field, Zero};
use lattirust_ring::{Cyclotomic, PolyRing};

use crate::SuitableRing;

pub fn rot_sum<R: SuitableRing>(
    a: R::CoefficientRepresentation,
    b: &[R::BaseRing],
) -> Vec<R::BaseRing> {
    assert_eq!(b.len(), R::CoefficientRepresentation::dimension());

    b.iter()
        .zip(a.into_rot_iter())
        .map(|(b_i, x_i_a)| {
            x_i_a
                .into_coeffs()
                .into_iter()
                .map(|x| <R::BaseRing as Field>::from_base_prime_field(x) * b_i)
                .collect::<Vec<R::BaseRing>>()
        })
        .fold(
            vec![R::BaseRing::zero(); R::CoefficientRepresentation::dimension()],
            |mut acc, b_i_times_coeff_x_i_a| {
                acc.iter_mut()
                    .zip(b_i_times_coeff_x_i_a)
                    .for_each(|(acc_j, b_i_times_coeff_x_i_a_j)| *acc_j += b_i_times_coeff_x_i_a_j);

                acc
            },
        )
}

#[cfg(test)]
mod tests {
    use ark_ff::UniformRand;
    use lattirust_ring::cyclotomic_ring::models::goldilocks::Fq3;
    use rand::thread_rng;

    use super::*;
    use crate::{GoldilocksRingNTT, GoldilocksRingPoly};

    #[test]
    fn test_rot_sum_with_coeffs() {
        let mut rng = thread_rng();
        let a = GoldilocksRingPoly::rand(&mut rng);
        let b = GoldilocksRingPoly::rand(&mut rng);

        // RotSum(a, coeff(b)) = coeff(a * b)
        assert_eq!(
            rot_sum::<GoldilocksRingNTT>(
                a,
                &b.coeffs()
                    .iter()
                    .map(|&x| Fq3::from_base_prime_field(x))
                    .collect::<Vec<Fq3>>()
            ),
            (a * b)
                .into_coeffs()
                .into_iter()
                .map(Fq3::from_base_prime_field)
                .collect::<Vec<Fq3>>()
        );
    }

    #[test]
    fn test_rot_sum_with_ring_elems() {
        let mut rng = thread_rng();

        let a = GoldilocksRingPoly::rand(&mut rng);
        let b: Vec<GoldilocksRingPoly> = (0..Fq3::extension_degree())
            .map(|_| GoldilocksRingPoly::rand(&mut rng))
            .collect();

        let capital_b: Vec<Fq3> = (0..GoldilocksRingPoly::dimension())
            .map(|i| {
                Fq3::from_base_prime_field_elems(&[
                    b[0].coeffs()[i],
                    b[1].coeffs()[i],
                    b[2].coeffs()[i],
                ])
                .unwrap()
            })
            .collect();

        // RotSum(a, B) = \sum coeff(a * b_i) * Y^{i-1}
        assert_eq!(
            rot_sum::<GoldilocksRingNTT>(a, &capital_b),
            (0..GoldilocksRingPoly::dimension())
                .map(|i| {
                    Fq3::from_base_prime_field_elems(&[
                        (a * b[0]).coeffs()[i],
                        (a * b[1]).coeffs()[i],
                        (a * b[2]).coeffs()[i],
                    ])
                    .unwrap()
                })
                .collect::<Vec<Fq3>>()
        );
    }
}
