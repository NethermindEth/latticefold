use std::sync::Arc;

use super::{SumCheckError, SumCheckIP, SumCheckProof};
use crate::transcript::Transcript;
use ark_crypto_primitives::sponge::Absorb;
use ark_ff::PrimeField;
use lattirust_arithmetic::{
    challenge_set::latticefold_challenge_set::{LatticefoldChallengeSet, OverField},
    mle::DenseMultilinearExtension,
    polynomials::fix_variables,
    polynomials::VirtualPolynomial,
};

pub struct SumCheckProver<F: PrimeField, R: OverField<F>, CS: LatticefoldChallengeSet<F, R>>
where
    F: Absorb,
{
    pub _marker: std::marker::PhantomData<(F, CS)>,
    pub polynomial: VirtualPolynomial<R>,
    pub claimed_sum: R,
}

impl<F: PrimeField, R: OverField<F>, CS: LatticefoldChallengeSet<F, R>> SumCheckProver<F, R, CS>
where
    F: Absorb,
{
    pub fn prove(
        &self,
        transcript: &mut impl Transcript<F, R, ChallengeSet = CS>,
    ) -> Result<(SumCheckIP<F, R>, SumCheckProof<F, R>), SumCheckError<R>> {
        let num_vars = self.polynomial.aux_info.num_variables;
        let mut poly = self.polynomial.clone();
        let mut sum_check_proof = SumCheckProof::<F, R>::new(num_vars);

        let protocol = SumCheckIP::new(self.claimed_sum, self.polynomial.aux_info.clone());

        for j in 0..num_vars {
            let mut flattened_ml_extensions: Vec<DenseMultilinearExtension<R>> = poly
                .flattened_ml_extensions
                .iter()
                .map(|x| x.as_ref().clone())
                .collect();

            let challenge = transcript.get_big_challenge();

            flattened_ml_extensions.iter_mut().for_each(|mle| {
                let eval0 = (0..1 << (num_vars - j - 1))
                    .fold(R::zero(), |acc, index| acc + mle.evaluations[index]);
                let eval1 = (1 << (num_vars - j - 1)..1 << (num_vars - j))
                    .fold(R::zero(), |acc, index| acc + mle.evaluations[index]);
                *mle = DenseMultilinearExtension::<R>::from_evaluations_slice(1, &[eval0, eval1]);
            });

            let mut uni: VirtualPolynomial<R> = VirtualPolynomial::new(1);

            for (coeff, indices) in &poly.products {
                let mut new_poly = VirtualPolynomial::new(1);
                new_poly.add_mle_list(
                    indices
                        .iter()
                        .map(|i| Arc::from(flattened_ml_extensions[*i].clone())),
                    *coeff,
                )?;

                uni = &uni + &new_poly;
            }

            sum_check_proof.add_round(transcript, uni);
            poly.flattened_ml_extensions.iter_mut().for_each(|mle| {
                *mle = Arc::from(fix_variables(mle, &[challenge.into()]));
            });
        }

        Ok((protocol, sum_check_proof))
    }
}
