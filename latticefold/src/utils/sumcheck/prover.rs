use std::sync::Arc;

use super::{univ_poly::UnivPoly, SumCheckError, SumCheckIP, SumCheckProof};
use crate::transcript::Transcript;
use lattirust_arithmetic::{
    challenge_set::latticefold_challenge_set::{LatticefoldChallengeSet, OverField},
    mle::DenseMultilinearExtension,
    polynomials::fix_variables,
    polynomials::VirtualPolynomial,
};

pub struct SumCheckProver<R: OverField, CS: LatticefoldChallengeSet<R>> {
    pub _marker: std::marker::PhantomData<CS>,
    pub polynomial: VirtualPolynomial<R>,
    pub claimed_sum: R,
}

impl<R: OverField, CS: LatticefoldChallengeSet<R>> SumCheckProver<R, CS> {
    pub fn prove(
        &self,
        transcript: &mut impl Transcript<R, ChallengeSet = CS>,
    ) -> Result<(SumCheckIP<R>, SumCheckProof<R>), SumCheckError<R>> {
        let num_vars = self.polynomial.aux_info.num_variables;
        let mut poly = self.polynomial.clone();
        let mut sum_check_proof = SumCheckProof::<R>::new(num_vars);

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

            sum_check_proof.add_round(transcript, UnivPoly::try_from(uni)?);
            poly.flattened_ml_extensions.iter_mut().for_each(|mle| {
                *mle = Arc::from(fix_variables(mle, &[challenge.into()]));
            });
        }
        Ok((protocol, sum_check_proof))
    }
}
