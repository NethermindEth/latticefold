use crate::{
    poly_utils::MultiPoly,
    sum_check_transcript::SumCheckTranscript,
    transcript::Transcript,
};
use ark_crypto_primitives::sponge::Absorb;
use ark_ff::PrimeField;
use lattirust_arithmetic::challenge_set::latticefold_challenge_set::{
    LatticefoldChallengeSet,
    OverField,
};

pub struct SumCheckProver<F: PrimeField, R: OverField<F>, CS: LatticefoldChallengeSet<F, R>>
    where F: Absorb {
    pub _marker: std::marker::PhantomData<(F, CS)>,
    pub polynomial: MultiPoly<R>,
    pub claimed_sum: R,
}
impl<F: PrimeField, R: OverField<F>, CS: LatticefoldChallengeSet<F, R>> SumCheckProver<F, R, CS>
    where F: Absorb
{
    pub fn prove(&self) -> SumCheckTranscript<F, R, CS> {
        let num_vars = self.polynomial.num_vars();
        let mut poly = self.polynomial.clone();
        let mut transcript = SumCheckTranscript::<F, R, CS>::new(
            self.claimed_sum,
            self.polynomial.clone().simplify(),
            self.polynomial.num_vars()
        );
        for j in (0..num_vars).rev() {
            let vals = (0..num_vars)
                .map(|n| {
                    if n < j { Some(vec![R::zero(), R::one()]) } else { None }
                })
                .collect();

            let partial = poly.partial_summation(&vals).simplify();

            let uni = partial.to_univariate(j);
            let challenge = transcript.hasher.get_big_challenge();
            transcript.add_round(challenge.into(), j, uni);
            poly = poly.partial_eval(challenge.into(), j).simplify().clone();
        }
        transcript
    }
}
