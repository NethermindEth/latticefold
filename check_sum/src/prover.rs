use crate::{ poly_utils::MultiPoly, transcript::{ SumCheckRound, SumCheckTranscript } };
use lattirust_arithmetic::ring::Ring;
use std::sync::Arc;
use rand::Rng;

pub struct SumCheckProver<R: Ring> {
    pub polynomial: MultiPoly<R>,
    pub cyclotomic_ring_modulus: u128,
    pub claimed_sum: R,
}
impl<R: Ring> SumCheckProver<R> {
    pub fn prove(&self) -> SumCheckTranscript<R> {
        let mut rng = rand::thread_rng();
        let num_vars = self.polynomial.num_vars();
        let mut poly = self.polynomial.clone();
        let mut transcript = SumCheckTranscript::<R>::new(
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

            // send univariate restriction to verifier

            let challenge = R::from(rng.gen_range(0..self.cyclotomic_ring_modulus));
            // restrict according to random challenge
            let round: SumCheckRound<R> = SumCheckRound::new(poly.clone(), uni, j, challenge);

            transcript.add_round(round);
            poly = poly.partial_eval(challenge, j).simplify().clone();
        }
        transcript
    }
}
