// Adapted from https://github.com/bgillesp/pazk/blob/main/src/bin/sum_check.rs
// Copyright (c) 2022 Bryan Gillespie

use std::sync::Arc;

use rand::Rng;
use lattirust_arithmetic::ring::Ring;
use crate::poly_utils::MultiPoly;
use crate::ip::{ Channel, Data, IPResult, IP };

pub struct SumCheckVerifier<R: Ring> {
    pub cyclotomic_ring_modulus: u128,
    pub polynomial: Arc<MultiPoly<R>>,
    pub degrees: Vec<usize>,
    pub claimed_sum: R,
}

impl<R: Ring> IP<Data<R>> for SumCheckVerifier<R> {
    fn execute(&self, ch: Channel<Data<R>>) -> Result<IPResult, IPResult> {
        let mut rng = rand::thread_rng();

        let zero = R::zero();
        let one = R::one();

        let mut check_value = self.claimed_sum;
        let mut challenges: Vec<R> = Vec::with_capacity(self.polynomial.num_vars());

        for j in (0..self.polynomial.num_vars()).rev() {
            let uni = ch.receive().to_polynomial().unwrap();
            if uni.degree() > self.degrees[j] {
                let data = Data::Decision(false);
                ch.send(data);
                return Err(IPResult::REJECT);
            }

            if uni.eval(&zero) + uni.eval(&one) != check_value {
                let data = Data::Decision(false);

                ch.send(data);
                return Err(IPResult::REJECT);
            }

            let challenge = R::from(rng.gen_range(0..self.cyclotomic_ring_modulus));

            check_value = uni.eval(&challenge);

            // record challenge for later reference
            challenges.push(challenge);
            let data = Data::Scalar(challenge);
            if j > 0 {
                ch.send(data);
            } else {
                challenges.reverse();
                let oracle_evaluation = self.polynomial.eval_poly(&challenges);

                // accept if the oracle evaluation equals the final check value
                // otherwise reject
                let decision: Data<R> = Data::Decision(oracle_evaluation == check_value);

                ch.send(Data::Decision(oracle_evaluation == check_value));
            }
            // evaluate polynomial at vector of challenge points
        }
        Ok(IPResult::ACCEPT)
    }
}
