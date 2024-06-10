// Adapted from https://github.com/bgillesp/pazk/blob/main/src/bin/sum_check.rs
// Copyright (c) 2022 Bryan Gillespie

use crate::poly_utils::MultiPoly;
use lattirust_arithmetic::ring::Ring;
use std::sync::Arc;
use crate::ip::{ IP, Data, IPResult };

pub struct SumCheckProver<R: Ring> {
    pub polynomial: Arc<MultiPoly<R>>,
}
impl<R: Ring> IP<Data<R>> for SumCheckProver<R> {
    fn execute(&self, ch: crate::ip::Channel<Data<R>>) -> Result<IPResult, IPResult> {
        let num_vars = self.polynomial.num_vars();
        let mut poly = self.polynomial.clone();
        for j in (0..num_vars).rev() {
            let vals = (0..num_vars)
                .map(|n| {
                    if n < j { Some(vec![R::zero(), R::one()]) } else { None }
                })
                .collect();
            let partial = self.polynomial.partial_summation(&vals).simplify();
            let uni = partial.to_univariate(j);
            // send univariate restriction to verifier
            let data = Data::Polynomial(uni);

            ch.send(data);
            if j > 0 {
                let data = ch.receive();
                if let Data::Decision(false) = data {
                    return Err(IPResult::REJECT);
                }
                let challenge = data.to_scalar().unwrap();
                print!("{:?}", challenge);
                // restrict according to random challenge

                poly = Arc::new(self.polynomial.partial_eval(challenge, j));
            }
        }
        Ok(IPResult::ACCEPT)
    }
}
