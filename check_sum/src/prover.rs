use crate::poly_utils::{ MultiPoly, UnivPoly };
use lattirust_arithmetic::ring::Ring;
use std::{ iter::Sum, sync::{ Arc, Mutex } };
use crate::ip::{ IP, Data };

use std::fmt;
struct SumCheckProver<R: Ring> {
    polynomial: Arc<MultiPoly<R>>,
}
impl<R: Ring> IP<Data<R>> for SumCheckProver<R> {
    fn execute(&self, ch: crate::ip::Channel<Data<R>>, log: crate::ip::Log) {
        let num_vars = self.polynomial.num_vars();
        let mut poly = self.polynomial.clone();

        for j in (0..num_vars).rev() {
            log.write(format!("P computes univariate polynomial g_{}", j));
            let vals = (0..num_vars)
                .map(|n| {
                    if n < j { Some(vec![R::zero(), R::one()]) } else { None }
                })
                .collect();
            let partial = self.polynomial.partial_summation(&vals).simplify();
            let uni = partial.to_univariate(j);
            // send univariate restriction to verifier
            let data = Data::Polynomial(uni);
            log.write(format!("P --> (g_{} = {})", j, data));
            ch.send(data);
            if j > 0 {
                let data = ch.receive();
                if let Data::Decision(false) = data {
                    return;
                }
                let challenge = data.to_scalar().unwrap();

                // restrict according to random challenge
                log.write(format!("P computes partial evaluation at x_{} = r_{}", j, j));
                poly = Arc::new(self.polynomial.partial_eval(challenge, j));
            }
        }
    }
}
