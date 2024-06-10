// ! Provides generic structure for an interactive proof
// Adapted from https://github.com/bgillesp/pazk/blob/main/src/ip.rs
// Copyright (c) 2022 Bryan Gillespie

use std::thread;
use std::fmt::{ self, Display };

use std::sync::mpsc;

use lattirust_arithmetic::ring::Ring;
use crate::poly_utils::UnivPoly;
#[derive(Debug)]
pub enum IPResult {
    ACCEPT,
    REJECT,
}

// 2-party interactive protocol
pub trait IP<T: Clone> {
    fn execute(&self, ch: Channel<T>) -> Result<IPResult, IPResult>;
}

pub fn execute<T: Clone + Send + 'static + Display>(
    prover: impl IP<T> + Send + 'static,
    verifier: impl IP<T> + Send + 'static
) -> Result<IPResult, IPResult> {
    let (ch1, ch2) = Channel::<T>::gen();

    let prover_handle = thread::spawn(move || { prover.execute(ch1) });
    let verifier_handle = thread::spawn(move || {
        verifier.execute(ch2);
    });
    verifier_handle.join().unwrap();
    prover_handle.join().unwrap()
}

// bidirectional channel
pub struct Channel<T: Clone> {
    tx: mpsc::Sender<T>,
    rx: mpsc::Receiver<T>,
}

impl<T: Clone> Channel<T> {
    pub fn gen() -> (Channel<T>, Channel<T>) {
        let (tx1, rx1) = mpsc::channel();
        let (tx2, rx2) = mpsc::channel();
        (
            Channel {
                tx: tx1,
                rx: rx2,
            },
            Channel {
                tx: tx2,
                rx: rx1,
            },
        )
    }

    pub fn send(&self, data: T) {
        // send data down channel, ignore error if send fails
        self.tx.send(data.clone()).ok();
    }

    pub fn receive(&self) -> T {
        self.rx.recv().unwrap()
    }
}

#[derive(Clone)]
pub enum Data<R: Ring> {
    Scalar(R),
    Polynomial(UnivPoly<R>),
    Decision(bool),
}

impl<R: Ring> Data<R> {
    pub fn to_scalar(self) -> Option<R> {
        if let Data::Scalar(x) = self { Some(x) } else { None }
    }

    pub fn to_polynomial(self) -> Option<UnivPoly<R>> {
        if let Data::Polynomial(p) = self { Some(p) } else { None }
    }

    pub fn _to_decision(self) -> Option<bool> {
        if let Data::Decision(d) = self { Some(d) } else { None }
    }
}
impl<R: Ring> fmt::Display for Data<R> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Data::Scalar(x) => { write!(f, "{}", x) }
            Data::Polynomial(p) => { write!(f, "{}", p.format_univ_poly("x")) }
            Data::Decision(b) => {
                if *b { write!(f, "Accept") } else { write!(f, "Reject") }
            }
        }
    }
}

#[cfg(test)]
mod test {
    use crate::poly_utils::MultiPoly;
    use crate::prover::SumCheckProver;
    use crate::verifier::SumCheckVerifier;

    use lattirust_arithmetic::ring::Z2_64;
    use crate::ip;
    #[test]
    fn test_sum_check() {
        // Define a polynomial 3 + 2x1^2 + x1*x2 + 5x2 in Z2_64
        let poly = MultiPoly {
            terms: vec![
                (Z2_64::from(3 as u8), vec![]),
                (Z2_64::from(2 as u8), vec![(0, 2)]),
                (Z2_64::from(1 as u8), vec![(0, 1), (1, 1)]),
                (Z2_64::from(5 as u8), vec![(1, 1)])
            ],
        };
        let claimed_sum = Z2_64::from(27u8);
        let degrees = poly.multi_degree();

        let prover = SumCheckProver::<Z2_64> {
            polynomial: poly.clone().into(),
        };

        let verifier = SumCheckVerifier::<Z2_64> {
            cyclotomic_ring_modulus: 2,
            polynomial: poly.clone().into(),
            degrees,
            claimed_sum,
        };

        let result = ip::execute(prover, verifier);

        assert!(result.is_ok())
    }
}
