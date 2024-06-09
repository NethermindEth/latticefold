// Copied from pazk
// TODO add a propoer attribution here

use std::thread;
use std::fmt::{ self, Display };

use std::sync::mpsc;
use std::sync::{ Arc, Mutex };

use lattirust_arithmetic::ring::Ring;
use crate::poly_utils::UnivPoly;

// 2-party interactive protocol
pub trait IP<T: Clone> {
    fn execute(&self, ch: Channel<T>, log: Log);
}

pub fn execute<T: Clone + Send + 'static + Display>(
    prover: impl IP<T> + Send + 'static,
    verifier: impl IP<T> + Send + 'static
) {
    let (ch1, ch2) = Channel::<T>::gen();
    let log = Log::new();
    let (lg1, lg2) = (log.clone(), log.clone());

    let prover_handle = thread::spawn(move || {
        prover.execute(ch1, lg1);
    });
    let verifier_handle = thread::spawn(move || {
        verifier.execute(ch2, lg2);
    });

    prover_handle.join().unwrap();
    verifier_handle.join().unwrap();

    let log = &*log.get_log().lock().unwrap();
    for message in log {
        println!("{}", message);
    }
}

#[derive(Clone)]
pub struct Log {
    log: Arc<Mutex<Vec<String>>>,
}

impl Log {
    pub fn new() -> Log {
        let vec: Vec<String> = Vec::new();
        Log {
            log: Arc::new(Mutex::new(vec)),
        }
    }

    pub fn write(&self, message: String) {
        let mut log = self.log.lock().unwrap();
        log.push(message);
    }

    pub fn get_log(&self) -> &Arc<Mutex<Vec<String>>> {
        &self.log
    }
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
