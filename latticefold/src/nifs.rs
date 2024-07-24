pub mod decomposition;
pub mod folding;
pub mod linearization;

use std::marker::PhantomData;

use ark_ff::Field;
use lattirust_arithmetic::challenge_set::latticefold_challenge_set::OverField;

use crate::transcript::Transcript;

pub struct NIFSProver<F: Field, R: OverField<F>, T: Transcript<F, R>> {
    _f: PhantomData<F>,
    _r: PhantomData<R>,
    _t: PhantomData<T>,
}

pub struct NIFSVerifier<F: Field, R: OverField<F>, T: Transcript<F, R>> {
    _f: PhantomData<F>,
    _r: PhantomData<R>,
    _t: PhantomData<T>,
}
