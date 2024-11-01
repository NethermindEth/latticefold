#![allow(non_snake_case)]
use ark_std::{fmt::Debug, UniformRand};
use cyclotomic_rings::SuitableRing;
use latticefold::arith::{r1cs::get_test_dummy_r1cs, CCS};

pub fn get_test_dummy_ccs<
    R: Clone + UniformRand + Debug + SuitableRing + for<'a> std::ops::AddAssign<&'a R>,
    const X_LEN: usize,
    const WIT_LEN: usize,
    const W: usize,
>(
    r1cs_rows: usize,
) -> CCS<R> {
    let r1cs = get_test_dummy_r1cs::<R, X_LEN, WIT_LEN>(r1cs_rows);
    CCS::<R>::from_r1cs(r1cs, W)
}
