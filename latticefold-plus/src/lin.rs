use ark_std::log2;
use latticefold::{
    arith::r1cs::R1CS,
    transcript::Transcript,
    utils::sumcheck::{
        utils::{build_eq_x_r, eq_eval},
        MLSumcheck,
    },
};
use stark_rings::{OverField, Ring};
use stark_rings_poly::mle::DenseMultilinearExtension;

pub trait Linearize<R: OverField> {
    type Proof: Verify<R>;
    fn linearize(&self, transcript: &mut impl Transcript<R>) -> Self::Proof;
}

pub trait Verify<R: OverField> {
    fn verify(&self, transcript: &mut impl Transcript<R>) -> bool;
}

#[derive(Debug)]
pub struct Lin<R: Ring> {
    pub C_Mf: Vec<R>,
    pub cm_mtau: Vec<R>,
}
