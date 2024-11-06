use cyclotomic_rings::SuitableRing;
use lattirust_poly::mle::DenseMultilinearExtension;

use crate::arith::LCCCS;

use super::{
    decomposition::DecompositionProof, folding::FoldingProof, linearization::LinearizationProof,
};

pub struct LatticefoldState<const C: usize, R: SuitableRing> {
    pub powers_of_b: Vec<R>,
    pub r_0: Vec<R>,
    pub mz_mles: Vec<DenseMultilinearExtension<R>>,
    pub lcccs: LCCCS<C, R>,
    pub decomposed_lcccs_s: Vec<LCCCS<C, R>>,
}
