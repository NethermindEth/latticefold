use std::marker::PhantomData;

use crate::arith::{Witness, LCCCS};
use crate::commitment::Commitment;
use lattirust_poly::mle::DenseMultilinearExtension;
use lattirust_ring::OverField;
use num_traits::Zero;

use super::decomposition::DecompositionProof;
use super::folding::FoldingProof;
use super::linearization::LinearizationProof;

pub struct LatticefoldState<const C: usize, R: OverField> {
    pub powers_of_b: Vec<R>,
    pub mz_mles: Vec<DenseMultilinearExtension<R>>,
    pub lcccs: LCCCS<C, R>,
    pub decomposed_lcccs_s: Vec<LCCCS<C, R>>,
    pub wit_s: Vec<Witness<R>>,
}

impl<const C: usize, R: OverField + Default> Default for LatticefoldState<C, R> {
    fn default() -> Self {
        LatticefoldState {
            powers_of_b: Vec::new(),

            mz_mles: Vec::new(),
            lcccs: LCCCS {
                r: vec![],
                v: R::default(),
                cm: Commitment::zero(),
                u: vec![],
                x_w: vec![],
                h: R::default(),
            },
            decomposed_lcccs_s: Vec::new(),
            wit_s: Vec::new(),
        }
    }
}

/// `C` is the length of Ajtai commitment vectors.
/// `NTT` is a cyclotomic ring in the NTT form.
#[derive(Clone)]
pub struct LFProof<const C: usize, NTT: OverField> {
    pub linearization_proof: LinearizationProof<NTT>,
    pub decomposition_proof_l: DecompositionProof<C, NTT>,
    pub decomposition_proof_r: DecompositionProof<C, NTT>,
    pub folding_proof: FoldingProof<NTT>,
}

/// `C` is the length of commitment vectors or, equivalently, the number of rows of the Ajtai matrix.
/// `W` is the length of witness vectors or, equivalently, the number of columns of the Ajtai matrix.
/// `NTT` is a suitable cyclotomic ring.
/// `P` is the decomposition parameters.
/// `T` is the FS-transform transcript.
pub struct NIFSProver<const C: usize, const W: usize, NTT, P, T> {
    _r: PhantomData<NTT>,
    _p: PhantomData<P>,
    _t: PhantomData<T>,
}