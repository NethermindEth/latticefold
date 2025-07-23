use ark_std::{
    iter::once,
    log2,
    ops::{Mul, Sub},
    One, Zero,
};
use latticefold::{
    transcript::Transcript,
    utils::sumcheck::{
        utils::{build_eq_x_r, eq_eval},
        MLSumcheck, Proof, SumCheckError,
    },
};
use stark_rings::{
    balanced_decomposition::{
        convertible_ring::ConvertibleRing, Decompose, DecomposeToVec, GadgetDecompose,
    },
    exp, psi, psi_range_check, unit_monomial, CoeffRing, OverField, PolyRing, Ring, Zq,
};
use stark_rings_linalg::{ops::Transpose, Matrix, SparseMatrix};
use stark_rings_poly::mle::{DenseMultilinearExtension, SparseMultilinearExtension};
use thiserror::Error;

use crate::{
    cm::{Cm, CmComs, CmProof},
    lin::{Lin, Parameters},
    mlin::Mlin,
    rgchk::{Dcom, Rg, RgInstance},
    utils::{tensor, tensor_product},
};

#[derive(Debug)]
pub struct Plus<R> {
    instances: Mlin<R>,
}

impl<R: CoeffRing> Plus<R>
where
    R::BaseRing: ConvertibleRing + Decompose + Zq,
    R: Decompose,
{
    /// Fold
    pub fn fold(&self, transcript: &mut impl Transcript<R>) {
        let (linb2, cmproof) = self.instances.mlin(transcript);
    }
}
