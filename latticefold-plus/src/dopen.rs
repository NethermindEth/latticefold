use stark_rings::{OverField, PolyRing, Ring};
use stark_rings_linalg::SparseMatrix;

pub struct Dopen<R: PolyRing> {
    cM: Vec<R>,            // kappa
    tau: Vec<R::BaseRing>, // n
    M: SparseMatrix<R>,    // n x m
}

impl<R: OverField> Dopen<R> {}
