use crate::arith::{Witness, LCCCS};
use crate::commitment::Commitment;
use lattirust_poly::mle::DenseMultilinearExtension;
use lattirust_ring::OverField;
use num_traits::Zero;

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
