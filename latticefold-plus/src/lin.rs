use latticefold::transcript::Transcript;
use stark_rings::{
    balanced_decomposition::{convertible_ring::ConvertibleRing, Decompose},
    CoeffRing, OverField, Zq,
};
use stark_rings_linalg::{Matrix, SparseMatrix};

use crate::{
    cm::CmProof,
    mlin::{LinB2, Mlin},
    rgchk::DecompParameters,
};

pub trait Linearize<R: OverField> {
    type Proof: Verify<R>;
    fn linearize(&self, transcript: &mut impl Transcript<R>) -> (LinB<R>, Self::Proof);
}

pub trait Verify<R: OverField> {
    fn verify(&self, transcript: &mut impl Transcript<R>) -> bool;
}

#[derive(Clone, Debug)]
pub struct LinParameters {
    pub kappa: usize,
    pub decomp: DecompParameters,
}

#[derive(Clone, Debug)]
pub struct LinBX<R> {
    pub cm_f: Vec<R>,
    pub r: Vec<(R, R)>,
    pub v: Vec<(R, R)>,
}

#[derive(Clone, Debug)]
pub struct LinB<R> {
    pub f: Vec<R>,
    pub x: LinBX<R>,
}

impl<R: CoeffRing> LinB<R>
where
    R::BaseRing: ConvertibleRing + Decompose + Zq,
    R: Decompose,
{
    /// Πlin protocol
    ///
    /// Runs the Πmlin protocol with only L=1 instance
    pub fn lin(
        &self,
        M: &[SparseMatrix<R>],
        A: &Matrix<R>,
        params: &LinParameters,
        transcript: &mut impl Transcript<R>,
    ) -> (LinB2<R>, CmProof<R>) {
        let mlin = Mlin {
            lins: vec![self.clone()],
            A: A.clone(),
            params: params.clone(),
        };

        mlin.mlin(M, transcript)
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::PrimeField;
    use ark_std::One;
    use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
    use latticefold::{arith::r1cs::R1CS, transcript::poseidon::PoseidonTS};
    use stark_rings::{
        balanced_decomposition::GadgetDecompose, cyclotomic_ring::models::frog_ring::RqPoly as R,
        PolyRing,
    };
    use stark_rings_linalg::SparseMatrix;

    use super::*;
    use crate::r1cs::CommittedR1CS;

    #[test]
    fn test_lin() {
        let n = 1 << 15;
        let k = 2;
        let m = n / k;
        let z = vec![R::one(); m];

        let mut r1cs = R1CS::<R> {
            l: 1,
            A: SparseMatrix::identity(m),
            B: SparseMatrix::identity(m),
            C: SparseMatrix::identity(m),
        };

        r1cs.A.coeffs[0][0].0 = 2u128.into();
        r1cs.C.coeffs[0][0].0 = 2u128.into();

        r1cs.A = r1cs.A.gadget_decompose(2, k);
        r1cs.B = r1cs.B.gadget_decompose(2, k);
        r1cs.C = r1cs.C.gadget_decompose(2, k);

        let f = z.gadget_decompose(2, k);
        r1cs.check_relation(&f).unwrap();

        let cr1cs = CommittedR1CS {
            r1cs,
            f,
            x: vec![1u128.into()],
            cm: vec![],
        };

        let kappa = 2;
        let b = (R::dimension() / 2) as u128;
        // log_d' (q)
        let l = ((<<R as PolyRing>::BaseRing>::MODULUS.0[0] as f64).ln()
            / ((R::dimension() / 2) as f64).ln())
        .ceil() as usize;
        let params = LinParameters {
            kappa,
            decomp: DecompParameters { b, k, l },
        };
        let A = Matrix::<R>::rand(&mut rand::thread_rng(), params.kappa, n);
        let M = vec![
            cr1cs.r1cs.A.clone(),
            cr1cs.r1cs.B.clone(),
            cr1cs.r1cs.C.clone(),
        ];

        let mut ts = PoseidonTS::default::<PC>();
        let (linb, lproof) = cr1cs.linearize(&mut ts);
        let (_linb2, cmproof) = linb.lin(&M, &A, &params, &mut ts);

        let mut ts = PoseidonTS::default::<PC>();
        lproof.verify(&mut ts);
        cmproof.verify(&mut ts).unwrap();
    }
}
