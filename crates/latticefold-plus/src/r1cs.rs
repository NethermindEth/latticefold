use ark_std::log2;
use latticefold::{
    arith::r1cs::R1CS,
    transcript::Transcript,
    utils::sumcheck::{
        utils::{build_eq_x_r, eq_eval},
        MLSumcheck, Proof,
    },
};
use stark_rings::{
    balanced_decomposition::{Decompose, GadgetDecompose},
    OverField, Ring,
};
use stark_rings_linalg::{Matrix, SparseMatrix};
use stark_rings_poly::mle::DenseMultilinearExtension;

use crate::lin::{LinB, LinBX, Linearize, LinearizedVerify};

/// Committed R1CS
///
/// Assume $n=m*\hat{l}$.
#[derive(Clone, Debug)]
pub struct ComR1CS<R: Ring> {
    pub x: ComR1CSX<R>,
    pub f: Vec<R>, // n
}

#[derive(Clone, Debug)]
pub struct ComR1CSX<R: Ring> {
    pub r1cs: R1CS<R>,
    pub z: Vec<R>,    // m
    pub cm_f: Vec<R>, // kappa
    /// Public input length
    pub l_in: usize,
}

#[derive(Clone, Debug)]
pub struct ComR1CSProof<R: Ring> {
    pub sumcheck_proof: Proof<R>,
    pub nvars: usize,
    pub r: Vec<R>,
    pub v: R,
    pub va: R,
    pub vb: R,
    pub vc: R,
}

impl<R: Decompose + Ring> ComR1CS<R> {
    pub fn new(r1cs: R1CS<R>, z: Vec<R>, l_in: usize, b: u128, k: usize, A: &Matrix<R>) -> Self {
        let f = z.gadget_decompose(b, k);
        let cm_f = A.try_mul_vec(&f).unwrap();
        let x = ComR1CSX {
            r1cs,
            z,
            cm_f,
            l_in,
        };
        Self { x, f }
    }
}

impl<R: Ring> ComR1CSX<R> {
    pub fn matrices(&self) -> Vec<SparseMatrix<R>> {
        vec![
            self.r1cs.A.clone(),
            self.r1cs.B.clone(),
            self.r1cs.C.clone(),
        ]
    }
}

impl<R: OverField> Linearize<R> for ComR1CS<R> {
    type Proof = ComR1CSProof<R>;
    fn linearize(&self, transcript: &mut impl Transcript<R>) -> (LinB<R>, Self::Proof) {
        let nvars = log2(self.f.len().next_power_of_two()) as usize;
        let ga = self.x.r1cs.A.try_mul_vec(&self.f).unwrap();
        let gb = self.x.r1cs.B.try_mul_vec(&self.f).unwrap();
        let gc = self.x.r1cs.C.try_mul_vec(&self.f).unwrap();
        let r: Vec<R> = transcript
            .get_challenges(nvars)
            .into_iter()
            .map(|x| x.into())
            .collect();
        let eq = build_eq_x_r(&r).unwrap();
        let mle_f = DenseMultilinearExtension::from_evaluations_vec(nvars, self.f.clone());
        let mle_ga = DenseMultilinearExtension::from_evaluations_vec(nvars, ga);
        let mle_gb = DenseMultilinearExtension::from_evaluations_vec(nvars, gb);
        let mle_gc = DenseMultilinearExtension::from_evaluations_vec(nvars, gc);

        let mles = vec![eq, mle_ga.clone(), mle_gb.clone(), mle_gc.clone()];

        let comb_fn = |vals: &[R]| -> R { vals[0] * (vals[1] * vals[2] - vals[3]) };

        let (sumcheck_proof, prover_state) =
            MLSumcheck::prove_as_subprotocol(transcript, mles, nvars, 3, comb_fn);
        let ro = prover_state
            .randomness
            .into_iter()
            .map(|x| x.into())
            .collect::<Vec<R>>();

        let v = mle_f.evaluate(&ro).unwrap();
        let va = mle_ga.evaluate(&ro).unwrap();
        let vb = mle_gb.evaluate(&ro).unwrap();
        let vc = mle_gc.evaluate(&ro).unwrap();

        absorb_evaluations(&[v, va, vb, vc], transcript);

        let proof = Self::Proof {
            sumcheck_proof,
            nvars,
            r: ro.clone(),
            v,
            va,
            vb,
            vc,
        };

        let r = ro.iter().map(|&r| (r, r)).collect::<Vec<_>>();
        let v = vec![(v, v), (va, va), (vb, vb), (vc, vc)];

        let x = LinBX {
            cm_f: self.x.cm_f.clone(),
            r,
            v,
        };
        let linb = LinB {
            f: self.f.clone(),
            x,
        };

        (linb, proof)
    }
}

impl<R: OverField> LinearizedVerify<R> for ComR1CSProof<R> {
    fn verify(&self, transcript: &mut impl Transcript<R>) -> bool {
        let r: Vec<R> = transcript
            .get_challenges(self.nvars)
            .into_iter()
            .map(|x| x.into())
            .collect();
        let subclaim = MLSumcheck::verify_as_subprotocol(
            transcript,
            self.nvars,
            3,
            R::zero(),
            &self.sumcheck_proof,
        )
        .unwrap();

        let ro: Vec<R> = subclaim.point.into_iter().map(|x| x.into()).collect();
        let s = subclaim.expected_evaluation;

        absorb_evaluations(&[self.v, self.va, self.vb, self.vc], transcript);

        let e = eq_eval(&r, &ro).unwrap();

        assert_eq!(e * (self.va * self.vb - self.vc), s);

        true
    }
}

fn absorb_evaluations<R: OverField>(evals: &[R; 4], transcript: &mut impl Transcript<R>) {
    transcript.absorb_slice(evals);
}

/// Decomposes and squares a R1CS
///
/// n x m -> n x n, where m * k = n
pub fn r1cs_decomposed_square<R: Decompose + Ring>(
    mut r1cs: R1CS<R>,
    n: usize,
    b: u128,
    k: usize,
) -> R1CS<R> {
    r1cs.A = r1cs.A.gadget_decompose(b, k);
    r1cs.B = r1cs.B.gadget_decompose(b, k);
    r1cs.C = r1cs.C.gadget_decompose(b, k);
    r1cs.A.pad_rows(n);
    r1cs.B.pad_rows(n);
    r1cs.C.pad_rows(n);
    r1cs
}

#[cfg(test)]
mod tests {
    use ark_std::One;
    use latticefold_rings::rings::GoldilocksPoseidonConfig as PC;
    use stark_rings::{
        balanced_decomposition::GadgetDecompose, cyclotomic_ring::models::goldilocks::RqPoly as R,
    };
    use stark_rings_linalg::SparseMatrix;

    use super::*;
    use crate::transcript::PoseidonTranscript;

    fn identity_cs(n: usize) -> (R1CS<R>, Vec<R>) {
        let r1cs = R1CS::<R> {
            l: 1,
            A: SparseMatrix::identity(n),
            B: SparseMatrix::identity(n),
            C: SparseMatrix::identity(n),
        };
        let z = vec![R::one(); n];
        (r1cs, z)
    }

    #[test]
    fn test_linearization() {
        let n = 1 << 7;
        let k = 4;
        let m = n / k;
        let b = 2;
        let kappa = 2;
        let (mut r1cs, z) = identity_cs(m);

        r1cs.A = r1cs.A.gadget_decompose(b, k);
        r1cs.B = r1cs.B.gadget_decompose(b, k);
        r1cs.C = r1cs.C.gadget_decompose(b, k);

        let A = Matrix::<R>::rand(&mut ark_std::test_rng(), kappa, n);
        let cr1cs = ComR1CS::new(r1cs, z, 1, b, k, &A);
        cr1cs.x.r1cs.check_relation(&cr1cs.f).unwrap();

        let mut ts = PoseidonTranscript::empty::<PC>();
        let (_linb, lproof) = cr1cs.linearize(&mut ts);

        let mut ts = PoseidonTranscript::empty::<PC>();
        lproof.verify(&mut ts);
    }
}
