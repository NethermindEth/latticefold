use ark_std::log2;
use latticefold::{
    arith::r1cs::R1CS,
    transcript::Transcript,
    utils::sumcheck::{
        utils::{build_eq_x_r, eq_eval},
        MLSumcheck, Proof,
    },
};
use stark_rings::{OverField, Ring};
use stark_rings_poly::mle::DenseMultilinearExtension;

/// Assume $n=m*\hat{l}$.
#[derive(Debug)]
pub struct CommittedR1CS<R: Ring> {
    r1cs: R1CS<R>,
    cm: Vec<R>, // kappa
    x: Vec<R>,  // l_in
    f: Vec<R>,  // n
}

#[derive(Debug)]
pub struct Lin<R: Ring> {
    sumcheck_proof: Proof<R>,
    nvars: usize,
    degree: usize,
    v: R,
    va: R,
    vb: R,
    vc: R,
}

impl<R: OverField> CommittedR1CS<R> {
    pub fn linearize(&self, transcript: &mut impl Transcript<R>) -> Lin<R> {
        let nvars = log2(self.f.len().next_power_of_two()) as usize;
        let ga = self.r1cs.A.try_mul_vec(&self.f).unwrap();
        let gb = self.r1cs.B.try_mul_vec(&self.f).unwrap();
        let gc = self.r1cs.C.try_mul_vec(&self.f).unwrap();
        let r: Vec<R> = transcript
            .get_challenges(nvars)
            .into_iter()
            .map(|x| x.into())
            .collect();
        let eq = build_eq_x_r(&r).unwrap();
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

        let mle_f = DenseMultilinearExtension::from_evaluations_slice(nvars, &self.f);
        let v = mle_f.evaluate(&ro).unwrap();
        let va = mle_ga.evaluate(&ro).unwrap();
        let vb = mle_gb.evaluate(&ro).unwrap();
        let vc = mle_gc.evaluate(&ro).unwrap();
        Lin {
            sumcheck_proof,
            nvars,
            degree: 3,
            v,
            va,
            vb,
            vc,
        }
    }
}

impl<R: OverField> Lin<R> {
    pub fn verify(&self, transcript: &mut impl Transcript<R>) {
        let r: Vec<R> = transcript
            .get_challenges(self.nvars)
            .into_iter()
            .map(|x| x.into())
            .collect();
        let subclaim = MLSumcheck::verify_as_subprotocol(
            transcript,
            self.nvars,
            self.degree,
            R::zero(),
            &self.sumcheck_proof,
        )
        .unwrap();

        let ro: Vec<R> = subclaim.point.into_iter().map(|x| x.into()).collect();
        let s = subclaim.expected_evaluation;
        let e = eq_eval(&r, &ro).unwrap();

        assert_eq!(e * (self.va * self.vb - self.vc), s)
    }
}

#[cfg(test)]
mod tests {
    use ark_std::One;
    use cyclotomic_rings::rings::GoldilocksPoseidonConfig as PC;
    use latticefold::transcript::poseidon::PoseidonTS;
    use stark_rings::{
        balanced_decomposition::GadgetDecompose, cyclotomic_ring::models::goldilocks::RqPoly as R,
    };
    use stark_rings_linalg::SparseMatrix;

    use super::*;

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
        let (mut r1cs, z) = identity_cs(1 << 5);

        r1cs.A.coeffs[0][0].0 = 9u128.into();
        r1cs.C.coeffs[0][0].0 = 9u128.into();

        r1cs.A = r1cs.A.gadget_decompose(2, 4);
        r1cs.B = r1cs.B.gadget_decompose(2, 4);
        r1cs.C = r1cs.C.gadget_decompose(2, 4);

        let f = z.gadget_decompose(2, 4);
        r1cs.check_relation(&f).unwrap();

        let cr1cs = CommittedR1CS {
            r1cs,
            f,
            x: vec![1u128.into()],
            cm: vec![],
        };

        let mut ts = PoseidonTS::default::<PC>();
        let lin = cr1cs.linearize(&mut ts);

        let mut ts = PoseidonTS::default::<PC>();
        lin.verify(&mut ts);
    }
}
