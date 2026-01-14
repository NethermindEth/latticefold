use ark_std::log2;
use latticefold::{
    arith::r1cs::R1CS,
    transcript::Transcript,
    utils::sumcheck::{
        utils::eq_eval,
        MLSumcheck, Proof,
    },
};
use latticefold::commitment::AjtaiCommitmentScheme;
use stark_rings::{
    balanced_decomposition::{Decompose, GadgetDecompose},
    OverField, PolyRing, Ring,
};
use stark_rings_linalg::{Matrix, SparseMatrix};
use std::sync::Arc;

use crate::lin::{LinB, LinBX, Linearize, LinearizedVerify};
use crate::rgchk::WitnessVec;
use crate::streaming_sumcheck::{StreamingMleEnum, StreamingSumcheck};

/// Committed R1CS
///
/// Assume $n=m*\hat{l}$.
#[derive(Clone, Debug)]
pub struct ComR1CS<R: Ring + PolyRing> {
    pub x: ComR1CSX<R>,
    pub f: WitnessVec<R>, // n (prover-only representation; verifier unaffected)
}

#[derive(Clone, Debug)]
pub struct ComR1CSX<R: Ring> {
    // Store matrices as Arc to avoid catastrophic deep clones.
    pub a: Arc<SparseMatrix<R>>,
    pub b: Arc<SparseMatrix<R>>,
    pub c: Arc<SparseMatrix<R>>,
    /// Public input length
    pub l_in: usize,
    /// Number of public inputs in the underlying R1CS (kept for tests/debug).
    pub l: usize,
    pub z: Vec<R>,    // m (unused by verifier; retained for legacy/tests)
    pub cm_f: Vec<R>, // kappa
}

/// Committed R1CS (const-coeff/SP1 regime) where A/B/C are stored over the **base ring**.
///
/// This avoids storing `SparseMatrix<R>` (catastrophic for large ring dimension) while keeping the
/// verifier-side protocol unchanged: all proof objects and transcript flow remain identical.
#[derive(Clone, Debug)]
pub struct ComR1CSBase<R: Ring + PolyRing> {
    pub x: ComR1CSXBase<R>,
    pub f: WitnessVec<R>, // prover-only; should be `ConstCoeffBase` in this construction
}

#[derive(Clone, Debug)]
pub struct ComR1CSXBase<R: Ring + PolyRing> {
    pub a: Arc<SparseMatrix<R::BaseRing>>,
    pub b: Arc<SparseMatrix<R::BaseRing>>,
    pub c: Arc<SparseMatrix<R::BaseRing>>,
    pub l_in: usize,
    pub l: usize,
    pub cm_f: Vec<R>, // kappa
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

impl<R: Decompose + Ring + PolyRing> ComR1CS<R> {
    pub fn new(r1cs: R1CS<R>, z: Vec<R>, l_in: usize, b: u128, k: usize, A: &Matrix<R>) -> Self {
        let f = z.gadget_decompose(b, k);
        let cm_f = A.try_mul_vec(&f).unwrap();
        let l = r1cs.l;
        let x = ComR1CSX {
            a: Arc::new(r1cs.A),
            b: Arc::new(r1cs.B),
            c: Arc::new(r1cs.C),
            z,
            cm_f,
            l_in,
            l,
        };
        Self {
            x,
            f: WitnessVec::Ring(Arc::new(f)),
        }
    }

    /// Construct a committed R1CS instance from a **fully materialized witness vector** `f`.
    ///
    /// This bypasses `gadget_decompose` on an external `z` and is intended for integrations where the
    /// witness is already represented in the ring domain expected by the protocol (e.g. SP1 witness
    /// values embedded as constant-coeff ring elements).
    ///
    /// NOTE: `z` is left empty in this constructor (it is not used by the verifier-side relation).
    pub fn from_f(r1cs: R1CS<R>, f: Vec<R>, l_in: usize, A: &Matrix<R>) -> Self {
        let cm_f = A.try_mul_vec(&f).unwrap();
        let l = r1cs.l;
        let x = ComR1CSX {
            a: Arc::new(r1cs.A),
            b: Arc::new(r1cs.B),
            c: Arc::new(r1cs.C),
            z: Vec::new(),
            cm_f,
            l_in,
            l,
        };
        Self {
            x,
            f: WitnessVec::Ring(Arc::new(f)),
        }
    }

    /// Construct a committed R1CS instance from `f`, using a seeded implicit Ajtai matrix.
    pub fn from_f_seeded(r1cs: R1CS<R>, f: Vec<R>, l_in: usize, scheme: &AjtaiCommitmentScheme<R>) -> Self
    where
        R: stark_rings::PolyRing,
        R::BaseRing: stark_rings::Ring,
        R: core::ops::Mul<R::BaseRing, Output = R>,
    {
        let cm_f = scheme
            .commit_const_coeff_fast(&f)
            .expect("commit_const_coeff_fast")
            .as_ref()
            .to_vec();
        let l = r1cs.l;
        let x = ComR1CSX {
            a: Arc::new(r1cs.A),
            b: Arc::new(r1cs.B),
            c: Arc::new(r1cs.C),
            z: Vec::new(),
            cm_f,
            l_in,
            l,
        };
        Self {
            x,
            f: WitnessVec::Ring(Arc::new(f)),
        }
    }

    /// Construct a committed R1CS instance from a **constant-coefficient** witness (base scalars),
    /// using a seeded implicit Ajtai matrix.
    ///
    /// This avoids allocating `Vec<R>` for the witness in the common SP1 regime.
    pub fn from_f0_seeded(
        r1cs: R1CS<R>,
        f0: Arc<Vec<R::BaseRing>>,
        l_in: usize,
        scheme: &AjtaiCommitmentScheme<R>,
    ) -> Self
    where
        R::BaseRing: stark_rings::Ring,
        R: core::ops::Mul<R::BaseRing, Output = R>,
    {
        let n = scheme.width();
        // Commit as if `f0` were zero-padded to length `n`.
        let cm_f = scheme
            .commit_many_const_coeff_base_fast(n, 1, {
                let f0 = f0.clone();
                move |j, out| {
                    out[0] = f0.get(j).copied().unwrap_or(R::BaseRing::ZERO);
                }
            })
            .expect("commit_many_const_coeff_base_fast (f0 padded)")
            .into_iter()
            .next()
            .expect("t=1")
            .as_ref()
            .to_vec();
        let l = r1cs.l;
        let x = ComR1CSX {
            a: Arc::new(r1cs.A),
            b: Arc::new(r1cs.B),
            c: Arc::new(r1cs.C),
            z: Vec::new(),
            cm_f,
            l_in,
            l,
        };
        Self {
            x,
            f: WitnessVec::ConstCoeffBase {
                values: f0,
                domain_len: n,
            },
        }
    }
}

impl<R: Ring + PolyRing> ComR1CSXBase<R> {
    pub fn matrices_arc_base(&self) -> Vec<Arc<SparseMatrix<R::BaseRing>>> {
        vec![self.a.clone(), self.b.clone(), self.c.clone()]
    }
}

impl<R: Decompose + Ring + PolyRing> ComR1CSBase<R> {
    /// Construct a committed R1CS instance from **base-ring matrices** (const-coeff) and a
    /// constant-coefficient witness (base scalars), using a seeded implicit Ajtai matrix.
    pub fn from_f0_seeded_base(
        r1cs: R1CS<R::BaseRing>,
        f0: Arc<Vec<R::BaseRing>>,
        l_in: usize,
        scheme: &AjtaiCommitmentScheme<R>,
    ) -> Self
    where
        R::BaseRing: stark_rings::Ring,
        R: core::ops::Mul<R::BaseRing, Output = R>,
    {
        let n = scheme.width();
        let cm_f = scheme
            .commit_many_const_coeff_base_fast(n, 1, {
                let f0 = f0.clone();
                move |j, out| {
                    out[0] = f0.get(j).copied().unwrap_or(R::BaseRing::ZERO);
                }
            })
            .expect("commit_many_const_coeff_base_fast (f0 padded)")
            .into_iter()
            .next()
            .expect("t=1")
            .as_ref()
            .to_vec();
        let l = r1cs.l;
        let x = ComR1CSXBase {
            a: Arc::new(r1cs.A),
            b: Arc::new(r1cs.B),
            c: Arc::new(r1cs.C),
            cm_f,
            l_in,
            l,
        };
        Self {
            x,
            f: WitnessVec::ConstCoeffBase {
                values: f0,
                domain_len: n,
            },
        }
    }
}

impl<R: Ring> ComR1CSX<R> {
    pub fn matrices_arc(&self) -> Vec<Arc<SparseMatrix<R>>> {
        vec![self.a.clone(), self.b.clone(), self.c.clone()]
    }

    /// Legacy helper for small tests/debugging: materialize a temporary `R1CS<R>` by cloning matrices.
    pub fn r1cs_cloned(&self) -> R1CS<R> {
        R1CS {
            l: self.l,
            A: (*self.a).clone(),
            B: (*self.b).clone(),
            C: (*self.c).clone(),
        }
    }
}

impl<R: OverField + PolyRing> Linearize<R> for ComR1CS<R> {
    type Proof = ComR1CSProof<R>;
    fn linearize(&self, transcript: &mut impl Transcript<R>) -> (LinB<R>, Self::Proof) {
        let n = self.f.len().next_power_of_two();
        let nvars = log2(n) as usize;

        // Streaming sumcheck (memory-friendly) producing the same `Proof<R>` type.
        let r0 = transcript.get_challenges(nvars);
        let one_minus_r0 = r0.iter().copied().map(|x| R::BaseRing::ONE - x).collect();

        let mles = match &self.f {
            WitnessVec::Ring(f) => {
                let ga = self.x.a.try_mul_vec(f.as_ref()).unwrap();
                let gb = self.x.b.try_mul_vec(f.as_ref()).unwrap();
                let gc = self.x.c.try_mul_vec(f.as_ref()).unwrap();
                vec![
                    // eq(x, r0) (constant-coeff)
                    StreamingMleEnum::EqBase {
                        scale: R::BaseRing::ONE,
                        r: r0,
                        one_minus_r: one_minus_r0,
                    },
                    StreamingMleEnum::DenseArc {
                        evals: Arc::new(ga),
                        num_vars: nvars,
                    },
                    StreamingMleEnum::DenseArc {
                        evals: Arc::new(gb),
                        num_vars: nvars,
                    },
                    StreamingMleEnum::DenseArc {
                        evals: Arc::new(gc),
                        num_vars: nvars,
                    },
                    StreamingMleEnum::DenseArc {
                        evals: f.clone(),
                        num_vars: nvars,
                    },
                ]
            }
            WitnessVec::ConstCoeffBase { values: f0, .. } => {
                // Constant-coefficient witness: represent f as base scalars and avoid allocating `Vec<R>`.
                //
                // NOTE: This does not change the verifier relation; it only changes how the prover
                // evaluates MLEs during sumcheck.
                vec![
                    StreamingMleEnum::EqBase {
                        scale: R::BaseRing::ONE,
                        r: r0,
                        one_minus_r: one_minus_r0,
                    },
                    // For now we keep ga/gb/gc as dense vectors (computed on demand by downstream callers).
                    // Large-scale SP1 proving uses chunking; this path is intended for const-coeff chunks.
                    StreamingMleEnum::SparseMatVecConstCoeff {
                        matrix: self.x.a.clone(),
                        witness0: f0.clone(),
                        num_vars: nvars,
                    },
                    StreamingMleEnum::SparseMatVecConstCoeff {
                        matrix: self.x.b.clone(),
                        witness0: f0.clone(),
                        num_vars: nvars,
                    },
                    StreamingMleEnum::SparseMatVecConstCoeff {
                        matrix: self.x.c.clone(),
                        witness0: f0.clone(),
                        num_vars: nvars,
                    },
                    StreamingMleEnum::BaseScalarArc {
                        evals: f0.clone(),
                        num_vars: nvars,
                        square: false,
                    },
                ]
            }
        };

        let comb_fn = |vals: &[R]| -> R { vals[0] * (vals[1] * vals[2] - vals[3]) };

        let (sumcheck_proof, randomness, final_vals) =
            StreamingSumcheck::prove_as_subprotocol(transcript, mles, nvars, 3, comb_fn);

        let ro = randomness.into_iter().map(|x| x.into()).collect::<Vec<R>>();
        let va = final_vals[1];
        let vb = final_vals[2];
        let vc = final_vals[3];
        let v = final_vals[4];

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

impl<R: OverField + PolyRing> Linearize<R> for ComR1CSBase<R> {
    type Proof = ComR1CSProof<R>;
    fn linearize(&self, transcript: &mut impl Transcript<R>) -> (LinB<R>, Self::Proof) {
        let n = self.f.len().next_power_of_two();
        let nvars = log2(n) as usize;

        let r0 = transcript.get_challenges(nvars);
        let one_minus_r0 = r0.iter().copied().map(|x| R::BaseRing::ONE - x).collect();

        let f0 = match &self.f {
            WitnessVec::ConstCoeffBase { values, .. } => values.clone(),
            _ => panic!("ComR1CSBase requires ConstCoeffBase witness"),
        };

        let mles = vec![
            StreamingMleEnum::EqBase {
                scale: R::BaseRing::ONE,
                r: r0,
                one_minus_r: one_minus_r0,
            },
            StreamingMleEnum::SparseMatVecConstCoeffBase {
                matrix: self.x.a.clone(),
                witness0: f0.clone(),
                num_vars: nvars,
            },
            StreamingMleEnum::SparseMatVecConstCoeffBase {
                matrix: self.x.b.clone(),
                witness0: f0.clone(),
                num_vars: nvars,
            },
            StreamingMleEnum::SparseMatVecConstCoeffBase {
                matrix: self.x.c.clone(),
                witness0: f0.clone(),
                num_vars: nvars,
            },
            StreamingMleEnum::BaseScalarArc {
                evals: f0.clone(),
                num_vars: nvars,
                square: false,
            },
        ];

        let comb_fn = |vals: &[R]| -> R { vals[0] * (vals[1] * vals[2] - vals[3]) };

        let (sumcheck_proof, randomness, final_vals) =
            StreamingSumcheck::prove_as_subprotocol(transcript, mles, nvars, 3, comb_fn);

        let ro = randomness.into_iter().map(|x| x.into()).collect::<Vec<R>>();
        let va = final_vals[1];
        let vb = final_vals[2];
        let vc = final_vals[3];
        let v = final_vals[4];

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
        let v_pairs = vec![(v, v), (va, va), (vb, vb), (vc, vc)];

        let x = LinBX {
            cm_f: self.x.cm_f.clone(),
            r,
            v: v_pairs,
        };
        let linb = LinB { f: self.f.clone(), x };

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
    use cyclotomic_rings::rings::GoldilocksPoseidonConfig as PC;
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
        let f_ring = cr1cs.f.as_ring_arc().expect("test uses ring witness");
        cr1cs.x.r1cs_cloned().check_relation(f_ring.as_ref()).unwrap();

        let mut ts = PoseidonTranscript::empty::<PC>();
        let (_linb, lproof) = cr1cs.linearize(&mut ts);

        let mut ts = PoseidonTranscript::empty::<PC>();
        lproof.verify(&mut ts);
    }
}
