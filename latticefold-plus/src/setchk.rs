use ark_std::log2;
use latticefold::{
    transcript::Transcript,
    utils::sumcheck::{
        utils::{build_eq_x_r, eq_eval},
        MLSumcheck, Proof, SumCheckError,
    },
};
use stark_rings::{OverField, PolyRing, Ring};
use stark_rings_linalg::{ops::Transpose, SparseMatrix};
use stark_rings_poly::mle::{DenseMultilinearExtension, SparseMultilinearExtension};
use thiserror::Error;

// cM: double commitment, commitment to M
// M: witness matrix of monomials

#[derive(Debug)]
pub enum Arity<T> {
    Single(T),
    Batch(Vec<T>),
}

impl<T> Arity<T> {
    pub fn as_slice(&self) -> &[T] {
        match self {
            Self::Single(a) => ark_std::slice::from_ref(a),
            Self::Batch(a) => a.as_slice(),
        }
    }
}

impl<T> From<Vec<T>> for Arity<T> {
    fn from(v: Vec<T>) -> Self {
        if v.len() == 1 {
            Self::Single(v.into_iter().next().unwrap())
        } else {
            Self::Batch(v)
        }
    }
}

#[derive(Debug)]
pub struct In<R> {
    pub nvars: usize,
    pub M: Arity<SparseMatrix<R>>, // n x m
}

#[derive(Debug)]
pub struct Out<R: Ring> {
    pub nvars: usize,
    pub r: Vec<R>, // log n
    pub sumcheck_proof: Proof<R>,
    pub e: Arity<Vec<R>>, // m
}

#[derive(Debug, Error)]
pub enum SetCheckError<R: Ring> {
    #[error("Sumcheck failed: {0}")]
    Sumcheck(#[from] SumCheckError<R>),
    #[error("Recomputed claim `v` mismatch")]
    ExpectedClaim,
}

fn ev<R: PolyRing>(r: &R, x: R::BaseRing) -> R::BaseRing {
    r.coeffs()
        .iter()
        .enumerate()
        .map(|(i, c)| *c * x.pow([i as u64]))
        .sum()
}

impl<R: OverField> In<R> {
    /// Monomial set check
    pub fn set_check(&self, transcript: &mut impl Transcript<R>) -> Out<R> {
        let Ms = self.M.as_slice();
        let ncols = Ms[0].ncols;
        let MTs = Ms.iter().map(|M| M.transpose()).collect::<Vec<_>>();
        let tnvars = log2(Ms[0].nrows.next_power_of_two()) as usize;

        let mut mles = Vec::with_capacity(Ms.len() * (ncols * 2 + 1));
        let mut alphas = Vec::with_capacity(Ms.len());

        // loop for batch support
        for M in Ms {
            // Step 1
            let c: Vec<R> = transcript
                .get_challenges(self.nvars)
                .into_iter()
                .map(|x| x.into())
                .collect();
            let beta = transcript.get_challenge();

            // Step 2
            let MT = M.transpose();

            // explore sMLE
            for row in MT.coeffs.iter() {
                let mut m_j = vec![R::zero(); M.nrows];
                row.iter().for_each(|(r, i)| m_j[*i] = R::from(ev(r, beta)));
                // ev(x^2) = ev(x)^2, iif monomial
                let m_prime_j = m_j.iter().map(|z| *z * z).collect::<Vec<_>>();

                let mle_m_j = DenseMultilinearExtension::from_evaluations_vec(tnvars, m_j);
                let mle_m_prime_j =
                    DenseMultilinearExtension::from_evaluations_vec(tnvars, m_prime_j);

                mles.push(mle_m_j);
                mles.push(mle_m_prime_j);
            }

            let eq = build_eq_x_r(&c).unwrap();
            mles.push(eq);

            let alpha: R = transcript.get_challenge().into();
            alphas.push(alpha);
        }

        // random linear combinator, for batching
        let rc: Option<R::BaseRing> =
            matches!(self.M, Arity::Batch(_)).then(|| transcript.get_challenge().into());

        let comb_fn = |vals: &[R]| -> R {
            let mut lc = R::zero();
            for i in 0..Ms.len() {
                // 2 * ncols for (m_j, m_prime_j), +1 for eq
                let s = i * (2 * ncols + 1);
                let mut res = R::zero();
                for j in 0..ncols {
                    res += alphas[i].pow([j as u64])
                        * (vals[s + j * 2] * vals[s + j * 2] - vals[s + j * 2 + 1])
                }
                res *= vals[s + 2 * ncols]; // eq
                lc += if let Some(rc) = &rc {
                    res * rc.pow([i as u64])
                } else {
                    return res;
                };
            }
            lc
        };

        let (sumcheck_proof, prover_state) =
            MLSumcheck::prove_as_subprotocol(transcript, mles, self.nvars, 3, comb_fn);

        let r = prover_state
            .randomness
            .into_iter()
            .map(|x| x.into())
            .collect::<Vec<R>>();

        let e: Arity<Vec<R>> = MTs
            .iter()
            .map(|MT| {
                MT.coeffs
                    .iter()
                    .map(|row| {
                        let evals: Vec<(usize, R)> = row.iter().map(|&(r, i)| (i, r)).collect();
                        let mle = SparseMultilinearExtension::from_evaluations(tnvars, &evals);
                        mle.evaluate(&r)
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<Vec<_>>>() // batch
            .into();

        Out {
            nvars: self.nvars,
            e,
            r,
            sumcheck_proof,
        }
    }
}

impl<R: OverField> Out<R> {
    pub fn verify(&self, transcript: &mut impl Transcript<R>) -> Result<(), SetCheckError<R>> {
        let es = self.e.as_slice();

        let cba: Vec<(Vec<R>, R::BaseRing, R)> = (0..es.len())
            .map(|_| {
                let c: Vec<R> = transcript
                    .get_challenges(self.nvars)
                    .into_iter()
                    .map(|x| x.into())
                    .collect();
                let beta = transcript.get_challenge();
                let alpha: R = transcript.get_challenge().into();
                (c, beta, alpha)
            })
            .collect();

        let rc: Option<R::BaseRing> =
            matches!(self.e, Arity::Batch(_)).then(|| transcript.get_challenge().into());

        let subclaim = MLSumcheck::verify_as_subprotocol(
            transcript,
            self.nvars,
            3,
            R::zero(),
            &self.sumcheck_proof,
        )?;

        let r: Vec<R> = subclaim.point.into_iter().map(|x| x.into()).collect();

        let v = subclaim.expected_evaluation;

        use ark_std::One;
        let mut ver = R::zero();
        for (i, e) in es.iter().enumerate() {
            let c = &cba[i].0;
            let beta = &cba[i].1;
            let alpha = &cba[i].2;
            let eq = eq_eval(&c, &r).unwrap();
            let e_sum = e
                .iter()
                .enumerate()
                .map(|(j, e_j)| {
                    let ev1 = R::from(ev(e_j, *beta));
                    let ev2 = R::from(ev(e_j, *beta * beta));
                    alpha.pow([j as u64]) * (ev1 * ev1 - ev2)
                })
                .sum::<R>();
            ver += eq * e_sum * rc.as_ref().unwrap_or(&R::BaseRing::one()).pow([i as u64]);
        }

        (ver == v).then(|| ()).ok_or(SetCheckError::ExpectedClaim)?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use ark_std::One;
    use cyclotomic_rings::rings::GoldilocksPoseidonConfig as PC;
    use latticefold::transcript::poseidon::PoseidonTS;
    use stark_rings::cyclotomic_ring::models::goldilocks::RqPoly as R;
    use stark_rings_linalg::SparseMatrix;

    use super::*;

    #[test]
    fn test_set_check() {
        let n = 4;
        let M = SparseMatrix::<R>::identity(n);

        let scin = In {
            M: Arity::Single(M),
            nvars: log2(n) as usize,
        };

        let mut ts = PoseidonTS::default::<PC>();
        let out = scin.set_check(&mut ts);

        let mut ts = PoseidonTS::default::<PC>();
        out.verify(&mut ts).unwrap();
    }

    #[test]
    fn test_set_check_bad() {
        let n = 4;
        let mut M = SparseMatrix::<R>::identity(n);
        // 1 + X, not a monomial
        let mut onepx = R::one();
        onepx.coeffs_mut()[1] = 1u128.into();
        M.coeffs[0][0].0 = onepx;

        let scin = In {
            M: Arity::Single(M),
            nvars: log2(n) as usize,
        };

        let mut ts = PoseidonTS::default::<PC>();
        let out = scin.set_check(&mut ts);

        let mut ts = PoseidonTS::default::<PC>();
        assert!(out.verify(&mut ts).is_err());
    }

    #[test]
    fn test_set_check_batched() {
        let n = 4;
        let M0 = SparseMatrix::<R>::identity(n);
        let M1 = SparseMatrix::<R>::identity(n);

        let scin = In {
            M: Arity::Batch(vec![M0, M1]),
            nvars: log2(n) as usize,
        };

        let mut ts = PoseidonTS::default::<PC>();
        let out = scin.set_check(&mut ts);

        let mut ts = PoseidonTS::default::<PC>();
        out.verify(&mut ts).unwrap();
    }

    #[test]
    fn test_set_check_batched_bad() {
        let n = 4;
        let M0 = SparseMatrix::<R>::identity(n);
        let mut M1 = SparseMatrix::<R>::identity(n);
        // 1 + X, not a monomial
        let mut onepx = R::one();
        onepx.coeffs_mut()[1] = 1u128.into();
        M1.coeffs[0][0].0 = onepx;

        let scin = In {
            M: Arity::Batch(vec![M0, M1]),
            nvars: log2(n) as usize,
        };

        let mut ts = PoseidonTS::default::<PC>();
        let out = scin.set_check(&mut ts);

        let mut ts = PoseidonTS::default::<PC>();
        assert!(out.verify(&mut ts).is_err());
    }
}
