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

#[derive(Clone, Debug)]
pub enum MonomialSet<R> {
    Matrix(SparseMatrix<R>),
    Vector(Vec<R>),
}

#[derive(Clone, Debug)]
pub struct In<R> {
    pub nvars: usize,
    pub sets: Vec<MonomialSet<R>>, // Ms and ms: n x m, or n
}

#[derive(Clone, Debug)]
pub struct Out<R: PolyRing> {
    pub nvars: usize,
    pub r: Vec<R::BaseRing>, // log n
    pub sumcheck_proof: Proof<R>,
    pub e: Vec<Vec<Vec<R>>>, // m, matrices outputs
    pub b: Vec<R>,           // vectors outputs
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
    ///
    /// Proves sets rings are all unit monomials.
    /// Currently requires k >= 1 monomial matrices sets. TODO support other scenarios.
    /// If k > 1, sumcheck batching is employed.
    pub fn set_check(&self, M: &[SparseMatrix<R>], transcript: &mut impl Transcript<R>) -> Out<R> {
        let Ms: Vec<&SparseMatrix<R>> = self
            .sets
            .iter()
            .filter_map(|set| match set {
                MonomialSet::Matrix(m) => Some(m),
                _ => None,
            })
            .collect();
        let ms: Vec<&Vec<R>> = self
            .sets
            .iter()
            .filter_map(|set| match set {
                MonomialSet::Vector(v) => Some(v),
                _ => None,
            })
            .collect();

        let ncols = Ms[0].ncols;
        let MTs = Ms.iter().map(|M| M.transpose()).collect::<Vec<_>>();
        let tnvars = log2(Ms[0].nrows.next_power_of_two()) as usize;

        let mut mles = Vec::with_capacity((Ms.len() + ms.len()) * (ncols * 2 + 1));
        let mut alphas = Vec::with_capacity(Ms.len());

        // matrix sets
        for M in Ms.iter() {
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

        // vector sets
        for m in ms.iter() {
            // Step 1
            let c: Vec<R> = transcript
                .get_challenges(self.nvars)
                .into_iter()
                .map(|x| x.into())
                .collect();
            let beta = transcript.get_challenge();

            let m_j = m.iter().map(|r| R::from(ev(r, beta))).collect::<Vec<_>>();
            // ev(x^2) = ev(x)^2, iif monomial
            let m_prime_j = m_j.iter().map(|z| *z * z).collect::<Vec<_>>();

            let mle_m_j = DenseMultilinearExtension::from_evaluations_vec(tnvars, m_j);
            let mle_m_prime_j = DenseMultilinearExtension::from_evaluations_vec(tnvars, m_prime_j);

            mles.push(mle_m_j);
            mles.push(mle_m_prime_j);

            let eq = build_eq_x_r(&c).unwrap();
            mles.push(eq);

            let alpha: R = transcript.get_challenge().into();
            alphas.push(alpha);
        }

        // random linear combinator, for batching
        let rc: Option<R::BaseRing> = (Ms.len() > 1).then(|| transcript.get_challenge());

        let comb_fn = |vals: &[R]| -> R {
            let mut lc = R::zero();
            for (i, alpha) in alphas.iter().enumerate().take(Ms.len()) {
                // 2 * ncols for (m_j, m_prime_j), +1 for eq
                let s = i * (2 * ncols + 1);
                let mut res = R::zero();
                for j in 0..ncols {
                    res += alpha.pow([j as u64])
                        * (vals[s + j * 2] * vals[s + j * 2] - vals[s + j * 2 + 1])
                }
                res *= vals[s + 2 * ncols]; // eq
                lc += if let Some(rc) = &rc {
                    res * rc.pow([i as u64])
                } else {
                    return res;
                };
            }
            for i in 0..ms.len() {
                let s_base = Ms.len() * (2 * ncols + 1);
                let s = s_base + i * 3;
                let mut res = R::zero();
                let alpha_idx = Ms.len() + i;
                res += alphas[alpha_idx] * (vals[s] * vals[s] - vals[s + 1]);
                res *= vals[s + 2]; // eq
                lc += if let Some(rc) = &rc {
                    res * rc.pow([alpha_idx as u64])
                } else {
                    return res;
                };
            }
            lc
        };

        let (sumcheck_proof, prover_state) =
            MLSumcheck::prove_as_subprotocol(transcript, mles, self.nvars, 3, comb_fn);

        let r = prover_state.randomness.clone();
        let r_poly = prover_state
            .randomness
            .into_iter()
            .map(|x| x.into())
            .collect::<Vec<R>>();

        // Step 3
        let e: Vec<Vec<Vec<R>>> = {
            let mut e = Vec::with_capacity(1 + M.len());

            let e0 = MTs
                .iter()
                .map(|MT| {
                    MT.coeffs
                        .iter()
                        .map(|row| {
                            let evals: Vec<(usize, R)> = row.iter().map(|&(r, i)| (i, r)).collect();
                            let mle = SparseMultilinearExtension::from_evaluations(tnvars, &evals);
                            mle.evaluate(&r_poly)
                        })
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<Vec<R>>>();
            e.push(e0);

            // Mf
            for Mi in M {
                let ei = MTs
                    .iter()
                    .map(|MT| {
                        MT.coeffs
                            .iter()
                            .map(|row| {
                                let mut drow = vec![R::zero(); MT.ncols];
                                row.iter().for_each(|&(r, i)| {
                                    drow[i] = r;
                                });
                                let evals = Mi.try_mul_vec(&drow).unwrap();
                                let mle =
                                    DenseMultilinearExtension::from_evaluations_vec(tnvars, evals);
                                mle.evaluate(&r_poly).unwrap()
                            })
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<Vec<R>>>();
                e.push(ei);
            }
            e
        };

        let b: Vec<R> = ms
            .iter()
            .map(|m| {
                let mle = DenseMultilinearExtension::from_evaluations_slice(tnvars, m);
                mle.evaluate(&r_poly).unwrap()
            })
            .collect::<Vec<_>>();

        // Prover to Verifier messages
        absorb_evaluations(&e, &b, transcript);

        Out {
            nvars: self.nvars,
            e,
            b,
            r,
            sumcheck_proof,
        }
    }
}

impl<R: OverField> Out<R> {
    pub fn verify(&self, transcript: &mut impl Transcript<R>) -> Result<(), SetCheckError<R>> {
        let nclaims = self.e[0].len() + self.b.len();

        let cba: Vec<(Vec<R>, R::BaseRing, R)> = (0..nclaims)
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

        let rc: Option<R::BaseRing> = (self.e[0].len() > 1).then(|| transcript.get_challenge());

        let subclaim = MLSumcheck::verify_as_subprotocol(
            transcript,
            self.nvars,
            3,
            R::zero(),
            &self.sumcheck_proof,
        )?;

        let r: Vec<R> = subclaim.point.into_iter().map(|x| x.into()).collect();

        let v = subclaim.expected_evaluation;

        // Prover to Verifier messages
        absorb_evaluations(&self.e, &self.b, transcript);

        use ark_std::One;
        let mut ver = R::zero();
        for (i, e) in self.e[0].iter().enumerate() {
            let c = &cba[i].0;
            let beta = &cba[i].1;
            let alpha = &cba[i].2;
            let eq = eq_eval(c, &r).unwrap();
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
        for (i, b) in self.b.iter().enumerate() {
            let offset = self.e[0].len();
            let c = &cba[i + offset].0;
            let beta = &cba[i + offset].1;
            let alpha = &cba[i + offset].2;
            let eq = eq_eval(c, &r).unwrap();
            let b_claim = {
                let ev1 = R::from(ev(b, *beta));
                let ev2 = R::from(ev(b, *beta * *beta));
                ev1 * ev1 - ev2
            };
            ver += eq
                * *alpha
                * b_claim
                * rc.as_ref()
                    .unwrap_or(&R::BaseRing::one())
                    .pow([(i + offset) as u64]);
        }

        (ver == v)
            .then_some(())
            .ok_or(SetCheckError::ExpectedClaim)?;

        Ok(())
    }
}

fn absorb_evaluations<R: OverField>(
    e: &[Vec<Vec<R>>],
    b: &[R],
    transcript: &mut impl Transcript<R>,
) {
    for ek in e {
        for ej in ek {
            transcript.absorb_slice(ej);
        }
    }
    transcript.absorb_slice(b);
}

#[cfg(test)]
mod tests {
    use ark_std::One;
    use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
    use stark_rings::{cyclotomic_ring::models::frog_ring::RqPoly as R, unit_monomial};
    use stark_rings_linalg::SparseMatrix;

    use super::*;
    use crate::transcript::PoseidonTranscript;

    #[test]
    fn test_set_check() {
        let n = 4;
        let M = SparseMatrix::<R>::identity(n);

        let scin = In {
            sets: vec![MonomialSet::Matrix(M)],
            nvars: log2(n) as usize,
        };

        let mut ts = PoseidonTranscript::empty::<PC>();
        let out = scin.set_check(&[], &mut ts);

        let mut ts = PoseidonTranscript::empty::<PC>();
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
            sets: vec![MonomialSet::Matrix(M)],
            nvars: log2(n) as usize,
        };

        let mut ts = PoseidonTranscript::empty::<PC>();
        let out = scin.set_check(&[], &mut ts);

        let mut ts = PoseidonTranscript::empty::<PC>();
        assert!(out.verify(&mut ts).is_err());
    }

    #[test]
    fn test_set_check_batched() {
        let n = 4;
        let M0 = SparseMatrix::<R>::identity(n);
        let M1 = SparseMatrix::<R>::identity(n);

        let scin = In {
            sets: vec![MonomialSet::Matrix(M0), MonomialSet::Matrix(M1)],
            nvars: log2(n) as usize,
        };

        let mut ts = PoseidonTranscript::empty::<PC>();
        let out = scin.set_check(&[], &mut ts);

        let mut ts = PoseidonTranscript::empty::<PC>();
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
            sets: vec![MonomialSet::Matrix(M0), MonomialSet::Matrix(M1)],
            nvars: log2(n) as usize,
        };

        let mut ts = PoseidonTranscript::empty::<PC>();
        let out = scin.set_check(&[], &mut ts);

        let mut ts = PoseidonTranscript::empty::<PC>();
        assert!(out.verify(&mut ts).is_err());
    }

    #[test]
    fn test_set_check_mix() {
        let n = 4;
        let M0 = SparseMatrix::<R>::identity(n);
        let M1 = SparseMatrix::<R>::identity(n);
        let m0 = vec![R::one(); n];
        let m1 = vec![unit_monomial(2); n];

        let scin = In {
            sets: vec![
                MonomialSet::Matrix(M0),
                MonomialSet::Matrix(M1),
                MonomialSet::Vector(m0),
                MonomialSet::Vector(m1),
            ],
            nvars: log2(n) as usize,
        };

        let mut ts = PoseidonTranscript::empty::<PC>();
        let out = scin.set_check(&[], &mut ts);

        let mut ts = PoseidonTranscript::empty::<PC>();
        out.verify(&mut ts).unwrap();
    }

    #[test]
    fn test_set_check_mix_bad() {
        let n = 4;
        let M0 = SparseMatrix::<R>::identity(n);
        let M1 = SparseMatrix::<R>::identity(n);
        let mut m0 = vec![R::one(); n];
        let mut onepx = R::one();
        onepx.coeffs_mut()[1] = 1u128.into();
        m0[0] = onepx;

        let scin = In {
            sets: vec![
                MonomialSet::Matrix(M0),
                MonomialSet::Matrix(M1),
                MonomialSet::Vector(m0),
            ],
            nvars: log2(n) as usize,
        };

        let mut ts = PoseidonTranscript::empty::<PC>();
        let out = scin.set_check(&[], &mut ts);

        let mut ts = PoseidonTranscript::empty::<PC>();
        assert!(out.verify(&mut ts).is_err());
    }
}
