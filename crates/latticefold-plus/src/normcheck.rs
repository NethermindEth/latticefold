//! ell-2 norm-check (ePrint 2026/721, section 4), per-instance (unbatched).
//!
//! This module provides a [`NormCheck`] strategy abstraction with two
//! implementations:
//!
//! - [`LinfNormCheck`]: the existing ell-infinity monomial range check
//!   ([`crate::rgchk::Rg::range_check`]), wrapped for a fair comparison baseline.
//! - [`L2NormCheck`]: the ell-2 check built from two reductions of knowledge,
//!   `Pi^exact` (exact shortening via complex conjugation, section 4.4) followed
//!   by `Pi^proj` (JL random projection, section 4.3 "Attempt II").
//!
//! Forward (proving) order follows the composition diagram on page 10:
//! `Pi^exact` then `Pi^proj`. Extraction (reverse) upgrades a relaxed witness to
//! shortish (`Pi^proj`) then to exactly short (`Pi^exact`).
//!
//! Scope and honesty notes:
//!
//! - This is the per-instance version (paper Remark 1); batching across the `L`
//!   instances (shared `Pi`, batched inner-product sumcheck) is a follow-up.
//! - SOUNDNESS: 2026/721 provides no formal proofs; it states security "follows
//!   verbatim" from BC25b / KLOT25 / KLNO25. Every load-bearing constant or step
//!   below carries a `// SOUNDNESS:` cite for human review.
//! - The standalone gadget verifier checks the protocol's local acceptance
//!   predicates (norm bounds, projection/folding consistency, sumcheck validity).
//!   Binding the opened sumcheck evaluations to the Ajtai commitments `cm_w` /
//!   `cm_v` (and hence full knowledge soundness / extraction) is provided by the
//!   committed-relation flow in [`crate::cm`] when this is integrated (Phase C):
//!   here those evaluations are carried in the proof as the claims that the
//!   homogenization step would forward.

use ark_std::log2;
use latticefold::{
    transcript::Transcript,
    utils::sumcheck::{MLSumcheck, Proof},
};
use stark_rings::{CoeffRing, OverField, PolyRing, Zq};
use stark_rings_linalg::{Matrix, SparseMatrix};
use stark_rings_poly::mle::DenseMultilinearExtension;
use thiserror::Error;

use crate::{
    conjugation::conjugate,
    rgchk::{Dcom, DecompParameters, Rg, RgInstance},
    utils::{short_challenge, tensor_product},
};

/// Number of JL projection rows.
///
/// SOUNDNESS: fixed by Lemma 1 of 2026/721 (Corollary 3.2 [GHL22] with
/// Corollary 4.2 [BS23]). Do not tune to pass tests.
pub const JL_ROWS: usize = 256;

/// JL squared lower factor (`sqrt(30)^2`).
///
/// SOUNDNESS: fixed by Lemma 1 of 2026/721. `sqrt(30) ||w|| < ||Pi w|| < sqrt(337) ||w||`
/// except w.p. `2^-128` over `Pi <- C^{256 x d}`, for `||w|| <= q/125`. The
/// union bound over fold steps (caveat 1, account per 2025/1220) is not folded
/// into these constants and must be handled by parameter selection.
pub const JL_LOWER_SQ: u128 = 30;
/// JL squared upper factor (`sqrt(337)^2`).
///
/// SOUNDNESS: fixed by Lemma 1 of 2026/721 (see [`JL_LOWER_SQ`]).
pub const JL_UPPER_SQ: u128 = 337;

#[derive(Debug, Error)]
pub enum NormCheckError {
    #[error("sumcheck verification failed")]
    Sumcheck,
    #[error("ell-2 norm exceeds bound: ct(u) = {got} > beta^2 = {bound}")]
    NormExceeded { got: u128, bound: u128 },
    #[error("projected witness norm outside JL band")]
    ProjectionNorm,
    #[error("folding consistency failed: sum_j c0_j t_j != s")]
    FoldingConsistency,
    #[error("exact-shortening claim mismatch: u' * conj(u'') != u*")]
    ExactClaim,
    #[error("projection claim mismatch: MLE eval inconsistent with sumcheck")]
    ProjectionClaim,
    #[error("ell-infinity range check failed")]
    RangeCheck,
}

/// Strategy abstraction over the two norm-check paths.
///
/// `prove` consumes the `L` witnesses `w` (each in `R^m`), the Ajtai matrix `A`,
/// public `params`, and a transcript, producing a proof. `verify` re-derives all
/// Fiat-Shamir challenges from the transcript and checks the proof.
pub trait NormCheck<R: CoeffRing>
where
    R::BaseRing: Zq,
{
    type Params;
    type Proof;

    fn prove(
        &self,
        w: &[Vec<R>],
        A: &Matrix<R>,
        params: &Self::Params,
        transcript: &mut impl Transcript<R>,
    ) -> Self::Proof;

    fn verify(
        &self,
        proof: &Self::Proof,
        params: &Self::Params,
        transcript: &mut impl Transcript<R>,
    ) -> Result<(), NormCheckError>;
}

// ---------------------------------------------------------------------------
// ell-infinity baseline (wraps the existing monomial range check)
// ---------------------------------------------------------------------------

/// ell-infinity baseline strategy: a thin wrapper over [`Rg::range_check`].
pub struct LinfNormCheck;

pub struct LinfParams {
    pub kappa: usize,
    pub decomp: DecompParameters,
}

pub struct LinfProof<R: PolyRing> {
    pub dcom: Dcom<R>,
}

impl<R: CoeffRing> NormCheck<R> for LinfNormCheck
where
    R::BaseRing: stark_rings::balanced_decomposition::convertible_ring::ConvertibleRing
        + stark_rings::balanced_decomposition::Decompose
        + Zq,
    R: stark_rings::balanced_decomposition::Decompose,
{
    type Params = LinfParams;
    type Proof = LinfProof<R>;

    fn prove(
        &self,
        w: &[Vec<R>],
        A: &Matrix<R>,
        params: &Self::Params,
        transcript: &mut impl Transcript<R>,
    ) -> Self::Proof {
        let nvars = log2(w[0].len()) as usize;
        let instances = w
            .iter()
            .map(|w_j| RgInstance::from_f(w_j.clone(), A, &params.decomp))
            .collect::<Vec<_>>();
        let rg = Rg {
            nvars,
            instances,
            dparams: params.decomp.clone(),
        };
        // No constraint matrices: this baseline measures the norm-check path only.
        let dcom = rg.range_check(&[] as &[SparseMatrix<R>], transcript);
        LinfProof { dcom }
    }

    fn verify(
        &self,
        proof: &Self::Proof,
        _params: &Self::Params,
        transcript: &mut impl Transcript<R>,
    ) -> Result<(), NormCheckError> {
        proof
            .dcom
            .verify(transcript)
            .map_err(|_| NormCheckError::RangeCheck)
    }
}

// ---------------------------------------------------------------------------
// ell-2 norm check
// ---------------------------------------------------------------------------

/// ell-2 norm-check strategy (`Pi^exact` then `Pi^proj`).
pub struct L2NormCheck;

#[derive(Clone, Debug)]
pub struct L2Params {
    /// Witness length per instance (must be a power of two).
    pub m: usize,
    /// Target squared ell-2 bound `beta^2` (one instance).
    pub beta_sq: u128,
}

#[derive(Debug, Error)]
pub enum L2ParamError {
    #[error("modulus too small for exact shortening: alpha^2 * beta^2 >= q = {q}")]
    ModulusTooSmall { q: u128 },
    #[error("parameter overflow computing alpha^2 * beta^2")]
    Overflow,
}

impl L2Params {
    /// Build ell-2 params from the target ell-2 bound `beta` and the norm growth
    /// factor `alpha` (the factor by which a fold inflates the norm before the
    /// exact-shortening step), validating the modulus.
    ///
    /// SOUNDNESS gate 6: the identity `||w||_2^2 = ct(<w, conj(w)>)` holds over
    /// the integers (not just mod q) only if `alpha^2 * beta^2 < q`. We fail
    /// loudly otherwise rather than silently proving a wrong statement.
    pub fn from_security(m: usize, q: u128, alpha: u128, beta: u128) -> Result<Self, L2ParamError> {
        let beta_sq = beta.checked_mul(beta).ok_or(L2ParamError::Overflow)?;
        let need = alpha
            .checked_mul(alpha)
            .and_then(|a2| a2.checked_mul(beta_sq))
            .ok_or(L2ParamError::Overflow)?;
        if need >= q {
            return Err(L2ParamError::ModulusTooSmall { q });
        }
        Ok(L2Params { m, beta_sq })
    }
}

#[derive(Clone, Debug)]
pub struct L2Proof<R: OverField> {
    /// Ajtai commitments of the `L` witnesses.
    pub cm_w: Vec<Vec<R>>,
    // Pi^exact
    /// `u_j = <w_j, conj(w_j)>`, with `ct(u_j) = ||w_j||_2^2`.
    pub u: Vec<R>,
    pub exact_sc: Vec<Proof<R>>,
    /// `u'_j = MLE[w_j](r*)`.
    pub u_prime: Vec<R>,
    /// `u''_j = MLE[w_j](conj r*)`.
    pub u_dprime: Vec<R>,
    // Pi^proj
    /// Ajtai commitment of the concatenated projection `v` (length `m`).
    pub cm_v: Vec<R>,
    /// `t_j = <c1, v_j> = c1^T (I_{m/b} (x) Pi) w_j`.
    pub t: Vec<R>,
    /// `s = (c0 (x) c1) v`.
    pub s: R,
    pub projw_sc: Vec<Proof<R>>,
    /// `MLE[w_j](ro_j)` opened for the `t_j` sumcheck.
    pub w_eval: Vec<R>,
    pub projv_sc: Proof<R>,
    /// `MLE[v](ro)` opened for the `s` sumcheck.
    pub v_eval: R,
}

/// Sample the JL projection `Pi <- C^{rows x cols}` from the transcript.
///
/// `C` puts mass `1/2` on `0` and `1/4` on each of `+1`, `-1`. We draw two bits
/// per entry: `{00, 01} -> 0`, `10 -> +1`, `11 -> -1`.
///
/// SOUNDNESS: the distribution `C` and the row count are fixed by Lemma 1 of
/// 2026/721. The caller MUST squeeze this only after the witness commitments are
/// absorbed (gate 1): otherwise a prover can choose `w` to pass a known `Pi`.
fn sample_pi<R: OverField>(
    rows: usize,
    cols: usize,
    transcript: &mut impl Transcript<R>,
) -> Vec<Vec<i8>> {
    let n = rows * cols;
    let nbytes = n.div_ceil(4);
    let bytes = transcript.squeeze_bytes(nbytes);
    let mut pi = vec![vec![0i8; cols]; rows];
    let mut idx = 0usize;
    for row in pi.iter_mut() {
        for entry in row.iter_mut() {
            let two_bits = (bytes[idx / 4] >> (2 * (idx % 4))) & 0b11;
            *entry = match two_bits {
                0 | 1 => 0,
                2 => 1,
                _ => -1,
            };
            idx += 1;
        }
    }
    pi
}

/// Project a single witness `w_j in R^m` block-wise: `v_j[i*256 + r] =
/// sum_p Pi[r][p] w_j[i*b + p]` over blocks `i in 0..m/b`. Result length `m/L`.
fn project<R: OverField>(w_j: &[R], pi: &[Vec<i8>], b: usize) -> Vec<R> {
    let rows = pi.len();
    let blocks = w_j.len() / b;
    let mut v_j = Vec::with_capacity(blocks * rows);
    for i in 0..blocks {
        let block = &w_j[i * b..(i + 1) * b];
        for pi_row in pi.iter() {
            let mut acc = R::zero();
            for (p, &coeff) in pi_row.iter().enumerate() {
                match coeff {
                    1 => acc += block[p],
                    -1 => acc -= block[p],
                    _ => {}
                }
            }
            v_j.push(acc);
        }
    }
    v_j
}

/// Build `a = (I_{m/b} (x) Pi)^T c1`, the length-`m` left vector of the `t_j`
/// sumcheck (same for all instances). `a[i*b + p] = sum_r Pi[r][p] c1[i*256 + r]`.
fn proj_left<R: OverField>(pi: &[Vec<i8>], c1: &[R], b: usize, m: usize) -> Vec<R> {
    let rows = pi.len();
    let blocks = m / b;
    let mut a = vec![R::zero(); m];
    for i in 0..blocks {
        for (r, pi_row) in pi.iter().enumerate() {
            let c1_ir = c1[i * rows + r];
            for (p, &coeff) in pi_row.iter().enumerate() {
                match coeff {
                    1 => a[i * b + p] += c1_ir,
                    -1 => a[i * b + p] -= c1_ir,
                    _ => {}
                }
            }
        }
    }
    a
}

fn lift_point<R: OverField>(point: &[R::BaseRing]) -> Vec<R> {
    point.iter().map(|x| R::from(*x)).collect()
}

impl<R: CoeffRing> NormCheck<R> for L2NormCheck
where
    R::BaseRing: Zq,
{
    type Params = L2Params;
    type Proof = L2Proof<R>;

    fn prove(
        &self,
        w: &[Vec<R>],
        A: &Matrix<R>,
        params: &Self::Params,
        transcript: &mut impl Transcript<R>,
    ) -> Self::Proof {
        let L = w.len();
        let m = params.m;
        let mu = log2(m) as usize;
        assert!(m.is_power_of_two(), "m must be a power of two");
        assert!(L.is_power_of_two(), "L must be a power of two");
        let b = JL_ROWS * L;
        assert!(m % b == 0, "m must be a multiple of b = 256 * L");

        // Commit each witness and absorb (gate 1: before any challenge).
        let cm_w = w
            .iter()
            .map(|w_j| A.try_mul_vec(w_j).unwrap())
            .collect::<Vec<Vec<R>>>();
        cm_w.iter().for_each(|c| transcript.absorb_slice(c));

        // ---- Pi^exact (section 4.4) ----
        let conj_w = w
            .iter()
            .map(|w_j| w_j.iter().map(|x| conjugate(x)).collect::<Vec<R>>())
            .collect::<Vec<_>>();
        let u = w
            .iter()
            .zip(&conj_w)
            .map(|(w_j, cw_j)| {
                w_j.iter()
                    .zip(cw_j)
                    .fold(R::zero(), |acc, (a, b)| acc + *a * *b)
            })
            .collect::<Vec<R>>();
        transcript.absorb_slice(&u);

        let mut exact_sc = Vec::with_capacity(L);
        let mut u_prime = Vec::with_capacity(L);
        let mut u_dprime = Vec::with_capacity(L);
        for j in 0..L {
            let w_mle = DenseMultilinearExtension::from_evaluations_vec(mu, w[j].clone());
            let cw_mle = DenseMultilinearExtension::from_evaluations_vec(mu, conj_w[j].clone());
            let comb = |vals: &[R]| -> R { vals[0] * vals[1] };
            let (proof, state) =
                MLSumcheck::prove_as_subprotocol(transcript, vec![w_mle, cw_mle], mu, 2, comb);
            let r_star: Vec<R> = state.randomness.iter().map(|x| R::from(*x)).collect();
            let conj_r_star: Vec<R> = r_star.iter().map(conjugate).collect();
            let w_mle = DenseMultilinearExtension::from_evaluations_vec(mu, w[j].clone());
            u_prime.push(w_mle.evaluate(&r_star).unwrap());
            u_dprime.push(w_mle.evaluate(&conj_r_star).unwrap());
            exact_sc.push(proof);
        }

        // ---- Pi^proj (section 4.3 Attempt II) ----
        // SOUNDNESS gate 1: Pi is drawn only now, after cm_w was absorbed.
        let pi = sample_pi(JL_ROWS, b, transcript);
        let v_blocks = w
            .iter()
            .map(|w_j| project(w_j, &pi, b))
            .collect::<Vec<Vec<R>>>();
        // v = [v_0 || ... || v_{L-1}], length L * (m/L) = m.
        let v: Vec<R> = v_blocks.iter().flatten().copied().collect();
        let cm_v = A.try_mul_vec(&v).unwrap();
        transcript.absorb_slice(&cm_v);

        // c = c0 (x) c1 ; c0 in R^L (instances), c1 in R^{m/L}.
        let c0 = (0..L)
            .map(|_| short_challenge(128, transcript))
            .collect::<Vec<R>>();
        let c1 = (0..m / L)
            .map(|_| short_challenge(128, transcript))
            .collect::<Vec<R>>();
        let c = tensor_product(&c0, &c1); // length m

        // t_j = <c1, v_j> ; s = <c, v>.
        let t = v_blocks
            .iter()
            .map(|v_j| {
                v_j.iter()
                    .zip(&c1)
                    .fold(R::zero(), |acc, (vi, ci)| acc + *vi * *ci)
            })
            .collect::<Vec<R>>();
        let s = v
            .iter()
            .zip(&c)
            .fold(R::zero(), |acc, (vi, ci)| acc + *vi * *ci);

        // Sumcheck for each t_j: sum_b MLE[a](b) MLE[w_j](b) = t_j.
        let a = proj_left(&pi, &c1, b, m);
        let mut projw_sc = Vec::with_capacity(L);
        let mut w_eval = Vec::with_capacity(L);
        for w_j in w.iter() {
            let a_mle = DenseMultilinearExtension::from_evaluations_vec(mu, a.clone());
            let w_mle = DenseMultilinearExtension::from_evaluations_vec(mu, w_j.clone());
            let comb = |vals: &[R]| -> R { vals[0] * vals[1] };
            let (proof, state) =
                MLSumcheck::prove_as_subprotocol(transcript, vec![a_mle, w_mle], mu, 2, comb);
            let ro: Vec<R> = state.randomness.iter().map(|x| R::from(*x)).collect();
            let w_mle = DenseMultilinearExtension::from_evaluations_vec(mu, w_j.clone());
            w_eval.push(w_mle.evaluate(&ro).unwrap());
            projw_sc.push(proof);
        }

        // Sumcheck for s: sum_b MLE[c](b) MLE[v](b) = s.
        let c_mle = DenseMultilinearExtension::from_evaluations_vec(mu, c.clone());
        let v_mle = DenseMultilinearExtension::from_evaluations_vec(mu, v.clone());
        let comb = |vals: &[R]| -> R { vals[0] * vals[1] };
        let (projv_sc, state) =
            MLSumcheck::prove_as_subprotocol(transcript, vec![c_mle, v_mle], mu, 2, comb);
        let ro: Vec<R> = state.randomness.iter().map(|x| R::from(*x)).collect();
        let v_mle = DenseMultilinearExtension::from_evaluations_vec(mu, v.clone());
        let v_eval = v_mle.evaluate(&ro).unwrap();

        L2Proof {
            cm_w,
            u,
            exact_sc,
            u_prime,
            u_dprime,
            cm_v,
            t,
            s,
            projw_sc,
            w_eval,
            projv_sc,
            v_eval,
        }
    }

    fn verify(
        &self,
        proof: &Self::Proof,
        params: &Self::Params,
        transcript: &mut impl Transcript<R>,
    ) -> Result<(), NormCheckError> {
        let L = proof.cm_w.len();
        let m = params.m;
        let mu = log2(m) as usize;
        let b = JL_ROWS * L;

        proof.cm_w.iter().for_each(|c| transcript.absorb_slice(c));

        // ---- Pi^exact ----
        transcript.absorb_slice(&proof.u);
        for j in 0..L {
            // Norm bound: ct(u_j) = ||w_j||_2^2 <= beta^2.
            // SOUNDNESS: bound fixed by the protocol; ct is taken over Z (modulus
            // must satisfy ||w||^2 < q, gate 6).
            let got = proof.u[j]
                .ct()
                .to_u64()
                .map_err(|_| NormCheckError::NormExceeded {
                    got: u128::MAX,
                    bound: params.beta_sq,
                })? as u128;
            if got > params.beta_sq {
                return Err(NormCheckError::NormExceeded {
                    got,
                    bound: params.beta_sq,
                });
            }
            let subclaim = MLSumcheck::verify_as_subprotocol(
                transcript,
                mu,
                2,
                proof.u[j],
                &proof.exact_sc[j],
            )
            .map_err(|_| NormCheckError::Sumcheck)?;
            // expected = MLE[w_j](r*) MLE[conj w_j](r*) = u'_j conj(u''_j).
            let expected = proof.u_prime[j] * conjugate(&proof.u_dprime[j]);
            if subclaim.expected_evaluation != expected {
                return Err(NormCheckError::ExactClaim);
            }
        }

        // ---- Pi^proj ----
        let pi = sample_pi(JL_ROWS, b, transcript);
        transcript.absorb_slice(&proof.cm_v);
        let c0 = (0..L)
            .map(|_| short_challenge(128, transcript))
            .collect::<Vec<R>>();
        let c1 = (0..m / L)
            .map(|_| short_challenge(128, transcript))
            .collect::<Vec<R>>();
        let c = tensor_product(&c0, &c1);

        // Folding consistency: sum_j c0_j t_j == s.
        let lhs = proof
            .t
            .iter()
            .zip(&c0)
            .fold(R::zero(), |acc, (t_j, c0_j)| acc + *t_j * *c0_j);
        if lhs != proof.s {
            return Err(NormCheckError::FoldingConsistency);
        }

        // t_j sumchecks: expected = MLE[a](ro) w_eval_j.
        let a = proj_left(&pi, &c1, b, m);
        let a_mle = DenseMultilinearExtension::from_evaluations_vec(mu, a);
        for j in 0..L {
            let subclaim = MLSumcheck::verify_as_subprotocol(
                transcript,
                mu,
                2,
                proof.t[j],
                &proof.projw_sc[j],
            )
            .map_err(|_| NormCheckError::Sumcheck)?;
            let ro = lift_point::<R>(&subclaim.point);
            let expected = a_mle.evaluate(&ro).unwrap() * proof.w_eval[j];
            if subclaim.expected_evaluation != expected {
                return Err(NormCheckError::ProjectionClaim);
            }
        }

        // s sumcheck: expected = MLE[c](ro) v_eval.
        let subclaim =
            MLSumcheck::verify_as_subprotocol(transcript, mu, 2, proof.s, &proof.projv_sc)
                .map_err(|_| NormCheckError::Sumcheck)?;
        let ro = lift_point::<R>(&subclaim.point);
        let c_mle = DenseMultilinearExtension::from_evaluations_vec(mu, c);
        let expected = c_mle.evaluate(&ro).unwrap() * proof.v_eval;
        if subclaim.expected_evaluation != expected {
            return Err(NormCheckError::ProjectionClaim);
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::PrimeField;
    use ark_std::Zero;
    use cyclotomic_rings::rings::FrogPoseidonConfig as PC;
    use rand::Rng;
    use stark_rings::{cyclotomic_ring::models::frog_ring::RqPoly as R, PolyRing};
    use stark_rings_linalg::Matrix;

    use super::*;
    use crate::{conjugation::l2_norm_sq, transcript::PoseidonTranscript};

    /// L witnesses in R^m with coefficients in {-1, 0, 1}.
    fn small_witness(rng: &mut impl Rng, l: usize, m: usize) -> Vec<Vec<R>> {
        (0..l)
            .map(|_| {
                (0..m)
                    .map(|_| {
                        let mut c = R::zero();
                        for k in 0..R::dimension() {
                            let v: i64 = rng.gen_range(-1..=1);
                            c.coeffs_mut()[k] = <R as PolyRing>::BaseRing::from(v.unsigned_abs())
                                * if v < 0 {
                                    -<R as PolyRing>::BaseRing::from(1u128)
                                } else {
                                    <R as PolyRing>::BaseRing::from(1u128)
                                };
                        }
                        c
                    })
                    .collect()
            })
            .collect()
    }

    fn norm_sq_u128(w_j: &[R]) -> u128 {
        l2_norm_sq(w_j).into_bigint().0[0] as u128
    }

    fn setup(l: usize, m: usize) -> (Vec<Vec<R>>, Matrix<R>, u128) {
        let mut rng = ark_std::test_rng();
        let w = small_witness(&mut rng, l, m);
        let kappa = 2;
        let A = Matrix::<R>::rand(&mut ark_std::test_rng(), kappa, m);
        let beta_sq = w.iter().map(|w_j| norm_sq_u128(w_j)).max().unwrap();
        (w, A, beta_sq)
    }

    #[test]
    fn test_l2_completeness() {
        for &l in &[2usize, 4] {
            let m = JL_ROWS * l; // smallest valid: one block per instance
            let (w, A, beta_sq) = setup(l, m);
            let params = L2Params { m, beta_sq };

            let mut ts = PoseidonTranscript::empty::<PC>();
            let proof = L2NormCheck.prove(&w, &A, &params, &mut ts);

            let mut ts = PoseidonTranscript::empty::<PC>();
            L2NormCheck.verify(&proof, &params, &mut ts).unwrap();
        }
    }

    #[test]
    fn test_l2_over_norm_rejected() {
        let l = 2;
        let m = JL_ROWS * l;
        let (w, A, beta_sq) = setup(l, m);
        // Tighten the bound just below the true norm: honest prover, failing check.
        let params = L2Params {
            m,
            beta_sq: beta_sq - 1,
        };

        let mut ts = PoseidonTranscript::empty::<PC>();
        let proof = L2NormCheck.prove(&w, &A, &params, &mut ts);

        let mut ts = PoseidonTranscript::empty::<PC>();
        let err = L2NormCheck.verify(&proof, &params, &mut ts).unwrap_err();
        assert!(matches!(err, NormCheckError::NormExceeded { .. }));
    }

    #[test]
    fn test_l2_inconsistent_s_rejected() {
        let l = 2;
        let m = JL_ROWS * l;
        let (w, A, beta_sq) = setup(l, m);
        let params = L2Params { m, beta_sq };

        let mut ts = PoseidonTranscript::empty::<PC>();
        let mut proof = L2NormCheck.prove(&w, &A, &params, &mut ts);
        // s no longer equals sum_j c0_j t_j.
        proof.s += R::from(1u128);

        let mut ts = PoseidonTranscript::empty::<PC>();
        let err = L2NormCheck.verify(&proof, &params, &mut ts).unwrap_err();
        assert!(matches!(err, NormCheckError::FoldingConsistency));
    }

    #[test]
    fn test_l2_tampered_cm_v_rejected() {
        // A cm_v inconsistent with the projection v shifts the c0/c1 challenges,
        // so the folding/projection checks fail (commitment binding of v).
        let l = 2;
        let m = JL_ROWS * l;
        let (w, A, beta_sq) = setup(l, m);
        let params = L2Params { m, beta_sq };

        let mut ts = PoseidonTranscript::empty::<PC>();
        let mut proof = L2NormCheck.prove(&w, &A, &params, &mut ts);
        proof.cm_v[0] += R::from(1u128);

        let mut ts = PoseidonTranscript::empty::<PC>();
        assert!(L2NormCheck.verify(&proof, &params, &mut ts).is_err());
    }

    #[test]
    fn test_l2_transcript_binding() {
        // Verifying against a transcript in a different Fiat-Shamir state (an
        // extra absorb up front) diverges all challenges and is rejected. This is
        // the structural defense behind soundness gate 1 (challenges, including
        // Pi, are bound to the absorbed commitments).
        let l = 2;
        let m = JL_ROWS * l;
        let (w, A, beta_sq) = setup(l, m);
        let params = L2Params { m, beta_sq };

        let mut ts = PoseidonTranscript::empty::<PC>();
        let proof = L2NormCheck.prove(&w, &A, &params, &mut ts);

        let mut ts = PoseidonTranscript::empty::<PC>();
        ts.absorb(&R::from(42u128));
        assert!(L2NormCheck.verify(&proof, &params, &mut ts).is_err());
    }

    #[test]
    fn test_l2params_modulus_guard() {
        // alpha^2 * beta^2 = 4 * 100 = 400 < q = 1000 -> ok.
        let ok = L2Params::from_security(512, 1000, 2, 10).unwrap();
        assert_eq!(ok.beta_sq, 100);
        // 4 * 100 = 400 >= q = 200 -> rejected (gate 6).
        assert!(matches!(
            L2Params::from_security(512, 200, 2, 10),
            Err(L2ParamError::ModulusTooSmall { .. })
        ));
    }

    #[test]
    fn test_calculate_bound_l2_is_public() {
        // Smoke test that the now-public bound is reachable for beta derivation.
        let q_log2 = (<<R as PolyRing>::BaseRing>::MODULUS.0[0] as f64).log2();
        let beta =
            latticefold::utils::security_check::calculate_bound_l2(R::dimension(), 2, q_log2);
        assert!(beta > num_bigint::BigUint::from(0u32));
    }

    #[test]
    fn test_linf_baseline_completeness() {
        // The ell-infinity path commits the gadget-decomposed witness, so it
        // needs a large n (here n = m): with kappa=1, k=2 the unpadded tau is
        // ~11K, hence m = 2^14.
        let l = 1;
        let m = 1 << 14;
        let mut rng = ark_std::test_rng();
        let w = small_witness(&mut rng, l, m);
        let kappa = 1;
        let A = Matrix::<R>::rand(&mut ark_std::test_rng(), kappa, m);

        let b = (R::dimension() / 2) as u128;
        let k = 2;
        let l_decomp = ((<<R as PolyRing>::BaseRing>::MODULUS.0[0] as f64).ln()
            / ((R::dimension() / 2) as f64).ln())
        .ceil() as usize;
        let params = LinfParams {
            kappa,
            decomp: DecompParameters { b, k, l: l_decomp },
        };

        let mut ts = PoseidonTranscript::empty::<PC>();
        let proof = LinfNormCheck.prove(&w, &A, &params, &mut ts);

        let mut ts = PoseidonTranscript::empty::<PC>();
        LinfNormCheck.verify(&proof, &params, &mut ts).unwrap();
    }
}
