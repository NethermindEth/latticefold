//! Tiny-field “lockable DPP” (A = {0,1}) for WE-style arming.
//!
//! This module provides a **non-enumerating** variant suitable for “arm-before-proof” demos:
//! - The armer samples a **single sparse query** `q`, but treats it as **toxic waste** (never published).
//! - The armer publishes only **public coins** plus an **accepting set** `A={0,1}` and uses `q` only
//!   internally to form a lock artifact (in production: LWE hints / AEAD KDF inputs).
//! - The prover/decapper later produces a proof **for the public coins**; decap uses the lock artifact
//!   without ever revealing `q`.
//!
//! This avoids any “enumerate all verifier randomness strings” blowup. It is intentionally
//! **not** a standard “proof-first” DPP in the strict sense; it matches a WE gate setting where
//! the query is hidden before the proof exists.
//!
//! Internally it uses the tiny-field gadgets from TR24-114 rev2 (Section 4.3.2–4.3.3):
//! - UV → powering via Vandermonde (gives Sq accepting set {0,1})
//! - Multiplication relation checked via (ρ,σ) linearization and the Sq test
//! - Outer FLPCP for NP dR1CS provides (α,β,γ) = (E(Az)[i], E(Bz)[i], ...), with proof π₀=(z_w||w)

use ark_ff::{BigInteger, PrimeField};
use rand::RngCore;
use std::collections::BTreeMap;

use crate::dr1cs_flpcp::RsDr1csNpFlpcpSparse;
use crate::packing::BoundedFlpcpSparse;
use crate::sparse::SparseVec;
use crate::toy_poseidon2::{estimate_cost, ToyPoseidon2Cost, ToyPoseidon2Params, ToyPoseidon2Sponge};

/// Coins defining the single lockable query.
#[derive(Clone, Debug)]
pub struct Theorem43Coins<F: PrimeField> {
    pub idx: usize,
    pub lambda: F,
    pub rho: F,
    pub sigma: F,
}

/// Public arming artifact for the hidden-query lockable DPP.
///
/// IMPORTANT:
/// - The query vector `q` is *not* part of the public API (it is stored privately here only as a
///   stand-in for “LWE hints” in this toy implementation).
/// - In the real PVUGC flow, this object would contain noisy LWE hints + AEAD material, not `q`.
#[derive(Clone, Debug)]
pub struct Theorem43LockArtifact<F: PrimeField> {
    /// Public statement-binding commitment (domain-separated).
    ///
    /// This is witness/proof-independent and can be agreed by all armers.
    pub c_stmt: Vec<F>,
    /// Accepting set (always `{0,1}`).
    pub accepting_set: [F; 2],
    /// Total length of `(x || π)` for sanity checks.
    pub len: usize,
    /// Public coins used by the prover to generate the matching proof.
    pub coins: Theorem43Coins<F>,
    /// Lightweight “full gate” stats (FS cost proxy + query sparsity).
    pub stats: Theorem43ArmingStats,
    /// Hidden query vector over `(x || π)` (TOXIC WASTE; do not publish).
    q: SparseVec<F>,
}

#[derive(Clone, Debug, Default)]
pub struct Theorem43ArmingStats {
    /// Number of toy permutation calls used for FS derivation.
    pub fs_permutes: u64,
    /// Estimated cost per permutation call (field adds/muls).
    pub fs_cost_permute: ToyPoseidon2Cost,
    /// Estimated total cost (fs_permutes × fs_cost_permute).
    pub fs_cost_total: ToyPoseidon2Cost,
    /// Number of nonzero terms in the hidden sparse query.
    pub q_nnz: usize,
}

impl<F: PrimeField> Theorem43LockArtifact<F> {
    /// Compute the DPP “answer” \(a = \langle q, (x || \pi)\rangle\) using the hidden query.
    ///
    /// In production, this would be recovered implicitly via the LWE/AEAD decapsulation path,
    /// without ever exposing `q`.
    pub fn answer_for(&self, x: &[F], pi: &[F]) -> Result<F, String> {
        if x.len() + pi.len() != self.len {
            return Err("bad (x||pi) length".to_string());
        }
        let mut v = Vec::with_capacity(self.len);
        v.extend_from_slice(x);
        v.extend_from_slice(pi);
        Ok(self.q.dot(&v))
    }

    /// Split the hidden query `q` into `(q_x, q_pi)` at the boundary `x_len`.
    ///
    /// - `q_x` is indexed over `x[0..x_len)`.
    /// - `q_pi` is indexed over `pi[0..pi_len)` (indices are **shifted down** by `x_len`).
    ///
    /// This matches PVUGC’s “shift away the public part” interface, where locks are built only
    /// for the proof vector `pi` using `q_pi`.
    pub fn split_query(&self, x_len: usize, pi_len: usize) -> Result<(SparseVec<F>, SparseVec<F>), String> {
        if x_len + pi_len != self.len {
            return Err("split_query: bad (x_len + pi_len)".to_string());
        }
        let mut qx = SparseVecBuilder::<F>::new();
        let mut qpi = SparseVecBuilder::<F>::new();
        for (c, idx) in self.q.terms.iter().copied() {
            if idx < x_len {
                qx.add_term(idx, c);
            } else {
                let j = idx - x_len;
                if j >= pi_len {
                    return Err("split_query: q index out of range".to_string());
                }
                qpi.add_term(j, c);
            }
        }
        Ok((qx.build(), qpi.build()))
    }
}

/// Non-enumerating Theorem-4.3-style lockable query generator for the NP dR1CS FLPCP.
#[derive(Clone, Debug)]
pub struct Theorem43Dpp<F: PrimeField> {
    flpcp: RsDr1csNpFlpcpSparse<F>,
    p: u64,
    q_minus_3: usize, // p-3
    /// Toy permutation params used for FS coin derivation (cost proxy).
    pub fs_params: ToyPoseidon2Params,
}

impl<F: PrimeField> Theorem43Dpp<F> {
    pub fn new(flpcp: RsDr1csNpFlpcpSparse<F>) -> Result<Self, String> {
        let p = field_modulus_u64::<F>().ok_or("field modulus does not fit u64")?;
        if p < 3 || p % 2 == 0 {
            return Err("field must have odd characteristic >= 3".to_string());
        }
        if p > 257 {
            return Err("field too large for tiny-field theorem43 demo".to_string());
        }
        Ok(Self {
            flpcp,
            p,
            q_minus_3: (p as usize).saturating_sub(3),
            // Default “cost proxy” parameters. These are NOT standardized Poseidon2 params.
            fs_params: ToyPoseidon2Params {
                width: 12,
                full_rounds: 8,
                partial_rounds: 22,
            },
        })
    }

    /// Length of the produced proof `π` (field elements).
    ///
    /// Layout:
    /// - π0 = FLPCP NP proof (z_w || w) of length `flpcp.m()`
    /// - μ, ν (2 elems)
    /// - u^3..u^{p-1} (p-3 elems)
    pub fn proof_len(&self) -> usize {
        self.flpcp.m() + 2 + self.q_minus_3
    }

    /// Produce the proof `π` matching a particular sampled query.
    pub fn prove_for_query(&self, x: &[F], z_w: &[F], coins: &Theorem43Coins<F>) -> Result<Vec<F>, String> {
        if x.len() != self.flpcp.n() {
            return Err("bad public input length".to_string());
        }
        if coins.idx >= self.flpcp.ell {
            return Err("bad idx coin".to_string());
        }

        // Base FLPCP proof π0 = (z_w || w).
        let pi0 = self.flpcp.prove(x, z_w);

        // Evaluate outer query answers (α,β,γ) on v0=(x||π0).
        let mut v0 = Vec::with_capacity(self.flpcp.n() + self.flpcp.m());
        v0.extend_from_slice(x);
        v0.extend_from_slice(&pi0);

        let (qs, _pred) = self
            .flpcp
            .queries_for_coins_sparse(coins.idx, coins.lambda, x)
            .map_err(|e| format!("outer coins->queries failed: {e}"))?;
        debug_assert_eq!(qs.len(), 3);

        let alpha = qs[0].dot(&v0);
        let beta = qs[1].dot(&v0);
        let gamma = qs[2].dot(&v0);

        // Multiplication gadget proof components.
        let mu = alpha * alpha;
        let nu = beta * beta;
        let u = coins.rho * alpha + coins.sigma * beta;

        // Output π = π0 || μ || ν || (u^3..u^{p-1}).
        let mut pi = Vec::with_capacity(self.proof_len());
        pi.extend_from_slice(&pi0);
        pi.push(mu);
        pi.push(nu);

        // u^3..u^{p-1}. (Empty if p=3.)
        if self.q_minus_3 > 0 {
            let mut cur = (u * u) * u; // u^3
            for _ in 0..self.q_minus_3 {
                pi.push(cur);
                cur *= u;
            }
        }
        debug_assert_eq!(pi.len(), self.proof_len());

        // NOTE: `gamma` is not placed in π; it is already determined by π0 and is used via q3.
        let _ = gamma;
        Ok(pi)
    }

    /// Arm a hidden-query lock artifact using fresh randomness for the UV secret.
    ///
    /// This returns a public artifact (coins + `C_stmt` + stats) and keeps `q` hidden inside.
    pub fn arm(&self, rng: &mut dyn RngCore, c_stmt: Vec<F>, x: &[F]) -> Result<Theorem43LockArtifact<F>, String> {
        if x.len() != self.flpcp.n() {
            return Err("bad public input length".to_string());
        }

        // Public coins.
        let idx = (rng.next_u64() as usize) % self.flpcp.ell;
        let lambda = self.random_field_elem(rng);
        let rho = self.random_field_elem(rng);
        let sigma = self.random_field_elem(rng);

        // UV randomness (secret) -> Sq coefficients for powers u^1..u^{p-1}.
        let coeffs = self.sample_sq_coeffs(rng)?;
        let (q, len) = self.build_query_from_sq_coeffs(
            x,
            &Theorem43Coins { idx, lambda, rho, sigma },
            &coeffs,
        )?;

        Ok(Theorem43LockArtifact {
            c_stmt,
            accepting_set: [F::ZERO, F::ONE],
            len,
            coins: Theorem43Coins { idx, lambda, rho, sigma },
            stats: Theorem43ArmingStats {
                fs_permutes: 0,
                fs_cost_permute: ToyPoseidon2Cost::default(),
                fs_cost_total: ToyPoseidon2Cost::default(),
                q_nnz: q.terms.len(),
            },
            q,
        })
    }

    /// FS variant for the **full-gate cost shape**, while keeping `q` hidden.
    ///
    /// - Public coins are derived from `(domain_sep, C_stmt, x)`.
    /// - The hidden query randomness (UV bits) is derived from the same transcript **plus**
    ///   an armer-private secret salt `armer_secret`.
    ///
    /// This is the version you’d arithmetize inside the “full gate” to de-risk constraint cost.
    /// The permutation is a **toy cost proxy**, not a standardized hash.
    pub fn arm_fs(&self, c_stmt: &[F], x: &[F], armer_secret: &[F]) -> Result<Theorem43LockArtifact<F>, String> {
        if x.len() != self.flpcp.n() {
            return Err("bad public input length".to_string());
        }

        // Sponge: absorb domain-sep + statement-binding commitment + statement x.
        let mut sp = ToyPoseidon2Sponge::<F>::new(self.fs_params.clone());
        sp.absorb(&[F::from(43u64)]); // domain sep for theorem43
        sp.absorb(c_stmt);
        sp.absorb(x);

        // Squeeze PUBLIC coins.
        let idx_f = sp.squeeze(1)[0];
        let idx_u = field_elem_to_u64(idx_f).ok_or("fs idx did not fit u64")? as usize;
        let idx = idx_u % self.flpcp.ell;

        let lambda = sp.squeeze(1)[0];
        let rho = sp.squeeze(1)[0];
        let sigma = sp.squeeze(1)[0];

        // Mix in armer-private secret so the query `q` is not reconstructible from public data.
        sp.absorb(armer_secret);

        // Use remaining squeezes as UV bits (0/1), by taking LSB of the field encoding.
        // In-circuit you would constrain these to be boolean; for now we just derive them.
        let q_minus_1 = (self.p as usize) - 1;
        let mut q_bits: Vec<u8> = Vec::with_capacity(q_minus_1);
        let bits_src = sp.squeeze(q_minus_1);
        for (i, b) in bits_src.into_iter().enumerate() {
            let mut u = field_elem_to_u64(b).ok_or("fs bit elem did not fit u64")?;
            // Tiny-field toy: de-bias the single-bit extraction with the index.
            // (In-circuit we still constrain bits to be boolean; this just avoids degenerate all-zero patterns.)
            u ^= i as u64;
            q_bits.push((u & 1) as u8);
        }

        let coins = Theorem43Coins { idx, lambda, rho, sigma };
        let coeffs = self.sq_coeffs_from_uv_bits(&q_bits)?;
        let (q, len) = self.build_query_from_sq_coeffs(x, &coins, &coeffs)?;

        let fs_permutes = sp.permute_count();
        let per = estimate_cost(&self.fs_params);
        let tot = ToyPoseidon2Cost {
            muls: per.muls.saturating_mul(fs_permutes),
            adds: per.adds.saturating_mul(fs_permutes),
        };

        Ok(Theorem43LockArtifact {
            c_stmt: c_stmt.to_vec(),
            accepting_set: [F::ZERO, F::ONE],
            len,
            coins,
            stats: Theorem43ArmingStats {
                fs_permutes,
                fs_cost_permute: per,
                fs_cost_total: tot,
                q_nnz: q.terms.len(),
            },
            q,
        })
    }

    fn build_query_from_sq_coeffs(
        &self,
        x: &[F],
        coins: &Theorem43Coins<F>,
        coeffs: &[F],
    ) -> Result<(SparseVec<F>, usize), String> {
        // Outer FLPCP queries (q1,q2,q3) on v0=(x||π0).
        let (qs, _pred) = self
            .flpcp
            .queries_for_coins_sparse(coins.idx, coins.lambda, x)
            .map_err(|e| format!("outer coins->queries failed: {e}"))?;
        debug_assert_eq!(qs.len(), 3);

        if coeffs.len() != (self.p as usize) - 1 {
            return Err("bad Sq coeff length".to_string());
        }
        let c1 = coeffs[0];
        let c2 = coeffs[1];

        let coeff_alpha = c1 * coins.rho;
        let coeff_beta = c1 * coins.sigma;
        let coeff_gamma = {
            let two = F::from(2u64);
            c2 * (two * coins.rho * coins.sigma)
        };

        let mut builder = SparseVecBuilder::<F>::new();
        builder.add_scaled_sparse(&qs[0], coeff_alpha);
        builder.add_scaled_sparse(&qs[1], coeff_beta);
        builder.add_scaled_sparse(&qs[2], coeff_gamma);

        let m0 = self.flpcp.m();
        let base_pi = x.len();
        let base_tail = base_pi + m0;

        let coeff_mu = c2 * (coins.rho * coins.rho);
        let coeff_nu = c2 * (coins.sigma * coins.sigma);
        builder.add_term(base_tail + 0, coeff_mu);
        builder.add_term(base_tail + 1, coeff_nu);
        for (t, c) in coeffs.iter().copied().skip(2).enumerate() {
            if c.is_zero() {
                continue;
            }
            builder.add_term(base_tail + 2 + t, c);
        }

        let len = x.len() + self.proof_len();
        Ok((builder.build(), len))
    }

    #[inline]
    pub fn accept_answer(&self, a: &F) -> bool {
        *a == F::ZERO || *a == F::ONE
    }

    fn random_field_elem(&self, rng: &mut dyn RngCore) -> F {
        let r = rng.next_u64() % self.p;
        F::from(r)
    }

    /// Sample UV randomness and compute Sq coefficients `c_i = -Σ_{λ∈F*} q_λ * λ^i`.
    fn sample_sq_coeffs(&self, rng: &mut dyn RngCore) -> Result<Vec<F>, String> {
        let q_minus_1 = (self.p as usize) - 1;
        let mut q_bits = Vec::with_capacity(q_minus_1);
        for _ in 0..q_minus_1 {
            q_bits.push((rng.next_u64() & 1) as u8);
        }
        self.sq_coeffs_from_uv_bits(&q_bits)
    }

    fn sq_coeffs_from_uv_bits(&self, q_bits: &[u8]) -> Result<Vec<F>, String> {
        let q_minus_1 = (self.p as usize) - 1;
        if q_bits.len() != q_minus_1 {
            return Err("bad uv bit length".to_string());
        }

        let mut coeffs = vec![F::ZERO; q_minus_1];
        for i in 1..=q_minus_1 {
            let mut acc = F::ZERO;
            for lam_u in 1..self.p {
                let bit = q_bits[(lam_u as usize) - 1];
                if bit == 0 {
                    continue;
                }
                let lam = F::from(lam_u);
                acc += lam.pow([i as u64]);
            }
            coeffs[i - 1] = -acc;
        }

        if coeffs.len() < 2 {
            return Err("field too small for Sq coefficients".to_string());
        }
        Ok(coeffs)
    }
}

/// Accumulate sparse terms with coefficient merging (by index).
#[derive(Clone, Debug)]
struct SparseVecBuilder<F: PrimeField> {
    acc: BTreeMap<usize, F>,
}

impl<F: PrimeField> SparseVecBuilder<F> {
    fn new() -> Self {
        Self {
            acc: BTreeMap::new(),
        }
    }

    fn add_term(&mut self, idx: usize, coeff: F) {
        if coeff.is_zero() {
            return;
        }
        self.acc
            .entry(idx)
            .and_modify(|c| *c += coeff)
            .or_insert(coeff);
    }

    fn add_scaled_sparse(&mut self, v: &SparseVec<F>, scale: F) {
    if scale.is_zero() {
        return;
    }
        for (c, idx) in v.terms.iter().copied() {
            self.add_term(idx, c * scale);
        }
    }

    fn build(mut self) -> SparseVec<F> {
        self.acc.retain(|_, c| !c.is_zero());
        SparseVec::new(self.acc.into_iter().map(|(i, c)| (c, i)).collect())
    }
}

fn field_modulus_u64<F: PrimeField>() -> Option<u64> {
    let bytes = F::MODULUS.to_bytes_le();
    if bytes.len() > 8 {
        return None;
    }
    let mut buf = [0u8; 8];
    buf[..bytes.len()].copy_from_slice(&bytes);
    Some(u64::from_le_bytes(buf))
}

fn field_elem_to_u64<F: PrimeField>(a: F) -> Option<u64> {
    let bytes = a.into_bigint().to_bytes_le();
    if bytes.len() > 8 {
        return None;
    }
    let mut buf = [0u8; 8];
    buf[..bytes.len()].copy_from_slice(&bytes);
    Some(u64::from_le_bytes(buf))
}

