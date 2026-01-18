//! Symphony-style high-arity *linear* folding helpers (Figure 4, Step 4–6).
//!
//! This module does **not** implement the full `Π_fold` protocol (which includes running many
//! `Π_gr1cs` instances with shared randomness and batched sumchecks). Instead, it provides the
//! deterministic linear-combination step used by `Π_fold` once the verifier supplies the
//! low-norm challenge vector `β`.
//!
//! We keep this small and explicit because downstream Witness Encryption / arming wants access
//! to the folded accumulator’s public `(r, v)` values.

use stark_rings::PolyRing;

use crate::rp_rgchk::PiRgVerifiedOutput;

/// A minimal public instance in the natural Symphony `(c, r, v)` shape.
///
/// - `c`: commitment (Ajtai commitment vector)
/// - `r`: public point (in K^{log m} or K^{log n}, depending on the subprotocol)
/// - `v`: public value (in K^d, i.e. coefficient-space form)
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SymphonyInstance<R: PolyRing> {
    pub c: Vec<R>,
    /// Shared verifier randomness point (e.g. `r̄||s̄` in Figure 4).
    pub r: Vec<R::BaseRing>,
    /// Value in the tensor ring `E`; in our concrete instantiations we represent it as an `R_q` element.
    pub v: R,
}

/// The batched-linear output (Eq. (28)): shared `r'` and per-digit `u(i)`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SymphonyBatchLin<R: PolyRing> {
    pub r_prime: Vec<R::BaseRing>,
    /// Per-digit monomial commitments `c(i)` (Figure 2 Step 3 / Eq. (28)).
    ///
    /// Shape: `[k_g][kappa]`.
    pub c_g: Vec<Vec<R>>,
    pub u: Vec<R>,
}

/// Fold `ℓ` instances linearly with respect to low-norm coefficients `β ∈ K^ℓ`.
///
/// This corresponds to Figure 4 (Step 5), where `(c*, x*_in, v*) := Σ β_ℓ (c_ℓ, x_{in,ℓ}, v_ℓ)`.
pub fn fold_instances<R>(
    beta: &[R],
    instances: &[SymphonyInstance<R>],
) -> SymphonyInstance<R>
where
    R: PolyRing,
{
    assert!(!instances.is_empty());
    assert_eq!(beta.len(), instances.len());

    let c_len = instances[0].c.len();
    let r_len = instances[0].r.len();
    for inst in instances {
        assert_eq!(inst.c.len(), c_len);
        assert_eq!(inst.r.len(), r_len);
        // r is shared randomness in Π_fold, so it must match across instances.
        assert_eq!(inst.r, instances[0].r, "r must be identical across instances");
    }

    let mut c = vec![R::ZERO; c_len];
    let r = instances[0].r.clone();
    let mut v = R::ZERO;

    for (b, inst) in beta.iter().zip(instances.iter()) {
        for (acc, x) in c.iter_mut().zip(inst.c.iter()) {
            *acc += *b * *x;
        }
        v += *b * inst.v;
    }

    SymphonyInstance { c, r, v }
}

/// Fold the batch-linear outputs (Eq. (28)) linearly with respect to low-norm `β`.
///
/// This corresponds to Figure 4 (Step 5), where for each digit `i`,
/// `(c(i), u(i)) := Σ β_ℓ (c(i)_ℓ, u(i)_ℓ)`.
pub fn fold_batchlin<R>(
    beta: &[R],
    batch: &[SymphonyBatchLin<R>],
) -> SymphonyBatchLin<R>
where
    R: PolyRing,
{
    assert!(!batch.is_empty());
    assert_eq!(beta.len(), batch.len());

    let r_len = batch[0].r_prime.len();
    let kg = batch[0].u.len();
    assert_eq!(batch[0].c_g.len(), kg, "c_g k_g mismatch");
    let kappa = batch[0].c_g[0].len();
    for b in batch {
        assert_eq!(b.r_prime.len(), r_len, "r' must be shared");
        assert_eq!(b.u.len(), kg, "k_g mismatch");
        assert_eq!(b.r_prime, batch[0].r_prime, "r' must be identical across instances");
        assert_eq!(b.c_g.len(), kg, "c_g k_g mismatch");
        for ci in &b.c_g {
            assert_eq!(ci.len(), kappa, "c_g kappa mismatch");
        }
    }

    let mut c_g = vec![vec![R::ZERO; kappa]; kg];
    let mut u = vec![R::ZERO; kg];
    for (b, inst) in beta.iter().zip(batch.iter()) {
        for dig in 0..kg {
            for j in 0..kappa {
                c_g[dig][j] += *b * inst.c_g[dig][j];
            }
        }
        for (acc, x) in u.iter_mut().zip(inst.u.iter()) {
            *acc += *b * *x;
        }
    }

    SymphonyBatchLin {
        r_prime: batch[0].r_prime.clone(),
        c_g,
        u,
    }
}

/// Convenience: turn a verified Π_rg output into the `(c,r,v)` + batchlin shape used by folding.
pub fn pi_rg_to_fold_shapes<R: PolyRing>(
    cm_f: Vec<R>,
    c_g: Vec<Vec<R>>,
    out: &PiRgVerifiedOutput<R>,
) -> (SymphonyInstance<R>, SymphonyBatchLin<R>) {
    // Convert the coefficient vector `v ∈ K^d` into an `R_q` element by placing it in the
    // coefficient representation. This matches the paper's view of `E` as an `R_q`-module.
    let mut v_rq = R::ZERO;
    for (i, c) in out.v.iter().enumerate() {
        v_rq.coeffs_mut()[i] = *c;
    }
    (
        SymphonyInstance {
            c: cm_f,
            r: out.r.clone(),
            v: v_rq,
        },
        SymphonyBatchLin {
            r_prime: out.r_prime.clone(),
            c_g,
            u: out.u.clone(),
        },
    )
}

