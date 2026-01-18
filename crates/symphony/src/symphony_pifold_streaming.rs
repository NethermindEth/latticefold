//! Streaming version of Symphony Π_fold prover.
//!
//! This module provides a memory-efficient version of the batched Π_fold prover
//! that uses streaming MLEs instead of materializing O(2^n) evaluation vectors.
//!
//! ## Key Optimizations
//!
//! 1. **Hadamard MLEs**: sparse mat-vec is evaluated on-demand, then fixed rounds materialize dense
//! 2. **Eq MLEs**: evaluated in the base ring (constant-coefficient) and lifted to `R`
//! 3. **Monomial MLEs**: base-field scalars (m_j/m') + base-field eq(c), both streaming/compact
//!
//! ## Memory Comparison (1M constraints, d=16)
//!
//! | Component | Dense | Streaming |
//! |-----------|-------|-----------|
//! | y = M*w MLEs | 3 × 16MB | O(1) per eval |
//! | eq MLE | 16MB | O(1) per eval |
//! | m_j/m' MLEs | 2 × k_g × 256MB | Same (complex) |

use ark_std::log2;
use ark_ff::PrimeField;
use std::collections::HashMap;
use std::sync::Arc;
use std::time::Instant;
use rayon::prelude::*;
use stark_rings::{
    balanced_decomposition::{Decompose, DecomposeToVec},
    exp, CoeffRing, OverField, Ring, Zq,
};
use stark_rings_linalg::SparseMatrix;
use cyclotomic_rings::rings::GetPoseidonParams;
use latticefold::{
    commitment::AjtaiCommitmentScheme,
    transcript::Transcript,
};

use crate::recording_transcript::RecordingTranscriptRef;
use crate::rp_rgchk::RPParams;
use crate::streaming_sumcheck::{StreamingMleEnum, StreamingProof, StreamingSumcheck};
use crate::symphony_cm::SymphonyCoins;
use crate::symphony_coins::{derive_beta_chi, derive_J, ev, ts_weights};
use crate::symphony_pifold_batched::{PiFoldAuxWitness, PiFoldBatchedProof, PiFoldProverOutput};
use crate::pcs::batchlin_pcs::BATCHLIN_PCS_DOMAIN_SEP;
use crate::pcs::cmf_pcs::{cmf_pcs_coin_bytes_len, cmf_pcs_params_for_flat_len, CMF_PCS_DOMAIN_SEP};
use crate::pcs::batchlin_pcs::{batchlin_scalar_pcs_params, BatchlinPcsConfig, mle_point_to_tensor};
use crate::pcs::folding_pcs_l2::{commit as pcs_commit, open as pcs_open, BinMatrix};

/// Fixed seed for the cm_g aggregate commitment scheme.
/// This derives a dedicated Ajtai matrix for binding all cm_g commitments into a single short value.
pub const CM_G_AGG_SEED: [u8; 32] = *b"SYMPHONY_CM_G_AGG_V1_0000000000_";

/// Fixed seed for the mon_b aggregate commitment scheme.
/// This derives a dedicated Ajtai matrix for binding all mon_b values into a single short value.
pub const MON_B_AGG_SEED: [u8; 32] = *b"SYMPHONY_MON_B_AGG_V1_000000000_";

/// Fixed seed for the v_digits_folded aggregate commitment scheme.
/// This derives a dedicated Ajtai matrix for binding the hook message `v_digits_folded` into a short value.
pub const V_DIGITS_AGG_SEED: [u8; 32] = *b"SYMPHONY_VDIGITS_AGG_V1_00000000";

/// Compute the aggregate commitment for all cm_g values.
///
/// This flattens `cm_g[inst][dig][*]` into a single vector and commits via a dedicated
/// Ajtai scheme. The result is a short commitment that binds all cm_g values, reducing
/// Poseidon absorption from O(ell * k_g * kappa) field elements to O(kappa).
pub fn compute_cm_g_aggregate<R: CoeffRing>(
    cm_g: &[Vec<Vec<R>>],
    kappa: usize,
) -> Result<Vec<R>, String>
where
    R::BaseRing: Zq + stark_rings::balanced_decomposition::Decompose,
{
    if cm_g.is_empty() {
        return Err("cm_g_aggregate: empty cm_g".to_string());
    }
    if kappa == 0 {
        return Err("cm_g_aggregate: kappa=0 is invalid".to_string());
    }
    let k_g = cm_g[0].len();
    if k_g == 0 {
        return Err("cm_g_aggregate: k_g=0 (no digits) is invalid".to_string());
    }
    // Flatten: concat_{inst, dig, j} cm_g[inst][dig][j]
    let mut cm_g_flat: Vec<R> = Vec::new();
    for (inst_idx, inst) in cm_g.iter().enumerate() {
        if inst.is_empty() {
            return Err(format!("cm_g_aggregate: empty cm_g[{inst_idx}]"));
        }
        if inst.len() != k_g {
            return Err(format!(
                "cm_g_aggregate: cm_g[{inst_idx}] digit count mismatch (got {}, expected {k_g})",
                inst.len()
            ));
        }
        for (dig_idx, dig) in inst.iter().enumerate() {
            if dig.len() != kappa {
                return Err(format!(
                    "cm_g_aggregate: cm_g[{inst_idx}][{dig_idx}] len mismatch (got {}, expected {kappa})",
                    dig.len()
                ));
            }
            cm_g_flat.extend(dig.iter().copied());
        }
    }
    let n = cm_g_flat.len();
    if n == 0 {
        return Err("cm_g_aggregate: flattened cm_g is empty".to_string());
    }
    
    // Create the aggregate scheme with the dedicated domain and seed.
    let agg_scheme = AjtaiCommitmentScheme::<R>::seeded(b"cm_g_agg", CM_G_AGG_SEED, kappa, n);
    let c_agg = agg_scheme
        .commit(&cm_g_flat)
        .map_err(|e| format!("cm_g_aggregate: commit failed: {e:?}"))?;
    Ok(c_agg.as_ref().to_vec())
}

/// Compute the aggregate commitment for all mon_b values.
///
/// This flattens `mon_b[inst][*]` into a single vector and commits via a dedicated
/// Ajtai scheme. The result is a short commitment that binds all mon_b values, reducing
/// Poseidon absorption from O(ell * k_g) ring elements to O(kappa).
pub fn compute_mon_b_aggregate<R: CoeffRing>(
    mon_b: &[Vec<R>],
    kappa: usize,
) -> Result<Vec<R>, String>
where
    R::BaseRing: Zq + stark_rings::balanced_decomposition::Decompose,
{
    if mon_b.is_empty() {
        return Err("mon_b_aggregate: empty mon_b".to_string());
    }
    if kappa == 0 {
        return Err("mon_b_aggregate: kappa=0 is invalid".to_string());
    }
    let k_g = mon_b[0].len();
    if k_g == 0 {
        return Err("mon_b_aggregate: k_g=0 (no digits) is invalid".to_string());
    }
    // Flatten: concat_{inst, dig} mon_b[inst][dig]
    let mut mon_b_flat: Vec<R> = Vec::new();
    for (inst_idx, inst) in mon_b.iter().enumerate() {
        if inst.len() != k_g {
            return Err(format!(
                "mon_b_aggregate: mon_b[{inst_idx}] digit count mismatch (got {}, expected {k_g})",
                inst.len()
            ));
        }
        mon_b_flat.extend(inst.iter().copied());
    }
    let n = mon_b_flat.len();
    if n == 0 {
        return Err("mon_b_aggregate: flattened mon_b is empty".to_string());
    }
    
    // Create the aggregate scheme with the dedicated domain and seed.
    let agg_scheme = AjtaiCommitmentScheme::<R>::seeded(b"mon_b_agg", MON_B_AGG_SEED, kappa, n);
    let c_agg = agg_scheme
        .commit(&mon_b_flat)
        .map_err(|e| format!("mon_b_aggregate: commit failed: {e:?}"))?;
    Ok(c_agg.as_ref().to_vec())
}

/// Compute an aggregate commitment for `v_digits_folded` (shape k_g × d, base-ring scalars).
///
/// The scalars are embedded into the ring as constant-coefficient elements and committed via a
/// dedicated Ajtai scheme. This is used as a **transcript compression**: absorb only `kappa` ring
/// elements instead of `k_g*d` scalars.
pub fn compute_v_digits_folded_aggregate<R: CoeffRing>(
    v_digits_folded: &[Vec<R::BaseRing>],
    kappa: usize,
) -> Result<Vec<R>, String>
where
    R::BaseRing: Zq + stark_rings::balanced_decomposition::Decompose,
    R: From<R::BaseRing>,
{
    if v_digits_folded.is_empty() {
        return Err("v_digits_aggregate: empty v_digits_folded".to_string());
    }
    if kappa == 0 {
        return Err("v_digits_aggregate: kappa=0 is invalid".to_string());
    }
    let d = v_digits_folded[0].len();
    if d == 0 {
        return Err("v_digits_aggregate: d=0 is invalid".to_string());
    }
    for (dig, v) in v_digits_folded.iter().enumerate() {
        if v.len() != d {
            return Err(format!(
                "v_digits_aggregate: v_digits_folded[{dig}] len mismatch (got {}, expected {d})",
                v.len()
            ));
        }
    }

    // Flatten as constant-coefficient ring elements.
    let mut flat: Vec<R> = Vec::with_capacity(v_digits_folded.len() * d);
    for v in v_digits_folded {
        for &x in v {
            flat.push(R::from(x));
        }
    }
    let n = flat.len();
    if n == 0 {
        return Err("v_digits_aggregate: empty flattened vector".to_string());
    }
    let agg_scheme = AjtaiCommitmentScheme::<R>::seeded(b"v_digits_agg", V_DIGITS_AGG_SEED, kappa, n);
    let c_agg = agg_scheme
        .commit_const_coeff_fast(&flat)
        .map_err(|e| format!("v_digits_aggregate: commit failed: {e:?}"))?;
    Ok(c_agg.as_ref().to_vec())
}

// Re-export streaming types
pub use crate::streaming_sumcheck::{StreamingProverMsg, StreamingSumcheckState};

/// Configuration for the streaming Π_fold prover.
///
/// Note: the symphony library does **not** read environment variables. If you want
/// runtime configuration, parse env vars in your binary/example and pass them through this struct.
#[derive(Clone, Debug)]
pub struct PiFoldStreamingConfig {
    /// Request the constant-coefficient fast path when inputs appear to be constant-coeff.
    ///
    /// Even when `true`, we still run a deterministic sampling check and fall back to the generic
    /// ring path if the check fails.
    pub request_const_coeff_fastpath: bool,
    /// Witness sample count used by the deterministic constant-coeff check.
    pub const_coeff_witness_samples: usize,
    /// Matrix nonzero-entry sample budget used by the deterministic constant-coeff check.
    pub const_coeff_matrix_entries: usize,
    /// Print internal phase timings (intended for benchmarks/profiling).
    pub profile: bool,
}

impl Default for PiFoldStreamingConfig {
    fn default() -> Self {
        Self {
            request_const_coeff_fastpath: true,
            const_coeff_witness_samples: 1024,
            const_coeff_matrix_entries: 4096,
            profile: false,
        }
    }
}

fn is_const_coeff_ring_elem<R: stark_rings::PolyRing>(x: &R) -> bool {
    x.coeffs()
        .iter()
        .skip(1)
        .all(|c| *c == <R as stark_rings::PolyRing>::BaseRing::ZERO)
}

fn seems_const_coeff_inputs<R: stark_rings::PolyRing>(
    witnesses: &[Arc<Vec<R>>],
    mats_to_check: &[Arc<stark_rings_linalg::SparseMatrix<R>>],
    w_samples: usize,
    m_entries: usize,
) -> bool {
    // Cheap, deterministic sampling to avoid unsound use of the fast path.
    // Defaults are tiny vs proving time; callers can adjust via `PiFoldStreamingConfig`.

    // Witness sample
    if w_samples > 0 {
        for w in witnesses.iter().take(1) {
            let step = (w.len() / w_samples.max(1)).max(1);
            let mut checked = 0usize;
            for i in (0..w.len()).step_by(step) {
                if !is_const_coeff_ring_elem(&w[i]) {
                    return false;
                }
                checked += 1;
                if checked >= w_samples {
                    break;
                }
            }
        }
    }

    // Matrix coefficient sample (scan nonzeros)
    if m_entries > 0 {
        let mut checked = 0usize;
        for m in mats_to_check.iter().take(3) {
            for row in &m.coeffs {
                for (coeff, _col) in row {
                    if !is_const_coeff_ring_elem(coeff) {
                        return false;
                    }
                    checked += 1;
                    if checked >= m_entries {
                        return true;
                    }
                }
            }
        }
    }

    true
}

fn absorb_public_inputs<R: OverField>(ts: &mut impl Transcript<R>, public_inputs: &[R::BaseRing])
where
    R::BaseRing: Zq,
{
    // Keep identical statement-binding schedule to symphony_pifold_batched.
    ts.absorb_field_element(&R::BaseRing::from(0x4c465053_5055424cu128)); // "LFPS_PUBL"
    ts.absorb_field_element(&R::BaseRing::from(public_inputs.len() as u128));
    for x in public_inputs {
        ts.absorb_field_element(x);
    }
}

fn encode_had_u_instance<R: OverField>(u: &[Vec<R::BaseRing>; 3]) -> Vec<R>
where
    R::BaseRing: Zq,
{
    let d = R::dimension();
    debug_assert_eq!(u[0].len(), d);
    debug_assert_eq!(u[1].len(), d);
    debug_assert_eq!(u[2].len(), d);
    let mut out = Vec::with_capacity(3 * d);
    out.extend(u[0].iter().copied().map(R::from));
    out.extend(u[1].iter().copied().map(R::from));
    out.extend(u[2].iter().copied().map(R::from));
    out
}

/// Canonical Poseidon-FS prover for Π_fold (streaming sumcheck).
///
/// This is the canonical prover entrypoint. It batches
/// a list of instances where each instance can have its own matrices `(M1,M2,M3)` (hetero-M).
///
/// (If you have a shared matrix across all instances, pass the same `(M1,M2,M3)` for each
/// instance by cloning the `Arc`s.)
pub fn prove_pi_fold_poseidon_fs<R: CoeffRing, PC>(
    Ms: &[[Arc<SparseMatrix<R>>; 3]],
    cms: &[Vec<R>],
    witnesses: &[Arc<Vec<R>>],
    public_inputs: &[R::BaseRing],
    cfs_had_u_scheme: Option<&AjtaiCommitmentScheme<R>>,
    cfs_mon_b_scheme: Option<&AjtaiCommitmentScheme<R>>,
    cm_g_scheme: &AjtaiCommitmentScheme<R>,
    rg_params: RPParams,
    config: &PiFoldStreamingConfig,
) -> Result<PiFoldProverOutput<R>, String>
where
    R::BaseRing: Zq + Decompose + PrimeField,
    PC: GetPoseidonParams<<<R>::BaseRing as ark_ff::Field>::BasePrimeField>,
{
    if cms.is_empty() {
        return Err("PiFold: empty batch".to_string());
    }
    if Ms.len() != cms.len() || cms.len() != witnesses.len() {
        return Err("PiFold: Ms/cms/witnesses length mismatch".to_string());
    }

    let mut ts = crate::transcript::PoseidonTranscript::<R>::empty::<PC>();
    let mut rts = RecordingTranscriptRef::<R, _>::new(&mut ts);
    rts.absorb_field_element(&R::BaseRing::from(0x4c465053_50494250u128)); // "LFPS_PIBP"
    absorb_public_inputs::<R>(&mut rts, public_inputs);

    let mut out = prove_pi_fold_streaming_hetero_m_with_config(
        &mut rts,
        Ms,
        cms,
        witnesses,
        rg_params,
        config,
        cm_g_scheme,
    )
        .map_err(|e| format!("PiFold: poseidon-fs prove failed: {e}"))?;

    out.proof.coins = SymphonyCoins::<R> {
        challenges: rts.coins_challenges,
        bytes: rts.coins_bytes,
        events: rts.events,
    };

    // Optionally commit to CP transcript witness messages (WE/DPP-facing path).
    out.cfs_had_u.clear();
    out.cfs_mon_b.clear();
    match (cfs_had_u_scheme, cfs_mon_b_scheme) {
        (Some(had_s), Some(mon_s)) => {
            out.cfs_had_u = out
                .aux
                .had_u
                .iter()
                .map(|u| {
                    had_s
                        .commit_const_coeff_fast(&encode_had_u_instance::<R>(u))
                        .map_err(|e| format!("PiFold: cfs_had_u commit failed: {e:?}"))
                        .map(|c| c.as_ref().to_vec())
                })
                .collect::<Result<Vec<_>, _>>()?;
            out.cfs_mon_b = out
                .aux
                .mon_b
                .iter()
                .map(|b| {
                    mon_s
                        .commit_const_coeff_fast(b)
                        .map_err(|e| format!("PiFold: cfs_mon_b commit failed: {e:?}"))
                        .map(|c| c.as_ref().to_vec())
                })
                .collect::<Result<Vec<_>, _>>()?;
        }
        (None, None) => {}
        _ => {
            return Err(
                "PiFold: must provide both cfs_had_u_scheme and cfs_mon_b_scheme or neither".to_string(),
            );
        }
    }

    Ok(out)
}


pub fn prove_pi_fold_streaming_hetero_m_with_config<R: CoeffRing>(
    transcript: &mut impl Transcript<R>,
    Ms: &[[Arc<SparseMatrix<R>>; 3]],
    cms: &[Vec<R>],
    witnesses: &[Arc<Vec<R>>],
    rg_params: RPParams,
    config: &PiFoldStreamingConfig,
    cm_g_scheme: &AjtaiCommitmentScheme<R>,
) -> Result<PiFoldProverOutput<R>, String>
where
    R::BaseRing: Zq + Decompose + PrimeField,
{
    // Cheap sample: first witness and first instance matrices.
    let const_coeff_fastpath = config.request_const_coeff_fastpath
        && seems_const_coeff_inputs::<R>(
            witnesses,
            &[Ms[0][0].clone(), Ms[0][1].clone(), Ms[0][2].clone()],
            config.const_coeff_witness_samples,
            config.const_coeff_matrix_entries,
        );
    // Keep this function bit-for-bit aligned with `prove_pi_fold_streaming_impl` except that each
    // instance pulls y = M*w from its own matrices.
    let prof = config.profile;
    if prof {
        eprintln!(
            "[PiFold streaming hetero] pool: rayon_threads={}, const_coeff_fastpath={}",
            rayon::current_num_threads(),
            const_coeff_fastpath
        );
    }
    let t_all = Instant::now();
    let ell = cms.len();
    let beta_cts = derive_beta_chi::<R>(transcript, ell);
    if ell == 0 || witnesses.len() != ell || Ms.len() != ell {
        return Err("PiFold: length mismatch".to_string());
    }

    // -----------------
    // Shared Π_had coins
    // -----------------
    let m = Ms[0][0].nrows;
    if !m.is_power_of_two() {
        return Err("PiFold: m must be power-of-two".to_string());
    }
    for inst_idx in 0..ell {
        for i in 0..3 {
            if Ms[inst_idx][i].nrows != m {
                return Err("PiFold: inconsistent m across instances".to_string());
            }
        }
    }
    let log_m = log2(m) as usize;
    let d = R::dimension();

    transcript.absorb_field_element(&R::BaseRing::from(0x4c465053_50494841u128));
    let s_base = transcript.get_challenges(log_m);
    let alpha_base = transcript.get_challenge();

    let mut alpha_pows = Vec::with_capacity(d);
    let mut pow = R::BaseRing::ONE;
    for _ in 0..d {
        alpha_pows.push(R::from(pow));
        pow *= alpha_base;
    }

    // -----------------
    // Per-instance prover messages + coins
    // -----------------
    //
    // Ordering (Figure 2 / Figure 4):
    // - absorb cm_f, derive J
    // - compute g^(i) and commit c(i)=Commit(g^(i)); absorb c(i)
    // - derive Π_mon coins (c,beta,alpha,...) that depend on the commitments
    let mut cba_all: Vec<Vec<(Vec<R>, R::BaseRing, R::BaseRing)>> = Vec::with_capacity(ell);
    let mut rc_all: Vec<Option<R::BaseRing>> = Vec::with_capacity(ell);
    let mut Js: Vec<Vec<Vec<R::BaseRing>>> = Vec::with_capacity(ell);
    let mut proj_digits_by_inst: Vec<Arc<Vec<R::BaseRing>>> = Vec::with_capacity(ell);
    let mut cm_g: Vec<Vec<Vec<R>>> = Vec::with_capacity(ell);
    let t_coins = Instant::now();

    let n_f = witnesses[0].len();
    if n_f == 0 || n_f % rg_params.l_h != 0 {
        return Err("PiFold: invalid witness length".to_string());
    }
    for w in witnesses.iter() {
        if w.len() != n_f {
            return Err("PiFold: inconsistent witness lengths".to_string());
        }
    }
    for inst_idx in 0..ell {
        for i in 0..3 {
            if Ms[inst_idx][i].ncols != n_f {
                return Err("PiFold: matrix/witness width mismatch".to_string());
            }
        }
    }

    let blocks = n_f / rg_params.l_h;
    let m_j = blocks
        .checked_mul(rg_params.lambda_pj)
        .ok_or_else(|| "PiFold: m_J overflow".to_string())?;
    if m < m_j || m % m_j != 0 {
        return Err("PiFold: require m_J <= m and m multiple of m_J".to_string());
    }

    let g_len = m
        .checked_mul(d)
        .ok_or_else(|| "PiFold: g_len overflow".to_string())?;
    if !g_len.is_power_of_two() {
        return Err("PiFold: require m*d power-of-two".to_string());
    }
    let g_nvars = log2(g_len) as usize;

    // Const-coeff cm_g optimization caches A·1 once (it is independent of the instance).
    let mut a1_cached: Option<Vec<R>> = None;

    for (inst_idx, cm_f) in cms.iter().enumerate() {
        let f: &[R] = &witnesses[inst_idx];

        transcript.absorb_slice(cm_f);
        // PCS#1 (cm_f PCS) coin splice: bind to cm_f surface and derive C1/C2 bytes.
        //
        // We can deterministically derive the PCS params from:
        // - the witness length (n_f ring elements -> flat_len = n_f*d base-field elems)
        // - the public commitment surface length (kappa_commit = cm_f.len()*d base-field elems)
        //
        // This keeps the protocol self-contained: no additional params need to be carried in the proof.
        transcript.absorb_field_element(&R::BaseRing::from(CMF_PCS_DOMAIN_SEP));
        let flat_len = n_f * d;
        let kappa_commit = cm_f.len() * d;
        let pcs_params_cmf = cmf_pcs_params_for_flat_len::<R::BaseRing>(flat_len, kappa_commit)?;
        let n_bytes_cmf = cmf_pcs_coin_bytes_len(&pcs_params_cmf);
        let _ = transcript.squeeze_bytes(n_bytes_cmf);
        let J = derive_J::<R>(transcript, rg_params.lambda_pj, rg_params.l_h);
        Js.push(J.clone());

        // Precompute projected digits once per instance.
        // Layout: digits[(row * d + col) * k_g + dig] where
        // - row ∈ [0, m_j)
        // - col ∈ [0, d)
        // - dig ∈ [0, k_g)
        let k_g = rg_params.k_g;
        let lambda_pj = rg_params.lambda_pj;
        let l_h = rg_params.l_h;
        let d_prime = rg_params.d_prime;

        let mut digits_flat = vec![R::BaseRing::ZERO; m_j * d * k_g];
            let out_ptr = digits_flat.as_mut_ptr() as usize;
        (0..m_j).into_par_iter().for_each(|out_row| {
                let i = out_row % lambda_pj;
                let b = out_row / lambda_pj;
                if const_coeff_fastpath {
                    // Constant-coeff specialization: only coeff[0] can be nonzero.
                    let mut h0 = R::BaseRing::ZERO;
                    for t in 0..l_h {
                        let in_row = b * l_h + t;
                        if in_row >= f.len() {
                            continue;
                        }
                        let coef = J[i][t];
                        h0 += coef * f[in_row].coeffs()[0];
                    }
                    let mut h_row = vec![R::BaseRing::ZERO; d];
                    h_row[0] = h0;
                    let digits = h_row.decompose_to_vec(d_prime, k_g);
                    for dig in 0..k_g {
                        let idx = (out_row * d) * k_g + dig; // col=0
                        unsafe {
                            *(out_ptr as *mut R::BaseRing).add(idx) = digits[0][dig];
                        }
                    }
                } else {
                    let mut h_row = vec![R::BaseRing::ZERO; d];
                    for t in 0..l_h {
                        let in_row = b * l_h + t;
                        if in_row >= f.len() {
                            continue;
                        }
                        let coef = J[i][t];
                        let coeffs = f[in_row].coeffs();
                        for col in 0..d {
                            h_row[col] += coef * coeffs[col];
                        }
                    }
                    let digits = h_row.decompose_to_vec(d_prime, k_g);
                    for col in 0..d {
                        for dig in 0..k_g {
                            let idx = (out_row * d + col) * k_g + dig;
                            unsafe {
                                *(out_ptr as *mut R::BaseRing).add(idx) = digits[col][dig];
                            }
                        }
                    }
                }
            });

        let digits_flat = Arc::new(digits_flat);
        proj_digits_by_inst.push(digits_flat.clone());

        // Commit to g^(i) and bind commitments into the transcript before Π_mon challenges.
        //
        // Const-coeff fast path optimization:
        // - for col != 0, digits are 0 so g^(dig)[j] = exp(0) is a constant g0
        // - write g = g0 * 1 + Δ where Δ is supported only on col=0
        // - then Commit(g) = g0*Commit(1) + Commit(Δ)
        //
        // This reduces the expensive Ajtai matvec from n=m*d positions to just m positions.
        let k_g = rg_params.k_g;
        let n = m * d;
        let g0 = exp::<R>(R::BaseRing::ZERO).expect("Exp failed");
        let cm_g_inst: Vec<Vec<R>> = if const_coeff_fastpath && d > 1 {
            // Compute A·1 once (shared across instances). Then scale by g0 per proof.
            if a1_cached.is_none() {
                let a1 = cm_g_scheme
                    .commit_many_with(n, 1, |_j, out: &mut [R]| out[0] = R::ONE)
                    .map_err(|e| format!("PiFold: cm_g commit(A·1) failed: {e:?}"))?;
                a1_cached = Some(a1[0].as_ref().to_vec());
            }
            let a1_vec = a1_cached.as_ref().expect("a1_cached present");
            let base = a1_vec.iter().map(|x| *x * g0).collect::<Vec<_>>();

            // Commit the correction Δ over only the col=0 range j ∈ [0, m).
            let digits_flat_for_commit = digits_flat.clone();
            let delta = cm_g_scheme
                .commit_many_with_range(n, k_g, 0, m, move |j, out: &mut [R]| {
                    let out_row = j % m_j; // since col=0, row=j
                    for dig in 0..k_g {
                        let digit0 = digits_flat_for_commit[(out_row * d) * k_g + dig]; // col=0
                        out[dig] = exp::<R>(digit0).expect("Exp failed") - g0;
                    }
                })
                .map_err(|e| format!("PiFold: cm_g commit(Δ) failed: {e:?}"))?;

            // Add base + delta.
            delta
                .into_iter()
                .map(|c| {
                    let mut v = base.clone();
                    for (acc, x) in v.iter_mut().zip(c.as_ref().iter()) {
                        *acc += *x;
                    }
                    v
                })
                .collect::<Vec<_>>()
        } else {
            // Generic path: commit the full g^(i) vectors.
            let digits_flat_for_commit = digits_flat.clone();
            cm_g_scheme
                .commit_many_with(n, k_g, move |j, out: &mut [R]| {
                    let col = j / m;
                    let row = j - col * m;
                    let out_row = row % m_j;
                    for dig in 0..k_g {
                        let digit = digits_flat_for_commit[(out_row * d + col) * k_g + dig];
                        out[dig] = exp::<R>(digit).expect("Exp failed");
                    }
                })
                .map_err(|e| format!("PiFold: cm_g commit failed: {e:?}"))?
                .into_iter()
                .map(|c| c.as_ref().to_vec())
                .collect::<Vec<_>>()
        };

        // Store cm_g for this instance (absorption happens after all instances computed).
        cm_g.push(cm_g_inst);
    }

    // -----------------------------------------------------------------------
    // Aggregate cm_g absorption: bind all cm_g commitments into transcript via a single
    // short Ajtai commitment (reduces Poseidon permutations from O(ell*k_g*kappa) to O(kappa)).
    // -----------------------------------------------------------------------
    let kappa = cm_g_scheme.kappa();
    let cm_g_agg = compute_cm_g_aggregate(&cm_g, kappa)?;
    transcript.absorb_slice(&cm_g_agg);

    // -----------------------------------------------------------------------
    // Phase 2: Derive Π_mon coins for all instances (now bound to the aggregate).
    // -----------------------------------------------------------------------
    for _inst_idx in 0..ell {
        let mut cba: Vec<(Vec<R>, R::BaseRing, R::BaseRing)> = Vec::with_capacity(rg_params.k_g);
        for _ in 0..rg_params.k_g {
            let c: Vec<R> = transcript
                .get_challenges(g_nvars)
                .into_iter()
                .map(|x| x.into())
                .collect();
            let beta = transcript.get_challenge();
            let alpha = transcript.get_challenge();
            cba.push((c, beta, alpha));
        }
        let rc: Option<R::BaseRing> = (rg_params.k_g > 1).then(|| transcript.get_challenge());
        cba_all.push(cba);
        rc_all.push(rc);
    }

    if prof {
        eprintln!(
            "[PiFold streaming hetero] coins+J+agg: {:?} (ell={}, k_g={}, g_nvars={})",
            t_coins.elapsed(),
            ell,
            rg_params.k_g,
            g_nvars
        );
    }

    let rhos = transcript
        .get_challenges(ell)
        .into_iter()
        .map(R::from)
        .collect::<Vec<R>>();

    // -----------------
    // STREAMING Hadamard MLEs
    // -----------------
    let t_had_build = Instant::now();
    let mut t_w0 = std::time::Duration::from_secs(0);
    // `eq(s, ·)` is shared across all instances; include it once.
    let mut mles_had: Vec<StreamingMleEnum<R>> = Vec::with_capacity(1 + ell * 3);
    mles_had.push(StreamingMleEnum::eq_base(s_base.clone()));
    let mut w0_cache: HashMap<usize, Arc<Vec<R::BaseRing>>> = HashMap::new();
    for inst_idx in 0..ell {
        let w0 = if const_coeff_fastpath {
            let key = Arc::as_ptr(&witnesses[inst_idx]) as usize;
            if let Some(v) = w0_cache.get(&key) {
                Some(v.clone())
            } else {
                let t0 = Instant::now();
                let v = Arc::new(
                    witnesses[inst_idx]
                        .iter()
                        .map(|x| x.coeffs()[0])
                        .collect::<Vec<_>>(),
                );
                t_w0 += t0.elapsed();
                w0_cache.insert(key, v.clone());
                Some(v)
            }
        } else {
            None
        };
        for i in 0..3 {
            if const_coeff_fastpath {
                mles_had.push(StreamingMleEnum::sparse_mat_vec_const_coeff(
                    Ms[inst_idx][i].clone(),
                    w0.as_ref().expect("w0 present").clone(),
                ));
            } else {
                mles_had.push(StreamingMleEnum::sparse_mat_vec(
                    Ms[inst_idx][i].clone(),
                    witnesses[inst_idx].clone(),
                ));
            }
        }
    }
    if prof {
        eprintln!(
            "[PiFold streaming hetero] had MLE build: {:?} (w0={:?}, const_coeff_fastpath={})",
            t_had_build.elapsed(),
            t_w0,
            const_coeff_fastpath
        );
    }

    let rhos_had = rhos.clone();
    let comb_had0: Option<Box<dyn Fn(&[R::BaseRing]) -> R::BaseRing + Sync + Send>> =
        if const_coeff_fastpath {
            let rhos0 = rhos_had.iter().map(|x| x.coeffs()[0]).collect::<Vec<_>>();
            let alpha0 = alpha_pows[0].coeffs()[0];
            Some(Box::new(move |vals: &[R::BaseRing]| -> R::BaseRing {
                let eq0 = vals[0];
                let mut acc_all0 = R::BaseRing::ZERO;
                for inst_idx in 0..ell {
                    let base = 1 + inst_idx * 3;
                    let y10 = vals[base];
                    let y20 = vals[base + 1];
                    let y30 = vals[base + 2];
                    let term0 = y10 * y20 - y30;
                    acc_all0 += rhos0[inst_idx] * (alpha0 * eq0 * term0);
                }
                acc_all0
            }))
        } else {
            None
        };

    let comb_had: Box<dyn Fn(&[R]) -> R + Sync + Send> = Box::new(move |vals: &[R]| -> R {
        let eq = vals[0];
        let mut acc_all = R::ZERO;
        for inst_idx in 0..ell {
            let base = 1 + inst_idx * 3;
            let y1 = &vals[base];
            let y2 = &vals[base + 1];
            let y3 = &vals[base + 2];
            let mut acc = R::ZERO;
            for j in 0..d {
                let term = y1.coeffs()[j] * y2.coeffs()[j] - y3.coeffs()[j];
                acc += alpha_pows[j] * eq * R::from(term);
            }
            acc_all += rhos_had[inst_idx] * acc;
        }
        acc_all
    });

    // -----------------
    // Build monomial side MLEs (identical to shared-M path)
    // -----------------
    let mon_mles_per = 3 * rg_params.k_g;
    let mut mles_mon: Vec<StreamingMleEnum<R>> = Vec::with_capacity(ell * mon_mles_per);

    let alphas_ring = cba_all
        .iter()
        .map(|cba| cba.iter().map(|(_, _, a)| R::from(*a)).collect::<Vec<_>>())
        .collect::<Vec<_>>();

    // NOTE: projected digits and cm_g commitments were computed above, and cm_g commitments were
    // absorbed into the transcript before Π_mon coins were derived.

    let t_mon_build = Instant::now();
    for inst_idx in 0..ell {
        let digits_flat = &proj_digits_by_inst[inst_idx];
        for dig in 0..rg_params.k_g {
            let (_c, beta_i, _alpha_i) = &cba_all[inst_idx][dig];

            // mj_compact[col*m_j + out_row] := ev(exp(digit), beta_i)
            let mut mj_compact = vec![R::BaseRing::ZERO; m_j * d];
                let out_ptr = mj_compact.as_mut_ptr() as usize;
                let digits_flat = digits_flat.clone();
                let beta_i = *beta_i;
                let dig_local = dig;
                if const_coeff_fastpath && d > 1 {
                    let g0 = exp::<R>(R::BaseRing::ZERO).expect("Exp failed");
                    let mjv_zero = ev(&g0, beta_i);
                (0..m_j).into_par_iter().for_each(|out_row| {
                        let digit0 = digits_flat[(out_row * d) * rg_params.k_g + dig_local]; // col=0
                        let g = exp::<R>(digit0).expect("Exp failed");
                        let mjv0 = ev(&g, beta_i);
                        unsafe {
                            let out = (out_ptr as *mut R::BaseRing).add(out_row);
                            *out = mjv0;
                        }
                    });
                    for col in 1..d {
                        mj_compact[col * m_j..(col + 1) * m_j].fill(mjv_zero);
                    }
                } else {
                (0..m_j).into_par_iter().for_each(|out_row| {
                        for col in 0..d {
                            let digit = digits_flat[(out_row * d + col) * rg_params.k_g + dig_local];
                            let g = exp::<R>(digit).expect("Exp failed");
                            let mjv = ev(&g, beta_i);
                            let idx = col * m_j + out_row;
                            unsafe {
                                let out = (out_ptr as *mut R::BaseRing).add(idx);
                                *out = mjv;
                            }
                        }
                    });
            }
            let mj_compact = Arc::new(mj_compact);
            mles_mon.push(StreamingMleEnum::periodic_base_scalar_vec(
                g_nvars,
                m,
                m_j,
                d,
                mj_compact.clone(),
                false,
            ));
            mles_mon.push(StreamingMleEnum::periodic_base_scalar_vec(
                g_nvars,
                m,
                m_j,
                d,
                mj_compact,
                true,
            ));

            // eq(c) as streaming MLE.
            let (c, _beta_i2, _alpha_i2) = &cba_all[inst_idx][dig];
            let c_base = c.iter().map(|x| x.coeffs()[0]).collect::<Vec<_>>();
            mles_mon.push(StreamingMleEnum::eq_base(c_base));
        }
    }
    if prof {
        eprintln!(
            "[PiFold streaming hetero] built mon MLEs in {:?} (g_len={}, g_nvars={}, k_g={}, ell={})",
            t_mon_build.elapsed(),
            g_len,
            g_nvars,
            rg_params.k_g,
            ell
        );
    }

    let rc_all_clone = rc_all.clone();
    let comb_mon = move |vals: &[R]| -> R {
        let mut acc_all = R::ZERO;
        for inst_idx in 0..ell {
            let mut lc = R::ZERO;
            let mut rc_pow = R::BaseRing::ONE;
            for dig in 0..rg_params.k_g {
                let base = inst_idx * mon_mles_per + dig * 3;
                let b_claim = vals[base] * vals[base] - vals[base + 1];
                let mut res = b_claim * alphas_ring[inst_idx][dig];
                res *= vals[base + 2];
                lc += if let Some(rc) = &rc_all_clone[inst_idx] {
                    let scaled = res * R::from(rc_pow);
                    rc_pow *= *rc;
                    scaled
                } else {
                    res
                };
            }
            acc_all += rhos[inst_idx] * lc;
        }
        acc_all
    };

    // -----------------
    // Run sumchecks with SHARED challenges (must match verifier schedule)
    // -----------------
    let hook_round = log2(m_j.next_power_of_two()) as usize;
    let mut v_digits_folded: Vec<Vec<R::BaseRing>> =
        vec![vec![R::BaseRing::ZERO; d]; rg_params.k_g];

    // Use absorb_field_element for scalars (saves ~(d-1) base-field absorbs each).
    transcript.absorb_field_element(&R::BaseRing::from(log_m as u128));
    transcript.absorb_field_element(&R::BaseRing::from(3u128));
    transcript.absorb_field_element(&R::BaseRing::from(g_nvars as u128));
    transcript.absorb_field_element(&R::BaseRing::from(3u128));

    let mut had_state = StreamingSumcheck::prover_init(mles_had, log_m, 3);
    let mut mon_state = StreamingSumcheck::prover_init(mles_mon, g_nvars, 3);

    let mut had_msgs = Vec::with_capacity(log_m);
    let mut mon_msgs = Vec::with_capacity(g_nvars);
    let mut sampled: Vec<R::BaseRing> = Vec::with_capacity(log_m.max(g_nvars));
    let mut v_msg_had: Option<R::BaseRing> = None;
    let mut v_msg_mon: Option<R::BaseRing> = None;

    let rounds = log_m.max(g_nvars);
    let t_sumcheck = Instant::now();
    for round_idx in 0..rounds {
        if round_idx < log_m {
            let pm_had = if const_coeff_fastpath {
                StreamingSumcheck::prove_round_base(
                    &mut had_state,
                    v_msg_had,
                    comb_had0.as_ref().expect("comb_had0").as_ref(),
                )
            } else {
                StreamingSumcheck::prove_round(&mut had_state, v_msg_had, comb_had.as_ref())
            };
            transcript.absorb_slice(&pm_had.evaluations);
            had_msgs.push(pm_had);
        }
        if round_idx < g_nvars {
            let pm_mon = StreamingSumcheck::prove_round(&mut mon_state, v_msg_mon, &comb_mon);
            transcript.absorb_slice(&pm_mon.evaluations);
            mon_msgs.push(pm_mon);
        }

        let r = transcript.get_challenge();
        transcript.absorb_field_element(&r);  // scalar absorb (saves ~(d-1) per round)
        sampled.push(r);

        if hook_round != 0 && sampled.len() == hook_round {
            let ts_r_full = ts_weights(&sampled);
            let ts_r = &ts_r_full[..m_j];
            for dig in 0..rg_params.k_g {
                for col in 0..d {
                    v_digits_folded[dig][col] = R::BaseRing::ZERO;
                }
            }
            for inst_idx in 0..ell {
                let b = beta_cts[inst_idx];
                let digits_flat = &proj_digits_by_inst[inst_idx];
                for row in 0..m_j {
                    let w = ts_r[row];
                    for col in 0..d {
                        for dig in 0..rg_params.k_g {
                            let digit = digits_flat[(row * d + col) * rg_params.k_g + dig];
                            v_digits_folded[dig][col] += b * digit * w;
                        }
                    }
                }
            }
            // Transcript compression: absorb a short aggregate commitment to v_digits_folded
            // (instead of k_g*d scalar absorbs).
            let v_digits_agg =
                compute_v_digits_folded_aggregate::<R>(&v_digits_folded, kappa)?;
            transcript.absorb_slice(&v_digits_agg);
        }

        v_msg_had = Some(r);
        v_msg_mon = Some(r);
    }
    if prof {
        eprintln!(
            "[PiFold streaming hetero] sumcheck: {:?} (rayon_threads={}, const_coeff_fastpath={})",
            t_sumcheck.elapsed(),
            rayon::current_num_threads(),
            const_coeff_fastpath
        );
        eprintln!("[PiFold streaming hetero] total: {:?}", t_all.elapsed());
    }

    let had_sumcheck = convert_streaming_proof(&StreamingProof(had_msgs));
    let mon_sumcheck = convert_streaming_proof(&StreamingProof(mon_msgs));

    // -----------------
    // Compute aux witness (match prove_pi_fold_streaming)
    // -----------------
    // Final randomness vectors for each sumcheck.
    let had_rand = sampled[..log_m].to_vec();
    let mon_rand = sampled[..g_nvars].to_vec();
    had_state.randomness = had_rand.clone();
    mon_state.randomness = mon_rand.clone();
    if log_m > 0 {
        had_state.fix_last_variable(*had_rand.last().expect("had_rand non-empty"));
    }
    if g_nvars > 0 {
        mon_state.fix_last_variable(*mon_rand.last().expect("mon_rand non-empty"));
    }

    let had_evals = had_state.final_evals();
    let mut had_u: Vec<[Vec<R::BaseRing>; 3]> = Vec::with_capacity(ell);
    for inst_idx in 0..ell {
        // MLE order: [eq(s,r)] then per instance: [y1, y2, y3]
        let base = 1 + inst_idx * 3;
        let y1 = had_evals[base];
        let y2 = had_evals[base + 1];
        let y3 = had_evals[base + 2];
        let mut U: [Vec<R::BaseRing>; 3] =
            [Vec::with_capacity(d), Vec::with_capacity(d), Vec::with_capacity(d)];
        U[0].extend_from_slice(y1.coeffs());
        U[1].extend_from_slice(y2.coeffs());
        U[2].extend_from_slice(y3.coeffs());
        had_u.push(U);
    }

    // IMPORTANT: transcript schedule alignment with verifier.
    //
    // The verifier absorbs the hadamard-evaluated U vectors (or aux-provided U) into the transcript
    // before proceeding to Π_mon / batchlin PCS. Without this, the transcript diverges and the
    // batchlin PCS coins will mismatch (making PCS#2 UNSAT in-gate).
    for U in &had_u {
        for x in &U[0] {
            transcript.absorb_field_element(x);
        }
        for x in &U[1] {
            transcript.absorb_field_element(x);
        }
        for x in &U[2] {
            transcript.absorb_field_element(x);
        }
    }

    let reps = m / m_j;
    let ts_r_mon = ts_weights(&mon_rand);
    let mut mon_b: Vec<Vec<R>> = Vec::with_capacity(ell);
    for inst_idx in 0..ell {
        let mut b_inst = Vec::with_capacity(rg_params.k_g);
        let digits_flat = proj_digits_by_inst[inst_idx].clone();
        for dig in 0..rg_params.k_g {
                let dig_local = dig;
                let digits_flat = digits_flat.clone();
            let acc = (0..m_j)
                .into_par_iter()
                    .map(|out_row| {
                        let mut local = R::ZERO;
                        for col in 0..d {
                            let digit = digits_flat[(out_row * d + col) * rg_params.k_g + dig_local];
                            let g = exp::<R>(digit).expect("Exp failed");
                            for rep in 0..reps {
                                let r = out_row + rep * m_j;
                                let idx = col * m + r;
                                local += R::from(ts_r_mon[idx]) * g;
                            }
                        }
                        local
                    })
                .reduce(|| R::ZERO, |a, b| a + b);
            b_inst.push(acc);
        }
        mon_b.push(b_inst);
    }

    // -----------------------------------------------------------------------
    // Aggregate mon_b absorption: bind all mon_b values into transcript via a single
    // short Ajtai commitment (reduces Poseidon permutations from O(ell*k_g) to O(kappa)).
    // -----------------------------------------------------------------------
    let mon_b_agg = compute_mon_b_aggregate(&mon_b, kappa)?;
    transcript.absorb_slice(&mon_b_agg);

    // -----------------------------------------------------------------------
    // Batchlin PCS transcript splice (domain-separated, REQUIRED):
    //
    // This is part of the canonical transcript schedule. It binds the batchlin PCS commitment
    // surface `batchlin_pcs_t` and derives PCS coins (C1/C2 bits) from Poseidon via SqueezeBytes.
    // The WE gate verifies this PCS (PCS#2) in-gate, so this phase is not optional.
    // -----------------------------------------------------------------------
    // Batched scalar PCS:
    // - We commit to a single scalar-valued MLE f[j] = Σ_dig γ^dig * ct(psi * g_folded^{(dig)}[j]).
    // - The gate binds the PCS evaluation to the Π_fold Step-5 scalars via batching with the same γ.
    transcript.absorb_field_element(&R::BaseRing::from(BATCHLIN_PCS_DOMAIN_SEP));
    let gamma = transcript.get_challenge();

    // Build x0,x1,x2 from r' (monomial sumcheck point).
    let log_n = g_nvars;
    let cfg_pcs = BatchlinPcsConfig::for_sp1(log_n, 1, rg_params.k_g)?;
    let (x0, x1, x2) = mle_point_to_tensor::<R::BaseRing>(&mon_rand, &cfg_pcs)?;

    // Compute f_batch (length 2^{log_n}) over the base field.
    // Use the same column-major layout as cm_g commit: idx = col*m + row.
    let psi_poly = stark_rings::psi::<R>();
    let psi_coeffs = psi_poly.coeffs().to_vec();
    let ctpsi = |g: &R| -> R::BaseRing {
        let coeffs = g.coeffs();
        let mut acc = psi_coeffs[0] * coeffs[0];
        for i_c in 1..d {
            acc -= psi_coeffs[i_c] * coeffs[d - i_c];
            }
                acc
            };

    let beta_sum = beta_cts.iter().copied().fold(R::BaseRing::ZERO, |a, b| a + b);
    let k_g = rg_params.k_g;

    // Precompute Σ_{dig=0..k_g-1} gamma^dig
    let mut gamma_pow = R::BaseRing::ONE;
    let mut gamma_series = R::BaseRing::ZERO;
    for _ in 0..k_g {
        gamma_series += gamma_pow;
        gamma_pow *= gamma;
    }

    let g0 = exp::<R>(R::BaseRing::ZERO).expect("Exp failed");
    let ctpsi_g0 = ctpsi(&g0);

    let mut f_batch = vec![R::BaseRing::ZERO; g_len];
    f_batch
        .par_iter_mut()
        .enumerate()
        .for_each(|(j, slot)| {
            let col = j / m;
            if const_coeff_fastpath && d > 1 && col != 0 {
                // All digits are 0 => g is constant exp(0).
                *slot = beta_sum * ctpsi_g0 * gamma_series;
                return;
            }
            let row = j - col * m;
            let out_row = row % m_j;
            let mut acc = R::BaseRing::ZERO;
            let mut gamma_pow = R::BaseRing::ONE;
            for dig in 0..k_g {
                let mut dig_acc = R::BaseRing::ZERO;
                for inst_idx in 0..ell {
                    let digits_flat = &proj_digits_by_inst[inst_idx];
                    let digit = digits_flat[(out_row * d + col) * k_g + dig];
                    let g = exp::<R>(digit).expect("Exp failed");
                    let v = ctpsi(&g);
                    dig_acc += beta_cts[inst_idx] * v;
                }
                acc += gamma_pow * dig_acc;
                gamma_pow *= gamma;
            }
            *slot = acc;
        });

    // PCS params and commitment/opening.
    //
    // NOTE: We may pad the batchlin PCS domain to the next cube power (r^3 >= 2^{log_n}).
    // `mle_point_to_tensor` already extends the evaluation point with zeros accordingly.
    let pcs_params = batchlin_scalar_pcs_params::<R::BaseRing>(log_n)?;
    if pcs_params.f_len() < g_len {
        return Err("PiFold: batchlin scalar pcs params mismatch (f_len < g_len)".to_string());
    }
    if pcs_params.f_len() != g_len {
        // Zero-pad f_batch on the larger domain.
        f_batch.resize(pcs_params.f_len(), R::BaseRing::ZERO);
    }
    let (t_vec, s_open) = pcs_commit::<R::BaseRing>(&pcs_params, &f_batch)?;

    // Bind commitment surface into transcript (required by verifier schedule).
    let batchlin_pcs_t: Vec<Vec<R::BaseRing>> = vec![t_vec.clone()];
    transcript.absorb_field_element(&R::BaseRing::from(batchlin_pcs_t.len() as u128));
    transcript.absorb_field_element(&R::BaseRing::from(batchlin_pcs_t[0].len() as u128));
    for x in &batchlin_pcs_t[0] {
        transcript.absorb_field_element(x);
    }

    // Derive PCS coins from transcript bytes (C1/C2 bits).
    let c_bytes = transcript.squeeze_bytes(64);
    let mut bits = Vec::with_capacity(c_bytes.len() * 8);
    for &b in &c_bytes {
        for i in 0..8 {
            bits.push(((b >> i) & 1) == 1);
        }
    }
    let r_pcs = pcs_params.r;
    let c1 = BinMatrix {
        rows: r_pcs * pcs_params.kappa,
        cols: pcs_params.kappa,
        data: (0..(r_pcs * pcs_params.kappa * pcs_params.kappa))
            .map(|i| if bits[i] { R::BaseRing::ONE } else { R::BaseRing::ZERO })
            .collect(),
    };
    let c2 = BinMatrix {
        rows: r_pcs * pcs_params.kappa,
        cols: pcs_params.kappa,
        data: (0..(r_pcs * pcs_params.kappa * pcs_params.kappa))
            .map(|i| if bits[(r_pcs * pcs_params.kappa * pcs_params.kappa) + i] { R::BaseRing::ONE } else { R::BaseRing::ZERO })
            .collect(),
    };
    let (_u_scalar, batchlin_pcs_core) = pcs_open::<R::BaseRing>(&pcs_params, &f_batch, &s_open, &x0, &x1, &x2, &c1, &c2)?;

    // Sanity: native verify should pass.
    crate::pcs::folding_pcs_l2::verify_folding_pcs_l2_with_c_matrices(
        &pcs_params,
        &t_vec,
        &x0,
        &x1,
        &x2,
        &vec![crate::pcs::folding_pcs_l2::kron_ikn_xt_mul(&x2, pcs_params.kappa, pcs_params.n, &batchlin_pcs_core.v0)[0]],
        &batchlin_pcs_core,
        &c1,
        &c2,
    ).map_err(|e| format!("PiFold: batchlin PCS native verify failed: {e}"))?;

    Ok(PiFoldProverOutput {
        proof: PiFoldBatchedProof {
            coins: SymphonyCoins {
                challenges: vec![],
                bytes: vec![],
                events: vec![],
            },
            rg_params,
            m_j,
            m,
            cm_g,
            batchlin_pcs_t,
            batchlin_pcs_x0: x0,
            batchlin_pcs_x1: x1,
            batchlin_pcs_x2: x2,
            batchlin_pcs_core,
            v_digits_folded,
            had_sumcheck,
            mon_sumcheck,
        },
        aux: PiFoldAuxWitness { had_u, mon_b },
        cfs_had_u: vec![],
        cfs_mon_b: vec![],
    })
}

/// Convert streaming proof to latticefold Proof format.
pub fn convert_streaming_proof<R: OverField>(
    proof: &StreamingProof<R>,
) -> latticefold::utils::sumcheck::Proof<R> {
    use latticefold::utils::sumcheck::prover::ProverMsg;
    latticefold::utils::sumcheck::Proof::new(
        proof
            .0
            .iter()
            .map(|m| ProverMsg::new(m.evaluations.clone()))
            .collect(),
    )
}