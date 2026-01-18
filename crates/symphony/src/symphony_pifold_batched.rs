//! Symphony Π_fold (Figure 4) — **verifier** functions and data types.
//!
//! This module provides:
//! - Proof and output data types used by both prover and verifier
//! - Verifier functions for batched Π_fold proofs (shared-M and hetero-M variants)
//!
//! For the prover, use `prove_pi_fold_poseidon_fs` from `symphony_pifold_streaming`.
//!
//! Notes:
//! - We still recompute some per-instance values during verification from witness openings
//!   (correctness-first; not ZK yet).
//! - Full paper compression replaces these with commitment+opening checks (`VfyOpen`).

use ark_std::log2;
use ark_ff::PrimeField;
use cyclotomic_rings::rings::GetPoseidonParams;
use latticefold::{
    transcript::Transcript,
    utils::sumcheck::{
        utils::eq_eval,
        MLSumcheck, Proof,
    },
};
use stark_rings::{
    balanced_decomposition::{Decompose, DecomposeToVec},
    exp,
    psi,
    CoeffRing,
    OverField,
    Ring,
    Zq,
};
use stark_rings_linalg::SparseMatrix;
use stark_rings_poly::mle::DenseMultilinearExtension;

use crate::{
    symphony_coins::{derive_beta_chi, derive_J, ev, ts_weights},
    symphony_open::{NoOpen, VfyOpen},
    symphony_pifold_streaming::{compute_cm_g_aggregate, compute_mon_b_aggregate},
    pcs::batchlin_pcs::BATCHLIN_PCS_DOMAIN_SEP,
    pcs::cmf_pcs::{cmf_pcs_coin_bytes_len, cmf_pcs_params_for_flat_len, CMF_PCS_DOMAIN_SEP},
    pcs::folding_pcs_l2::FoldingPcsL2ProofCore,
    rp_rgchk::{compose_v_digits, RPParams},
    symphony_cm::SymphonyCoins,
    symphony_fold::{SymphonyBatchLin, SymphonyInstance},
    transcript::PoseidonTranscriptMetrics,
};

/// Domain separator for Π_fold transcript schedule (Poseidon-FS and FS-replay variants).
///
/// This tag is absorbed at the start of the Π_fold verifier schedule. Outer protocols
/// may add their own domain separators *before* invoking Π_fold in the same transcript instance.
const DS_PI_FOLD: u128 = 0x4c465053_50494250u128; // "LFPS_PIBP"

/// Matrix input for Π_fold verification.
///
/// - `Shared` is the classic case where all instances share the same (M1,M2,M3).
/// - `Hetero` is the SP1/chunking case where each instance i has its own (M1_i,M2_i,M3_i).
#[derive(Clone, Copy, Debug)]
pub enum PiFoldMatrices<'a, R: CoeffRing>
where
    R::BaseRing: Zq,
{
    Shared([&'a SparseMatrix<R>; 3]),
    Hetero(&'a [[&'a SparseMatrix<R>; 3]]),
}

/// Output of CP-style Π_fold verification under Poseidon-FS.
///
/// This is the single canonical “relation verifier” output:
/// - `result`: the folded accumulator output (or a verifier error string), and
/// - transcript metrics + Poseidon trace (for algebraic/DPP frontends, even on failure).
#[derive(Clone, Debug)]
pub struct PiFoldCpPoseidonFsAttempt<R: CoeffRing>
where
    R::BaseRing: Zq,
{
    pub result: Result<(SymphonyInstance<R>, SymphonyBatchLin<R>), String>,
    pub metrics: PoseidonTranscriptMetrics,
    pub trace: crate::transcript::PoseidonTranscriptTrace<<R::BaseRing as ark_ff::Field>::BasePrimeField>,
}



#[derive(Clone, Debug)]
pub struct PiFoldBatchedProof<R: OverField>
where
    R::BaseRing: Zq + PrimeField,
{
    pub coins: SymphonyCoins<R>,
    pub rg_params: RPParams,

    /// m_J for the projection (n*lambda_pj/l_h).
    pub m_j: usize,
    /// m used by Π_mon (row-domain size).
    pub m: usize,

    /// Per-instance monomial commitment vectors `c(i)` (Figure 2 Step 3), i.e. Ajtai commitments to
    /// the large monomial vectors `g^(i)`.
    ///
    /// Shape: `[ell][k_g][kappa]`.
    ///
    /// Notes:
    /// - These are prover messages that must be fixed **before** the Π_mon challenges are derived.
    /// - We do not open them in `R_cp`; they are part of the reduced/output relation `R_o`
    ///   (paper `R_batchlin`) and will be proved succinctly by `π_lin`.
    pub cm_g: Vec<Vec<Vec<R>>>,

    /// Batchlin PCS: PCS commitment surface for the **batched scalar** folded batchlin object.
    ///
    /// This commitment is transcript-bound and verified in the WE gate (PCS#2) after `r'` and
    /// `u_*` are fixed.
    ///
    /// Shape: `[k_g][t_len]` where `t_len` is the PCS commitment length for a single digit.
    pub batchlin_pcs_t: Vec<Vec<R::BaseRing>>,

    /// Batchlin PCS: tensor eval vectors derived from `r'` (length r each).
    ///
    /// For the batched (single-proof) variant we commit to a *single* scalar-valued MLE, so these
    /// are shared across all digits and are public transcript-bound values.
    pub batchlin_pcs_x0: Vec<R::BaseRing>,
    pub batchlin_pcs_x1: Vec<R::BaseRing>,
    pub batchlin_pcs_x2: Vec<R::BaseRing>,

    /// Batchlin PCS: PCS proof core for the **batched scalar** commitment (verified in the WE gate as PCS#2).
    pub batchlin_pcs_core: FoldingPcsL2ProofCore<R::BaseRing>,

    /// Folded Π_rg hook message v_digits* (k_g × d), using scalar β (constants in the base field).
    pub v_digits_folded: Vec<Vec<R::BaseRing>>,

    /// Batched sumcheck proofs.
    pub had_sumcheck: Proof<R>,
    pub mon_sumcheck: Proof<R>,
}

/// Canonical prover output for Π_fold in this crate.
///
/// - `proof`: the public proof object (what a verifier receives)
/// - `aux`: the prover-side transcript witness messages (`m_i`) that a **CP-style relation**
///   binds via commitments (`cfs_*`) and opens during verification.
#[derive(Clone, Debug)]
pub struct PiFoldProverOutput<R: OverField>
where
    R::BaseRing: Zq + PrimeField,
{
    pub proof: PiFoldBatchedProof<R>,
    pub aux: PiFoldAuxWitness<R>,
    /// Optional CP transcript-message commitments (paper `c_fs,i = Commit(m_i)`), computed from `aux`.
    ///
    /// These are only populated when the prover is given commitment schemes for these messages.
    pub cfs_had_u: Vec<Vec<R>>,
    pub cfs_mon_b: Vec<Vec<R>>,
}

/// Auxiliary witness messages for CP-style relation checks.
///
/// In the paper’s CP relation, these are part of the hidden folding transcript messages `m_i`
/// that are committed to; here we plumb them explicitly so the verifier (and thus a DPP target)
/// does **not** need to recompute them from the full witness openings.
#[derive(Clone, Debug)]
pub struct PiFoldAuxWitness<R: OverField>
where
    R::BaseRing: Zq + PrimeField,
{
    /// Per-instance Π_had “U” message at the had-sumcheck evaluation point:
    /// for each instance i: U[i] = [U1(d), U2(d), U3(d)] where each is length `d = dim(R)`.
    pub had_u: Vec<[Vec<R::BaseRing>; 3]>,
    /// Per-instance Π_rg monomial message `b` at the mon-sumcheck evaluation point:
    /// for each instance i: b[i][dig] is a ring element for dig ∈ [0, k_g).
    pub mon_b: Vec<Vec<R>>,
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

fn absorb_public_inputs<R: OverField>(ts: &mut impl Transcript<R>, public_inputs: &[R::BaseRing])
where
    R::BaseRing: Zq,
{
    // Domain separate the *statement binding* portion from the rest of the folding transcript.
    ts.absorb_field_element(&R::BaseRing::from(0x4c465053_5055424cu128)); // "LFPS_PUBL"
    ts.absorb_field_element(&R::BaseRing::from(public_inputs.len() as u128));
    for x in public_inputs {
        ts.absorb_field_element(x);
    }
}

/// CP-style Π_fold verifier, but **run inside an existing transcript instance**.
///
/// This is the key plumbing for a *single transcript with domain-separated phases*:
/// an outer protocol can run Π_fold (R_cp) first, then continue with another phase (e.g. π_lin / R_o)
/// without resetting the transcript.
///
/// Notes:
/// - This absorbs Π_fold's own domain separators and public inputs.
/// - `cm_f` is treated as part of the statement (so we do not open it here).
/// - The CP transcript-message commitments (`cfs_*`) are verified (opened to `aux`) here.
fn verify_pi_fold_cp_in_transcript_hetero_m<R: CoeffRing>(
    ts: &mut impl Transcript<R>,
    Ms: &[[&SparseMatrix<R>; 3]],
    cm_f: &[Vec<R>],
    proof: &PiFoldBatchedProof<R>,
    open: &impl VfyOpen<R>,
    cfs_had_u: &[Vec<R>],
    cfs_mon_b: &[Vec<R>],
    aux: &PiFoldAuxWitness<R>,
    public_inputs: &[R::BaseRing],
) -> Result<(SymphonyInstance<R>, SymphonyBatchLin<R>), String>
where
    R::BaseRing: Zq + Decompose + PrimeField,
{
    ts.absorb_field_element(&R::BaseRing::from(DS_PI_FOLD));
    absorb_public_inputs::<R>(ts, public_inputs);

    if Ms.len() != cm_f.len() {
        return Err("PiFoldCP: Ms length mismatch".to_string());
    }
    if cfs_had_u.len() != cm_f.len() || cfs_mon_b.len() != cm_f.len() {
        return Err("PiFoldCP: cfs length mismatch".to_string());
    }
    let ell = cm_f.len();
    if ell == 0 {
        return Err("PiFoldCP: empty batch".to_string());
    }
    if aux.had_u.len() != ell || aux.mon_b.len() != ell {
        return Err("PiFoldCP: aux length mismatch".to_string());
    }

    // Core Π_fold verifier (aux supplies U/b, so no witness openings are required).
    let out = verify_pi_fold_batched_and_fold_outputs_with_openings_and_aux_hetero_m(
        ts,
        Ms,
        cm_f,
        proof,
        &NoOpen, // cm_f is not opened in CP relation
        &[],
        Some(aux),
    )?;

    // Verify CP commitment openings (domain-separated by `VfyOpen`).
    for i in 0..ell {
        let had_u_enc = encode_had_u_instance::<R>(&aux.had_u[i]);
        open.verify_opening(ts, "cfs_had_u", &cfs_had_u[i], &[], &had_u_enc, &[])?;
        open.verify_opening(ts, "cfs_mon_b", &cfs_mon_b[i], &[], &aux.mon_b[i], &[])?;
    }

    Ok(out)
}

/// CP-style Π_fold verifier, but **run inside an existing transcript instance** (shared-M).
///
/// This is the shared-matrix counterpart of `verify_pi_fold_cp_in_transcript_hetero_m`.
fn verify_pi_fold_cp_in_transcript<R: CoeffRing>(
    ts: &mut impl Transcript<R>,
    M: [&SparseMatrix<R>; 3],
    cm_f: &[Vec<R>],
    proof: &PiFoldBatchedProof<R>,
    open: &impl VfyOpen<R>,
    cfs_had_u: &[Vec<R>],
    cfs_mon_b: &[Vec<R>],
    aux: &PiFoldAuxWitness<R>,
    public_inputs: &[R::BaseRing],
) -> Result<(SymphonyInstance<R>, SymphonyBatchLin<R>), String>
where
    R::BaseRing: Zq + Decompose + PrimeField,
{
    ts.absorb_field_element(&R::BaseRing::from(DS_PI_FOLD));
    absorb_public_inputs::<R>(ts, public_inputs);

    if cfs_had_u.len() != cm_f.len() || cfs_mon_b.len() != cm_f.len() {
        return Err("PiFoldCP: cfs length mismatch".to_string());
    }
    let ell = cm_f.len();
    if ell == 0 {
        return Err("PiFoldCP: empty batch".to_string());
    }
    if aux.had_u.len() != ell || aux.mon_b.len() != ell {
        return Err("PiFoldCP: aux length mismatch".to_string());
    }

    // Core Π_fold verifier (aux supplies U/b, so no witness openings are required).
    let out = verify_pi_fold_batched_and_fold_outputs_with_openings_and_aux(
        ts,
        M,
        cm_f,
        proof,
        &NoOpen, // cm_f is not opened in CP relation
        &[],
        Some(aux),
    )?;

    // Verify CP commitment openings (domain-separated by `VfyOpen`).
    for i in 0..ell {
        let had_u_enc = encode_had_u_instance::<R>(&aux.had_u[i]);
        open.verify_opening(ts, "cfs_had_u", &cfs_had_u[i], &[], &had_u_enc, &[])?;
        open.verify_opening(ts, "cfs_mon_b", &cfs_mon_b[i], &[], &aux.mon_b[i], &[])?;
    }

    Ok(out)
}

/// Canonical CP-style Π_fold verifier under Poseidon-FS.
///
/// This is the **only** public verification API we keep:
/// - Poseidon-FS transcript binding is performed internally.
/// - CP transcript-message commitments (`cfs_*`) are opened to `aux` and verified.
/// - Returns the folded output and the full Poseidon trace (for algebraic/DPP frontends).
pub fn verify_pi_fold_cp_poseidon_fs<R: CoeffRing, PC>(
    matrices: PiFoldMatrices<'_, R>,
    cm_f: &[Vec<R>],
    proof: &PiFoldBatchedProof<R>,
    open: &impl VfyOpen<R>,
    cfs_had_u: &[Vec<R>],
    cfs_mon_b: &[Vec<R>],
    aux: &PiFoldAuxWitness<R>,
    public_inputs: &[R::BaseRing],
) -> PiFoldCpPoseidonFsAttempt<R>
where
    R::BaseRing: Zq + Decompose + PrimeField,
    PC: GetPoseidonParams<<<R>::BaseRing as ark_ff::Field>::BasePrimeField>,
{
    let mut ts = crate::transcript::TracePoseidonTranscript::<R>::empty::<PC>();

    let result = match matrices {
        PiFoldMatrices::Shared(M) => verify_pi_fold_cp_in_transcript::<R>(
            &mut ts,
            M,
            cm_f,
            proof,
            open,
            cfs_had_u,
            cfs_mon_b,
            aux,
            public_inputs,
        ),
        PiFoldMatrices::Hetero(Ms) => verify_pi_fold_cp_in_transcript_hetero_m::<R>(
            &mut ts,
            Ms,
            cm_f,
            proof,
            open,
            cfs_had_u,
            cfs_mon_b,
            aux,
            public_inputs,
        ),
    };

    PiFoldCpPoseidonFsAttempt {
        result,
        metrics: ts.metrics(),
        trace: ts.trace().clone(),
    }
}

/// Heterogeneous-matrix verifier with optional auxiliary witness messages (`had_u`, `mon_b`).
fn verify_pi_fold_batched_and_fold_outputs_with_openings_and_aux_hetero_m<R: CoeffRing>(
    transcript: &mut impl Transcript<R>,
    Ms: &[[&SparseMatrix<R>; 3]],
    cms: &[Vec<R>],
    proof: &PiFoldBatchedProof<R>,
    open: &impl VfyOpen<R>,
    cms_openings: &[Vec<R>],
    aux: Option<&PiFoldAuxWitness<R>>,
) -> Result<(SymphonyInstance<R>, SymphonyBatchLin<R>), String>
where
    R::BaseRing: Zq + Decompose + PrimeField,
{
    if cms.len() != Ms.len() {
        return Err("PiFold: Ms/cms length mismatch".to_string());
    }
    if !cms_openings.is_empty() {
        if cms.len() != cms_openings.len() {
            return Err("PiFold: cms/openings length mismatch".to_string());
        }
        for (cm, open_val) in cms.iter().zip(cms_openings.iter()) {
            open.verify_opening(transcript, "cm_witness", cm, &[], open_val, &[])?;
        }
    } else if aux.is_none() {
        return Err("PiFold: cms_openings required when aux is not provided".to_string());
    }

    let ell = cms.len();
    if ell == 0 {
        return Err("PiFold: length mismatch".to_string());
    }
    let rg_params = &proof.rg_params;

    let beta_cts = derive_beta_chi::<R>(transcript, ell);
    let beta_ring = beta_cts.iter().copied().map(R::from).collect::<Vec<R>>();

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

    transcript.absorb_field_element(&R::BaseRing::from(0x4c465053_50494841u128)); // "LFPS_PIHA"
    let s_base = transcript.get_challenges(log_m);
    let alpha_base = transcript.get_challenge();

    let s_poly: Vec<R> = s_base.iter().copied().map(R::from).collect();
    let m_j = proof.m_j;
    if m != proof.m {
        return Err("PiFold: proof m mismatch".to_string());
    }
    if m < m_j || m % m_j != 0 {
        return Err("PiFold: require m_J<=m".to_string());
    }

    let g_len = m * d;
    if !g_len.is_power_of_two() {
        return Err("PiFold: require m*d power-of-two".to_string());
    }
    let g_nvars = log2(g_len) as usize;

    let mut cba_all: Vec<Vec<(Vec<R>, R::BaseRing, R::BaseRing)>> = Vec::with_capacity(ell);
    let mut rc_all: Vec<Option<R::BaseRing>> = Vec::with_capacity(ell);
    let mut Js: Vec<Vec<Vec<R::BaseRing>>> = Vec::with_capacity(ell);

    let (_n_f, blocks) = if aux.is_none() {
        let n_f = cms_openings[0].len();
        if n_f == 0 || n_f % rg_params.l_h != 0 {
            return Err("PiFold: invalid witness length".to_string());
        }
        for w in cms_openings.iter() {
            if w.len() != n_f {
                return Err("PiFold: inconsistent witness lengths".to_string());
            }
        }
        let blocks = n_f / rg_params.l_h;
        let m_j_expected = blocks * rg_params.lambda_pj;
        if m_j_expected != m_j {
            return Err("PiFold: m_J mismatch".to_string());
        }
        // Sanity-check witness width matches matrices.
        for inst_idx in 0..ell {
            for i in 0..3 {
                if Ms[inst_idx][i].ncols != n_f {
                    return Err("PiFold: matrix/witness width mismatch".to_string());
                }
            }
        }
        (n_f, blocks)
    } else {
        (0usize, 0usize)
    };

    if proof.cm_g.len() != ell {
        return Err("PiFold: cm_g length mismatch".to_string());
    }

    // Phase 1: Absorb cm_f and derive J for each instance; validate cm_g structure.
    for (inst_idx, cm_f) in cms.iter().enumerate() {
        transcript.absorb_slice(cm_f);
        // PCS#1 (cm_f PCS) coin splice: bind to cm_f surface and derive C1/C2 bytes.
        // Deterministically parameterized by (n_f, cm_f surface length).
        transcript.absorb_field_element(&R::BaseRing::from(CMF_PCS_DOMAIN_SEP));
        let n_f = Ms[inst_idx][0].ncols;
        let flat_len = n_f * R::dimension();
        let kappa_commit = cm_f.len() * R::dimension();
        let pcs_params_cmf = cmf_pcs_params_for_flat_len::<R::BaseRing>(flat_len, kappa_commit)?;
        let n_bytes_cmf = cmf_pcs_coin_bytes_len(&pcs_params_cmf);
        let _ = transcript.squeeze_bytes(n_bytes_cmf);
        let J = derive_J::<R>(transcript, rg_params.lambda_pj, rg_params.l_h);
        Js.push(J);

        if proof.cm_g[inst_idx].len() != rg_params.k_g {
            return Err("PiFold: cm_g k_g mismatch".to_string());
        }
    }

    // Aggregate cm_g absorption: bind all cm_g commitments via a single short Ajtai commit.
    // Determine kappa from the first commitment (all should have same length).
    let kappa = if !proof.cm_g.is_empty() && !proof.cm_g[0].is_empty() {
        proof.cm_g[0][0].len()
    } else {
        return Err("PiFold: empty cm_g".to_string());
    };
    let cm_g_agg = compute_cm_g_aggregate(&proof.cm_g, kappa)?;
    transcript.absorb_slice(&cm_g_agg);

    // Phase 2: Derive Π_mon coins for all instances (now bound to the aggregate).
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

    let rhos = transcript
        .get_challenges(ell)
        .into_iter()
        .map(R::from)
        .collect::<Vec<R>>();

    // Transcript compression: bind v_digits_folded via a short aggregate commitment.
    let v_digits_agg =
        crate::symphony_pifold_streaming::compute_v_digits_folded_aggregate::<R>(
            &proof.v_digits_folded,
            kappa,
        )?;

    let hook_round = log2(m_j.next_power_of_two()) as usize;
    let (had_sc, mon_sc) = MLSumcheck::<R, _>::verify_two_as_subprotocol_shared_with_hook(
        transcript,
        log_m,
        3,
        R::ZERO,
        &proof.had_sumcheck,
        g_nvars,
        3,
        R::ZERO,
        &proof.mon_sumcheck,
        hook_round,
        |t, _sampled| {
            t.absorb_slice(&v_digits_agg);
        },
    )
    .map_err(|e| format!("PiFold: sumcheck verify failed: {e}"))?;

    // Recompute the expected hadamard sumcheck evaluation at the had point.
    let r_had = had_sc.point.clone();
    let eq_sr = eq_eval(&s_poly, &r_had.iter().copied().map(R::from).collect::<Vec<_>>())
        .map_err(|e| format!("PiFold: eq_eval failed: {e}"))?;

    let mut pow = R::BaseRing::ONE;
    let mut alpha_pows = Vec::with_capacity(d);
    for _ in 0..d {
        alpha_pows.push(R::from(pow));
        pow *= alpha_base;
    }

    let mut lhs = R::ZERO;
    for inst_idx in 0..ell {
        let U: [Vec<R::BaseRing>; 3] = if let Some(auxw) = aux {
            if auxw.had_u.len() != ell {
                return Err("PiFold: aux.had_u length mismatch".to_string());
            }
            auxw.had_u[inst_idx].clone()
        } else {
            let y = Ms[inst_idx]
                .iter()
                .map(|Mi| Mi.try_mul_vec(&cms_openings[inst_idx]).expect("mat-vec mul failed"))
                .collect::<Vec<Vec<R>>>();

            let mut U: [Vec<R::BaseRing>; 3] =
                [Vec::with_capacity(d), Vec::with_capacity(d), Vec::with_capacity(d)];
            for i in 0..3 {
                for j in 0..d {
                    let evals = (0..m)
                        .map(|row| R::from(y[i][row].coeffs()[j]))
                        .collect::<Vec<_>>();
                    let mle = DenseMultilinearExtension::from_evaluations_vec(log_m, evals);
                    let v = mle
                        .evaluate(&r_had.iter().copied().map(R::from).collect::<Vec<_>>())
                        .expect("MLE evaluate returned None");
                    U[i].push(v.ct());
                }
            }
            U
        };

        if U[0].len() != d || U[1].len() != d || U[2].len() != d {
            return Err("PiFold: aux had_U has wrong dimension".to_string());
        }

        let mut acc = R::ZERO;
        for j in 0..d {
            let u1 = R::from(U[0][j]);
            let u2 = R::from(U[1][j]);
            let u3 = R::from(U[2][j]);
            acc += alpha_pows[j] * (u1 * u2 - u3);
        }
        lhs += rhos[inst_idx] * (eq_sr * acc);

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
    if lhs != had_sc.expected_evaluation {
        return Err("PiFold: hadamard recomputation mismatch".to_string());
    }

    // -----------------
    // Verify the batched Π_mon linkage and folded Π_rg Step-5 check (does not depend on M entries).
    // -----------------
    let r_mon_r: Vec<R> = mon_sc.point.iter().copied().map(R::from).collect();

    let mon_b: Vec<Vec<R>> = if let Some(auxw) = aux {
        if auxw.mon_b.len() != ell {
            return Err("PiFold: aux.mon_b length mismatch".to_string());
        }
        for b_inst in auxw.mon_b.iter() {
            if b_inst.len() != rg_params.k_g {
                return Err("PiFold: aux.mon_b wrong k_g".to_string());
            }
        }
        auxw.mon_b.clone()
    } else {
        let mut mon_b: Vec<Vec<R>> = Vec::with_capacity(ell);
        let expand_row = |row: usize| -> usize { row % m_j };
        for inst_idx in 0..ell {
            // Build H_digits for this instance (on m_J rows).
            let f = &cms_openings[inst_idx];
            let mut H = vec![vec![R::BaseRing::ZERO; d]; m_j];
            let Jref = &Js[inst_idx];
            for b in 0..blocks {
                for i in 0..rg_params.lambda_pj {
                    let out_row = b * rg_params.lambda_pj + i;
                    for t in 0..rg_params.l_h {
                        let in_row = b * rg_params.l_h + t;
                        let coef = Jref[i][t];
                        for col in 0..d {
                            H[out_row][col] += coef * f[in_row].coeffs()[col];
                        }
                    }
                }
            }
            let mut H_digits: Vec<Vec<Vec<R::BaseRing>>> =
                vec![vec![vec![R::BaseRing::ZERO; d]; m_j]; rg_params.k_g];
            for r in 0..m_j {
                let row_digits = H[r].decompose_to_vec(rg_params.d_prime, rg_params.k_g);
                for c in 0..d {
                    for i in 0..rg_params.k_g {
                        H_digits[i][r][c] = row_digits[c][i];
                    }
                }
            }

            let mut b_inst = Vec::with_capacity(rg_params.k_g);
            for dig in 0..rg_params.k_g {
                let mut gi = Vec::with_capacity(g_len);
                for c in 0..d {
                    for r in 0..m {
                        gi.push(exp::<R>(H_digits[dig][expand_row(r)][c]).expect("Exp failed"));
                    }
                }
                let mle = DenseMultilinearExtension::from_evaluations_vec(g_nvars, gi);
                b_inst.push(mle.evaluate(&r_mon_r).unwrap());
            }
            mon_b.push(b_inst);
        }
        mon_b
    };

    // Aggregate mon_b absorption: bind all mon_b values into transcript via a single
    // short Ajtai commitment (reduces Poseidon permutations from O(ell*k_g) to O(kappa)).
    let mon_b_agg = compute_mon_b_aggregate(&mon_b, kappa)?;
    transcript.absorb_slice(&mon_b_agg);

    let v_expected = mon_sc.expected_evaluation;

    let mut ver = R::ZERO;
    for inst_idx in 0..ell {
        let mut inst_acc = R::ZERO;
        // Precompute rc^dig iteratively to avoid repeated pow().
        let mut rc_pow = R::BaseRing::ONE;
        for dig in 0..rg_params.k_g {
            let (c, beta_i, alpha_i) = &cba_all[inst_idx][dig];
            let eq = eq_eval(c, &r_mon_r).unwrap();
            let b_i = mon_b[inst_idx][dig];
            let ev1 = R::from(ev(&b_i, *beta_i));
            let ev2 = R::from(ev(&b_i, *beta_i * *beta_i));
            let b_claim = ev1 * ev1 - ev2;
            let mut term = eq * R::from(*alpha_i) * b_claim;
            if let Some(rc) = &rc_all[inst_idx] {
                term *= R::from(rc_pow);
                rc_pow *= *rc;
            }
            inst_acc += term;
        }
        ver += rhos[inst_idx] * inst_acc;
    }

    if ver != v_expected {
        return Err("PiFold: batched monomial recomputation mismatch".to_string());
    }

    // Π_rg Step-5 check on the **folded** values.
    let log_m = log2(m.next_power_of_two()) as usize;
    let s_chals = mon_sc.point[log_m..].to_vec();
    let ts_s_full = ts_weights(&s_chals);
    let ts_s = &ts_s_full[..d];

    // Fold u*: for each digit, u*(dig) := Σ_i beta_i * u_i(dig)
    let mut u_folded = vec![R::ZERO; rg_params.k_g];
    for inst_idx in 0..ell {
        let b = beta_ring[inst_idx];
        for dig in 0..rg_params.k_g {
            u_folded[dig] += b * mon_b[inst_idx][dig];
        }
    }

    for dig in 0..rg_params.k_g {
        let lhs = (psi::<R>() * u_folded[dig]).ct();
        let rhs = proof.v_digits_folded[dig]
            .iter()
            .zip(ts_s.iter())
            .fold(R::BaseRing::ZERO, |acc, (&vij, &t)| acc + vij * t);
        if lhs != rhs {
            return Err(format!("PiFold: folded Step5 mismatch at dig={dig}"));
        }
    }

    // Construct folded outputs directly.
    let mut c_star = vec![R::ZERO; cms[0].len()];
    for (inst_idx, cm) in cms.iter().enumerate() {
        for (acc, x) in c_star.iter_mut().zip(cm.iter()) {
            *acc += beta_ring[inst_idx] * *x;
        }
    }

    // r and r' are shared across the batch.
    let log_mj = log2(m_j.next_power_of_two()) as usize;
    let r_prime = mon_sc.point.clone();
    let r_star = r_prime[..log_mj].to_vec();

    // u* is the folded u_folded computed above; v* is composed from folded v_digits.
    let v_star = compose_v_digits::<R>(&proof.v_digits_folded, rg_params.d_prime);

    let mut v_rq = R::ZERO;
    for (i, c) in v_star.iter().enumerate() {
        v_rq.coeffs_mut()[i] = *c;
    }

    let folded_inst = SymphonyInstance {
        c: c_star,
        r: r_star,
        v: v_rq,
    };
    // Fold commitments c*(i) := Σ beta_j * c_j(i).
    let kappa = proof.cm_g[0][0].len();
    let mut c_g_folded = vec![vec![R::ZERO; kappa]; rg_params.k_g];
    for inst_idx in 0..ell {
        let b = beta_ring[inst_idx];
        for dig in 0..rg_params.k_g {
            for j in 0..kappa {
                c_g_folded[dig][j] += b * proof.cm_g[inst_idx][dig][j];
            }
        }
    }

    let folded_bat = SymphonyBatchLin { r_prime, c_g: c_g_folded, u: u_folded };

    // -----------------------------------------------------------------------
    // Batchlin PCS transcript splice (domain-separated, REQUIRED):
    //
    // If the proof includes a PCS commitment surface for the folded batchlin object, bind it
    // into the Poseidon transcript *after* `r'` and `u_*` are fixed, and then squeeze bytes to
    // provide transcript-derived coins for the in-gate PCS verifier.
    //

    if proof.batchlin_pcs_t.is_empty() {
        return Err("PiFold: missing batchlin_pcs_t (batchlin PCS commitment surface)".to_string());
    }
    // Batched batchlin PCS: single commitment (len=1).
    if proof.batchlin_pcs_t.len() != 1 {
        return Err("PiFold: batchlin_pcs_t expected len=1 (batched scalar PCS)".to_string());
    }
    transcript.absorb_field_element(&R::BaseRing::from(BATCHLIN_PCS_DOMAIN_SEP));
    // Batch scalar γ used for digit batching in PCS#2. This is transcript-derived and must be
    // sampled at a fixed point in the schedule (after `r'` and `u_*` are fixed, before PCS coins).
    transcript.get_challenge();
    // NOTE: `get_challenge()` already re-absorbs internally (Poseidon transcript); no extra absorb.
    transcript.absorb_field_element(&R::BaseRing::from(proof.batchlin_pcs_t.len() as u128));
    for dig in 0..proof.batchlin_pcs_t.len() {
        if proof.batchlin_pcs_t[dig].is_empty() {
            return Err(format!("PiFold: empty batchlin_pcs_t[{dig}]"));
        }
        transcript.absorb_field_element(&R::BaseRing::from(proof.batchlin_pcs_t[dig].len() as u128));
        for x in &proof.batchlin_pcs_t[dig] {
            transcript.absorb_field_element(x);
        }
    }
    let _ = transcript.squeeze_bytes(64);

    Ok((folded_inst, folded_bat))
}



/// Verifier with optional auxiliary witness messages (`had_u`, `mon_b`) to avoid
/// recomputing them from full witness openings.
fn verify_pi_fold_batched_and_fold_outputs_with_openings_and_aux<R: CoeffRing>(
    transcript: &mut impl Transcript<R>,
    M: [&SparseMatrix<R>; 3],
    cms: &[Vec<R>],
    proof: &PiFoldBatchedProof<R>,
    open: &impl VfyOpen<R>,
    cms_openings: &[Vec<R>],
    aux: Option<&PiFoldAuxWitness<R>>,
) -> Result<(SymphonyInstance<R>, SymphonyBatchLin<R>), String>
where
    R::BaseRing: Zq + Decompose + PrimeField,
{
    // Commitment-opening layer is verified **outside** the transcript schedule used for the
    // folding/sumcheck subprotocols.
    //
    // When `aux` is provided (CP/WE-facing path), the verifier does not need the full witness
    // openings; callers may pass `cms_openings = &[]` to skip opening checks for `cm_witness`.
    // (In that mode, `cm_f` is treated as part of the public statement.)
    if !cms_openings.is_empty() {
        if cms.len() != cms_openings.len() {
            return Err("PiFold: cms/openings length mismatch".to_string());
        }
        for (cm, open_val) in cms.iter().zip(cms_openings.iter()) {
            open.verify_opening(transcript, "cm_witness", cm, &[], open_val, &[])?;
        }
    } else if aux.is_none() {
        return Err("PiFold: cms_openings required when aux is not provided".to_string());
    }

    let ell = cms.len();
    if ell == 0 {
        return Err("PiFold: length mismatch".to_string());
    }
    let rg_params = &proof.rg_params;

    let beta_cts = derive_beta_chi::<R>(transcript, ell);
    let beta_ring = beta_cts.iter().copied().map(R::from).collect::<Vec<R>>();
    let m = M[0].nrows;
    if !m.is_power_of_two() {
        return Err("PiFold: m must be power-of-two".to_string());
    }
    let log_m = log2(m) as usize;
    let d = R::dimension();

    transcript.absorb_field_element(&R::BaseRing::from(0x4c465053_50494841u128)); // "LFPS_PIHA"
    let s_base = transcript.get_challenges(log_m);
    let alpha_base = transcript.get_challenge();

    let s_poly: Vec<R> = s_base.iter().copied().map(R::from).collect();
    let m_j = proof.m_j;
    if m != proof.m {
        return Err("PiFold: proof m mismatch".to_string());
    }
    if m < m_j || m % m_j != 0 {
        return Err("PiFold: require m_J<=m".to_string());
    }

    let g_len = m * d;
    if !g_len.is_power_of_two() {
        return Err("PiFold: require m*d power-of-two".to_string());
    }
    let g_nvars = log2(g_len) as usize;

    let mut cba_all: Vec<Vec<(Vec<R>, R::BaseRing, R::BaseRing)>> = Vec::with_capacity(ell);
    let mut rc_all: Vec<Option<R::BaseRing>> = Vec::with_capacity(ell);
    let mut Js: Vec<Vec<Vec<R::BaseRing>>> = Vec::with_capacity(ell);

    // If we don't have auxiliary transcript messages, we will recompute them from witness openings,
    // which requires a consistent witness length to define m_J.
    let (_n_f, blocks) = if aux.is_none() {
        let n_f = cms_openings[0].len();
        if n_f == 0 || n_f % rg_params.l_h != 0 {
            return Err("PiFold: invalid witness length".to_string());
        }
        for w in cms_openings.iter() {
            if w.len() != n_f {
                return Err("PiFold: inconsistent witness lengths".to_string());
            }
        }
        let blocks = n_f / rg_params.l_h;
        let m_j_expected = blocks * rg_params.lambda_pj;
        if m_j_expected != m_j {
            return Err("PiFold: m_J mismatch".to_string());
        }
        (n_f, blocks)
    } else {
        (0usize, 0usize)
    };

    if proof.cm_g.len() != ell {
        return Err("PiFold: cm_g length mismatch".to_string());
    }

    // Phase 1: Absorb cm_f and derive J for each instance; validate cm_g structure.
    for (inst_idx, cm_f) in cms.iter().enumerate() {
        transcript.absorb_slice(cm_f);
        transcript.absorb_field_element(&R::BaseRing::from(CMF_PCS_DOMAIN_SEP));
        let n_f = M[0].ncols;
        let flat_len = n_f * R::dimension();
        let kappa_commit = cm_f.len() * R::dimension();
        let pcs_params_cmf = cmf_pcs_params_for_flat_len::<R::BaseRing>(flat_len, kappa_commit)?;
        let n_bytes_cmf = cmf_pcs_coin_bytes_len(&pcs_params_cmf);
        let _ = transcript.squeeze_bytes(n_bytes_cmf);
        let J = derive_J::<R>(transcript, rg_params.lambda_pj, rg_params.l_h);
        Js.push(J);

        if proof.cm_g[inst_idx].len() != rg_params.k_g {
            return Err("PiFold: cm_g k_g mismatch".to_string());
        }
    }

    // Aggregate cm_g absorption: bind all cm_g commitments via a single short Ajtai commit.
    let kappa = if !proof.cm_g.is_empty() && !proof.cm_g[0].is_empty() {
        proof.cm_g[0][0].len()
    } else {
        return Err("PiFold: empty cm_g".to_string());
    };
    let cm_g_agg = compute_cm_g_aggregate(&proof.cm_g, kappa)?;
    transcript.absorb_slice(&cm_g_agg);

    // Phase 2: Derive Π_mon coins for all instances (now bound to the aggregate).
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

    let rhos = transcript
        .get_challenges(ell)
        .into_iter()
        .map(R::from)
        .collect::<Vec<R>>();

    // Transcript compression: bind v_digits_folded via a short aggregate commitment.
    let v_digits_agg =
        crate::symphony_pifold_streaming::compute_v_digits_folded_aggregate::<R>(
            &proof.v_digits_folded,
            kappa,
        )?;

    // Verify the two batched sumchecks with shared challenges.
    let hook_round = log2(m_j.next_power_of_two()) as usize;
    let (had_sc, mon_sc) = MLSumcheck::<R, _>::verify_two_as_subprotocol_shared_with_hook(
        transcript,
        log_m,
        3,
        R::ZERO,
        &proof.had_sumcheck,
        g_nvars,
        3,
        R::ZERO,
        &proof.mon_sumcheck,
        hook_round,
        |t, _sampled| {
            t.absorb_slice(&v_digits_agg);
        },
    )
    .map_err(|e| format!("PiFold: sumcheck verify failed: {e}"))?;

    // Batched Eq.(26) check for Π_had.
    // If `aux` is present, use it and avoid recomputing U from the full witness opening.
    let r_poly_had: Vec<R> = had_sc.point.iter().copied().map(R::from).collect();
    let eq_sr = eq_eval(&s_poly, &r_poly_had).map_err(|e| format!("PiFold: eq_eval failed: {e}"))?;

    let mut pow = R::BaseRing::ONE;
    let mut alpha_pows = Vec::with_capacity(d);
    for _ in 0..d {
        alpha_pows.push(R::from(pow));
        pow *= alpha_base;
    }

    let mut lhs = R::ZERO;
    for inst_idx in 0..ell {
        let U: [Vec<R::BaseRing>; 3] = if let Some(auxw) = aux {
            if auxw.had_u.len() != ell {
                return Err("PiFold: aux.had_u length mismatch".to_string());
            }
            auxw.had_u[inst_idx].clone()
        } else {
            // Recompute y = (M1 f, M2 f, M3 f).
            let y = M
                .iter()
                .map(|Mi| Mi.try_mul_vec(&cms_openings[inst_idx]).expect("mat-vec mul failed"))
                .collect::<Vec<Vec<R>>>();

            let mut U: [Vec<R::BaseRing>; 3] =
                [Vec::with_capacity(d), Vec::with_capacity(d), Vec::with_capacity(d)];
            for i in 0..3 {
                for j in 0..d {
                    let evals = (0..m)
                        .map(|row| R::from(y[i][row].coeffs()[j]))
                        .collect::<Vec<_>>();
                    let mle = DenseMultilinearExtension::from_evaluations_vec(log_m, evals);
                    let v = mle
                        .evaluate(&r_poly_had)
                        .expect("MLE evaluate returned None");
                    U[i].push(v.ct());
                }
            }
            U
        };

        if U[0].len() != d || U[1].len() != d || U[2].len() != d {
            return Err("PiFold: aux had_U has wrong dimension".to_string());
        }

        let mut acc = R::ZERO;
        for j in 0..d {
            let u1 = R::from(U[0][j]);
            let u2 = R::from(U[1][j]);
            let u3 = R::from(U[2][j]);
            acc += alpha_pows[j] * (u1 * u2 - u3);
        }
        lhs += rhos[inst_idx] * (eq_sr * acc);

        // absorb U (prover did this after sumcheck)
        for x in &U[0] { transcript.absorb_field_element(x); }
        for x in &U[1] { transcript.absorb_field_element(x); }
        for x in &U[2] { transcript.absorb_field_element(x); }
    }

    if lhs != had_sc.expected_evaluation {
        return Err("PiFold: batched had Eq(26) mismatch".to_string());
    }

    // Per-instance monomial `b` values at the mon point, and absorb them (prover does this after sumchecks).
    // If `aux` is present, use it and avoid recomputing from witness openings.
    let r_mon_r: Vec<R> = mon_sc.point.iter().copied().map(R::from).collect();
    let mon_b: Vec<Vec<R>> = if let Some(auxw) = aux {
        if auxw.mon_b.len() != ell {
            return Err("PiFold: aux.mon_b length mismatch".to_string());
        }
        for b_inst in auxw.mon_b.iter() {
            if b_inst.len() != rg_params.k_g {
                return Err("PiFold: aux.mon_b wrong k_g".to_string());
            }
        }
        auxw.mon_b.clone()
    } else {
        let mut mon_b: Vec<Vec<R>> = Vec::with_capacity(ell);
        let expand_row = |row: usize| -> usize { row % m_j };
        for inst_idx in 0..ell {
            // Build H_digits for this instance (on m_J rows).
            let f = &cms_openings[inst_idx];
            let mut H = vec![vec![R::BaseRing::ZERO; d]; m_j];
            let Jref = &Js[inst_idx];
            for b in 0..blocks {
                for i in 0..rg_params.lambda_pj {
                    let out_row = b * rg_params.lambda_pj + i;
                    for t in 0..rg_params.l_h {
                        let in_row = b * rg_params.l_h + t;
                        let coef = Jref[i][t];
                        for col in 0..d {
                            H[out_row][col] += coef * f[in_row].coeffs()[col];
                        }
                    }
                }
            }
            let mut H_digits: Vec<Vec<Vec<R::BaseRing>>> =
                vec![vec![vec![R::BaseRing::ZERO; d]; m_j]; rg_params.k_g];
            for r in 0..m_j {
                let row_digits = H[r].decompose_to_vec(rg_params.d_prime, rg_params.k_g);
                for c in 0..d {
                    for i in 0..rg_params.k_g {
                        H_digits[i][r][c] = row_digits[c][i];
                    }
                }
            }

            let mut b_inst = Vec::with_capacity(rg_params.k_g);
            for dig in 0..rg_params.k_g {
                let mut gi = Vec::with_capacity(g_len);
                for c in 0..d {
                    for r in 0..m {
                        gi.push(exp::<R>(H_digits[dig][expand_row(r)][c]).expect("Exp failed"));
                    }
                }
                let mle = DenseMultilinearExtension::from_evaluations_vec(g_nvars, gi);
                b_inst.push(mle.evaluate(&r_mon_r).unwrap());
            }
            mon_b.push(b_inst);
        }
        mon_b
    };

    // Aggregate mon_b absorption: bind all mon_b values into transcript via a single
    // short Ajtai commitment (reduces Poseidon permutations from O(ell*k_g) to O(kappa)).
    let mon_b_agg = compute_mon_b_aggregate(&mon_b, kappa)?;
    transcript.absorb_slice(&mon_b_agg);

    let v_expected = mon_sc.expected_evaluation;

    let mut ver = R::ZERO;
    for inst_idx in 0..ell {
        let mut inst_acc = R::ZERO;
        // Precompute rc^dig iteratively to avoid repeated pow().
        let mut rc_pow = R::BaseRing::ONE;
        for dig in 0..rg_params.k_g {
            let (c, beta_i, alpha_i) = &cba_all[inst_idx][dig];
            let eq = eq_eval(c, &r_mon_r).unwrap();
            let b_i = mon_b[inst_idx][dig];
            let ev1 = R::from(ev(&b_i, *beta_i));
            let ev2 = R::from(ev(&b_i, *beta_i * *beta_i));
            let b_claim = ev1 * ev1 - ev2;
            let mut term = eq * R::from(*alpha_i) * b_claim;
            if let Some(rc) = &rc_all[inst_idx] {
                term *= R::from(rc_pow);
                rc_pow *= *rc;
            }
            inst_acc += term;
        }
        ver += rhos[inst_idx] * inst_acc;
    }

    if ver != v_expected {
        return Err("PiFold: batched monomial recomputation mismatch".to_string());
    }

    // Π_rg Step-5 check on the **folded** values.
    //
    // With β restricted to base-field scalars (embedded as constants), this is sound and
    // matches the folding algebra used by Figure 4.
    let log_m = log2(m.next_power_of_two()) as usize;
    let s_chals = mon_sc.point[log_m..].to_vec();
    let ts_s_full = ts_weights(&s_chals);
    let ts_s = &ts_s_full[..d];

    // Fold u*: for each digit, u*(dig) := Σ_i beta_i * u_i(dig)
    let mut u_folded = vec![R::ZERO; rg_params.k_g];
    for inst_idx in 0..ell {
        let b = beta_ring[inst_idx];
        for dig in 0..rg_params.k_g {
            u_folded[dig] += b * mon_b[inst_idx][dig];
        }
    }

    for dig in 0..rg_params.k_g {
        let lhs = (psi::<R>() * u_folded[dig]).ct();
        let rhs = proof.v_digits_folded[dig]
                .iter()
                .zip(ts_s.iter())
                .fold(R::BaseRing::ZERO, |acc, (&vij, &t)| acc + vij * t);
            if lhs != rhs {
            return Err(format!("PiFold: folded Step5 mismatch at dig={dig}"));
        }
    }

    // Construct folded outputs directly.
    // c* = Σ beta_i * c_i  (beta_i are constants, but we still use ring arithmetic)
    let mut c_star = vec![R::ZERO; cms[0].len()];
    for (inst_idx, cm) in cms.iter().enumerate() {
        for (acc, x) in c_star.iter_mut().zip(cm.iter()) {
            *acc += beta_ring[inst_idx] * *x;
        }
    }

    // r and r' are shared across the batch.
    let log_mj = log2(m_j.next_power_of_two()) as usize;
    let r_prime = mon_sc.point.clone();
    let r_star = r_prime[..log_mj].to_vec();

    // u* is the folded u_folded computed above; v* is composed from folded v_digits.
    let v_star = compose_v_digits::<R>(&proof.v_digits_folded, rg_params.d_prime);

    let mut v_rq = R::ZERO;
    for (i, c) in v_star.iter().enumerate() {
        v_rq.coeffs_mut()[i] = *c;
    }

    let folded_inst = SymphonyInstance { c: c_star, r: r_star, v: v_rq };
    // Fold commitments c*(i) := Σ beta_j * c_j(i).
    let kappa = proof.cm_g[0][0].len();
    let mut c_g_folded = vec![vec![R::ZERO; kappa]; rg_params.k_g];
    for inst_idx in 0..ell {
        let b = beta_ring[inst_idx];
        for dig in 0..rg_params.k_g {
            for j in 0..kappa {
                c_g_folded[dig][j] += b * proof.cm_g[inst_idx][dig][j];
            }
        }
    }
    let folded_bat = SymphonyBatchLin { r_prime, c_g: c_g_folded, u: u_folded };

    // -----------------------------------------------------------------------
    // Batchlin PCS transcript splice (domain-separated, REQUIRED). See the hetero-M variant.
    // -----------------------------------------------------------------------
    // Batchlin PCS transcript splice (REQUIRED). See hetero-M variant for rationale.
    if proof.batchlin_pcs_t.is_empty() {
        return Err("PiFold: missing batchlin_pcs_t (batchlin PCS commitment surface)".to_string());
    }
    if proof.batchlin_pcs_t.len() != 1 {
        return Err("PiFold: batchlin_pcs_t expected len=1 (batched scalar PCS)".to_string());
    }
    transcript.absorb_field_element(&R::BaseRing::from(BATCHLIN_PCS_DOMAIN_SEP));
    transcript.get_challenge();
    transcript.absorb_field_element(&R::BaseRing::from(proof.batchlin_pcs_t.len() as u128));
    for dig in 0..proof.batchlin_pcs_t.len() {
        if proof.batchlin_pcs_t[dig].is_empty() {
            return Err(format!("PiFold: empty batchlin_pcs_t[{dig}]"));
        }
        transcript.absorb_field_element(&R::BaseRing::from(proof.batchlin_pcs_t[dig].len() as u128));
        for x in &proof.batchlin_pcs_t[dig] {
            transcript.absorb_field_element(x);
        }
    }
    let _ = transcript.squeeze_bytes(64);

    Ok((folded_inst, folded_bat))
}
