//! WE/DPP gate arithmetization helpers (prime-field sparse dR1CS).
//!
//! This module is the *reusable wrapper* around the individual arithmetizers:
//! - Poseidon transcript trace -> sparse dR1CS (`dpp_poseidon`)
//! - AjtaiOpen(commitment, message) -> sparse dR1CS (`dpp_ajtai`)
//! - PCS evaluation verification -> sparse dR1CS (`pcs::dpp_folding_pcs_l2`)
//! - Merging multiple sparse dR1CS instances into one (`merge_sparse_dr1cs_share_one`)
//!
//! ## Architecture
//!
//! The WE gate (`R_WE`) is a single-phase design with the following components:
//!
//! - **Poseidon-FS transcript binding**: Ensures all FS challenges are bound to the protocol transcript.
//! - **Ajtai opening checks**: Verifies CP transcript message commitments (`cfs_*`).
//! - **Π_fold verifier math**: Sumcheck verification, monomial checks, and folding equations.
//! - **PCS evaluation proofs**: Binding aux values to the committed witness via lattice-based PCS.
//!
//! With PCS, there is no separate "π_lin" phase. The evaluation claims `aux = f(r)` are proven
//! directly by the PCS within the single WE gate verification.

use ark_crypto_primitives::sponge::poseidon::PoseidonConfig;
use ark_ff::{Field, PrimeField};
use ark_std::log2;
use latticefold::commitment::AjtaiCommitmentScheme;
use latticefold::transcript::Transcript;
use latticefold::utils::sumcheck::MLSumcheck;
use rayon::prelude::*;
use stark_rings::{balanced_decomposition::Decompose, OverField, PolyRing, Ring, Zq};

use crate::dpp_ajtai::{ajtai_open_dr1cs_from_scheme, ajtai_open_dr1cs_from_scheme_full};
use crate::dpp_poseidon::{
    merge_sparse_dr1cs_share_one,
    merge_sparse_dr1cs_share_one_with_glue,
    poseidon_sponge_dr1cs_from_trace,
    poseidon_sponge_dr1cs_from_trace_with_wiring_and_bytes,
    SparseDr1csInstance,
};
use crate::dpp_pifold_math::pifold_verifier_math_dr1cs;
use crate::pcs::dpp_folding_pcs_l2::folding_pcs_l2_verify_dr1cs_with_c_bytes;
use crate::pcs::folding_pcs_l2::{FoldingPcsL2Params, FoldingPcsL2ProofCore};
use crate::public_coin_transcript::FixedTranscript;
use crate::symphony_coins::{derive_beta_chi, derive_J};
use crate::symphony_pifold_batched::{PiFoldAuxWitness, PiFoldBatchedProof};
use crate::symphony_pifold_streaming::{CM_G_AGG_SEED, MON_B_AGG_SEED, V_DIGITS_AGG_SEED};
use crate::pcs::batchlin_pcs::BATCHLIN_PCS_DOMAIN_SEP;
use crate::pcs::cmf_pcs::{cmf_pcs_coin_bytes_len, cmf_pcs_params_for_flat_len, CMF_PCS_DOMAIN_SEP};
use crate::transcript::PoseidonTraceOp;

pub struct WeGateDr1csBuilder;

fn to_bf<R: OverField>(x: R::BaseRing) -> <R::BaseRing as Field>::BasePrimeField
where
    R::BaseRing: Field,
{
    x.to_base_prime_field_elements()
        .into_iter()
        .next()
        .expect("to_bf expects extension degree 1")
}

impl WeGateDr1csBuilder {
    fn squeeze_bytes_outputs<F: PrimeField>(ops: &[PoseidonTraceOp<F>]) -> Vec<Vec<u8>> {
        ops.iter()
            .filter_map(|op| match op {
                PoseidonTraceOp::SqueezeBytes { out, .. } => Some(out.clone()),
                _ => None,
            })
            .collect()
    }

    fn build_r_cp_poseidon_pifold_math_and_cfs_parts<R>(
        poseidon_cfg: &PoseidonConfig<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
        ops: &[PoseidonTraceOp<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>],
        cms: &[Vec<R>],
        proof: &PiFoldBatchedProof<R>,
        scheme_had: &AjtaiCommitmentScheme<R>,
        scheme_mon: &AjtaiCommitmentScheme<R>,
        aux: &PiFoldAuxWitness<R>,
        cfs_had_u: &[Vec<R>],
        cfs_mon_b: &[Vec<R>],
    ) -> Result<
        (
            Vec<(
            SparseDr1csInstance<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
            Vec<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
            )>,
            Vec<(usize, usize, usize, usize)>,
            crate::dpp_poseidon::PoseidonByteWiring,
            crate::dpp_poseidon::PoseidonDr1csWiring,
            crate::dpp_pifold_math::PiFoldMathWiring,
        ),
        String,
    >
    where
        R: OverField + Ring + PolyRing,
        R::BaseRing: Zq + Field + Decompose + PrimeField,
    {
        if cms.len() != aux.had_u.len()
            || cms.len() != aux.mon_b.len()
            || cms.len() != cfs_had_u.len()
            || cms.len() != cfs_mon_b.len()
        {
            return Err("WeGateDr1csBuilder: cms/aux/cfs length mismatch".to_string());
        }
        let ell = cms.len();
        if ell == 0 {
            return Err("WeGateDr1csBuilder: empty batch".to_string());
        }

        // Poseidon trace -> dR1CS (+ wiring with squeeze-field + squeeze-byte var indices).
        let (poseidon_inst, poseidon_asg, _replay2, _byte_wit, wiring, byte_wiring) =
            poseidon_sponge_dr1cs_from_trace_with_wiring_and_bytes::<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>(
                poseidon_cfg,
                ops,
            )
            .map_err(|e| format!("poseidon_sponge_dr1cs_from_trace_with_wiring_and_bytes: {e}"))?;

        // Ajtai-open checks for cfs_* - built in parallel per instance.
        let cfs_parts_nested: Vec<Vec<_>> = (0..ell)
            .into_par_iter()
            .map(|i| {
                let mut had_msg: Vec<R> = Vec::with_capacity(3 * R::dimension());
                for blk in 0..3 {
                    for &x in &aux.had_u[i][blk] {
                        had_msg.push(R::from(x));
                    }
                }
                let (inst_had, asg_had) =
                    ajtai_open_dr1cs_from_scheme::<R>(scheme_had, &had_msg, &cfs_had_u[i])
                        .expect("ajtai_open_dr1cs_from_scheme failed");
                let (inst_mon, asg_mon) =
                    ajtai_open_dr1cs_from_scheme_full::<R>(scheme_mon, &aux.mon_b[i], &cfs_mon_b[i])
                        .expect("ajtai_open_dr1cs_from_scheme_full failed");
                vec![(inst_had, asg_had), (inst_mon, asg_mon)]
            })
            .collect();

        let mut cfs_parts: Vec<(
            SparseDr1csInstance<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
            Vec<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
        )> = Vec::with_capacity(2 * ell);
        for inst_part in cfs_parts_nested {
            cfs_parts.extend(inst_part);
        }
        let (cfs_inst, cfs_asg) = merge_sparse_dr1cs_share_one(&cfs_parts)?;

        let mut parts: Vec<(
            SparseDr1csInstance<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
            Vec<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
        )> = Vec::new();
        parts.push((poseidon_inst, poseidon_asg));
        parts.push((cfs_inst, cfs_asg));

        // ---------------------------------------------------------------------
        // Ajtai-open check for cm_g aggregate: enforce c_agg = A_agg · cm_g_flat.
        // This binds the aggregate (absorbed into transcript) to the individual cm_g values.
        // ---------------------------------------------------------------------
        let kappa = if !proof.cm_g.is_empty() && !proof.cm_g[0].is_empty() {
            proof.cm_g[0][0].len()
        } else {
            return Err("WeGateDr1csBuilder: empty cm_g for aggregate".to_string());
        };

        // Flatten cm_g: concat_{inst, dig, j} cm_g[inst][dig][j]
        let mut cm_g_flat: Vec<R> = Vec::new();
        for inst in proof.cm_g.iter() {
            for dig in inst.iter() {
                cm_g_flat.extend(dig.iter().copied());
            }
        }
        let n_agg = cm_g_flat.len();

        // Create the aggregate scheme (same as used in compute_cm_g_aggregate).
        let cm_g_agg_scheme =
            AjtaiCommitmentScheme::<R>::seeded(b"cm_g_agg", CM_G_AGG_SEED, kappa, n_agg);

        // Compute the expected aggregate (verifier-side recomputation).
        let cm_g_agg = cm_g_agg_scheme
            .commit(&cm_g_flat)
            .map_err(|e| format!("WeGateDr1csBuilder: cm_g_agg commit failed: {e:?}"))?
            .as_ref()
            .to_vec();

        // Enforce: A_agg · cm_g_flat = cm_g_agg (linear constraints).
        let (cm_g_agg_inst, cm_g_agg_asg) =
            ajtai_open_dr1cs_from_scheme_full::<R>(&cm_g_agg_scheme, &cm_g_flat, &cm_g_agg)?;
        parts.push((cm_g_agg_inst, cm_g_agg_asg));

        // Flatten mon_b: concat_{inst, dig} mon_b[inst][dig]
        let mut mon_b_flat: Vec<R> = Vec::new();
        for inst in aux.mon_b.iter() {
            mon_b_flat.extend(inst.iter().copied());
        }
        let n_mon_b_agg = mon_b_flat.len();

        // Create the aggregate scheme for mon_b (same as used in compute_mon_b_aggregate).
        let mon_b_agg_scheme =
            AjtaiCommitmentScheme::<R>::seeded(b"mon_b_agg", MON_B_AGG_SEED, kappa, n_mon_b_agg);

        // Compute the expected aggregate (verifier-side recomputation).
        let mon_b_agg = mon_b_agg_scheme
            .commit(&mon_b_flat)
            .map_err(|e| format!("WeGateDr1csBuilder: mon_b_agg commit failed: {e:?}"))?
            .as_ref()
            .to_vec();

        // Enforce: A_agg · mon_b_flat = mon_b_agg (linear constraints).
        let (mon_b_agg_inst, mon_b_agg_asg) =
            ajtai_open_dr1cs_from_scheme_full::<R>(&mon_b_agg_scheme, &mon_b_flat, &mon_b_agg)?;
        parts.push((mon_b_agg_inst, mon_b_agg_asg));

        // ---------------------------------------------------------------------
        // Ajtai-open check for v_digits_folded aggregate: enforce c_agg = A_agg · v_digits_flat.
        // This binds the (compressed) hook message absorbed into transcript to the actual
        // v_digits_folded values used by Π_fold math.
        // ---------------------------------------------------------------------
        let k_g = proof.rg_params.k_g;
        if proof.v_digits_folded.len() != k_g {
            return Err(format!(
                "WeGateDr1csBuilder: v_digits_folded len mismatch (got {}, expected k_g={k_g})",
                proof.v_digits_folded.len()
            ));
        }
        let d = R::dimension();
        for (dig, v) in proof.v_digits_folded.iter().enumerate() {
            if v.len() != d {
                return Err(format!(
                    "WeGateDr1csBuilder: v_digits_folded[{dig}] len mismatch (got {}, expected d={d})",
                    v.len()
                ));
            }
        }

        // Flatten v_digits_folded as constant-coefficient ring elements.
        let mut v_digits_flat: Vec<R> = Vec::with_capacity(k_g * d);
        for v_i in &proof.v_digits_folded {
            for &x in v_i {
                v_digits_flat.push(R::from(x));
            }
        }
        let n_v_digits_agg = v_digits_flat.len();
        let v_digits_agg_scheme = AjtaiCommitmentScheme::<R>::seeded(
            b"v_digits_agg",
            V_DIGITS_AGG_SEED,
            kappa,
            n_v_digits_agg,
        );
        let v_digits_agg = v_digits_agg_scheme
            .commit_const_coeff_fast(&v_digits_flat)
            .map_err(|e| format!("WeGateDr1csBuilder: v_digits_agg commit failed: {e:?}"))?
            .as_ref()
            .to_vec();
        let (v_digits_agg_inst, v_digits_agg_asg) = ajtai_open_dr1cs_from_scheme_full::<R>(
            &v_digits_agg_scheme,
            &v_digits_flat,
            &v_digits_agg,
        )?;
        parts.push((v_digits_agg_inst, v_digits_agg_asg));

        // ---------------------------------------------------------------------
        // Extract the Π_fold verifier coin pieces from the proof's recorded coin stream
        // and compute rs_shared via sumcheck verification replay.
        // ---------------------------------------------------------------------
        let mut ts = FixedTranscript::<R>::new_with_coins_and_events(
            proof.coins.challenges.clone(),
            proof.coins.bytes.clone(),
            proof.coins.events.clone(),
        );

        // NOTE: FixedTranscript absorb is a no-op; this is purely to replay the same coin schedule.
        let beta_cts = derive_beta_chi::<R>(&mut ts, ell);
        let beta_cts_bf = beta_cts.iter().copied().map(to_bf::<R>).collect::<Vec<_>>();

        let m = proof.m;
        let log_m = log2(m) as usize;
        ts.absorb_field_element(&R::BaseRing::from(0x4c465053_50494841u128)); // "LFPS_PIHA"
        let s_base = ts.get_challenges(log_m);
        let s_base_bf = s_base.iter().copied().map(to_bf::<R>).collect::<Vec<_>>();
        let alpha_base_bf = to_bf::<R>(ts.get_challenge());

        let d = R::dimension();
        let g_len = m * d;
        let g_nvars = log2(g_len) as usize;
        // k_g already checked above for v_digits_folded aggregation.

        // Reconstruct witness length `n_f` from `m_j = (n_f/l_h) * lambda_pj`.
        // This is needed to deterministically parameterize the cm_f PCS coin splice.
        let lambda_pj = proof.rg_params.lambda_pj;
        let l_h = proof.rg_params.l_h;
        if lambda_pj == 0 || l_h == 0 {
            return Err(
                "WeGateDr1csBuilder: invalid rg params (lambda_pj/l_h must be > 0)".to_string(),
            );
        }
        if proof.m_j % lambda_pj != 0 {
            return Err("WeGateDr1csBuilder: m_j not divisible by lambda_pj".to_string());
        }
        let blocks = proof.m_j / lambda_pj;
        let n_f = blocks * l_h;

        let mut cba_all_bf: Vec<Vec<(Vec<<R::BaseRing as Field>::BasePrimeField>, _, _)>> =
            Vec::with_capacity(ell);
        let mut rc_all_bf: Vec<Option<<R::BaseRing as Field>::BasePrimeField>> =
            Vec::with_capacity(ell);

        // Phase 1: absorb cm_f, do the cm_f PCS coin splice, then derive J for each instance.
        for cm_f in cms {
            ts.absorb_slice(cm_f);
            // PCS#1 (cm_f PCS) coin splice: bind to cm_f surface and derive C1/C2 bytes.
            // Deterministically parameterized by (n_f, cm_f surface length).
            ts.absorb_field_element(&R::BaseRing::from(CMF_PCS_DOMAIN_SEP));
            let flat_len = n_f * d;
            let kappa_commit = cm_f.len() * d;
            let pcs_params_cmf =
                cmf_pcs_params_for_flat_len::<R::BaseRing>(flat_len, kappa_commit)?;
            let n_bytes_cmf = cmf_pcs_coin_bytes_len(&pcs_params_cmf);
            let _ = ts.squeeze_bytes(n_bytes_cmf);
            let _j = derive_J::<R>(&mut ts, proof.rg_params.lambda_pj, proof.rg_params.l_h);
        }

        // Aggregate cm_g absorption (reuse cm_g_agg computed above for Ajtai check).
        ts.absorb_slice(&cm_g_agg);

        // Phase 2: derive coins for each instance.
        for _ in 0..ell {
            let mut cba = Vec::with_capacity(k_g);
            for _dig in 0..k_g {
                let c = ts
                    .get_challenges(g_nvars)
                    .into_iter()
                    .map(to_bf::<R>)
                    .collect::<Vec<_>>();
                let beta = to_bf::<R>(ts.get_challenge());
                let alpha = to_bf::<R>(ts.get_challenge());
                cba.push((c, beta, alpha));
            }
            let rc = (k_g > 1).then(|| to_bf::<R>(ts.get_challenge()));
            cba_all_bf.push(cba);
            rc_all_bf.push(rc);
        }
        let rhos_bf = ts
            .get_challenges(ell)
            .into_iter()
            .map(to_bf::<R>)
            .collect::<Vec<_>>();

        // Sumcheck randomness (shared): replay verifier schedule to get `mon_sc.point`.
        let hook_round = log2(proof.m_j.next_power_of_two()) as usize;
        let (_had_sc, mon_sc) = MLSumcheck::<R, _>::verify_two_as_subprotocol_shared_with_hook(
            &mut ts,
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
                for v_i in &proof.v_digits_folded {
                    for x in v_i {
                        t.absorb_field_element(x);
                    }
                }
            },
        )
        .map_err(|e| format!("WeGateDr1csBuilder: sumcheck verify replay failed: {e}"))?;

        let rs_shared_bf = mon_sc
            .point
            .iter()
            .copied()
            .map(to_bf::<R>)
            .collect::<Vec<_>>();

        // ---------------------------------------------------------------------
        // Batchlin PCS transcript splice (domain-separated, REQUIRED):
        // ---------------------------------------------------------------------
        if proof.batchlin_pcs_t.is_empty() {
            return Err("WeGateDr1csBuilder: missing batchlin_pcs_t".to_string());
        }
        if proof.batchlin_pcs_t.len() != 1 {
            return Err("WeGateDr1csBuilder: batchlin_pcs_t expected len=1 (batched scalar PCS)".to_string());
        }
        if proof.batchlin_pcs_t[0].is_empty() {
            return Err("WeGateDr1csBuilder: empty batchlin_pcs_t[0]".to_string());
        }
        ts.absorb_field_element(&R::BaseRing::from(BATCHLIN_PCS_DOMAIN_SEP));
        let _gamma = ts.get_challenge();
        ts.absorb_field_element(&R::BaseRing::from(proof.batchlin_pcs_t.len() as u128));
        for dig in 0..proof.batchlin_pcs_t.len() {
            ts.absorb_field_element(&R::BaseRing::from(proof.batchlin_pcs_t[dig].len() as u128));
            for x in &proof.batchlin_pcs_t[dig] {
                ts.absorb_field_element(x);
            }
        }
        let _ = ts.squeeze_bytes(64);

        // Ensure we consumed the entire provided coin stream.
        if ts.remaining_challenges() != 0 || ts.remaining_bytes() != 0 || ts.remaining_events() != 0 {
            return Err(format!(
                "WeGateDr1csBuilder: leftover coins after replay: challenges={} bytes={} events={}",
                ts.remaining_challenges(),
                ts.remaining_bytes(),
                ts.remaining_events()
            ));
        }

        // Build Π_fold verifier math dR1CS.
        let (pifold_inst, pifold_asg, pifold_wiring) = pifold_verifier_math_dr1cs::<R>(
            proof,
            aux,
            &beta_cts_bf,
            &s_base_bf,
            alpha_base_bf,
            &cba_all_bf,
            &rc_all_bf,
            &rhos_bf,
            &rs_shared_bf,
        )?;
        parts.push((pifold_inst, pifold_asg));

        // ---------------------------------------------------------------------
        // Glue Π_fold challenge variables to Poseidon squeeze-field variables.
        // ---------------------------------------------------------------------
        // Expected number of `get_challenge` outputs for Π_fold verifier schedule.
        //
        // Note: The transcript may contain **additional** challenges after Π_fold completes
        // (e.g. batchlin PCS batching γ). We only require Poseidon provided at least the prefix
        // needed by Π_fold.
        let per_inst = k_g * (g_nvars + 2) + if k_g > 1 { 1 } else { 0 };
        let total_challenges = log_m + 1 + ell * per_inst + ell + g_nvars;
        if wiring.squeeze_field_vars.len() < total_challenges {
            return Err(format!(
                "WeGateDr1csBuilder: poseidon wiring squeeze_field_vars too short: got {} need_at_least {}",
                wiring.squeeze_field_vars.len(),
                total_challenges
            ));
        }

        // Build glue list in local indices: (part_a, var_a, part_b, var_b).
        // part 0 = poseidon, part 1 = cfs-openings, part 2 = cm_g_agg, part 3 = mon_b_agg,
        // part 4 = v_digits_agg, part 5 = pifold-math
        const PIFOLD_PART: usize = 5;
        let mut glue: Vec<(usize, usize, usize, usize)> = Vec::new();

        let mut ch = 0usize;
        // s_base
        for &v in &pifold_wiring.s_base {
            glue.push((PIFOLD_PART, v, 0, wiring.squeeze_field_vars[ch]));
            ch += 1;
        }
        // alpha_base
        glue.push((PIFOLD_PART, pifold_wiring.alpha_base, 0, wiring.squeeze_field_vars[ch]));
        ch += 1;
        // per-instance per-digit: c vector, beta, alpha, optional rc
        for inst_idx in 0..ell {
            for dig in 0..k_g {
                let idx_flat = inst_idx * k_g + dig;
                for &cv in &pifold_wiring.c_all[idx_flat] {
                    glue.push((PIFOLD_PART, cv, 0, wiring.squeeze_field_vars[ch]));
                    ch += 1;
                }
                glue.push((
                    PIFOLD_PART,
                    pifold_wiring.beta_i_all[idx_flat],
                    0,
                    wiring.squeeze_field_vars[ch],
                ));
                ch += 1;
                glue.push((
                    PIFOLD_PART,
                    pifold_wiring.alpha_i_all[idx_flat],
                    0,
                    wiring.squeeze_field_vars[ch],
                ));
                ch += 1;
            }
            if let Some(rcv) = pifold_wiring.rc_all[inst_idx] {
                glue.push((PIFOLD_PART, rcv, 0, wiring.squeeze_field_vars[ch]));
                ch += 1;
            }
        }
        // rhos
        for &rv in &pifold_wiring.rhos {
            glue.push((PIFOLD_PART, rv, 0, wiring.squeeze_field_vars[ch]));
            ch += 1;
        }
        // rs_shared
        for &rv in &pifold_wiring.rs_shared {
            glue.push((PIFOLD_PART, rv, 0, wiring.squeeze_field_vars[ch]));
            ch += 1;
        }
        debug_assert_eq!(ch, total_challenges);

        Ok((parts, glue, byte_wiring, wiring, pifold_wiring))
    }

    /// Build the current arithmetized **R_cp** fragment:
    /// - Poseidon transcript trace constraints (FS-in-gate)
    /// - AjtaiOpen for CP transcript-message commitments (`cfs_had_u`, `cfs_mon_b`)
    ///
    /// This is the part that enforces transcript binding and that “aux matches its CP commitments”.
    pub fn r_cp_poseidon_and_cfs_openings<R>(
        poseidon_cfg: &PoseidonConfig<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
        ops: &[PoseidonTraceOp<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>],
        scheme_had: &AjtaiCommitmentScheme<R>,
        scheme_mon: &AjtaiCommitmentScheme<R>,
        aux: &PiFoldAuxWitness<R>,
        cfs_had_u: &[Vec<R>],
        cfs_mon_b: &[Vec<R>],
    ) -> Result<
        (
            SparseDr1csInstance<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
            Vec<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
        ),
        String,
    >
    where
        R: OverField + Ring + PolyRing,
        R::BaseRing: Zq + Field + PrimeField,
    {
        // Poseidon trace -> dR1CS.
        let (poseidon_inst, poseidon_asg, _replay2, _byte_wit) =
            poseidon_sponge_dr1cs_from_trace::<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>(
                poseidon_cfg,
                ops,
            )
            .map_err(|e| format!("poseidon_sponge_dr1cs_from_trace: {e}"))?;

        if aux.had_u.len() != cfs_had_u.len() || aux.mon_b.len() != cfs_mon_b.len() {
            return Err("WeGateDr1csBuilder: aux/cfs length mismatch".to_string());
        }
        let ell = cfs_had_u.len();
        if ell == 0 {
            return Err("WeGateDr1csBuilder: empty batch".to_string());
        }

        // Ajtai-open checks for cfs_* - built in parallel per instance.
        let instance_parts: Vec<_> = (0..ell)
            .into_par_iter()
            .map(|i| {
                let mut had_msg: Vec<R> = Vec::with_capacity(3 * R::dimension());
                for blk in 0..3 {
                    for &x in &aux.had_u[i][blk] {
                        had_msg.push(R::from(x));
                    }
                }
                let (inst_had, asg_had) =
                    ajtai_open_dr1cs_from_scheme::<R>(scheme_had, &had_msg, &cfs_had_u[i])
                        .expect("ajtai_open_dr1cs_from_scheme failed");
                let (inst_mon, asg_mon) =
                    ajtai_open_dr1cs_from_scheme_full::<R>(scheme_mon, &aux.mon_b[i], &cfs_mon_b[i])
                        .expect("ajtai_open_dr1cs_from_scheme_full failed");
                vec![(inst_had, asg_had), (inst_mon, asg_mon)]
            })
            .collect();

        let mut parts: Vec<(
            SparseDr1csInstance<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
            Vec<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
        )> = Vec::with_capacity(1 + 2 * ell);
        parts.push((poseidon_inst, poseidon_asg));
        for inst_part in instance_parts {
            parts.extend(inst_part);
        }

        merge_sparse_dr1cs_share_one(&parts)
    }

    /// Build the arithmetized **R_cp** fragment including:
    /// - Poseidon transcript trace constraints (FS-in-gate)
    /// - Π_fold verifier arithmetic constraints (sumcheck verify + Eq(26) + monomial recomputation + Step-5)
    /// - AjtaiOpen for CP transcript-message commitments (`cfs_had_u`, `cfs_mon_b`)
    ///
    /// Additionally, this **glues** Π_fold's challenge variables to Poseidon's `SqueezeField` variables
    /// (so the Π_fold math is bound to the Poseidon transcript).
    ///
    /// Note: `derive_beta_chi` and `derive_J` both use `SqueezeBytes`.
    /// We now arithmetize `SqueezeBytes` in `dpp_poseidon`, but Π_fold's β/J derivation is still
    /// treated as an external witness value here (until we add an in-circuit byte->field parsing layer).
    pub fn r_cp_poseidon_pifold_math_and_cfs_openings<R>(
        poseidon_cfg: &PoseidonConfig<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
        ops: &[PoseidonTraceOp<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>],
        cms: &[Vec<R>],
        proof: &PiFoldBatchedProof<R>,
        scheme_had: &AjtaiCommitmentScheme<R>,
        scheme_mon: &AjtaiCommitmentScheme<R>,
        aux: &PiFoldAuxWitness<R>,
        cfs_had_u: &[Vec<R>],
        cfs_mon_b: &[Vec<R>],
    ) -> Result<
        (
            SparseDr1csInstance<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
            Vec<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
        ),
        String,
    >
    where
        R: OverField + Ring + PolyRing,
        R::BaseRing: Zq + Field + Decompose + PrimeField,
    {
        let (parts, glue, _byte_wiring, _pose_wiring, _pifold_wiring) = Self::build_r_cp_poseidon_pifold_math_and_cfs_parts::<R>(
            poseidon_cfg,
            ops,
            cms,
            proof,
            scheme_had,
            scheme_mon,
            aux,
            cfs_had_u,
            cfs_mon_b,
        )?;
        merge_sparse_dr1cs_share_one_with_glue(&parts, &glue)
    }

    /// Build and merge:
    /// - Poseidon sponge dR1CS for `ops`
    /// - Ajtai-open dR1CS for all `cfs_had_u[i]` and `cfs_mon_b[i]` openings to `aux`
    ///
    /// Returns a single merged `(inst, assignment)` sharing a single constant-1 slot.
    pub fn poseidon_plus_cfs_openings<R>(
        poseidon_cfg: &PoseidonConfig<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
        ops: &[PoseidonTraceOp<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>],
        scheme_had: &AjtaiCommitmentScheme<R>,
        scheme_mon: &AjtaiCommitmentScheme<R>,
        aux: &PiFoldAuxWitness<R>,
        cfs_had_u: &[Vec<R>],
        cfs_mon_b: &[Vec<R>],
    ) -> Result<
        (
            SparseDr1csInstance<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
            Vec<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
        ),
        String,
    >
    where
        R: OverField + Ring,
        R: PolyRing,
        R::BaseRing: Zq + Field + PrimeField,
    {
        Self::r_cp_poseidon_and_cfs_openings::<R>(
            poseidon_cfg, ops, scheme_had, scheme_mon, aux, cfs_had_u, cfs_mon_b,
        )
    }

    /// Build **R_cp × R_pifold × R_cfs × R_pcs_f × R_pcs_g**:
    /// - transcript binding (Poseidon trace)
    /// - Π_fold verifier math
    /// - CP message openings (`cfs_*`)
    /// - two PCS evaluation proof verifications (folding PCS ℓ=2), intended for
    ///   `R_o = R_auxcs_lin × R_batchlin`, with `C1/C2` coins derived from Poseidon `SqueezeBytes`
    ///   and glued by variable-equality constraints inside the merged system.
    pub fn poseidon_plus_pifold_plus_cfs_plus_pcs<R>(
        poseidon_cfg: &PoseidonConfig<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
        ops: &[PoseidonTraceOp<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>],
        cms: &[Vec<R>],
        proof: &PiFoldBatchedProof<R>,
        scheme_had: &AjtaiCommitmentScheme<R>,
        scheme_mon: &AjtaiCommitmentScheme<R>,
        aux: &PiFoldAuxWitness<R>,
        cfs_had_u: &[Vec<R>],
        cfs_mon_b: &[Vec<R>],
        pcs_f_params: &FoldingPcsL2Params<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
        pcs_f_t: &[<<R as PolyRing>::BaseRing as Field>::BasePrimeField],
        pcs_f_x0: &[<<R as PolyRing>::BaseRing as Field>::BasePrimeField],
        pcs_f_x1: &[<<R as PolyRing>::BaseRing as Field>::BasePrimeField],
        pcs_f_x2: &[<<R as PolyRing>::BaseRing as Field>::BasePrimeField],
        pcs_f_proof: &FoldingPcsL2ProofCore<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
        pcs_f_coin_squeeze_idx: usize,
        pcs_g_params: &FoldingPcsL2Params<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
        pcs_g_t: &[<<R as PolyRing>::BaseRing as Field>::BasePrimeField],
        pcs_g_x0: &[<<R as PolyRing>::BaseRing as Field>::BasePrimeField],
        pcs_g_x1: &[<<R as PolyRing>::BaseRing as Field>::BasePrimeField],
        pcs_g_x2: &[<<R as PolyRing>::BaseRing as Field>::BasePrimeField],
        pcs_g_proof: &FoldingPcsL2ProofCore<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
        pcs_g_coin_squeeze_idx: usize,
    ) -> Result<
        (
            SparseDr1csInstance<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
            Vec<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
        ),
        String,
    >
    where
        R: OverField + Ring + PolyRing,
        R::BaseRing: Zq + Field + Decompose + PrimeField,
    {
        let (mut parts, mut glue, byte_wiring, pose_wiring, pifold_wiring) = Self::build_r_cp_poseidon_pifold_math_and_cfs_parts::<R>(
            poseidon_cfg,
            ops,
            cms,
            proof,
            scheme_had,
            scheme_mon,
            aux,
            cfs_had_u,
            cfs_mon_b,
        )?;

        // ---------------------------------------------------------------------
        // PCS verifier dR1CS (folding PCS ℓ=2), with `C1/C2` derived from Poseidon `SqueezeBytes`.
        // We wire TWO PCS verifiers (auxcs_lin + batchlin).
        // ---------------------------------------------------------------------
        let squeeze_outs = Self::squeeze_bytes_outputs(ops);

        let mut add_pcs_part = |label: &str,
                                pcs_params: &FoldingPcsL2Params<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
                                pcs_t: &[<<R as PolyRing>::BaseRing as Field>::BasePrimeField],
                                pcs_x0: &[<<R as PolyRing>::BaseRing as Field>::BasePrimeField],
                                pcs_x1: &[<<R as PolyRing>::BaseRing as Field>::BasePrimeField],
                                pcs_x2: &[<<R as PolyRing>::BaseRing as Field>::BasePrimeField],
                                pcs_proof: &FoldingPcsL2ProofCore<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
                                pcs_coin_squeeze_idx: usize|
         -> Result<(usize, crate::pcs::dpp_folding_pcs_l2::FoldingPcsL2Wiring), String> {
            if pcs_coin_squeeze_idx >= squeeze_outs.len() {
                return Err(format!("{}: pcs_coin_squeeze_idx out of range: idx={} num_squeezes={}",
                    label, pcs_coin_squeeze_idx, squeeze_outs.len()));
            }
            if pcs_coin_squeeze_idx >= byte_wiring.squeeze_byte_ranges.len() {
                return Err(format!(
                    "{}: Poseidon byte wiring missing squeeze idx {} (have {})",
                    label,
                    pcs_coin_squeeze_idx,
                    byte_wiring.squeeze_byte_ranges.len()
                ));
            }
            let c_bytes = &squeeze_outs[pcs_coin_squeeze_idx];
            let (byte_start, byte_len) = byte_wiring.squeeze_byte_ranges[pcs_coin_squeeze_idx];
            if byte_len != c_bytes.len() {
                return Err(format!(
                    "{}: pcs coin byte length mismatch: poseidon_wiring_len={} ops_len={}",
                    label,
                    byte_len,
                    c_bytes.len()
                ));
            }
            let pose_byte_vars =
                &byte_wiring.squeeze_byte_vars[byte_start..byte_start + byte_len];

            let mut b = crate::dpp_sumcheck::Dr1csBuilder::<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>::new();
            let pcs_byte_vars: Vec<usize> = c_bytes
                .iter()
                .map(|&by| {
                    b.new_var(
                        <<<R as PolyRing>::BaseRing as Field>::BasePrimeField as From<u64>>::from(by as u64),
                    )
                })
                .collect();
            let pcs_wiring = folding_pcs_l2_verify_dr1cs_with_c_bytes(
                &mut b,
                pcs_params,
                pcs_t,
                pcs_x0,
                pcs_x1,
                pcs_x2,
                pcs_proof,
                &pcs_byte_vars,
            )?;
            let (pcs_inst, pcs_asg) = b.into_instance();
            let pcs_part_idx = parts.len();
            parts.push((pcs_inst, pcs_asg));

            for (&pv, &cv) in pose_byte_vars.iter().zip(pcs_byte_vars.iter()) {
                glue.push((0, pv, pcs_part_idx, cv));
            }
            Ok((pcs_part_idx, pcs_wiring))
        };

        let (pcs_f_part_idx, pcs_f_wiring) = add_pcs_part(
            "cm_f",
            pcs_f_params,
            pcs_f_t,
            pcs_f_x0,
            pcs_f_x1,
            pcs_f_x2,
            pcs_f_proof,
            pcs_f_coin_squeeze_idx,
        )?;
        let (pcs_g_part_idx, pcs_g_wiring) = add_pcs_part(
            "batchlin",
            pcs_g_params,
            pcs_g_t,
            pcs_g_x0,
            pcs_g_x1,
            pcs_g_x2,
            pcs_g_proof,
            pcs_g_coin_squeeze_idx,
        )?;

        // ---------------------------------------------------------------------
        // PCS_f commitment binding:
        //
        // `cm_f` is absorbed as ring elements, followed immediately by `CMF_PCS_DOMAIN_SEP`.
        //
        // We locate the **contiguous block** of `d`-sized Absorb ops immediately preceding the
        // first `Absorb([CMF_PCS_DOMAIN_SEP])`, and glue the first `t_len` base-field elements from
        // that block to the PCS_f commitment variables.
        // ---------------------------------------------------------------------
        let d = R::dimension();
        let marker_cmf = <<<R as PolyRing>::BaseRing as Field>::BasePrimeField as From<u128>>::from(
            CMF_PCS_DOMAIN_SEP,
        );
        let mut absorb_op_idx = 0usize;
        let mut cmf_block_vars: Vec<usize> = Vec::new();
        let mut found_marker = false;
        for op in ops {
            if let PoseidonTraceOp::Absorb(elems) = op {
                let (start, len) = pose_wiring
                    .absorb_ranges
                    .get(absorb_op_idx)
                    .copied()
                    .unwrap_or((0, 0));
                absorb_op_idx += 1;
                if len != elems.len() {
                    return Err("WeGateDr1csBuilder: poseidon absorb wiring length mismatch".to_string());
                }
                let vars = &pose_wiring.absorb_vars[start..start + len];
                if elems.len() == 1 && elems[0] == marker_cmf {
                    found_marker = true;
                    break;
                }
                if elems.len() == d {
                    cmf_block_vars.extend_from_slice(vars);
                } else {
                    // Break contiguity: the cm_f block must be immediately before the marker.
                    cmf_block_vars.clear();
                }
            }
        }
        if !found_marker {
            return Err("WeGateDr1csBuilder: failed to locate CMF_PCS_DOMAIN_SEP absorb op".to_string());
        }
        if cmf_block_vars.len() < pcs_f_wiring.t_vars.len() {
            return Err("WeGateDr1csBuilder: failed to locate enough cm_f absorb vars to glue pcs_f_t".to_string());
        }
        for (&pv, &tv) in cmf_block_vars[..pcs_f_wiring.t_vars.len()]
            .iter()
            .zip(pcs_f_wiring.t_vars.iter())
        {
            glue.push((0, pv, pcs_f_part_idx, tv));
        }

        // ---------------------------------------------------------------------
        // Batchlin PCS binding:
        //
        // Bind PCS_g to Π_fold Step-5 scalars and transcript absorption.
        // - Glue PCS_g commitment vector `t_vars` to the Poseidon `Absorb` vars for `batchlin_pcs_t`.
        // - Glue PCS_g computed evaluation `u_re_vars[0]` to a batched linear combination of
        //   Π_fold Step-5 LHS scalars `lhs_ct[dig] = ct(psi*u_folded[dig])` with batching scalar γ
        //   derived from the transcript after the batchlin PCS splice.
        // ---------------------------------------------------------------------
        // Locate the Poseidon absorb vars for `batchlin_pcs_t` by replaying the absorb element stream:
        // after the marker, skip 2 length elements, then take `t_len` elements.
        let marker_bf = <<<R as PolyRing>::BaseRing as Field>::BasePrimeField as From<u128>>::from(BATCHLIN_PCS_DOMAIN_SEP);
        let mut absorb_op_idx = 0usize;
        let mut after_marker = false;
        let mut skip_after_marker = 0usize;
        let mut t_absorb_vars: Vec<usize> = Vec::new();
        for op in ops {
            if let PoseidonTraceOp::Absorb(elems) = op {
                let (start, len) = pose_wiring
                    .absorb_ranges
                    .get(absorb_op_idx)
                    .copied()
                    .unwrap_or((0, 0));
                absorb_op_idx += 1;
                if len != elems.len() {
                    return Err("WeGateDr1csBuilder: poseidon absorb wiring length mismatch".to_string());
                }
                let vars = &pose_wiring.absorb_vars[start..start + len];

                // Marker absorb is always a single field element.
                if elems.len() == 1 && elems[0] == marker_bf {
                    after_marker = true;
                    // After the marker, we perform:
                    // - get_challenge(): SqueezeField(gamma) then Absorb(gamma)
                    // - absorb two length tags (outer len, inner len)
                    // before the t elements.
                    skip_after_marker = 3;
                    continue;
                }
                if after_marker {
                    for (e_idx, _e) in elems.iter().enumerate() {
                        if skip_after_marker > 0 {
                            skip_after_marker -= 1;
                            continue;
                        }
                        t_absorb_vars.push(vars[e_idx]);
                        if t_absorb_vars.len() == pcs_g_wiring.t_vars.len() {
                            break;
                        }
                    }
                    if t_absorb_vars.len() == pcs_g_wiring.t_vars.len() {
                        break;
                    }
                }
            }
        }
        if t_absorb_vars.len() != pcs_g_wiring.t_vars.len() {
            return Err("WeGateDr1csBuilder: failed to locate batchlin_pcs_t absorb vars for glue".to_string());
        }
        for (&pv, &tv) in t_absorb_vars.iter().zip(pcs_g_wiring.t_vars.iter()) {
            glue.push((0, pv, pcs_g_part_idx, tv));
        }

        // Glue gamma_local to the first SqueezeField output after the batchlin PCS marker absorb.
        let mut squeeze_idx = 0usize;
        let mut saw_marker2 = false;
        let mut gamma_pose_var: Option<usize> = None;
        for op in ops {
            match op {
                PoseidonTraceOp::Absorb(elems) => {
                    if elems.len() == 1 && elems[0] == marker_bf {
                        saw_marker2 = true;
                    }
                }
                PoseidonTraceOp::SqueezeField(out) => {
                    if saw_marker2 && !out.is_empty() {
                        gamma_pose_var = Some(pose_wiring.squeeze_field_vars[squeeze_idx]);
                        break;
                    }
                    squeeze_idx += out.len();
                }
                _ => {}
            }
        }
        let gamma_pose_var = gamma_pose_var.ok_or_else(|| {
            "WeGateDr1csBuilder: missing gamma SqueezeField after batchlin PCS splice".to_string()
        })?;
        let gamma_val = parts[0].1[gamma_pose_var];

        // Build a tiny binding sub-system: u_re == Σ_dig γ^dig * lhs_ct[dig]
        let mut bb = crate::dpp_sumcheck::Dr1csBuilder::<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>::new();
        let bf_zero = <<<R as PolyRing>::BaseRing as Field>::BasePrimeField as ark_ff::Field>::ZERO;
        let bf_one = <<<R as PolyRing>::BaseRing as Field>::BasePrimeField as ark_ff::Field>::ONE;
        // Allocate local vars with correct initial values so glue constraints are satisfied.
        let gamma_local = bb.new_var(gamma_val);
        let pcs_g_asg = &parts[pcs_g_part_idx].1;
        let u_val = pcs_g_asg[pcs_g_wiring.u_re_vars[0]];
        let u_local = bb.new_var(u_val);
        let pifold_asg = &parts[PIFOLD_PART].1;
        let mut lhs_locals: Vec<usize> = Vec::with_capacity(pifold_wiring.step5_lhs_ct.len());
        for &pv in &pifold_wiring.step5_lhs_ct {
            lhs_locals.push(bb.new_var(pifold_asg[pv]));
        }
        // Enforce u_local = Σ gamma^i * lhs_i
        let mut acc = bb.new_var(bf_zero);
        bb.enforce_var_eq_const(acc, bf_zero);
        let mut gamma_pow = bb.new_var(bf_one);
        bb.enforce_var_eq_const(gamma_pow, bf_one);
        for &lhs in &lhs_locals {
            let term = bb.new_var(bb.assignment[lhs] * bb.assignment[gamma_pow]);
            bb.enforce_mul(lhs, gamma_pow, term);
            let new_acc = bb.new_var(bb.assignment[acc] + bb.assignment[term]);
            bb.add_constraint(vec![(bf_one, acc), (bf_one, term)], vec![(bf_one, bb.one())], vec![(bf_one, new_acc)]);
            acc = new_acc;
            let new_pow = bb.new_var(bb.assignment[gamma_pow] * bb.assignment[gamma_local]);
            bb.enforce_mul(gamma_pow, gamma_local, new_pow);
            gamma_pow = new_pow;
        }
        bb.enforce_lc_times_one_eq_const(vec![(bf_one, acc), (-bf_one, u_local)]);
        let (bind_inst, bind_asg) = bb.into_instance();
        let bind_part_idx = parts.len();
        parts.push((bind_inst, bind_asg));

        // Glue u_local to PCS_g u_re (scalar).
        if pcs_g_wiring.u_re_vars.len() != 1 {
            return Err("WeGateDr1csBuilder: PCS_g expected scalar u_re (kappa*n=1)".to_string());
        }
        glue.push((pcs_g_part_idx, pcs_g_wiring.u_re_vars[0], bind_part_idx, u_local));
        // Glue lhs locals to Π_fold Step-5 lhs vars.
        const PIFOLD_PART: usize = 5;
        for (&pv, &lv) in pifold_wiring.step5_lhs_ct.iter().zip(lhs_locals.iter()) {
            glue.push((PIFOLD_PART, pv, bind_part_idx, lv));
        }
        glue.push((0, gamma_pose_var, bind_part_idx, gamma_local));

        merge_sparse_dr1cs_share_one_with_glue(&parts, &glue)
    }
}

