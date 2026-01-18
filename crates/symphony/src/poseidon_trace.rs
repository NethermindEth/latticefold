//! Poseidon sponge trace replay (base-prime-field level).
//!
//! This module is a **mechanical** replay implementation of arkworks'
//! `PoseidonSponge`, but with hooks to record permutation boundaries.
//!
//! It is used as a stepping stone toward algebraic/DPP frontends: the WE-facing verifier
//! can output a `PoseidonTranscriptTrace`, and this module can deterministically replay it,
//! count permutations, and (optionally) emit a detailed permutation trace.

use ark_crypto_primitives::sponge::{poseidon::PoseidonConfig, DuplexSpongeMode};
use ark_ff::{BigInteger, PrimeField};
use thiserror::Error;

use crate::transcript::{PoseidonTraceOp, PoseidonTranscriptTrace};

#[derive(Debug, Error)]
pub enum PoseidonReplayError {
    #[error("trace mismatch: {0}")]
    Mismatch(String),
    #[error("invalid trace: {0}")]
    Invalid(String),
}

#[derive(Clone, Debug)]
pub struct PoseidonPermutationTrace<F: PrimeField> {
    pub before: Vec<F>,
    pub after: Vec<F>,
}

/// Returns `(after, states_after_each_round)`, where `states_after_each_round.len() == full_rounds + partial_rounds`.
///
/// Each round state is recorded **after** applying ARK, S-box, and MDS (matching the permutation loop).
pub fn permute_with_round_trace<F: PrimeField>(
    cfg: &PoseidonConfig<F>,
    before: &[F],
) -> Result<(Vec<F>, Vec<Vec<F>>), PoseidonReplayError> {
    if before.len() != cfg.rate + cfg.capacity {
        return Err(PoseidonReplayError::Invalid(format!(
            "bad state length: expected {}, got {}",
            cfg.rate + cfg.capacity,
            before.len()
        )));
    }
    let mut state = before.to_vec();
    let mut round_states: Vec<Vec<F>> = Vec::with_capacity(cfg.full_rounds + cfg.partial_rounds);
    let full_rounds_over_2 = cfg.full_rounds / 2;

    for i in 0..full_rounds_over_2 {
        apply_ark(cfg, &mut state, i);
        apply_sbox(cfg, &mut state, true);
        apply_mds(cfg, &mut state);
        round_states.push(state.clone());
    }

    for i in full_rounds_over_2..(full_rounds_over_2 + cfg.partial_rounds) {
        apply_ark(cfg, &mut state, i);
        apply_sbox(cfg, &mut state, false);
        apply_mds(cfg, &mut state);
        round_states.push(state.clone());
    }

    for i in (full_rounds_over_2 + cfg.partial_rounds)..(cfg.partial_rounds + cfg.full_rounds) {
        apply_ark(cfg, &mut state, i);
        apply_sbox(cfg, &mut state, true);
        apply_mds(cfg, &mut state);
        round_states.push(state.clone());
    }

    Ok((state, round_states))
}

#[derive(Clone, Debug)]
pub struct PoseidonSpongeReplayResult<F: PrimeField> {
    pub final_state: Vec<F>,
    pub permutes: Vec<PoseidonPermutationTrace<F>>,
}

/// Replay a `PoseidonTranscriptTrace` and ensure its recorded outputs are consistent with Poseidon.
///
/// This uses the same squeeze-bytes conversion as arkworks:
/// take `(MODULUS_BIT_SIZE-1)/8` bytes from each squeezed field element (little-endian).
pub fn replay_poseidon_transcript_trace<F: PrimeField>(
    cfg: &PoseidonConfig<F>,
    trace: &PoseidonTranscriptTrace<F>,
) -> Result<PoseidonSpongeReplayResult<F>, PoseidonReplayError> {
    replay_ops(cfg, &trace.ops)
}

pub fn replay_ops<F: PrimeField>(
    cfg: &PoseidonConfig<F>,
    ops: &[PoseidonTraceOp<F>],
) -> Result<PoseidonSpongeReplayResult<F>, PoseidonReplayError> {
    let mut state = vec![F::zero(); cfg.rate + cfg.capacity];
    let mut mode = DuplexSpongeMode::Absorbing { next_absorb_index: 0 };
    let mut permutes: Vec<PoseidonPermutationTrace<F>> = Vec::new();

    for (op_idx, op) in ops.iter().enumerate() {
        match op {
            PoseidonTraceOp::Absorb(elems) => {
                if elems.is_empty() {
                    continue;
                }
                let mut absorb_index = match mode {
                    DuplexSpongeMode::Absorbing { next_absorb_index } => next_absorb_index,
                    DuplexSpongeMode::Squeezing { .. } => {
                        // switching from squeezing to absorbing permutes first.
                        permute_in_place_with_record(cfg, &mut state, &mut permutes);
                        0
                    }
                };
                if absorb_index == cfg.rate {
                    permute_in_place_with_record(cfg, &mut state, &mut permutes);
                    absorb_index = 0;
                }
                let next_absorb_index =
                    absorb_internal_with_record(cfg, &mut state, &mut permutes, absorb_index, elems);
                mode = DuplexSpongeMode::Absorbing {
                    next_absorb_index,
                };
            }
            PoseidonTraceOp::SqueezeField(out) => {
                if out.is_empty() {
                    continue;
                }
                let mut squeezed = vec![F::zero(); out.len()];

                match mode {
                    DuplexSpongeMode::Absorbing { .. } => {
                        permute_in_place_with_record(cfg, &mut state, &mut permutes);
                        let next = squeeze_internal_with_record(
                            cfg,
                            &mut state,
                            &mut permutes,
                            0,
                            &mut squeezed,
                        );
                        mode = DuplexSpongeMode::Squeezing {
                            next_squeeze_index: next,
                        };
                    }
                    DuplexSpongeMode::Squeezing { next_squeeze_index } => {
                        let mut squeeze_index = next_squeeze_index;
                        if squeeze_index == cfg.rate {
                            permute_in_place_with_record(cfg, &mut state, &mut permutes);
                            squeeze_index = 0;
                        }
                        let next = squeeze_internal_with_record(
                            cfg,
                            &mut state,
                            &mut permutes,
                            squeeze_index,
                            &mut squeezed,
                        );
                        mode = DuplexSpongeMode::Squeezing {
                            next_squeeze_index: next,
                        };
                    }
                }

                if &squeezed != out {
                    return Err(PoseidonReplayError::Mismatch(format!(
                        "SqueezeField mismatch at op #{op_idx}: expected {:?}, got {:?}",
                        out, squeezed
                    )));
                }
            }
            PoseidonTraceOp::SqueezeBytes { n, out } => {
                // Ark implementation:
                // - squeeze_native_field_elements(num_elements)
                // - for each element, take usable_bytes = (MODULUS_BIT_SIZE-1)/8 from LE bigint
                // - truncate
                let usable_bytes = ((F::MODULUS_BIT_SIZE - 1) / 8) as usize;
                if usable_bytes == 0 {
                    return Err(PoseidonReplayError::Invalid(
                        "usable_bytes computed as 0".to_string(),
                    ));
                }
                let num_elements = (*n + usable_bytes - 1) / usable_bytes;
                let mut src_elements = vec![F::zero(); num_elements];

                // squeeze_native_field_elements
                match mode {
                    DuplexSpongeMode::Absorbing { .. } => {
                        permute_in_place_with_record(cfg, &mut state, &mut permutes);
                        let next = squeeze_internal_with_record(
                            cfg,
                            &mut state,
                            &mut permutes,
                            0,
                            &mut src_elements,
                        );
                        mode = DuplexSpongeMode::Squeezing {
                            next_squeeze_index: next,
                        };
                    }
                    DuplexSpongeMode::Squeezing { next_squeeze_index } => {
                        let mut squeeze_index = next_squeeze_index;
                        if squeeze_index == cfg.rate {
                            permute_in_place_with_record(cfg, &mut state, &mut permutes);
                            squeeze_index = 0;
                        }
                        let next = squeeze_internal_with_record(
                            cfg,
                            &mut state,
                            &mut permutes,
                            squeeze_index,
                            &mut src_elements,
                        );
                        mode = DuplexSpongeMode::Squeezing {
                            next_squeeze_index: next,
                        };
                    }
                }

                // Convert to bytes as ark does.
                let mut bytes: Vec<u8> = Vec::with_capacity(usable_bytes * num_elements);
                for elem in &src_elements {
                    let elem_bytes = elem.into_bigint().to_bytes_le();
                    bytes.extend_from_slice(&elem_bytes[..usable_bytes]);
                }
                bytes.truncate(*n);

                if &bytes != out {
                    return Err(PoseidonReplayError::Mismatch(format!(
                        "SqueezeBytes mismatch at op #{op_idx}: expected {:?}, got {:?}",
                        out, bytes
                    )));
                }
            }
        }
    }

    Ok(PoseidonSpongeReplayResult { final_state: state, permutes })
}

/// Find the **SqueezeBytes index** of the first `PoseidonTraceOp::SqueezeBytes` that occurs
/// after observing an `PoseidonTraceOp::Absorb([marker])`.
///
/// Returns an index into the sequence of `SqueezeBytes` ops (0-based), which matches how
/// `WeGateDr1csBuilder` indexes `PoseidonByteWiring.squeeze_byte_ranges`.
pub fn find_squeeze_bytes_idx_after_absorb_marker<F: PrimeField>(
    ops: &[PoseidonTraceOp<F>],
    marker: F,
) -> Option<usize> {
    let mut saw_marker = false;
    let mut squeeze_idx = 0usize;
    for op in ops {
        match op {
            PoseidonTraceOp::Absorb(v) => {
                if v.len() == 1 && v[0] == marker {
                    saw_marker = true;
                }
            }
            PoseidonTraceOp::SqueezeBytes { .. } => {
                if saw_marker {
                    return Some(squeeze_idx);
                }
                squeeze_idx += 1;
            }
            _ => {}
        }
    }
    None
}

fn absorb_internal_with_record<F: PrimeField>(
    cfg: &PoseidonConfig<F>,
    state: &mut [F],
    permutes: &mut Vec<PoseidonPermutationTrace<F>>,
    mut rate_start_index: usize,
    mut remaining: &[F],
) -> usize {
    loop {
        if rate_start_index + remaining.len() <= cfg.rate {
            for (i, element) in remaining.iter().enumerate() {
                state[cfg.capacity + i + rate_start_index] += element;
            }
            // done
            return rate_start_index + remaining.len();
        }

        let num_absorb = cfg.rate - rate_start_index;
        for (i, element) in remaining.iter().enumerate().take(num_absorb) {
            state[cfg.capacity + i + rate_start_index] += element;
        }
        permute_in_place_with_record(cfg, state, permutes);
        remaining = &remaining[num_absorb..];
        rate_start_index = 0;
    }
}

/// Squeeze `output` field elements starting from `rate_start_index`, recording any intermediate
/// permutations into `permutes`. Returns the next squeeze index (in `0..=rate`).
fn squeeze_internal_with_record<F: PrimeField>(
    cfg: &PoseidonConfig<F>,
    state: &mut [F],
    permutes: &mut Vec<PoseidonPermutationTrace<F>>,
    mut rate_start_index: usize,
    output: &mut [F],
) -> usize {
    let mut out_pos = 0usize;
    while out_pos < output.len() {
        if rate_start_index == cfg.rate {
            permute_in_place_with_record(cfg, state, permutes);
            rate_start_index = 0;
        }
        let take = core::cmp::min(cfg.rate - rate_start_index, output.len() - out_pos);
        output[out_pos..out_pos + take].clone_from_slice(
            &state[cfg.capacity + rate_start_index..(cfg.capacity + rate_start_index + take)],
        );
        out_pos += take;
        rate_start_index += take;
        if out_pos < output.len() && rate_start_index == cfg.rate {
            permute_in_place_with_record(cfg, state, permutes);
            rate_start_index = 0;
        }
    }
    rate_start_index
}

fn permute_in_place_with_record<F: PrimeField>(
    cfg: &PoseidonConfig<F>,
    state: &mut [F],
    permutes: &mut Vec<PoseidonPermutationTrace<F>>,
) {
    let before = state.to_vec();
    permute_in_place(cfg, state);
    let after = state.to_vec();
    permutes.push(PoseidonPermutationTrace { before, after });
}

fn permute_in_place<F: PrimeField>(cfg: &PoseidonConfig<F>, state: &mut [F]) {
    let full_rounds_over_2 = cfg.full_rounds / 2;

    for i in 0..full_rounds_over_2 {
        apply_ark(cfg, state, i);
        apply_sbox(cfg, state, true);
        apply_mds(cfg, state);
    }

    for i in full_rounds_over_2..(full_rounds_over_2 + cfg.partial_rounds) {
        apply_ark(cfg, state, i);
        apply_sbox(cfg, state, false);
        apply_mds(cfg, state);
    }

    for i in (full_rounds_over_2 + cfg.partial_rounds)..(cfg.partial_rounds + cfg.full_rounds) {
        apply_ark(cfg, state, i);
        apply_sbox(cfg, state, true);
        apply_mds(cfg, state);
    }
}

fn apply_sbox<F: PrimeField>(cfg: &PoseidonConfig<F>, state: &mut [F], is_full_round: bool) {
    if is_full_round {
        for elem in state.iter_mut() {
            *elem = elem.pow(&[cfg.alpha]);
        }
    } else {
        state[0] = state[0].pow(&[cfg.alpha]);
    }
}

fn apply_ark<F: PrimeField>(cfg: &PoseidonConfig<F>, state: &mut [F], round_number: usize) {
    for (i, state_elem) in state.iter_mut().enumerate() {
        state_elem.add_assign(&cfg.ark[round_number][i]);
    }
}

fn apply_mds<F: PrimeField>(cfg: &PoseidonConfig<F>, state: &mut [F]) {
    let mut new_state = vec![F::zero(); state.len()];
    for i in 0..state.len() {
        let mut cur = F::zero();
        for (j, state_elem) in state.iter().enumerate() {
            cur.add_assign(&state_elem.mul(&cfg.mds[i][j]));
        }
        new_state[i] = cur;
    }
    state.clone_from_slice(&new_state[..state.len()]);
}

