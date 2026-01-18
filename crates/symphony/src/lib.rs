//! # Symphony - Lattice-based folding SNARK
//!
//! Implementation of the Symphony protocol for post-quantum witness encryption.

#![allow(non_snake_case)]

// Core modules
pub mod setchk;
pub mod transcript;
pub mod utils;
pub mod public_coin_transcript;
pub mod recording_transcript;

// Symphony protocol modules
pub mod symphony_cm;
pub mod symphony_fold;
pub mod symphony_coins;
pub mod symphony_open;
pub mod symphony_pifold_batched;
pub mod symphony_we_relation;
pub mod we_statement;
pub mod poseidon_trace;
pub mod dpp_poseidon;
pub mod dpp_ajtai;
pub mod we_gate_arith;
pub mod dpp_sumcheck;
pub mod dpp_pifold_math;

// SP1 R1CS integration
pub mod sp1_r1cs_loader;
pub mod symphony_sp1_r1cs;

// Symphony Π_rg implementation
pub mod rp_rgchk;

// Streaming sumcheck (memory-efficient MLE handling)
pub mod streaming_sumcheck;
pub mod symphony_pifold_streaming;

// Polynomial Commitment Scheme (PCS) + dR1CS arithmetization (for aux → cm_f binding)
pub mod pcs;