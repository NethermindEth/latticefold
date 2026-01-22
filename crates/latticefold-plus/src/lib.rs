//! # LatticeFold+

#![allow(non_snake_case)]

pub mod cm;
pub mod decomp;
pub mod lin;
pub mod mlin;
pub mod plus;
pub mod r1cs;
pub mod rgchk;
pub mod setchk;
pub mod streaming_sumcheck;
pub mod tensor_eval;
pub mod transcript;
pub mod utils;

// SP1 lift vacuity/soundness regression tests (WE-gate only).
#[cfg(all(test, feature = "we_gate"))]
mod sp1_lift_cheat_tests;

// WE/DPP arithmetization frontends (feature-gated; not needed in production proving path).
#[cfg(feature = "we_gate")]
pub mod recording_transcript;
#[cfg(feature = "we_gate")]
pub mod we_statement;
#[cfg(feature = "we_gate")]
pub mod we_gate_arith;

// SP1 shrink verifier R1LF loader helpers (feature-gated; research only).
// We gate these under `we_gate` so the WE/DPP benches can reuse them.
#[cfg(feature = "we_gate")]
pub mod sp1_r1lf;
#[cfg(feature = "we_gate")]
pub mod sp1_witness_io;

// API wrapper for the SP1 oneproof WE-gate harness (so downstream crates can call it in-process).
#[cfg(feature = "we_gate")]
pub mod sp1_oneproof_api;