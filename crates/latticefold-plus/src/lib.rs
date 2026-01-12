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

// WE/DPP arithmetization frontends (feature-gated; not needed in production proving path).
#[cfg(feature = "we_gate")]
pub mod recording_transcript;
#[cfg(feature = "we_gate")]
pub mod we_statement;
#[cfg(feature = "we_gate")]
pub mod we_gate_arith;

// SP1 shrink verifier R1CS loader helpers (feature-gated; research only).
// We gate these under `we_gate` so the WE/DPP benches can reuse them.
#[cfg(feature = "we_gate")]
pub mod sp1_r1cs;
#[cfg(feature = "we_gate")]
pub mod sp1_r1lf;
#[cfg(feature = "we_gate")]
pub mod sp1_lift_witness;