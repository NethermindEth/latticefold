//! PCS (Polynomial Commitment Scheme) modules for Symphony.
//!
//! This directory is the long-term home for:
//! - PCS prover/verifier (non-arithmetized) logic
//! - DPP/dR1CS arithmetization of PCS verification
//!
//! ## Modules
//!
//! - `folding_pcs_l2`: Core ℓ=2 folding PCS verifier (paper 2024-281, Figure 5)
//! - `dpp_folding_pcs_l2`: DPP/dR1CS arithmetization of PCS verification
//! - `batchlin_pcs`: Batchlin integration layer (point conversion, batching)

pub mod folding_pcs_l2;
pub mod dpp_folding_pcs_l2;
pub mod batchlin_pcs;
pub mod cmf_pcs;

