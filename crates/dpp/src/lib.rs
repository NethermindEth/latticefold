//! Dot-Product Proofs (DPP) — prototype implementation.
//!
//! This crate is a research-engineering implementation effort following:
//! - Bitansky, Harsha, Ishai, Rothblum, Wu, “Dot-Product Proofs and Their Applications”
//!   (FOCS 2024 / ePrint 2024/1138, Rev 2 on Dec 30 2025).
//!
//! Scope (initial):
//! - Implement the “query packing” transformation (Construction 5.21 / Theorem 5.7) that
//!   converts a bounded k-query (F)LPCP into a 1-query DPP.
//! - Implement a simple FLPCP for deterministic R1CS (dR1CS) based on systematic Reed–Solomon
//!   multiplication codes (Theorem 4.5/4.6), as a baseline outer FLPCP.
//!
//! Notes:
//! - This is NOT a production-ready cryptographic implementation.
//! - Soundness/parameter regimes require careful instantiation and measurement.

#![forbid(unsafe_code)]

pub mod subset_sum;
pub mod packing;
pub mod embedding;
pub mod boolean_proof;
pub mod pipeline;

pub mod rs;
pub mod dr1cs_flpcp;
pub mod sparse;

pub use packing::{
    BoundedFlpcpSparse, DppFromBoundedFlpcp, DppFromBoundedFlpcpSparse, PackedDppParams,
    PackedDppQuery, PackedDppQuerySparse,
};
pub use embedding::{EmbeddedFlpcpSparse, EmbeddingParams};
pub use boolean_proof::BooleanProofFlpcpSparse;
pub use subset_sum::{decode_bounded_subset_sum, SubsetSumError};
pub use sparse::SparseVec;
