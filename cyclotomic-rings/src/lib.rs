//!
//! A crate containing the trait definition of a ring suitable to be used in the LatticeFold protocol,
//! a few ready-to-use rings and short challenge set machinery.
//!

#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![forbid(unsafe_code)]

/// Short challenge set API.
pub mod challenge_set;
/// Cyclotomic ring API for the LatticeFold protocol.
pub mod rings;
mod rot_sum;

pub use rot_sum::*;
