//!
//! A crate containing the trait definition of a ring suitable to be used in the LatticeFold protocol,
//! a few ready-to-use rings and short challenge set machinery.
//!

#![cfg_attr(not(feature = "std"), no_std)]
#![allow(incomplete_features)]
#![forbid(unsafe_code)]

#[macro_use]
extern crate ark_std;

pub mod challenge_set;
pub mod rings;
pub mod rotation;

#[doc(hidden)]
mod ark_base {
    pub use ark_std::clone::Clone;
    pub use ark_std::convert::From;
    pub use ark_std::iter::Iterator;
    pub use ark_std::prelude::rust_2021::{derive, Debug};
    pub use ark_std::result::Result::{self, Err, Ok};
    pub use ark_std::vec::*;
}
