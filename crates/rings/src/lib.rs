//!
//! A crate containing the trait definition of a ring suitable to be used in the LatticeFold protocol,
//! a few ready-to-use rings and small challenge set machinery.
//!

#![cfg_attr(not(feature = "std"), no_std)]
#![forbid(unsafe_code)]

#[macro_use]
extern crate ark_std;

pub mod challenge_set;
pub mod rings;
pub mod rotation;

#[doc(hidden)]
mod ark_base {
    pub use ark_std::{
        clone::Clone,
        convert::From,
        iter::Iterator,
        prelude::rust_2021::{derive, Debug},
        result::Result::{self, Err, Ok},
        vec::*,
    };
}
