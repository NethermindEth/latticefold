#![cfg_attr(not(feature = "std"), no_std)]
#![allow(clippy::type_complexity)]
#![allow(non_snake_case)]
#![forbid(unsafe_code)]

#[macro_use]
extern crate ark_std;

extern crate alloc;

pub mod arith;
pub mod commitment;
pub mod decomposition_parameters;
pub mod nifs;
pub mod transcript;
pub mod utils;

#[doc(hidden)]
mod ark_base {
    pub use ark_std::{
        clone::Clone,
        convert::From,
        iter::Iterator,
        prelude::rust_2021::{derive, Debug},
        result::Result::{self, Err, Ok},
        string::{String, ToString},
        vec::*,
    };
}
