#![cfg_attr(not(feature = "std"), no_std)]
#![allow(clippy::type_complexity)]
#![allow(incomplete_features)]
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
    pub use ark_std::clone::Clone;
    pub use ark_std::convert::From;
    pub use ark_std::iter::Iterator;
    pub use ark_std::prelude::rust_2021::{derive, Debug};
    pub use ark_std::result::Result::{self, Err, Ok};
    pub use ark_std::string::{String, ToString};
    pub use ark_std::vec::*;
}
