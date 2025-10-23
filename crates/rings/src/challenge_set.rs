//!
//!  Small challenge set API.
//!

use error::ChallengeSetError;

use crate::{ark_base::*, rings::SuitableRing};

pub mod error;

/// A trait to specify small challenge set for use in the LatticeFold protocol.
pub trait ChallengeSet<R: SuitableRing> {
    /// Amount of bytes needed to obtain a single small challenge.
    const BYTES_NEEDED: usize;

    /// Given a slice of bytes `bs` returns the small challenge encode with these bytes
    /// in the coefficient form. Returns `TooFewBytes` error if there is not enough bytes
    /// to obtain a small challenge.
    fn small_challenge_from_random_bytes(
        bs: &[u8],
    ) -> Result<R::CoefficientRepresentation, ChallengeSetError>;
}
