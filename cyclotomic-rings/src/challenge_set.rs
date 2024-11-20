//!
//!  Short challenge set API.
//!

use crate::ark_base::*;
use error::ChallengeSetError;

use crate::rings::SuitableRing;

pub mod error;

/// A trait to specify short challenge set for use in the LatticeFold protocol.
pub trait LatticefoldChallengeSet<R: SuitableRing> {
    /// Amount of bytes needed to obtain a single short challenge.
    const BYTES_NEEDED: usize;

    /// Given a slice of bytes `bs` returns the short challenge encode with these bytes
    /// in the coefficient form. Returns `TooFewBytes` error if there is not enough bytes
    /// to obtain a short challenge.
    fn short_challenge_from_random_bytes(
        bs: &[u8],
    ) -> Result<R::CoefficientRepresentation, ChallengeSetError>;
}
