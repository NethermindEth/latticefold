//! Relation-layer glue types for consumers that want to treat Symphony as an NP relation.
//!
//! This module is intentionally **small**:
//! - It defines the folded output shape (`FoldedOutput`) and
//! - Traits for application-specific reduced relations (`ReducedRelation`).

use stark_rings::{CoeffRing, Zq};

use crate::symphony_fold::{SymphonyBatchLin, SymphonyInstance};

/// Public output of the folding layer (the reduced instance).
#[derive(Clone, Debug)]
pub struct FoldedOutput<R: CoeffRing>
where
    R::BaseRing: Zq,
{
    pub folded_inst: SymphonyInstance<R>,
    pub folded_bat: SymphonyBatchLin<R>,
}

/// Reduced relation `R_o` checker (application-specific).
///
/// In the Symphony paper, `R_o` is the reduced relation after folding.
/// For WE/DPP, we avoid verifying an external SNARK proof and instead check the
/// underlying reduced relation directly (with an explicit reduced witness).
pub trait ReducedRelation<R: CoeffRing>
where
    R::BaseRing: Zq,
{
    type Witness;

    fn check(public: &FoldedOutput<R>, witness: &Self::Witness) -> Result<(), String>;
}

// =============================================================================
// Test/scaffolding helpers
// =============================================================================

/// Trivial reduced relation for testing / scaffolding.
pub struct TrivialRo;

impl<R: CoeffRing> ReducedRelation<R> for TrivialRo
where
    R::BaseRing: Zq,
{
    type Witness = ();

    fn check(_public: &FoldedOutput<R>, _witness: &Self::Witness) -> Result<(), String> {
        Ok(())
    }
}

