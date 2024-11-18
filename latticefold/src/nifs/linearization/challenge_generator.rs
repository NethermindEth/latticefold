use crate::nifs::linearization::{BetaChallengeGenerator, ChallengeGenerator};
use crate::transcript::Transcript;
use ark_ff::Field;
use ark_ff::PrimeField;
use lattirust_ring::OverField;

impl<NTT: OverField> ChallengeGenerator<NTT> for BetaChallengeGenerator<NTT> {
    fn generate_challenges(transcript: &mut impl Transcript<NTT>, log_m: usize) -> Vec<NTT> {
        transcript.absorb_field_element(&<NTT::BaseRing as Field>::from_base_prime_field(
            <NTT::BaseRing as Field>::BasePrimeField::from_be_bytes_mod_order(b"beta_s"),
        ));

        transcript
            .get_challenges(log_m)
            .into_iter()
            .map(|x| x.into())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::transcript::poseidon::PoseidonTranscript;
    use cyclotomic_rings::challenge_set::LatticefoldChallengeSet;
    use cyclotomic_rings::rings::{
        BabyBearChallengeSet, FrogChallengeSet, GoldilocksChallengeSet, StarkChallengeSet,
        SuitableRing,
    };
    use lattirust_ring::cyclotomic_ring::models::{
        babybear::RqNTT as BabyBearRqNTT, frog_ring::RqNTT as FrogRqNTT,
        goldilocks::RqNTT as GoldilocksRqNTT, stark_prime::RqNTT as StarkRqNTT,
    };

    fn test_challenge_generator<R: SuitableRing, CS: LatticefoldChallengeSet<R>>() {
        let mut transcript = PoseidonTranscript::<R, CS>::default();
        let log_m = 3;

        let challenges = BetaChallengeGenerator::<R>::generate_challenges(&mut transcript, log_m);

        assert_eq!(challenges.len(), log_m);
        // Verify challenges are deterministic
        let mut transcript2 = PoseidonTranscript::<R, CS>::default();
        let challenges2 = BetaChallengeGenerator::<R>::generate_challenges(&mut transcript2, log_m);
        assert_eq!(challenges, challenges2);
    }

    #[test]
    fn test_stark_challenge_generator() {
        test_challenge_generator::<StarkRqNTT, StarkChallengeSet>();
    }

    #[test]
    fn test_goldilocks_challenge_generator() {
        test_challenge_generator::<GoldilocksRqNTT, GoldilocksChallengeSet>();
    }

    #[test]
    fn test_frog_challenge_generator() {
        test_challenge_generator::<FrogRqNTT, FrogChallengeSet>();
    }

    #[test]
    fn test_babybear_challenge_generator() {
        test_challenge_generator::<BabyBearRqNTT, BabyBearChallengeSet>();
    }
}
