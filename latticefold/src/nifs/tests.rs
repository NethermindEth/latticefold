use crate::arith::r1cs::get_test_z_split;
use crate::arith::tests::get_test_ccs;
use crate::arith::{Witness, CCCS, CCS, LCCCS};
use crate::commitment::AjtaiCommitmentScheme;
use crate::decomposition_parameters::DecompositionParams;
use crate::nifs::linearization::{LFLinearizationProver, LinearizationProver};
use crate::nifs::{NIFSProver, NIFSVerifier};
use crate::transcript::poseidon::PoseidonTranscript;
use crate::transcript::TranscriptWithShortChallenges;
use ark_std::test_rng;
use cyclotomic_rings::challenge_set::LatticefoldChallengeSet;
use cyclotomic_rings::rings::SuitableRing;
use rand::Rng;

fn setup_test_environment<
    const C: usize,
    RqNTT: SuitableRing,
    DP: DecompositionParams,
    const W: usize,
    const WIT_LEN: usize,
    CS: LatticefoldChallengeSet<RqNTT>,
>() -> (
    LCCCS<C, RqNTT>, // acc
    Witness<RqNTT>,  // w_acc
    CCCS<C, RqNTT>,  // cm_i
    Witness<RqNTT>,  // w_i
    CCS<RqNTT>,
    AjtaiCommitmentScheme<C, W, RqNTT>,
) {
    let ccs = get_test_ccs::<RqNTT>(W, DP::L);
    let mut rng = test_rng();
    let (_, x_ccs, w_ccs) = get_test_z_split::<RqNTT>(rng.gen_range(0..64));
    let scheme = AjtaiCommitmentScheme::rand(&mut rng);

    let wit_i = Witness::from_w_ccs::<DP>(w_ccs);
    let cm_i = CCCS {
        cm: wit_i.commit::<C, W, DP>(&scheme).unwrap(),
        x_ccs: x_ccs.clone(),
    };

    let rand_w_ccs: Vec<RqNTT> = (0..WIT_LEN).map(|i| RqNTT::from(i as u64)).collect();
    let wit_acc = Witness::from_w_ccs::<DP>(rand_w_ccs);

    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();

    let (acc, _) = LFLinearizationProver::<_, PoseidonTranscript<RqNTT, CS>>::prove(
        &cm_i,
        &wit_acc,
        &mut transcript,
        &ccs,
    )
    .unwrap();
    (acc, wit_acc, cm_i, wit_i, ccs, scheme)
}

fn test_nifs_prove<
    const C: usize,
    const W: usize,
    const WIT_LEN: usize,
    RqNTT: SuitableRing,
    CS: LatticefoldChallengeSet<RqNTT>,
    DP: DecompositionParams,
    T: TranscriptWithShortChallenges<RqNTT>,
>() {
    let (acc, w_acc, cm_i, w_i, ccs, scheme) =
        setup_test_environment::<C, RqNTT, DP, W, WIT_LEN, CS>();

    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();

    let result = NIFSProver::<C, W, RqNTT, DP, T>::prove(
        &acc,
        &w_acc,
        &cm_i,
        &w_i,
        &mut transcript,
        &ccs,
        &scheme,
    );

    assert!(result.is_ok());
}

fn test_nifs_verify<
    const C: usize,
    const W: usize,
    const WIT_LEN: usize,
    RqNTT: SuitableRing,
    CS: LatticefoldChallengeSet<RqNTT>,
    DP: DecompositionParams,
    T: TranscriptWithShortChallenges<RqNTT>,
>() {
    let (acc, w_acc, cm_i, w_i, ccs, scheme) =
        setup_test_environment::<C, RqNTT, DP, W, WIT_LEN, CS>();

    let mut prover_transcript = PoseidonTranscript::<RqNTT, CS>::default();
    let mut verifier_transcript = PoseidonTranscript::<RqNTT, CS>::default();

    let (_, _, proof) = NIFSProver::<C, W, RqNTT, DP, T>::prove(
        &acc,
        &w_acc,
        &cm_i,
        &w_i,
        &mut prover_transcript,
        &ccs,
        &scheme,
    )
    .unwrap();

    let result = NIFSVerifier::<C, RqNTT, DP, T>::verify(
        &acc,
        &cm_i,
        &proof,
        &mut verifier_transcript,
        &ccs,
    );

    assert!(result.is_ok());
}

mod stark {
    use super::*;
    use crate::decomposition_parameters::test_params::StarkDP;
    use cyclotomic_rings::rings::{StarkChallengeSet, StarkRingNTT};

    type RqNTT = StarkRingNTT;
    type CS = StarkChallengeSet;
    type DP = StarkDP;
    type T = PoseidonTranscript<RqNTT, CS>;

    const C: usize = 4;
    const WIT_LEN: usize = 4;
    const W: usize = WIT_LEN * DP::L;

    #[ignore]
    #[test]
    fn test_prove() {
        test_nifs_prove::<C, W, WIT_LEN, RqNTT, CS, DP, T>();
    }

    #[ignore]
    #[test]
    fn test_verify() {
        test_nifs_verify::<C, W, WIT_LEN, RqNTT, CS, DP, T>();
    }
}

mod goldilocks {
    use super::*;
    use crate::decomposition_parameters::test_params::GoldilocksDP;
    use cyclotomic_rings::rings::{GoldilocksChallengeSet, GoldilocksRingNTT};

    type RqNTT = GoldilocksRingNTT;
    type CS = GoldilocksChallengeSet;
    type DP = GoldilocksDP;
    type T = PoseidonTranscript<RqNTT, CS>;

    const C: usize = 4;
    const WIT_LEN: usize = 4;
    const W: usize = WIT_LEN * DP::L;

    #[test]
    fn test_prove() {
        test_nifs_prove::<C, W, WIT_LEN, RqNTT, CS, DP, T>();
    }

    #[test]
    fn test_verify() {
        test_nifs_verify::<C, W, WIT_LEN, RqNTT, CS, DP, T>();
    }
}

mod babybear {
    use super::*;
    use crate::decomposition_parameters::test_params::BabyBearDP;
    use cyclotomic_rings::rings::{BabyBearChallengeSet, BabyBearRingNTT};

    type RqNTT = BabyBearRingNTT;
    type CS = BabyBearChallengeSet;
    type DP = BabyBearDP;
    type T = PoseidonTranscript<RqNTT, CS>;

    const C: usize = 4;
    const WIT_LEN: usize = 4;
    const W: usize = WIT_LEN * DP::L;

    #[test]
    fn test_prove() {
        test_nifs_prove::<C, W, WIT_LEN, RqNTT, CS, DP, T>();
    }

    #[test]
    fn test_verify() {
        test_nifs_verify::<C, W, WIT_LEN, RqNTT, CS, DP, T>();
    }
}
