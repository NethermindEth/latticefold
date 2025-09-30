use ark_std::{test_rng, vec::Vec};
use cyclotomic_rings::{challenge_set::LatticefoldChallengeSet, rings::SuitableRing};
use rand::Rng;

use crate::{
    arith::{r1cs::get_test_z_split, tests::get_test_ccs, Witness, CCCS, CCS, LCCCS},
    commitment::AjtaiCommitmentScheme,
    decomposition_parameters::DecompositionParams,
    nifs::{
        linearization::{LFLinearizationProver, LinearizationProver},
        NIFSProver, NIFSVerifier,
    },
    transcript::{poseidon::PoseidonTranscript, TranscriptWithShortChallenges},
};

fn setup_test_environment<
    RqNTT: SuitableRing,
    DP: DecompositionParams,
    CS: LatticefoldChallengeSet<RqNTT>,
>(
    kappa: usize,
    n: usize,
    wit_len: usize,
) -> (
    LCCCS<RqNTT>,   // acc
    Witness<RqNTT>, // w_acc
    CCCS<RqNTT>,    // cm_i
    Witness<RqNTT>, // w_i
    CCS<RqNTT>,
    AjtaiCommitmentScheme<RqNTT>,
) {
    let ccs = get_test_ccs::<RqNTT>(n, DP::L);
    let mut rng = test_rng();
    let (_, x_ccs, w_ccs) = get_test_z_split::<RqNTT>(rng.gen_range(0..64));
    let scheme = AjtaiCommitmentScheme::rand(kappa, n, &mut rng);

    let wit_i = Witness::from_w_ccs::<DP>(w_ccs);
    let cm_i = CCCS {
        cm: wit_i.commit::<DP>(&scheme).unwrap(),
        x_ccs: x_ccs.clone(),
    };

    let rand_w_ccs: Vec<RqNTT> = (0..wit_len).map(|i| RqNTT::from(i as u64)).collect();
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
    RqNTT: SuitableRing,
    CS: LatticefoldChallengeSet<RqNTT>,
    DP: DecompositionParams,
    T: TranscriptWithShortChallenges<RqNTT>,
>(
    kappa: usize,
    n: usize,
    wit_len: usize,
) {
    let (acc, w_acc, cm_i, w_i, ccs, scheme) =
        setup_test_environment::<RqNTT, DP, CS>(kappa, n, wit_len);

    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();

    let result = NIFSProver::<RqNTT, DP, T>::prove(
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
    RqNTT: SuitableRing,
    CS: LatticefoldChallengeSet<RqNTT>,
    DP: DecompositionParams,
    T: TranscriptWithShortChallenges<RqNTT>,
>(
    kappa: usize,
    n: usize,
    wit_len: usize,
) {
    let (acc, w_acc, cm_i, w_i, ccs, scheme) =
        setup_test_environment::<RqNTT, DP, CS>(kappa, n, wit_len);

    let mut prover_transcript = PoseidonTranscript::<RqNTT, CS>::default();
    let mut verifier_transcript = PoseidonTranscript::<RqNTT, CS>::default();

    let (_, _, proof) = NIFSProver::<RqNTT, DP, T>::prove(
        &acc,
        &w_acc,
        &cm_i,
        &w_i,
        &mut prover_transcript,
        &ccs,
        &scheme,
    )
    .unwrap();

    let result =
        NIFSVerifier::<RqNTT, DP, T>::verify(&acc, &cm_i, &proof, &mut verifier_transcript, &ccs);

    assert!(result.is_ok());
}

mod e2e_tests {
    use super::*;
    mod stark {
        use cyclotomic_rings::rings::{StarkChallengeSet, StarkRingNTT};

        use crate::{
            decomposition_parameters::{test_params::StarkDP, DecompositionParams},
            nifs::tests::{test_nifs_prove, test_nifs_verify},
            transcript::poseidon::PoseidonTranscript,
        };

        type RqNTT = StarkRingNTT;
        type CS = StarkChallengeSet;
        type DP = StarkDP;
        type T = PoseidonTranscript<RqNTT, CS>;

        const KAPPA: usize = 4;
        const WIT_LEN: usize = 4;
        const N: usize = WIT_LEN * DP::L;

        #[ignore]
        #[test]
        fn test_prove() {
            test_nifs_prove::<RqNTT, CS, DP, T>(KAPPA, N, WIT_LEN);
        }

        #[ignore]
        #[test]
        fn test_verify() {
            test_nifs_verify::<RqNTT, CS, DP, T>(KAPPA, N, WIT_LEN);
        }
    }

    mod goldilocks {
        use cyclotomic_rings::rings::{GoldilocksChallengeSet, GoldilocksRingNTT};

        use super::*;
        use crate::decomposition_parameters::test_params::GoldilocksDP;

        type RqNTT = GoldilocksRingNTT;
        type CS = GoldilocksChallengeSet;
        type DP = GoldilocksDP;
        type T = PoseidonTranscript<RqNTT, CS>;

        const KAPPA: usize = 4;
        const WIT_LEN: usize = 4;
        const N: usize = WIT_LEN * DP::L;

        #[test]
        fn test_prove() {
            test_nifs_prove::<RqNTT, CS, DP, T>(KAPPA, N, WIT_LEN);
        }

        #[test]
        fn test_verify() {
            test_nifs_verify::<RqNTT, CS, DP, T>(KAPPA, N, WIT_LEN);
        }
    }

    mod babybear {
        use cyclotomic_rings::rings::{BabyBearChallengeSet, BabyBearRingNTT};

        use super::*;
        use crate::decomposition_parameters::test_params::BabyBearDP;

        type RqNTT = BabyBearRingNTT;
        type CS = BabyBearChallengeSet;
        type DP = BabyBearDP;
        type T = PoseidonTranscript<RqNTT, CS>;

        const KAPPA: usize = 4;
        const WIT_LEN: usize = 4;
        const N: usize = WIT_LEN * DP::L;

        #[test]
        fn test_prove() {
            test_nifs_prove::<RqNTT, CS, DP, T>(KAPPA, N, WIT_LEN);
        }

        #[test]
        fn test_verify() {
            test_nifs_verify::<RqNTT, CS, DP, T>(KAPPA, N, WIT_LEN);
        }
    }
}
