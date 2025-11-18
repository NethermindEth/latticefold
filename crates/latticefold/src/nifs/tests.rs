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
    transcript::poseidon::PoseidonTranscript,
};

fn setup_test_environment<RqNTT: SuitableRing, CS: LatticefoldChallengeSet<RqNTT>>(
    dparams: &DecompositionParams,
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
    let ccs = get_test_ccs::<RqNTT>(n, dparams.l);
    let mut rng = test_rng();
    let (_, x_ccs, w_ccs) = get_test_z_split::<RqNTT>(rng.gen_range(0..64));
    let scheme = AjtaiCommitmentScheme::rand(kappa, n, &mut rng);

    let wit_i = Witness::from_w_ccs(w_ccs, dparams.B, dparams.l);
    let cm_i = CCCS {
        cm: wit_i.commit(&scheme).unwrap(),
        x_ccs: x_ccs.clone(),
    };

    let rand_w_ccs: Vec<RqNTT> = (0..wit_len).map(|i| RqNTT::from(i as u64)).collect();
    let wit_acc = Witness::from_w_ccs(rand_w_ccs, dparams.B, dparams.l);

    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();

    let (acc, _) = LFLinearizationProver::new(&mut transcript)
        .prove(&cm_i, &wit_acc, &ccs)
        .unwrap();
    (acc, wit_acc, cm_i, wit_i, ccs, scheme)
}

fn test_nifs_prove<RqNTT: SuitableRing, CS: LatticefoldChallengeSet<RqNTT> + Sync>(
    dparams: DecompositionParams,
    kappa: usize,
    n: usize,
    wit_len: usize,
) {
    let (acc, w_acc, cm_i, w_i, ccs, scheme) =
        setup_test_environment::<RqNTT, CS>(&dparams, kappa, n, wit_len);

    let transcript = PoseidonTranscript::<RqNTT, CS>::default();

    let result =
        NIFSProver::new(dparams, transcript).prove(&acc, &w_acc, &cm_i, &w_i, &ccs, &scheme);

    assert!(result.is_ok());
}

fn test_nifs_verify<RqNTT: SuitableRing, CS: LatticefoldChallengeSet<RqNTT> + Sync>(
    dparams: DecompositionParams,
    kappa: usize,
    n: usize,
    wit_len: usize,
) {
    let (acc, w_acc, cm_i, w_i, ccs, scheme) =
        setup_test_environment::<RqNTT, CS>(&dparams, kappa, n, wit_len);

    let prover_transcript = PoseidonTranscript::<RqNTT, CS>::default();
    let verifier_transcript = PoseidonTranscript::<RqNTT, CS>::default();

    let (_, _, proof) = NIFSProver::new(dparams.clone(), prover_transcript)
        .prove(&acc, &w_acc, &cm_i, &w_i, &ccs, &scheme)
        .unwrap();

    let result = NIFSVerifier::new(dparams, verifier_transcript).verify(&acc, &cm_i, &proof, &ccs);

    assert!(result.is_ok());
}

mod e2e_tests {
    use super::*;
    mod stark {
        use cyclotomic_rings::rings::{StarkChallengeSet, StarkRingNTT};

        use crate::{
            decomposition_parameters::test_params::dp_stark,
            nifs::tests::{test_nifs_prove, test_nifs_verify},
        };

        type RqNTT = StarkRingNTT;
        type CS = StarkChallengeSet;

        const KAPPA: usize = 4;
        const WIT_LEN: usize = 4;

        #[ignore]
        #[test]
        fn test_prove() {
            let dparams = dp_stark();
            let n = WIT_LEN * dparams.l;
            test_nifs_prove::<RqNTT, CS>(dparams, KAPPA, n, WIT_LEN);
        }

        #[ignore]
        #[test]
        fn test_verify() {
            let dparams = dp_stark();
            let n = WIT_LEN * dparams.l;
            test_nifs_verify::<RqNTT, CS>(dparams, KAPPA, n, WIT_LEN);
        }
    }

    mod goldilocks {
        use cyclotomic_rings::rings::{GoldilocksChallengeSet, GoldilocksRingNTT};

        use super::*;
        use crate::decomposition_parameters::test_params::dp_goldilocks;

        type RqNTT = GoldilocksRingNTT;
        type CS = GoldilocksChallengeSet;

        const KAPPA: usize = 4;
        const WIT_LEN: usize = 4;

        #[test]
        fn test_prove() {
            let dparams = dp_goldilocks();
            let n = WIT_LEN * dparams.l;
            test_nifs_prove::<RqNTT, CS>(dparams, KAPPA, n, WIT_LEN);
        }

        #[test]
        fn test_verify() {
            let dparams = dp_goldilocks();
            let n = WIT_LEN * dparams.l;
            test_nifs_verify::<RqNTT, CS>(dparams, KAPPA, n, WIT_LEN);
        }
    }

    mod babybear {
        use cyclotomic_rings::rings::{BabyBearChallengeSet, BabyBearRingNTT};

        use super::*;
        use crate::decomposition_parameters::test_params::dp_babybear;

        type RqNTT = BabyBearRingNTT;
        type CS = BabyBearChallengeSet;

        const KAPPA: usize = 4;
        const WIT_LEN: usize = 4;

        #[test]
        fn test_prove() {
            let dparams = dp_babybear();
            let n = WIT_LEN * dparams.l;
            test_nifs_prove::<RqNTT, CS>(dparams, KAPPA, n, WIT_LEN);
        }

        #[test]
        fn test_verify() {
            let dparams = dp_babybear();
            let n = WIT_LEN * dparams.l;
            test_nifs_verify::<RqNTT, CS>(dparams, KAPPA, n, WIT_LEN);
        }
    }
}
