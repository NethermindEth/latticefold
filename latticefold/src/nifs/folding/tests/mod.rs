use ark_std::test_rng;
use rand::{Rng};
use cyclotomic_rings::challenge_set::LatticefoldChallengeSet;
use crate::{
    arith::{r1cs::get_test_z_split, tests::get_test_ccs, Witness, CCCS},
    commitment::AjtaiCommitmentScheme,
    decomposition_parameters::{DecompositionParams},
    nifs::{
        decomposition::{
            DecompositionProver, DecompositionVerifier, LFDecompositionProver,
            LFDecompositionVerifier,
        },
        linearization::{
            LFLinearizationProver, LFLinearizationVerifier, LinearizationProver,
            LinearizationVerifier,
        },
    },
    transcript::poseidon::PoseidonTranscript,
};
use cyclotomic_rings::rings::{ SuitableRing};
use crate::arith::{CCS, LCCCS};
use crate::decomposition_parameters::test_params::DP;
use lattirust_ring::cyclotomic_ring::models::{
    babybear::RqNTT as BabyBearRqNTT, frog_ring::RqNTT as FrogRqNTT,
    goldilocks::RqNTT as GoldilocksRqNTT, stark_prime::RqNTT as StarkRqNTT,
};
const C: usize = 4;
const WIT_LEN: usize = 4;
const W: usize = WIT_LEN * DP::L;
fn setup_test_environment<RqNTT, CS, DP>(
    input: Option<usize>,
) -> (
    Vec<LCCCS<C, RqNTT>>,
    Vec<Witness<RqNTT>>,
    PoseidonTranscript<RqNTT, CS>,
    CCS<RqNTT>
)
where
    RqNTT: SuitableRing,
    CS: LatticefoldChallengeSet<RqNTT>,
    DP: DecompositionParams,
{
    let ccs = get_test_ccs::<RqNTT>(W);
    let mut rng = test_rng();
    let (_, x_ccs, w_ccs) = get_test_z_split::<RqNTT>(input.unwrap_or(rng.gen_range(0..64)));
    let scheme = AjtaiCommitmentScheme::rand(&mut rng);

    let wit = Witness::from_w_ccs::<DP>(w_ccs);
    let cm_i = CCCS {
        cm: wit.commit::<C, W, DP>(&scheme).unwrap(),
        x_ccs,
    };
    let mut prover_transcript = PoseidonTranscript::<RqNTT, CS>::default();
    let mut verifier_transcript = PoseidonTranscript::<RqNTT, CS>::default();

    let (_, linearization_proof) = LFLinearizationProver::<_, PoseidonTranscript<RqNTT, CS>>::prove(
        &cm_i,
        &wit,
        &mut prover_transcript,
        &ccs,
    )
        .unwrap();

    let lcccs = LFLinearizationVerifier::<_, PoseidonTranscript<RqNTT, CS>>::verify(
        &cm_i,
        &linearization_proof,
        &mut verifier_transcript,
        &ccs,
    )
        .unwrap();

    let (_, wit_vec, decomposition_proof) =
        LFDecompositionProver::<_, PoseidonTranscript<RqNTT, CS>>::prove::<W, C, DP>(
            &lcccs,
            &wit,
            &mut prover_transcript,
            &ccs,
            &scheme,
        )
            .unwrap();
    
    let lcccs_vec = LFDecompositionVerifier::<_, PoseidonTranscript<RqNTT, CS>>::verify::<C, DP>(
        &lcccs,
        &decomposition_proof,
        &mut verifier_transcript,
        &ccs,
    )
        .unwrap();
    
    let (lcccs, wit_s) = {
        let mut lcccs = lcccs_vec.clone();
        let mut lcccs_r = lcccs_vec;
        lcccs.append(&mut lcccs_r);

        let mut wit_s = wit_vec.clone();
        let mut wit_s_r = wit_vec;
        wit_s.append(&mut wit_s_r);

        (lcccs, wit_s)
    };

    (lcccs, wit_s, prover_transcript, ccs)
}


#[test]
fn test_setup_f_hat_mles() {
    type RqNTT = StarkRqNTT;
    
    
}