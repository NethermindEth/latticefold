use crate::arith::{CCS, LCCCS};
use crate::decomposition_parameters::test_params::DP;
use crate::nifs::folding::{FoldingProver, LFFoldingProver, LFFoldingVerifier};
use crate::nifs::FoldingProof;
use crate::{
    arith::{r1cs::get_test_z_split, tests::get_test_ccs, Witness, CCCS},
    commitment::AjtaiCommitmentScheme,
    decomposition_parameters::DecompositionParams,
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
use ark_std::test_rng;
use cyclotomic_rings::challenge_set::LatticefoldChallengeSet;
use cyclotomic_rings::rings::{FrogChallengeSet, GoldilocksChallengeSet, SuitableRing};
use lattirust_ring::cyclotomic_ring::models::{
    babybear::RqNTT as BabyBearRqNTT, frog_ring::RqNTT as FrogRqNTT,
    goldilocks::RqNTT as GoldilocksRqNTT, stark_prime::RqNTT as StarkRqNTT,
};
use num_traits::{One, Zero};
use rand::Rng;

const C: usize = 4;
const WIT_LEN: usize = 4;
const W: usize = WIT_LEN * DP::L;
fn setup_test_environment<RqNTT, CS, DP>(
    input: Option<usize>,
    generate_proof: bool,
) -> (
    Vec<LCCCS<C, RqNTT>>,
    Vec<Witness<RqNTT>>,
    PoseidonTranscript<RqNTT, CS>,
    CCS<RqNTT>,
    Option<FoldingProof<RqNTT>>,
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

    let (_, linearization_proof) =
        LFLinearizationProver::<_, PoseidonTranscript<RqNTT, CS>>::prove(
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

    let folding_proof = if generate_proof {
        Some(generate_folding_proof(
            &ccs,
            &mut prover_transcript,
            &lcccs,
            &wit_s,
        ))
    } else {
        None
    };

    (lcccs, wit_s, verifier_transcript, ccs, folding_proof)
}

fn generate_folding_proof<RqNTT, CS>(
    ccs: &CCS<RqNTT>,
    mut prover_transcript: &mut PoseidonTranscript<RqNTT, CS>,
    lcccs: &Vec<LCCCS<4, RqNTT>>,
    wit_s: &Vec<Witness<RqNTT>>,
) -> FoldingProof<RqNTT>
where
    RqNTT: SuitableRing,
    CS: LatticefoldChallengeSet<RqNTT>,
{
    let (_, _, folding_proof) = LFFoldingProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::prove::<
        C,
        DP,
    >(&lcccs, &wit_s, prover_transcript, &ccs)
    .unwrap();
    folding_proof
}

#[test]
fn test_setup_f_hat_mles() {
    type RqNTT = StarkRqNTT;
}

#[test]
fn test_calculate_claims() {
    type RqNTT = GoldilocksRqNTT;
    type CS = GoldilocksChallengeSet;

    let (lcccs_vec, _, _, _, _) = setup_test_environment::<RqNTT, CS, DP>(None, false);

    let alpha_s = vec![RqNTT::one(); 2 * DP::K];
    let zeta_s = vec![RqNTT::one(); 2 * DP::K];

    let (claim_g1, claim_g3) =
        LFFoldingVerifier::<RqNTT, PoseidonTranscript<RqNTT, CS>>::calculate_claims::<C>(
            &alpha_s, &zeta_s, &lcccs_vec,
        );

    assert_ne!(claim_g1, RqNTT::zero());
    assert_ne!(claim_g3, RqNTT::zero());

    let zero_alphas = vec![RqNTT::zero(); 2 * DP::K];
    let zero_zetas = vec![RqNTT::zero(); 2 * DP::K];

    let (zero_claim_g1, zero_claim_g3) =
        LFFoldingVerifier::<RqNTT, PoseidonTranscript<RqNTT, CS>>::calculate_claims::<C>(
            &zero_alphas,
            &zero_zetas,
            &lcccs_vec,
        );

    assert_eq!(zero_claim_g1, RqNTT::zero());
    assert_eq!(zero_claim_g3, RqNTT::zero());
}

#[test]
fn test_verify_sumcheck_proof() {
    type RqNTT = FrogRqNTT;
    type CS = FrogChallengeSet;

    let (lcccs_vec, _, mut transcript, ccs, proof) =
        setup_test_environment::<RqNTT, CS, DP>(None, true);
    let proof = proof.unwrap();

    let (alpha_s, beta_s, zeta_s, mu_s, poly_info) =
        LFFoldingVerifier::<RqNTT, PoseidonTranscript<RqNTT, CS>>::get_alphas_betas_zetas_mus::<
            C,
            DP,
        >(&lcccs_vec, ccs.s, &mut transcript)
        .unwrap();

    let (claim_g1, claim_g3) =
        LFFoldingVerifier::<RqNTT, PoseidonTranscript<RqNTT, CS>>::calculate_claims::<C>(
            &alpha_s, &zeta_s, &lcccs_vec,
        );

    let result = LFFoldingVerifier::<RqNTT, PoseidonTranscript<RqNTT, CS>>::verify_sumcheck_proof(
        &mut transcript,
        &poly_info,
        claim_g1 + claim_g3,
        &proof,
    );

    match result {
        Ok((r_0, expected_eval)) => {
            assert_eq!(r_0.len(), ccs.s);
            // We can add more assertions here if needed
        }
        Err(e) => panic!("Sumcheck verification failed: {:?}", e),
    }
}
