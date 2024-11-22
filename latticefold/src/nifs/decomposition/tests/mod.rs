use cyclotomic_rings::{challenge_set::LatticefoldChallengeSet, rings::SuitableRing};
use lattirust_poly::mle::DenseMultilinearExtension;
use rand::Rng;

use crate::nifs::error::DecompositionError;
use crate::nifs::linearization::utils::compute_u;
use crate::nifs::mle_helpers::{evaluate_mles, to_mles};
use crate::{
    arith::{r1cs::get_test_z_split, tests::get_test_ccs, utils::mat_vec_mul, Witness, CCS, LCCCS},
    ark_base::*,
    commitment::{AjtaiCommitmentScheme, Commitment},
    decomposition_parameters::DecompositionParams,
    nifs::decomposition::{
        DecompositionProver, DecompositionVerifier, LFDecompositionProver, LFDecompositionVerifier,
    },
    transcript::poseidon::PoseidonTranscript,
};

fn generate_decomposition_args<RqNTT, CS, DP, const WIT_LEN: usize, const W: usize>() -> (
    LCCCS<4, RqNTT>,
    PoseidonTranscript<RqNTT, CS>,
    PoseidonTranscript<RqNTT, CS>,
    CCS<RqNTT>,
    Witness<RqNTT>,
    AjtaiCommitmentScheme<4, W, RqNTT>,
)
where
    RqNTT: SuitableRing,
    CS: LatticefoldChallengeSet<RqNTT>,
    DP: DecompositionParams,
{
    let mut rng = ark_std::test_rng();
    let input: usize = rng.gen_range(1..101);
    let ccs = get_test_ccs(W, DP::L);
    let log_m = ccs.s;

    let scheme = AjtaiCommitmentScheme::rand(&mut rng);
    let (_, x_ccs, _) = get_test_z_split::<RqNTT>(input);

    let wit = Witness::rand::<_, DP>(&mut rng, WIT_LEN);

    let mut z: Vec<RqNTT> = Vec::with_capacity(x_ccs.len() + WIT_LEN + 1);

    z.extend_from_slice(&x_ccs);
    z.push(RqNTT::one());
    z.extend_from_slice(&wit.w_ccs);

    let cm: Commitment<4, RqNTT> = scheme.commit_ntt(&wit.f).unwrap();

    let r: Vec<RqNTT> = (0..log_m).map(|_| RqNTT::rand(&mut rng)).collect();
    let Mz_mles: Vec<DenseMultilinearExtension<RqNTT>> = ccs
        .M
        .iter()
        .map(|M| DenseMultilinearExtension::from_slice(log_m, &mat_vec_mul(M, &z).unwrap()))
        .collect();

    let u = compute_u(&Mz_mles, &r).unwrap();

    let v = evaluate_mles::<RqNTT, &DenseMultilinearExtension<RqNTT>, _, DecompositionError>(
        &to_mles::<_, _, DecompositionError>(log_m, &wit.f_hat).unwrap(),
        &r,
    )
    .unwrap();

    let lcccs = LCCCS {
        r,
        v,
        cm,
        u,
        x_w: x_ccs,
        h: RqNTT::one(),
    };

    (
        lcccs,
        PoseidonTranscript::<RqNTT, CS>::default(),
        PoseidonTranscript::<RqNTT, CS>::default(),
        ccs,
        wit,
        scheme,
    )
}

fn test_decomposition<RqNTT, CS, DP, const WIT_LEN: usize, const W: usize>()
where
    RqNTT: SuitableRing,
    CS: LatticefoldChallengeSet<RqNTT>,
    DP: DecompositionParams,
{
    let (lcccs, mut verifier_transcript, mut prover_transcript, ccs, wit, scheme) =
        generate_decomposition_args::<RqNTT, CS, DP, WIT_LEN, W>();

    let (_, _, decomposition_proof) =
        LFDecompositionProver::<_, PoseidonTranscript<RqNTT, CS>>::prove::<W, 4, DP>(
            &lcccs,
            &wit,
            &mut prover_transcript,
            &ccs,
            &scheme,
        )
        .unwrap();

    let res = LFDecompositionVerifier::<_, PoseidonTranscript<RqNTT, CS>>::verify::<4, DP>(
        &lcccs,
        &decomposition_proof,
        &mut verifier_transcript,
        &ccs,
    );
    assert!(res.is_ok())
}

mod stark {
    use crate::arith::r1cs::get_test_dummy_z_split;
    use crate::arith::tests::get_test_dummy_ccs;
    use crate::arith::{Witness, CCCS};
    use crate::commitment::AjtaiCommitmentScheme;
    use crate::decomposition_parameters::{test_params::StarkDP as DP, DecompositionParams};
    use crate::nifs::linearization::{
        LFLinearizationProver, LFLinearizationVerifier, LinearizationProver, LinearizationVerifier,
    };
    use crate::transcript::poseidon::PoseidonTranscript;
    use cyclotomic_rings::rings::StarkChallengeSet;
    use lattirust_ring::cyclotomic_ring::models::stark_prime::RqNTT;

    type CS = StarkChallengeSet;
    const WIT_LEN: usize = 4;
    const W: usize = WIT_LEN * DP::L;
    #[test]
    fn test_decomposition() {
        super::test_decomposition::<RqNTT, CS, DP, WIT_LEN, W>();
    }

    #[test]
    fn test_dummy_decomposition() {
        type R = RqNTT;
        type CS = StarkChallengeSet;
        type T = PoseidonTranscript<R, CS>;

        const C: usize = 16;
        const X_LEN: usize = 1;
        const WIT_LEN: usize = 2048;
        const W: usize = WIT_LEN * DP::L; // the number of columns of the Ajtai matrix
        let r1cs_rows_size = X_LEN + WIT_LEN + 1; // Let's have a square matrix
        let ccs = get_test_dummy_ccs::<R, X_LEN, WIT_LEN, W>(r1cs_rows_size, DP::L);
        let (_, x_ccs, w_ccs) = get_test_dummy_z_split::<R, X_LEN, WIT_LEN>();
        let mut rng = ark_std::test_rng();
        let scheme = AjtaiCommitmentScheme::rand(&mut rng);
        let wit = Witness::from_w_ccs::<DP>(w_ccs);

        let cm_i = CCCS {
            cm: wit.commit::<C, W, DP>(&scheme).unwrap(),
            x_ccs,
        };
        let mut transcript = PoseidonTranscript::<R, CS>::default();
        let res = LFLinearizationProver::<_, T>::prove(&cm_i, &wit, &mut transcript, &ccs);
        let mut transcript = PoseidonTranscript::<R, CS>::default();
        let res = LFLinearizationVerifier::<_, PoseidonTranscript<R, CS>>::verify(
            &cm_i,
            &res.expect("Linearization proof generation error").1,
            &mut transcript,
            &ccs,
        );
        res.expect("Linearization Verification error");
    }
}

mod goldilocks {
    use crate::decomposition_parameters::{test_params::GoldilocksDP as DP, DecompositionParams};
    use cyclotomic_rings::rings::GoldilocksChallengeSet;
    use lattirust_ring::cyclotomic_ring::models::goldilocks::RqNTT;
    type CS = GoldilocksChallengeSet;

    const WIT_LEN: usize = 4;
    const W: usize = WIT_LEN * DP::L;
    #[test]
    fn test_decomposition() {
        super::test_decomposition::<RqNTT, CS, DP, WIT_LEN, W>();
    }
}

mod frog {
    use crate::decomposition_parameters::{test_params::GoldilocksDP as DP, DecompositionParams};
    use cyclotomic_rings::rings::FrogChallengeSet;
    use lattirust_ring::cyclotomic_ring::models::frog_ring::RqNTT;
    type CS = FrogChallengeSet;

    const WIT_LEN: usize = 4;
    const W: usize = WIT_LEN * DP::L;
    #[test]
    fn test_decomposition() {
        super::test_decomposition::<RqNTT, CS, DP, WIT_LEN, W>();
    }
}

mod babybear {
    use crate::decomposition_parameters::{test_params::BabyBearDP as DP, DecompositionParams};
    use cyclotomic_rings::rings::BabyBearChallengeSet;
    use lattirust_ring::cyclotomic_ring::models::babybear::RqNTT;
    type CS = BabyBearChallengeSet;

    const WIT_LEN: usize = 4;
    const W: usize = WIT_LEN * DP::L;
    #[test]
    fn test_decomposition() {
        super::test_decomposition::<RqNTT, CS, DP, WIT_LEN, W>();
    }
}
