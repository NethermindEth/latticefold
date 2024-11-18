use ark_std::marker::PhantomData;

use cyclotomic_rings::rings::SuitableRing;
use lattirust_ring::OverField;

use crate::{
    arith::{Witness, CCCS, CCS, LCCCS},
    commitment::AjtaiCommitmentScheme,
    decomposition_parameters::DecompositionParams,
    transcript::TranscriptWithShortChallenges,
};
use decomposition::*;
use error::LatticefoldError;
use folding::{FoldingProof, FoldingProver, FoldingVerifier, LFFoldingProver, LFFoldingVerifier};
use linearization::{
    LFLinearizationProver, LFLinearizationVerifier, LinearizationProof, LinearizationProver,
    LinearizationVerifier,
};
pub mod decomposition;
pub mod error;
pub mod folding;
pub mod linearization;

/// `C` is the length of Ajtai commitment vectors.
/// `NTT` is a cyclotomic ring in the NTT form.
#[derive(Clone)]
pub struct LFProof<const C: usize, NTT: OverField> {
    pub linearization_proof: LinearizationProof<NTT>,
    pub decomposition_proof_l: DecompositionProof<C, NTT>,
    pub decomposition_proof_r: DecompositionProof<C, NTT>,
    pub folding_proof: FoldingProof<NTT>,
}

/// `C` is the length of commitment vectors or, equivalently, the number of rows of the Ajtai matrix.
/// `W` is the length of witness vectors or, equivalently, the number of columns of the Ajtai matrix.
/// `NTT` is a suitable cyclotomic ring.
/// `P` is the decomposition parameters.
/// `T` is the FS-transform transcript.
pub struct NIFSProver<const C: usize, const W: usize, NTT, P, T> {
    _r: PhantomData<NTT>,
    _p: PhantomData<P>,
    _t: PhantomData<T>,
}

impl<
        const C: usize,
        const W: usize,
        NTT: SuitableRing,
        P: DecompositionParams,
        T: TranscriptWithShortChallenges<NTT>,
    > NIFSProver<C, W, NTT, P, T>
{
    pub fn prove(
        acc: &LCCCS<C, NTT>,
        w_acc: &Witness<NTT>,
        cm_i: &CCCS<C, NTT>,
        w_i: &Witness<NTT>,
        transcript: &mut impl TranscriptWithShortChallenges<NTT>,
        ccs: &CCS<NTT>,
        scheme: &AjtaiCommitmentScheme<C, W, NTT>,
    ) -> Result<(LCCCS<C, NTT>, Witness<NTT>, LFProof<C, NTT>), LatticefoldError<NTT>> {
        let (linearized_cm_i, linearization_proof) =
            LFLinearizationProver::<_, T>::prove(cm_i, w_i, transcript, ccs)?;
        let (decomposed_lcccs_l, decomposed_wit_l, decomposition_proof_l) =
            LFDecompositionProver::<_, T>::prove::<W, C, P>(acc, w_acc, transcript, ccs, scheme)?;
        let (decomposed_lcccs_r, decomposed_wit_r, decomposition_proof_r) =
            LFDecompositionProver::<_, T>::prove::<W, C, P>(
                &linearized_cm_i,
                w_i,
                transcript,
                ccs,
                scheme,
            )?;

        let (lcccs, wit_s) = {
            let mut lcccs = decomposed_lcccs_l;
            let mut lcccs_r = decomposed_lcccs_r;
            lcccs.append(&mut lcccs_r);

            let mut wit_s = decomposed_wit_l;
            let mut wit_s_r = decomposed_wit_r;
            wit_s.append(&mut wit_s_r);

            (lcccs, wit_s)
        };

        let (folded_lcccs, wit, folding_proof) =
            LFFoldingProver::<_, T>::prove::<C, P>(&lcccs, &wit_s, transcript, ccs)?;

        Ok((
            folded_lcccs,
            wit,
            LFProof {
                linearization_proof,
                decomposition_proof_l,
                decomposition_proof_r,
                folding_proof,
            },
        ))
    }
}

/// `C` is the length of commitment vectors or, equivalently, the number of rows of the Ajtai matrix.
/// `W` is the length of witness vectors or, equivalently, the number of columns of the Ajtai matrix.
/// `NTT` is a suitable cyclotomic ring.
/// `P` is the decomposition parameters.
/// `T` is the FS-transform transcript.
pub struct NIFSVerifier<const C: usize, NTT, P, T> {
    _r: PhantomData<NTT>,
    _p: PhantomData<P>,
    _t: PhantomData<T>,
}

impl<
        const C: usize,
        NTT: SuitableRing,
        P: DecompositionParams,
        T: TranscriptWithShortChallenges<NTT>,
    > NIFSVerifier<C, NTT, P, T>
{
    pub fn verify(
        acc: &LCCCS<C, NTT>,
        cm_i: &CCCS<C, NTT>,
        proof: &LFProof<C, NTT>,
        transcript: &mut impl TranscriptWithShortChallenges<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<LCCCS<C, NTT>, LatticefoldError<NTT>> {
        let linearized_cm_i = LFLinearizationVerifier::<_, T>::verify(
            cm_i,
            &proof.linearization_proof,
            transcript,
            ccs,
        )?;
        let decomposed_acc = LFDecompositionVerifier::<_, T>::verify::<C, P>(
            acc,
            &proof.decomposition_proof_l,
            transcript,
            ccs,
        )?;
        let decomposed_cm_i = LFDecompositionVerifier::<_, T>::verify::<C, P>(
            &linearized_cm_i,
            &proof.decomposition_proof_r,
            transcript,
            ccs,
        )?;

        let lcccs_s = {
            let mut decomposed_acc = decomposed_acc;
            let mut decomposed_cm_i = decomposed_cm_i;

            decomposed_acc.append(&mut decomposed_cm_i);

            decomposed_acc
        };

        Ok(LFFoldingVerifier::<NTT, T>::verify::<C, P>(
            &lcccs_s,
            &proof.folding_proof,
            transcript,
            ccs,
        )?)
    }
}

#[cfg(test)]
mod tests {
    use cyclotomic_rings::{challenge_set::LatticefoldChallengeSet, rings::SuitableRing};
    use rand::thread_rng;

    use crate::nifs::decomposition::{DecompositionProver, DecompositionVerifier};
    use crate::nifs::linearization::LinearizationProver;
    use crate::{
        arith::{r1cs::get_test_z_split, tests::get_test_ccs, Witness, CCCS, LCCCS},
        commitment::AjtaiCommitmentScheme,
        decomposition_parameters::{DecompositionParams, PP},
        nifs::{
            decomposition::{LFDecompositionProver, LFDecompositionVerifier},
            linearization::LFLinearizationProver,
            LFProof, NIFSProver, NIFSVerifier,
        },
        transcript::poseidon::PoseidonTranscript,
    };

    const WIT_LEN: usize = 4;
    const W: usize = WIT_LEN * PP::L;

    fn nifs_prove<RqNTT, CS>() -> (LFProof<4, RqNTT>, LCCCS<4, RqNTT>, CCCS<4, RqNTT>)
    where
        RqNTT: SuitableRing,
        CS: LatticefoldChallengeSet<RqNTT>,
    {
        let ccs = get_test_ccs::<RqNTT>(W);
        let (_, x_ccs, w_ccs) = get_test_z_split::<RqNTT>(3);
        let scheme = AjtaiCommitmentScheme::rand(&mut thread_rng());

        let wit: Witness<RqNTT> = Witness::from_w_ccs::<PP>(w_ccs);
        let cm_i: CCCS<4, RqNTT> = CCCS {
            cm: wit.commit::<4, W, PP>(&scheme).unwrap(),
            x_ccs,
        };

        let mut prover_transcript = PoseidonTranscript::<RqNTT, CS>::default();

        let (lcccs, _) = LFLinearizationProver::<_, PoseidonTranscript<RqNTT, CS>>::prove(
            &cm_i,
            &wit,
            &mut prover_transcript,
            &ccs,
        )
        .unwrap();

        let (lcccs_vec, wit_vec, decomposition_proof) =
            LFDecompositionProver::<_, PoseidonTranscript<RqNTT, CS>>::prove::<W, 4, PP>(
                &lcccs,
                &wit,
                &mut prover_transcript,
                &ccs,
                &scheme,
            )
            .unwrap();

        let verified_lcccs_vec =
            LFDecompositionVerifier::<_, PoseidonTranscript<RqNTT, CS>>::verify::<4, PP>(
                &lcccs,
                &decomposition_proof,
                &mut prover_transcript,
                &ccs,
            )
            .unwrap();

        let (_, _, proof) = NIFSProver::<4, W, RqNTT, PP, PoseidonTranscript<RqNTT, CS>>::prove(
            &verified_lcccs_vec[0],
            &wit_vec[0],
            &cm_i,
            &wit,
            &mut prover_transcript,
            &ccs,
            &scheme,
        )
        .unwrap();

        (proof, lcccs, cm_i)
    }

    fn nifs_verify<RqNTT, CS>(
        proof: LFProof<4, RqNTT>,
        lcccs: LCCCS<4, RqNTT>,
        cm_i: CCCS<4, RqNTT>,
    ) where
        RqNTT: SuitableRing,
        CS: LatticefoldChallengeSet<RqNTT>,
    {
        let ccs = get_test_ccs::<RqNTT>(W);
        let mut verifier_transcript = PoseidonTranscript::<RqNTT, CS>::default();

        let res = NIFSVerifier::<4, RqNTT, PP, PoseidonTranscript<RqNTT, CS>>::verify(
            &lcccs,
            &cm_i,
            &proof,
            &mut verifier_transcript,
            &ccs,
        );

        assert!(res.is_ok());
    }

    mod tests_pow2 {
        use super::*;
        use cyclotomic_rings::challenge_set::BinarySmallSet;
        use lattirust_ring::cyclotomic_ring::models::pow2_debug::Pow2CyclotomicPolyRingNTT;

        const Q: u64 = 17;
        const N: usize = 8;
        type RqNTT = Pow2CyclotomicPolyRingNTT<Q, N>;
        type CS = BinarySmallSet<Q, N>;

        #[test]
        fn test_nifs_eprove() {
            nifs_prove::<RqNTT, CS>();
        }

        #[test]
        fn test_nifs_verify() {
            let (proof, lcccs, cm_i) = nifs_prove::<RqNTT, CS>();
            nifs_verify::<RqNTT, CS>(proof, lcccs, cm_i);
        }
    }

    mod tests_stark {
        use super::*;
        use cyclotomic_rings::rings::{StarkChallengeSet, StarkRingNTT};

        type RqNTT = StarkRingNTT;
        type CS = StarkChallengeSet;

        #[test]
        fn test_nifs_prove() {
            nifs_prove::<RqNTT, CS>();
        }

        #[test]
        fn test_nifs_verify() {
            let (proof, lcccs, cm_i) = nifs_prove::<RqNTT, CS>();
            nifs_verify::<RqNTT, CS>(proof, lcccs, cm_i);
        }
    }

    mod tests_goldilocks {
        use super::*;
        use cyclotomic_rings::rings::{GoldilocksChallengeSet, GoldilocksRingNTT};

        type RqNTT = GoldilocksRingNTT;
        type CS = GoldilocksChallengeSet;

        #[test]
        fn test_nifs_prove() {
            nifs_prove::<RqNTT, CS>();
        }

        #[test]
        fn test_nifs_verify() {
            let (proof, lcccs, cm_i) = nifs_prove::<RqNTT, CS>();
            nifs_verify::<RqNTT, CS>(proof, lcccs, cm_i);
        }
    }

    mod tests_frog {
        use super::*;
        use cyclotomic_rings::rings::{FrogChallengeSet, FrogRingNTT};

        type RqNTT = FrogRingNTT;
        type CS = FrogChallengeSet;

        #[test]
        fn test_nifs_prove() {
            nifs_prove::<RqNTT, CS>();
        }

        #[test]
        fn test_nifs_verify() {
            let (proof, lcccs, cm_i) = nifs_prove::<RqNTT, CS>();
            nifs_verify::<RqNTT, CS>(proof, lcccs, cm_i);
        }
    }

    mod tests_babybear {
        use super::*;
        use cyclotomic_rings::rings::{BabyBearChallengeSet, BabyBearRingNTT};

        type RqNTT = BabyBearRingNTT;
        type CS = BabyBearChallengeSet;

        #[test]
        fn test_nifs_prove() {
            nifs_prove::<RqNTT, CS>();
        }

        #[test]
        fn test_nifs_verify() {
            let (proof, lcccs, cm_i) = nifs_prove::<RqNTT, CS>();
            nifs_verify::<RqNTT, CS>(proof, lcccs, cm_i);
        }
    }
}
