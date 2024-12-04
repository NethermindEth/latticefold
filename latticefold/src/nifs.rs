use ark_std::marker::PhantomData;

use cyclotomic_rings::rings::SuitableRing;
use lattirust_ring::OverField;

use crate::{
    arith::{error::CSError, Witness, CCCS, CCS, LCCCS},
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
        sanity_check::<NTT, P>(ccs)?;

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
            LFFoldingProver::<_, T>::prove::<C, P>(&lcccs, wit_s, transcript, ccs)?;

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
        sanity_check::<NTT, P>(ccs)?;

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

fn sanity_check<NTT: SuitableRing, DP: DecompositionParams>(
    ccs: &CCS<NTT>,
) -> Result<(), LatticefoldError<NTT>> {
    if ccs.m != usize::max((ccs.n - ccs.l - 1) * DP::L, ccs.m).next_power_of_two() {
        return Err(CSError::InvalidSizeBounds(ccs.m, ccs.n, DP::L).into());
    }

    Ok(())
}

#[cfg(test)]
#[macro_use]
mod tests {
    use core::fmt::Debug;

    use ark_ff::UniformRand;
    use cyclotomic_rings::{
        challenge_set::LatticefoldChallengeSet,
        rings::SuitableRing,
        rings::{GoldilocksChallengeSet, GoldilocksRingNTT},
    };

    use crate::{
        arith::{
            r1cs::{get_test_dummy_r1cs, get_test_dummy_z_split},
            Arith, Witness, CCCS, CCS,
        },
        commitment::AjtaiCommitmentScheme,
        decomposition_parameters::DecompositionParams,
        nifs::NIFSVerifier,
        transcript::poseidon::PoseidonTranscript,
    };

    use super::{
        linearization::{
            LFLinearizationProver, LFLinearizationVerifier, LinearizationProver,
            LinearizationVerifier,
        },
        NIFSProver,
    };

    #[macro_export]
    macro_rules! define_params {
        ($w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
            paste::paste! {
                #[derive(Clone)]
                struct [<DecompParamsWithB $b W $w b $b_small K $k>];

                impl DecompositionParams for [<DecompParamsWithB $b W $w b $b_small K $k>] {
                    const B: u128 = $b;
                    const L: usize = $l;
                    const B_SMALL: usize = $b_small;
                    const K: usize = $k;
                }
            }
        };
    }
    #[allow(unused_macros)]
    macro_rules! run_single_goldilocks_e2e {
        ($io:expr, $cw:expr, $w:expr, $b:expr, $l:expr, $b_small:expr, $k:expr) => {
            define_params!($w, $b, $l, $b_small, $k);
            paste::paste! {
                e2e::<$io, $cw, $w, {$w * $l}, GoldilocksChallengeSet, GoldilocksRingNTT, [<DecompParamsWithB $b W $w b $b_small K $k>]>();
            }
        };
    }

    #[test]
    fn test_e2e() {
        // X_LEN, C, WIT_LEN, B, L, B_SMALL, K
        run_single_goldilocks_e2e!(1, 39, 512, 4294967296, 2, 2, 32);
    }
    fn e2e<
        const X_LEN: usize,
        const C: usize,
        const WIT_LEN: usize,
        const W: usize,
        CS: LatticefoldChallengeSet<R> + Clone + 'static,
        R: SuitableRing,
        P: DecompositionParams,
    >() {
        let (cm_i, wit, ccs, scheme) =
            wit_and_ccs_gen::<X_LEN, C, WIT_LEN, W, P, R>(X_LEN + WIT_LEN + 1);
        verify_e2e::<C, W, P, R, CS>(&cm_i, &wit, &ccs, &scheme);
    }

    fn verify_e2e<
        const C: usize,
        const W: usize,
        P: DecompositionParams,
        R: SuitableRing,
        CS: LatticefoldChallengeSet<R> + Clone,
    >(
        cm_i: &CCCS<C, R>,
        wit: &Witness<R>,
        ccs: &CCS<R>,
        scheme: &AjtaiCommitmentScheme<C, W, R>,
    ) {
        let mut prover_transcript = PoseidonTranscript::<R, CS>::default();
        let mut verifier_transcript = PoseidonTranscript::<R, CS>::default();

        let (prover_lcccs_acc, acc_linearization_proof) = LFLinearizationProver::<
            _,
            PoseidonTranscript<R, CS>,
        >::prove(
            cm_i, wit, &mut prover_transcript, ccs
        )
        .expect("Failed to generate acc linearization proof");

        let verifier_lcccs_acc = LFLinearizationVerifier::<_, PoseidonTranscript<R, CS>>::verify(
            cm_i,
            &acc_linearization_proof,
            &mut verifier_transcript,
            ccs,
        )
        .expect("Failed to verify acc linearization");

        let (_, _, proof) = NIFSProver::<C, W, R, P, PoseidonTranscript<R, CS>>::prove(
            &prover_lcccs_acc,
            wit,
            cm_i,
            wit,
            &mut prover_transcript,
            ccs,
            scheme,
        )
        .expect("Failed to generate proof");

        NIFSVerifier::<C, R, P, PoseidonTranscript<R, CS>>::verify(
            &verifier_lcccs_acc,
            cm_i,
            &proof,
            &mut verifier_transcript,
            ccs,
        )
        .expect("Failed to verify NIFS proof");
    }

    #[allow(dead_code)]
    pub fn wit_and_ccs_gen<
        const X_LEN: usize,
        const C: usize,
        const WIT_LEN: usize,
        const W: usize,
        P: DecompositionParams,
        R: Clone + UniformRand + Debug + SuitableRing,
    >(
        r1cs_rows: usize,
    ) -> (
        CCCS<C, R>,
        Witness<R>,
        CCS<R>,
        AjtaiCommitmentScheme<C, W, R>,
    ) {
        let mut rng = ark_std::test_rng();

        let new_r1cs_rows = if P::L == 1 && (WIT_LEN > 0 && (WIT_LEN & (WIT_LEN - 1)) == 0) {
            r1cs_rows - 2
        } else {
            r1cs_rows // This makes a square matrix but is too much memory
        };
        let ccs: CCS<R> = get_test_dummy_ccs::<R, X_LEN, WIT_LEN, W>(new_r1cs_rows, P::L);
        let (one, x_ccs, w_ccs) = get_test_dummy_z_split::<R, X_LEN, WIT_LEN>();
        let mut z = vec![one];
        z.extend(&x_ccs);
        z.extend(&w_ccs);
        ccs.check_relation(&z).expect("R1CS invalid!");

        let scheme: AjtaiCommitmentScheme<C, W, R> = AjtaiCommitmentScheme::rand(&mut rng);
        let wit: Witness<R> = Witness::from_w_ccs::<P>(w_ccs);

        let cm_i: CCCS<C, R> = CCCS {
            cm: wit.commit::<C, W, P>(&scheme).unwrap(),
            x_ccs,
        };

        (cm_i, wit, ccs, scheme)
    }

    #[allow(dead_code)]
    pub fn get_test_dummy_ccs<
        R: Clone + UniformRand + Debug + SuitableRing,
        const X_LEN: usize,
        const WIT_LEN: usize,
        const W: usize,
    >(
        r1cs_rows: usize,
        L: usize,
    ) -> CCS<R> {
        let r1cs = get_test_dummy_r1cs::<R, X_LEN, WIT_LEN>(r1cs_rows);
        CCS::<R>::from_r1cs_padded(r1cs, W, L)
    }
}
