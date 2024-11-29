#![allow(non_snake_case, clippy::upper_case_acronyms)]
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::marker::PhantomData;
use lattirust_poly::mle::DenseMultilinearExtension;
use lattirust_ring::{OverField, Ring};

use crate::{
    arith::{Witness, CCS, LCCCS},
    ark_base::*,
    commitment::AjtaiCommitmentScheme,
    commitment::Commitment,
    decomposition_parameters::DecompositionParams,
    nifs::error::DecompositionError,
    transcript::Transcript,
};
use cyclotomic_rings::rings::SuitableRing;

/// The proof structure of the decomposition subprotocol.
#[derive(Clone, Debug, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct DecompositionProof<const C: usize, NTT: Ring> {
    /// CCS-linearization evaluation claims w.r.t. decomposed witnesses.
    ///
    /// After a run of the decomposition subprotocol prover this field contains
    /// `P::K` vectors of length `ccs.t`. `u_s[i][j]` is such that
    /// $$
    /// \text{u\\_s[i][j]} =  \sum\_{\vec{\mathbf{b}} \in \\{0,1\\}^{\log n\_c}}
    /// \text{mle}[\text{ccs.M[j]}](\vec{\mathbf{x}}, \vec{\mathbf{b}}) \cdot \text{mle}\[\mathbf{z}\_{i}\](\vec{\mathbf{b}})
    /// $$
    /// where $\mathbf{z}_i$ is obtained by concatenating `x_s[i] || w_s[i]` (`w_s` are the decomposed witnesses).
    pub u_s: Vec<Vec<NTT>>,
    /// Evaluation claims about rows of $\hat{f}$-matrices of decomposed witnesses.
    ///
    /// After a run of the decomposition subprotocol prover this field contains
    /// `P::K` vectors of length `NTT::CoefficientRepresentation::dimension() / NTT::dimension()`. `v_s[i][j]` is such that
    /// $$
    /// \text{v\\_s[i][j]}= \text{mle}[\text{w\\_s[i].f\\_hat[j]}] (\mathbf{r}).
    /// $$
    /// where $\mathbf{r}$ is the evaluation point from the LCCCS.
    pub v_s: Vec<Vec<NTT>>,
    /// Decomposed public parts of the statement-witness vectors `z_s`.
    ///
    /// It is expensive to compute them
    /// on the verifier's side thus the prover computes them by itself and sends to the verifier.
    pub x_s: Vec<Vec<NTT>>,
    /// Commitments to the decomposed witnesses.
    pub y_s: Vec<Commitment<C, NTT>>,
}

pub trait DecompositionProver<NTT: SuitableRing, T: Transcript<NTT>> {
    fn prove<const W: usize, const C: usize, P: DecompositionParams>(
        cm_i: &LCCCS<C, NTT>,
        wit: &Witness<NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
        scheme: &AjtaiCommitmentScheme<C, W, NTT>,
    ) -> Result<
        (
            Vec<Vec<DenseMultilinearExtension<NTT>>>,
            Vec<LCCCS<C, NTT>>,
            Vec<Witness<NTT>>,
            DecompositionProof<C, NTT>,
        ),
        DecompositionError,
    >;
}

pub trait DecompositionVerifier<NTT: OverField, T: Transcript<NTT>> {
    fn verify<const C: usize, P: DecompositionParams>(
        cm_i: &LCCCS<C, NTT>,
        proof: &DecompositionProof<C, NTT>,
        transcript: &mut impl Transcript<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<Vec<LCCCS<C, NTT>>, DecompositionError>;
}

pub struct LFDecompositionProver<NTT, T> {
    _ntt: PhantomData<NTT>,
    _t: PhantomData<T>,
}

pub struct LFDecompositionVerifier<NTT, T> {
    _ntt: PhantomData<NTT>,
    _t: PhantomData<T>,
}
