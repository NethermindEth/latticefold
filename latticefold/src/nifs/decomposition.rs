use lattirust_arithmetic::{
    balanced_decomposition::decompose_balanced_polyring,
    challenge_set::latticefold_challenge_set::OverField,
    mle::DenseMultilinearExtension,
};

use ark_std::iterable::Iterable;

use super::{error::{DecompositionError, NotSatisfiedError}, NIFSProver, NIFSVerifier};

use crate::{
    arith::{Witness, CCS, LCCCS},
    transcript::Transcript,
};

#[derive(Clone)]
pub struct DecompositionProof<R: OverField> {
    pub u_s: Vec<Vec<R>>,
    pub v_s: Vec<R>,
    pub x_s: Vec<Vec<R>>,
    pub y_s: Vec<Vec<R>>,
}

pub trait DecompositionProver<R: OverField, T: Transcript<R>> {
    type Proof: Clone;
    type Error: std::error::Error;

    fn prove(
        cm_i: &LCCCS<R>,
        wit: &Witness<R>,
        transcript: &mut impl Transcript<R>,
        ccs: &CCS<R>,
    ) -> Result<(Vec<LCCCS<R>>, Vec<Witness<R>>, Self::Proof), Self::Error>;
}

pub trait DecompositionVerifier<R: OverField, T: Transcript<R>> {
    type Prover: DecompositionProver<R, T>;
    type Error = <Self::Prover as DecompositionProver<R, T>>::Error;

    fn verify(
        cm_i: &LCCCS<R>,
        proof: &<Self::Prover as DecompositionProver<R, T>>::Proof,
        transcript: &mut impl Transcript<R>,
        ccs: &CCS<R>,
    ) -> Result<Vec<LCCCS<R>>, Self::Error>;
}

impl<R: OverField, T: Transcript<R>> DecompositionProver<R, T> for NIFSProver<R, T> {
    type Proof = DecompositionProof<R>;
    type Error = DecompositionError<R>;

    fn prove(
        cm_i: &LCCCS<R>,
        wit: &Witness<R>,
        transcript: &mut impl Transcript<R>,
        ccs: &CCS<R>,
    ) -> Result<(Vec<LCCCS<R>>, Vec<Witness<R>>, DecompositionProof<R>), DecompositionError<R>> {
        let mut u_s: Vec<Vec<R>> = Vec::new();
        let mut v_s: Vec<R> = Vec::new();
        let mut x_s: Vec<Vec<R>> = Vec::new();
        let mut y_s: Vec<Vec<R>> = Vec::new();

        for i in 0..cm_i.r_arr.len() {
            let f_i = compute_fi(wit.w_ccs.clone(), cm_i.r_arr[i], ccs.s);
            let l_f_i: Vec<R> = f_i.clone().iter().map(|f| l_function(f)).collect();
            let v_i = mle(&wit.f_arr, &cm_i.r_arr); // Use the correct field
            let u_i = inner(&f_i, &cm_i.r_arr);
            let x_w_i: Vec<R> = f_i.clone().iter().map(|f| l_w_function(f, ccs.s)).flatten().collect();

            y_s.push(l_f_i);
            v_s.push(v_i);
            u_s.push(u_i);
            x_s.push(x_w_i);
        }

        let proof = DecompositionProof {
            u_s,
            v_s,
            x_s,
            y_s,
        };

        // Placeholder for resulting LCCCS and Witnesses
        let lcccs_list: Vec<LCCCS<R>> = vec![];
        let witnesses: Vec<Witness<R>> = vec![];

        Ok((lcccs_list, witnesses, proof))
    }
}

impl<R: OverField, T: Transcript<R>> DecompositionVerifier<R, T> for NIFSVerifier<R, T> {
    type Prover = NIFSProver<R, T>;

    fn verify(
        cm_i: &LCCCS<R>,
        proof: &<Self::Prover as DecompositionProver<R, T>>::Proof,
        transcript: &mut impl Transcript<R>,
        ccs: &CCS<R>,
    ) -> Result<Vec<LCCCS<R>>, DecompositionError<R>> {
        for i in 0..proof.u_s.len() {
            let computed_v = mle(&proof.u_s[i], &cm_i.r_arr);
            if computed_v != proof.v_s[i] {
                return Err(DecompositionError::VerificationFailed(NotSatisfiedError));
            }

            let computed_x: Vec<R> = proof.x_s[i].iter().map(|x| l_w_function(x, ccs.s)).flatten().collect();
            if computed_x != proof.x_s[i] {
                return Err(DecompositionError::VerificationFailed(NotSatisfiedError));
            }

            let computed_y: Vec<R> = proof.y_s[i].iter().map(|y| l_function(y)).collect();
            if computed_y != proof.y_s[i] {
                return Err(DecompositionError::VerificationFailed(NotSatisfiedError));
            }
        }

        // Placeholder for resulting LCCCS
        let lcccs_list: Vec<LCCCS<R>> = vec![];

        Ok(lcccs_list)
    }
}

fn compute_fi<R: OverField>(f: Vec<R>, b: R, k: usize) -> Vec<R> {
    f.iter().map(|x| decompose_balanced_polyring(x, b.to_u128().unwrap(), Some(k)).iter().map(|e| e.clone()).collect::<Vec<R>>()).flatten().collect()
}

fn tensor<R: OverField>(r: &Vec<R>) -> Vec<R> {
    let log_m = r.len();
    let num_combinations = 1 << log_m; 

    let mut result = Vec::with_capacity(num_combinations);

    for i in 0..num_combinations {
        let mut product = R::one(); 
        for (j, &rj) in r.iter().enumerate() {
            if (i & (1 << j)) != 0 {
                product = product * rj;
            } else {
                product = product * (R::one() - rj);
            }
        }
        result.push(product);
    }

    result
}

fn mle<R: OverField>(f: &Vec<R>, r: &Vec<R>) -> R {
    let mle = DenseMultilinearExtension::from_evaluations_vec(r.len(), f.clone());
    mle.evaluate(r.as_slice()).unwrap()
}

fn l_function<R: OverField>(f: &R) -> R {
    // Implementation of L(f) (Ajtai Commitment)
    todo!()
}

fn l_w_function<R: OverField>(f: &R, n_initial: usize) -> Vec<R> {
    // Implementation of Lw(f) between 1 and n_initial
    todo!()
}

fn inner<R: OverField>(fi: &Vec<R>, r: &Vec<R>) -> Vec<R> {
    let lw = l_w_function(&fi[0], r.len());
    let tensor_r = tensor(r);
    let m_tensor_r: Vec<R> = lw.iter().zip(tensor_r.iter()).map(|(lw_i, t_r)| *lw_i * *t_r).collect();
    m_tensor_r
}
