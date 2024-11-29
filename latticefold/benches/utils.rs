#![allow(non_snake_case)]
use ark_std::{fmt::Debug, UniformRand};
use cyclotomic_rings::rings::SuitableRing;
use latticefold::arith::r1cs::get_test_dummy_z_split;
use latticefold::arith::{r1cs::get_test_dummy_r1cs, Arith, Witness, CCCS, CCS};
use latticefold::commitment::AjtaiCommitmentScheme;
use latticefold::decomposition_parameters::DecompositionParams;

#[allow(dead_code)]
pub fn wit_and_ccs_gen<
    const X_LEN: usize,
    const C: usize, // rows
    const WIT_LEN: usize,
    const W: usize, // columns
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
    match ccs.check_relation(&z) {
        Ok(_) => (),
        Err(e) => panic!("R1CS invalid: {:?}", e),
    }

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
