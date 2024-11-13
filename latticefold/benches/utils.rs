#![allow(non_snake_case)]
use ark_std::{fmt::Debug, UniformRand};
use cyclotomic_rings::rings::SuitableRing;
use latticefold::arith::r1cs::get_test_dummy_z_split;
use latticefold::arith::{r1cs::get_test_dummy_r1cs, Arith, Witness, CCCS, CCS};
use latticefold::commitment::AjtaiCommitmentScheme;
use latticefold::decomposition_parameters::DecompositionParams;
use rand::thread_rng;

pub fn wit_and_ccs_gen<
    const X_LEN: usize,
    const C: usize, // rows
    const WIT_LEN: usize,
    const W: usize, // columns
    P: DecompositionParams,
    R: Clone + UniformRand + Debug + SuitableRing + for<'a> std::ops::AddAssign<&'a R>,
>(
    r1cs_rows: usize,
) -> (
    CCCS<C, R>,
    Witness<R>,
    CCS<R>,
    AjtaiCommitmentScheme<C, W, R>,
) {
    println!("Step 1: Calculating number of R1CS rows");
    let new_r1cs_rows = if P::L == 1 && (WIT_LEN > 0 && (WIT_LEN & (WIT_LEN - 1)) == 0) {
        r1cs_rows - 2
    } else {
        r1cs_rows // This makes a square matrix but is too much memory
    };

    println!("Step 2: Generating test dummy CCS");
    let ccs: CCS<R> = get_test_dummy_ccs::<R, X_LEN, WIT_LEN, W>(new_r1cs_rows);

    println!("Step 3: Getting test dummy z split");
    let (one, x_ccs, w_ccs) = get_test_dummy_z_split::<R, X_LEN, WIT_LEN>();

    println!("Step 4: Building z vector");
    let mut z = vec![one];
    z.extend(&x_ccs);
    z.extend(&w_ccs);

    println!("Step 5: Checking R1CS relation");
    match ccs.check_relation(&z) {
        Ok(_) => println!("R1CS valid!"),
        Err(e) => println!("R1CS invalid: {:?}", e),
    }

    println!("Step 6: Generating random Ajtai commitment scheme");
    let scheme: AjtaiCommitmentScheme<C, W, R> = AjtaiCommitmentScheme::rand(&mut thread_rng());

    println!("Step 7: Creating witness from w_ccs");
    let wit: Witness<R> = Witness::from_w_ccs::<P>(&w_ccs);

    println!("Step 8: Creating CCCS with commitment and x_ccs");
    let cm_i: CCCS<C, R> = CCCS {
        cm: wit.commit::<C, W, P>(&scheme).unwrap(),
        x_ccs,
    };

    println!("Step 9: Returning generated values");
    (cm_i, wit, ccs, scheme)
}
pub fn get_test_dummy_ccs<
    R: Clone + UniformRand + Debug + SuitableRing + for<'a> std::ops::AddAssign<&'a R>,
    const X_LEN: usize,
    const WIT_LEN: usize,
    const W: usize,
>(
    r1cs_rows: usize,
) -> CCS<R> {
    let r1cs = get_test_dummy_r1cs::<R, X_LEN, WIT_LEN>(r1cs_rows);
    CCS::<R>::from_r1cs(r1cs, W)
}
