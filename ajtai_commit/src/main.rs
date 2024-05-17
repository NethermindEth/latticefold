use std::str::FromStr;

use ajtai_commit::*;
use qfall_math::{
    integer::Z,
    integer_mod_q::{ModulusPolynomialRingZq, Zq},
};
fn main() {
    // Test STARK prime
    let stark_prime_string =
        "3618502788666131213697322783095070105623107215331596699973092056135872020481";
    let two_in_stark_string = format!("2 mod {}", stark_prime_string);
    let two_in_stark = Zq::from_str(&two_in_stark_string).unwrap();
    let a_in_stark_string = format!(
        "1809251394333065606848661391547535052811553607665798349986546028067936010238 mod {}",
        stark_prime_string
    );
    let a_in_stark = Zq::from_str(&a_in_stark_string).unwrap();
    println!("Two: {}", &two_in_stark);
    println!("A: {}", &a_in_stark);
    let two_times_a = &two_in_stark * &a_in_stark;
    println!("Two time A:\n{}", two_times_a);
    let zero_string = format!("0 mod {}", stark_prime_string);
    let zero = Zq::from_str(&zero_string).unwrap();
    let five = &zero - &two_times_a;
    println!("Should be five: {}", five);
    //test Ajtai

    let po2_np = 2;
    let num_input_polys = 1 << po2_np as usize;
    let field_modulus = (15 * (1 << 27) + 1) as usize; // using M31 as a placeholder
                                                       // rou 27 for babybear
    let modulus_poly_degree = 2;

    let mut desc_mod_poly = format!("{}  1", modulus_poly_degree + 1);
    for _ in 0..modulus_poly_degree - 1 {
        desc_mod_poly.push_str(" 0");
    }
    desc_mod_poly.push_str(" 1");
    desc_mod_poly.push_str(&format!(" mod {}", field_modulus));

    let modulus_poly = ModulusPolynomialRingZq::from_str(&desc_mod_poly).unwrap();

    let ajtai_input = AjtaiVecRingElems::new(num_input_polys, field_modulus, modulus_poly.clone());
    println!("Naive Ajtai commitment");
    println!("Input");
    for (i, polyring) in ajtai_input.clone().polys.into_iter().enumerate() {
        println!("poly #{}: {}", i, polyring)
    }
    let ajtai_matrix = AjtaiMatrixRingElems::new(
        num_input_polys,
        num_input_polys,
        field_modulus,
        modulus_poly.clone(),
    ); // use a square matrix for the time being
    let commitment = ajtai_matrix.naive_commit(ajtai_input.clone());

    println!("Naive Commitment");
    for (i, polyring) in commitment.polys.into_iter().enumerate() {
        println!("poly #{}: {}", i, polyring)
    }

    print!("Ajtai commitment using FFT");
    let rou = Zq::from_str(format!("{} mod {}", 27, field_modulus).as_str()).unwrap();
    let ntt_domain = NTTDomain::new(
        rou,
        (modulus_poly_degree as u32).next_power_of_two() as usize,
    );
    let ajtai_matrix = AjtaiEvalsMatrix::sample_rand_mat_evals(
        num_input_polys,
        num_input_polys,
        field_modulus,
        (modulus_poly_degree + 1 as u32).next_power_of_two() as usize, //check that num of evals is the same as in domain
    );
    let ajtai_evals_input = ajtai_input.evaluate(&ntt_domain);
    let commitment = ajtai_matrix * ajtai_evals_input;
    let intt_domain = INTTDomain::new(&ntt_domain);
    let commitment = commitment.make_coeffs(&intt_domain, modulus_poly.clone());
    println!("Commitment");
    for (i, polyring) in commitment.polys.into_iter().enumerate() {
        println!("poly #{}: {}", i, polyring)
    }
}
