use std::str::FromStr;

use ajtai_commit::*;
use qfall_math::integer_mod_q::{ModulusPolynomialRingZq, Zq};
fn main() {
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
