use ajtai_commit::*;
fn main() {
    let po2_np = 2;
    let num_input_polys = 2_u32.pow(po2_np) as usize;
    let field_modulus = (2_u32.pow(31) - 1) as usize; // using M31 as a placeholder
    let modulus_poly_degree = 2;
    let ajtai_input = AjtaiVecRingElems::new(num_input_polys, field_modulus, modulus_poly_degree);
    println!("Naive Ajtai commitment");
    println!("Input");
    for (i, polyring) in ajtai_input.clone().polys.into_iter().enumerate() {
        println!("poly #{}: {}", i, polyring)
    }
    let ajtai_matrix = AjtaiMatrixRingElems::new(
        num_input_polys,
        num_input_polys,
        field_modulus,
        modulus_poly_degree,
    ); // use a square matrix for the time being
    let commitment = ajtai_matrix * ajtai_input.clone();

    println!("Commitment");
    for (i, polyring) in commitment.polys.into_iter().enumerate() {
        println!("poly #{}: {}", i, polyring)
    }

    print!("Ajtai commitment using FFT");
    let ntt_domain = NTTDomain::new();
    let ajtai_evals_input = ajtai_input.evaluate(&ntt_domain);
}
