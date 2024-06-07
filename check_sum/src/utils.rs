use lattirust_arithmetic::ring::PolyRing;

// Takes A ring elements and evaluates it
// Given an array of coefficients, smallest coefficient first

pub fn eval_poly<R: PolyRing>(r: &R, coeffs: &[R]) -> R {
    let mut result = R::zero();

    for (i, coeff) in coeffs.iter().enumerate() {
        result += r.pow(&[i as u64]) * coeff;
    }

    result
}
