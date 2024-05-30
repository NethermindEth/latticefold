use super::Field;

pub struct RqCRT<F: Field> {
    pub crt_coeffs: Vec<F>,
}
