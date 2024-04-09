use std::{collections::HashMap, fmt::Debug};

use super::Field;

use super::PrimeCyclotomicFieldElement;

#[derive(Clone)]
pub struct Number<F: Field> {
    order: usize,
    pub coeffs: HashMap<usize, F>,
}

pub fn print_cyclotomic_element<F: Field>(z: &Number<F>) -> String {
    let mut str_list: Vec<String> = vec![];
    let mut exp = 0;
    while &exp != &z.order {
        let zero = F::zero();
        let coeff = z.coeffs.get(&exp).unwrap_or(&zero);
        if !coeff.is_zero() {
            str_list.push(String::from(
                format!("{}*E({})^{}", coeff, z.order, exp).as_str(),
            ));
        }
        exp += 1;
    }
    "(".to_string() + &str_list.join(" + ") + ")"
}

impl<F: Field> Debug for Number<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "CyclotomicNumber ({})", print_cyclotomic_element(self))
    }
}

impl<F> Number<F>
where
    F: Field,
{
    pub fn new(order: usize, coeffs: &HashMap<usize, F>) -> Number<F> {
        // There should be a restriction on the coeffs keys depending on the order
        Number {
            order,
            coeffs: coeffs.clone(),
        }
    }

    pub fn increase_order_to(z: &mut Self, new_order: usize) {
        let mut new_coeffs = HashMap::default();
        for (exp, coeff) in &z.coeffs {
            new_coeffs.insert(new_order * exp / z.order, *coeff);
        }
        z.order = new_order;
        z.coeffs = new_coeffs;
    }
}

impl<F> PrimeCyclotomicFieldElement<F> for Number<F>
where
    F: Field,
{
    fn e(n: usize, k: usize) -> Self {
        let mut coeffs = HashMap::<usize, F>::default();
        coeffs.insert(k, F::one());
        Number::new(n, &coeffs)
    }

    fn scalar_mul(&mut self, scalar: &F) -> &mut Self {
        // TODO resolve clones
        let mut result = self.clone();
        for (_, coeff) in result.coeffs.iter_mut() {
            *coeff = coeff.clone() * scalar.clone();
        }
        *self = result;
        self
    }

    fn zero_order(n: usize) -> Self {
        todo!()
    }

    fn one_order(n: usize) -> Self {
        todo!()
    }
}

#[cfg(test)]
mod cyclotomic_prime_tests {
    use rand::random;

    use super::*;
    #[test]
    fn print_test() {
        let mut coeffs = HashMap::<usize, fields::M31>::default();
        for i in 0..5 {
            let random_u64: u64 = random();
            let m31_rand = fields::M31::reduce(random_u64);
            coeffs.insert(i, m31_rand);
            
        }
        let random_cyclotomic_element = Number::new(5, &coeffs);
        println!("{:?}", random_cyclotomic_element);
    }
}
