use std::{
    fmt::format,
    ops::{Add, Mul},
    str::FromStr,
};

use qfall_math::{
    integer::Z,
    integer_mod_q::{Modulus, ModulusPolynomialRingZq, PolyOverZq, PolynomialRingZq, Zq},
    traits::{GetCoefficient, Pow},
};

#[derive(Debug, Clone)]
pub struct AjtaiVecRingElems {
    pub polys: Vec<PolynomialRingZq>,
}

impl AjtaiVecRingElems {
    pub fn new(
        num_polys: usize,
        field_modulus: usize,
        modulus_poly: ModulusPolynomialRingZq,
    ) -> AjtaiVecRingElems {
        AjtaiVecRingElems {
            polys: sample_rand_vec_polys(num_polys, modulus_poly),
        }
    }

    pub fn evaluate(self, domain: &NTTDomain) -> AjtaiEvalsVec {
        let mut evals = Vec::new();
        let polys = self.polys;
        for poly in polys {
            evals.push(domain.ntt(poly));
        }
        AjtaiEvalsVec { vec: evals }
    }
}

/// This contains the root of unity and the twiddle factors
pub struct NTTDomain {
    omega_powers: Vec<Zq>,
}
impl NTTDomain {
    pub fn new(rou: Zq, domain_len: usize) -> Self {
        let q = rou.get_mod().to_string().parse::<u32>().unwrap();
        let resize_power = q / domain_len as u32;
        let new_rou = rou.pow(resize_power).unwrap();
        let mut omega_power = Zq::from_str(format!("1 mod {}", q).as_str()).unwrap();
        let mut omega_powers = Vec::new();
        for _ in 0..domain_len {
            omega_powers.push(omega_power.clone());
            omega_power = omega_power * new_rou.clone();
        }
        NTTDomain { omega_powers }
    }

    fn ntt(&self, poly: PolynomialRingZq) -> RqNTT {
        let mut coefficients = Vec::new();
        for i in 0..(poly.get_degree() + 1) {
            coefficients.push(poly.get_coeff(i).unwrap());
        }
        let modulus = poly.get_mod().get_q();
        let mut coeffs = coefficients
            .iter()
            .map(|coeff| Zq::from_z_modulus(&Z::from(coeff), modulus.clone()))
            .collect::<Vec<_>>(); //remove clones
        coeffs.resize(
            (coeffs.len() + 1).next_power_of_two(), // Currently Restricted to domain size of a power of two
            Zq::from_z_modulus(&Z::from(0), modulus.clone()),
        );
        let mut coeff_slice = coeffs.clone();
        radix2ntt(self.omega_powers.as_slice(), &mut coeff_slice); //remove clone
        return RqNTT { evals: coeff_slice };
    }
}
/// This contains the inverse of the NTTDomain
pub struct INTTDomain {
    inv_omega_powers: Vec<Zq>,
}
impl INTTDomain {
    pub fn new(domain: &NTTDomain) -> Self {
        let omega_inv_powers = domain
            .omega_powers
            .iter()
            .map(|w_power| w_power.inverse().unwrap())
            .collect::<Vec<_>>();
        INTTDomain {
            inv_omega_powers: omega_inv_powers,
        }
    }

    pub fn intt(&self, evals: RqNTT, modulus: ModulusPolynomialRingZq) -> PolynomialRingZq {
        let mut evals = evals.evals;

        radix2ntt(self.inv_omega_powers.as_slice(), &mut evals);
        let q = evals[0].get_mod();
        let n = Zq::from_z_modulus(&Z::from(evals.len() as u32), q.clone());
        let inv_n = n.inverse().unwrap(); // resolve this unwrap, tho it should exist

        for i in 0..evals.len() {
            evals[i] = evals[i].clone() * inv_n.clone();
        }

        let coeffs = evals.clone();
        let coeffs_string = coeffs
            .iter()
            .map(|c| c.get_value().to_string())
            .collect::<Vec<_>>()
            .join(" ");
        let poly_format = format!("{}  {} mod {}", evals.len(), coeffs_string, q.to_string());
        let poly = PolyOverZq::from_str(&poly_format.as_str()).unwrap();

        PolynomialRingZq::from((&poly, &modulus)) // from() reduces the polynomial automatically
    }
}

fn radix2ntt(omega_powers: &[Zq], coeffs: &mut [Zq]) {
    let n = coeffs.len();
    let q = omega_powers.first().unwrap().get_mod();
    if n == 1 {
        return;
    }

    let zero = Zq::from_z_modulus(&Z::from(0), q);
    let mut coeffs_even = vec![zero.clone(); n / 2];
    let mut coeffs_odd = vec![zero.clone(); n / 2];

    for i in 0..n / 2 {
        coeffs_even[i] = coeffs[2 * i].clone();
        coeffs_odd[i] = coeffs[2 * i + 1].clone(); // Remove clones
    }

    let half_omegas = omega_powers
        .iter()
        .step_by(2)
        .map(|z| z.clone())
        .collect::<Vec<_>>();
    radix2ntt(half_omegas.as_slice(), &mut coeffs_even);
    radix2ntt(half_omegas.as_slice(), &mut coeffs_odd);

    for i in 0..n / 2 {
        let t = coeffs_odd[i].clone() * omega_powers[i].clone();
        coeffs[i] = coeffs_even[i].clone() + t.clone();
        coeffs[i + n / 2] = coeffs_even[i].clone() - t;
    }
}

fn intt(domain: &INTTDomain, evals: RqNTT, modulus: ModulusPolynomialRingZq) -> PolynomialRingZq {
    todo!()
}

pub struct AjtaiMatrixRingElems {
    mat: Vec<Vec<PolynomialRingZq>>,
}

impl AjtaiMatrixRingElems {
    pub fn new(
        rows: usize,
        columns: usize,
        field_modulus: usize,
        modulus_poly: ModulusPolynomialRingZq,
    ) -> Self {
        AjtaiMatrixRingElems {
            mat: sample_rand_mat_polys(rows, columns, field_modulus, modulus_poly),
        }
    }
    pub fn naive_commit(self, to_commit: AjtaiVecRingElems) -> AjtaiVecRingElems {
        self * to_commit
    }
}
impl Mul<AjtaiVecRingElems> for AjtaiMatrixRingElems {
    type Output = AjtaiVecRingElems;

    fn mul(self, rhs: AjtaiVecRingElems) -> AjtaiVecRingElems {
        let mut output = Vec::with_capacity(self.mat.len());
        for row in self.mat {
            let to_sum = row
                .iter()
                .zip(rhs.polys.iter())
                .map(|(m, v)| m * v)
                .collect::<Vec<_>>(); // change to fold
            let (sum, rest) = to_sum.split_first().unwrap();
            let mut sum = sum.clone();
            for poly in rest {
                sum = sum + poly;
            }
            output.push(sum);
        }
        AjtaiVecRingElems { polys: output }
    }
}

fn sample_rand_vec_polys(
    num_polys: usize,
    modulus_poly: ModulusPolynomialRingZq,
) -> Vec<PolynomialRingZq> {
    let mut polys_vec = Vec::new();
    for _ in 0..num_polys {
        let rand_poly = PolynomialRingZq::sample_uniform(modulus_poly.clone());
        polys_vec.push(rand_poly);
    }
    polys_vec
}

fn sample_rand_mat_polys(
    rows: usize,
    columns: usize,
    field_modulus: usize,
    modulus_poly: ModulusPolynomialRingZq,
) -> Vec<Vec<PolynomialRingZq>> {
    let mut matrix = Vec::new();
    for _ in 0..rows {
        matrix.push(sample_rand_vec_polys(columns, modulus_poly.clone()));
    }
    matrix
}

/// Representation of an element of a ring where the coefficients are elements of a prime field
/// and \zeta_m powers maps to powers of the root of unity following usual NTT
#[derive(Clone)]
pub struct RqNTT {
    evals: Vec<Zq>,
}

impl Add for RqNTT {
    type Output = RqNTT;
    fn add(self, rhs: Self) -> Self::Output {
        assert_eq!(self.evals.len(), rhs.evals.len());
        let result = self
            .evals
            .iter()
            .zip(rhs.evals.iter())
            .map(|(l, r)| l + r)
            .collect::<Vec<_>>();
        RqNTT { evals: result }
    }
}

impl Mul for RqNTT {
    type Output = RqNTT;
    fn mul(self, rhs: Self) -> Self::Output {
        assert_eq!(self.evals.len(), rhs.evals.len());
        let result = self
            .evals
            .iter()
            .zip(rhs.evals.iter())
            .map(|(l, r)| l * r)
            .collect::<Vec<_>>();
        RqNTT { evals: result }
    }
}

pub struct AjtaiEvalsVec {
    vec: Vec<RqNTT>,
}

impl AjtaiEvalsVec {
    pub fn make_coeffs(
        self,
        domain: &INTTDomain,
        modulus: ModulusPolynomialRingZq,
    ) -> AjtaiVecRingElems {
        let result = self
            .vec
            .iter()
            .map(|poly_evals| domain.intt(poly_evals.clone(), modulus.clone()))
            .collect::<Vec<_>>();
        AjtaiVecRingElems { polys: result }
    }
}

pub struct AjtaiEvalsMatrix {
    mat: Vec<Vec<RqNTT>>,
}

impl AjtaiEvalsMatrix {
    pub fn sample_rand_mat_evals(
        rows: usize,
        columns: usize,
        field_modulus: usize,
        num_evals: usize,
    ) -> AjtaiEvalsMatrix {
        let modulus = Modulus::from_str(format!("{}", field_modulus).as_str()).unwrap();
        let mut mat = Vec::new();
        for _ in 0..rows {
            let mut row = Vec::new();
            for _ in 0..columns {
                let evals = (0..num_evals)
                    .map(|_| Zq::sample_uniform(modulus.clone()).unwrap())
                    .collect::<Vec<_>>();
                let poly_eval = RqNTT { evals };
                row.push(poly_eval);
            }
            mat.push(row);
        }
        AjtaiEvalsMatrix { mat }
    }
}

impl Mul<AjtaiEvalsVec> for AjtaiEvalsMatrix {
    type Output = AjtaiEvalsVec;

    fn mul(self, rhs: AjtaiEvalsVec) -> Self::Output {
        let lhs_mat = self.mat;
        let rhs_vec = rhs.vec;
        let modulus = rhs_vec.first().unwrap().evals.first().unwrap().get_mod();
        let zero = Zq::from_z_modulus(&Z::from(0), modulus);
        let zero_evals = RqNTT {
            evals: vec![zero; rhs_vec.first().unwrap().evals.len()],
        };
        let mut result = vec![zero_evals; lhs_mat.len()];

        for (i, row) in lhs_mat.iter().enumerate() {
            for (j, val) in row.into_iter().enumerate() {
                result[i] = result[i].clone() + val.clone() * rhs_vec[j].clone();
                // remove clones
            }
        }

        AjtaiEvalsVec { vec: result }
    }
}
