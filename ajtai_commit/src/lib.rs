use std::{
    ops::{Add, Mul},
    str::FromStr,
};

use qfall_math::{
    integer::Z,
    integer_mod_q::{Modulus, ModulusPolynomialRingZq, PolynomialRingZq, Zq},
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
            polys: sample_rand_vec_polys(num_polys, field_modulus, modulus_poly),
        }
    }

    pub fn evaluate(self, domain: &NTTDomain) -> AjtaiEvalsVec {
        let mut evals = Vec::new();
        let polys = self.polys;
        for poly in polys {
            evals.push(ntt(domain, poly));
        }
        AjtaiEvalsVec { vec: evals }
    }
}

/// This contains the root of unity and the twiddle factors
pub struct NTTDomain {}
impl NTTDomain {
    pub fn new() -> Self {
        todo!()
    }
}
/// This contains the inverse of the NTTDomain
pub struct INTTDomain {}
impl INTTDomain {
    pub fn new() -> Self {
        todo!()
    }
}

fn ntt(domain: &NTTDomain, poly: PolynomialRingZq) -> PolyEvaluation {
    todo!()
}

fn intt(
    domain: &INTTDomain,
    evals: PolyEvaluation,
    modulus: ModulusPolynomialRingZq,
) -> PolynomialRingZq {
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
    field_modulus: usize,
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
        matrix.push(sample_rand_vec_polys(
            columns,
            field_modulus,
            modulus_poly.clone(),
        ));
    }
    matrix
}

#[derive(Clone)]
pub struct PolyEvaluation {
    evals: Vec<Zq>,
}

impl Add for PolyEvaluation {
    type Output = PolyEvaluation;
    fn add(self, rhs: Self) -> Self::Output {
        assert_eq!(self.evals.len(), rhs.evals.len());
        let result = self
            .evals
            .iter()
            .zip(rhs.evals.iter())
            .map(|(l, r)| l + r)
            .collect::<Vec<_>>();
        PolyEvaluation { evals: result }
    }
}

impl Mul for PolyEvaluation {
    type Output = PolyEvaluation;
    fn mul(self, rhs: Self) -> Self::Output {
        assert_eq!(self.evals.len(), rhs.evals.len());
        let result = self
            .evals
            .iter()
            .zip(rhs.evals.iter())
            .map(|(l, r)| l * r)
            .collect::<Vec<_>>();
        PolyEvaluation { evals: result }
    }
}

pub struct AjtaiEvalsVec {
    vec: Vec<PolyEvaluation>,
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
            .map(|poly_evals| intt(domain, poly_evals.clone(), modulus.clone()))
            .collect::<Vec<_>>();
        AjtaiVecRingElems { polys: result }
    }
}

pub struct AjtaiEvalsMatrix {
    mat: Vec<Vec<PolyEvaluation>>,
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
                let poly_eval = PolyEvaluation { evals };
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
        let zero_evals = PolyEvaluation {
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
