use std::{
    ops::{Add, Mul},
    str::FromStr,
};

use qfall_math::{
    integer::Z,
    integer_mod_q::{ModulusPolynomialRingZq, PolynomialRingZq, Zq},
};

#[derive(Debug,Clone)]
pub struct AjtaiVecRingElems {
    pub polys: Vec<PolynomialRingZq>,
}

impl AjtaiVecRingElems {
    pub fn new(
        num_polys: usize,
        field_modulus: usize,
        modulus_poly_degree: usize,
    ) -> AjtaiVecRingElems {
        AjtaiVecRingElems {
            polys: sample_rand_vec_polys(num_polys, field_modulus, modulus_poly_degree),
        }
    }

    pub fn evaluate(self, root_of_unity: Zq, num_evals: usize) -> Vec<PolyEvaluation> {
        let mut evals = Vec::new();
        let polys = self.polys;
        for poly in polys {
            evals.push(ntt(poly));
        }
        evals
    }
}

fn ntt(rou: Zq, domain_size: usize, poly: PolynomialRingZq) -> PolyEvaluation {
    todo!()
}

fn intt(irou: Zq, evals: PolyEvaluation, modulus: ModulusPolynomialRingZq) -> PolynomialRingZq {
    todo!()
}

struct AjtaiMatrixRingElems {
    mat: Vec<Vec<PolynomialRingZq>>,
}

impl AjtaiMatrixRingElems {
    pub fn new(
        rows: usize,
        columns: usize,
        field_modulus: usize,
        modulus_poly_degree: usize,
    ) -> Self {
        AjtaiMatrixRingElems {
            mat: sample_rand_mat_polys(rows, columns, field_modulus, modulus_poly_degree),
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
    modulus_poly_degree: usize,
) -> Vec<PolynomialRingZq> {
    let mut desc_mod_poly = format!("{}  1", modulus_poly_degree + 1);
    for _ in 0..modulus_poly_degree - 1 {
        desc_mod_poly.push_str(" 0");
    }
    desc_mod_poly.push_str(" 1");
    desc_mod_poly.push_str(&format!(" mod {}", field_modulus));

    print!("desc: {}", desc_mod_poly);
    let modulus_poly = ModulusPolynomialRingZq::from_str(&desc_mod_poly).unwrap();

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
    modulus_poly_degree: usize,
) -> Vec<Vec<PolynomialRingZq>> {
    let mut matrix = Vec::new();
    for _ in 0..rows {
        matrix.push(sample_rand_vec_polys(
            columns,
            field_modulus,
            modulus_poly_degree,
        ));
    }
    matrix
}

struct PolyEvaluation {
    evals: Vec<Zq>,
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

struct EvalsVec {
    vec: Vec<PolyEvaluation>,
}

struct EvalsMatrix {
    mat: Vec<Vec<PolyEvaluation>>,
}
