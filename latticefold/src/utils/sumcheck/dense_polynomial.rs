// Copyright (c) 2023 Espresso Systems (espressosys.com)
// This file is part of the HyperPlonk library.

// Adapted for rings by Nethermind

//! This module defines our main mathematical object `DensePolynomial`; and
//! various functions associated with it.

use ark_serialize::CanonicalSerialize;
use ark_std::{
    cfg_iter_mut, end_timer,
    rand::{Rng, RngCore},
    start_timer,
    string::ToString,
    vec::*,
};
use ark_std::{cmp::max, marker::PhantomData};
use lattirust_poly::mle::DenseMultilinearExtension;
use lattirust_poly::polynomials::{random_mle_list, random_zero_mle_list, ArithErrors, RefCounter};

use lattirust_ring::Ring;
#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct DensePolynomial<R: Ring> {
    /// Aux information about the multilinear polynomial
    pub aux_info: DPAuxInfo<R>,
    /// Stores multilinear extensions in which product multiplicand can refer
    /// to.
    pub mles: Vec<RefCounter<DenseMultilinearExtension<R>>>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq, CanonicalSerialize)]
/// Auxiliary information about the multilinear polynomial
pub struct DPAuxInfo<R: Ring> {
    /// max number of multiplicands in each product
    pub max_degree: usize,
    /// number of variables of the polynomial
    pub num_variables: usize,
    /// Associated field
    #[doc(hidden)]
    pub phantom: PhantomData<R>,
}

// TODO: convert this into a trait
impl<R: Ring> DensePolynomial<R> {
    /// Creates an empty `DensePolynomial` with `num_variables`.
    pub fn new(num_variables: usize) -> Self {
        DensePolynomial {
            aux_info: DPAuxInfo {
                max_degree: 0,
                num_variables,
                phantom: PhantomData,
            },
            mles: Vec::new(),
        }
    }

    /// Creates an new DensePolynomial from a MLE.
    pub fn new_from_mle(mle: &RefCounter<DenseMultilinearExtension<R>>) -> Self {
        DensePolynomial {
            aux_info: DPAuxInfo {
                // The max degree is the max degree of any individual variable
                max_degree: 1,
                num_variables: mle.num_vars,
                phantom: PhantomData,
            },
            // here `0` points to the first polynomial of `mles`
            mles: vec![mle.clone()],
        }
    }

    /// Add MLEs, increasing the degree if the list is larger than the current degree.
    ///
    /// Returns an error if the list is empty, or the MLE has a different `num_vars` from self.
    pub fn add_mles(
        &mut self,
        mle_list: impl IntoIterator<Item = RefCounter<DenseMultilinearExtension<R>>>,
    ) -> Result<(), ArithErrors> {
        let mle_list: Vec<RefCounter<DenseMultilinearExtension<R>>> =
            mle_list.into_iter().collect();

        if mle_list.is_empty() {
            return Err(ArithErrors::InvalidParameters(
                "input mle_list is empty".to_string(),
            ));
        }

        self.aux_info.max_degree = max(self.aux_info.max_degree, mle_list.len());

        for mle in mle_list {
            if mle.num_vars != self.aux_info.num_variables {
                return Err(ArithErrors::InvalidParameters(format!(
                    "product has a multiplicand with wrong number of variables {} vs {}",
                    mle.num_vars, self.aux_info.num_variables
                )));
            }

            self.mles.push(mle);
        }

        Ok(())
    }

    /// Adds an MLE, incrementing also the degree.
    ///
    /// Returns an error if the MLE has a different `num_vars` from self.
    pub fn mul_mle(
        &mut self,
        mle: RefCounter<DenseMultilinearExtension<R>>,
    ) -> Result<(), ArithErrors> {
        let start = start_timer!(|| "mul by mle");

        if mle.num_vars != self.aux_info.num_variables {
            return Err(ArithErrors::InvalidParameters(format!(
                "product has a multiplicand with wrong number of variables {} vs {}",
                mle.num_vars, self.aux_info.num_variables
            )));
        }

        self.mles.push(mle);

        // increase the max degree by one as the MLE has degree 1.
        self.aux_info.max_degree += 1;
        end_timer!(start);
        Ok(())
    }

    /// Evaluate the `DensePolynomial` at point `point`.
    /// Returns an error is point.len() does not match `num_variables`.
    pub fn evaluate(&self, point: &[R]) -> Result<R, ArithErrors> {
        let start = start_timer!(|| "evaluation");

        if self.aux_info.num_variables != point.len() {
            return Err(ArithErrors::InvalidParameters(format!(
                "wrong number of variables {} vs {}",
                self.aux_info.num_variables,
                point.len()
            )));
        }

        let _evals: Vec<R> = self
            .mles
            .iter()
            .map(|x| {
                x.evaluate(point).unwrap() // safe unwrap here since we have
                                           // already checked that num_var
                                           // matches
            })
            .collect();
        let res = R::zero();
        //let res = self
        //    .products
        //    .iter()
        //    .map(|(c, p)| *c * p.iter().map(|&i| evals[i]).product::<R>())
        //    .sum();

        end_timer!(start);
        Ok(res)
    }

    /// Sample a random `DensePolynomial`, return the polynomial and its sum.
    pub fn rand<Rn: RngCore>(
        nv: usize,
        num_multiplicands_range: (usize, usize),
        num_products: usize,
        rng: &mut Rn,
    ) -> Result<(Self, Vec<(R, Vec<usize>)>, R), ArithErrors> {
        let start = start_timer!(|| "sample random dense polynomial");

        let mut sum = R::zero();
        let mut poly = DensePolynomial::new(nv);
        let mut products = Vec::with_capacity(num_products);
        let mut current_mle_index = 0;
        for _ in 0..num_products {
            let num_multiplicands =
                rng.gen_range(num_multiplicands_range.0..num_multiplicands_range.1);
            let (product, product_sum) = random_mle_list(nv, num_multiplicands, rng);

            let coefficient = R::rand(rng);
            poly.add_mles(product.into_iter())?;
            sum += product_sum * coefficient;

            let indices: Vec<usize> =
                (current_mle_index..current_mle_index + num_multiplicands).collect();
            products.push((coefficient, indices));
            current_mle_index += num_multiplicands;
        }

        end_timer!(start);
        Ok((poly, products, sum))
    }

    /// Sample a random dense polynomial that evaluates to zero everywhere
    /// over the boolean hypercube.
    pub fn rand_zero<Rn: RngCore>(
        nv: usize,
        num_multiplicands_range: (usize, usize),
        num_products: usize,
        rng: &mut Rn,
    ) -> Result<Self, ArithErrors> {
        let mut poly = DensePolynomial::new(nv);
        for _ in 0..num_products {
            let num_multiplicands =
                rng.gen_range(num_multiplicands_range.0..num_multiplicands_range.1);
            let product = random_zero_mle_list(nv, num_multiplicands, rng);
            poly.add_mles(product.into_iter())?;
        }

        Ok(poly)
    }

    // Input poly f(x) and a random vector r, output
    //      \hat f(x) = \sum_{x_i \in eval_x} f(x_i) eq(x, r)
    // where
    //      eq(x,y) = \prod_i=1^num_var (x_i * y_i + (1-x_i)*(1-y_i))
    //
    // This function is used in ZeroCheck.
    pub fn build_f_hat(&self, r: &[R]) -> Result<Self, ArithErrors> {
        let start = start_timer!(|| "zero check build hat f");

        if self.aux_info.num_variables != r.len() {
            return Err(ArithErrors::InvalidParameters(format!(
                "r.len() is different from number of variables: {} vs {}",
                r.len(),
                self.aux_info.num_variables
            )));
        }

        let eq_x_r = build_eq_x_r(r)?;
        let mut res = self.clone();
        res.mul_mle(eq_x_r)?;

        end_timer!(start);
        Ok(res)
    }
}

pub fn rand_poly<R: Ring>(
    nv: usize,
    num_multiplicands_range: (usize, usize),
    num_products: usize,
    rng: &mut impl RngCore,
) -> Result<
    (
        (Vec<RefCounter<DenseMultilinearExtension<R>>>, usize),
        Vec<(R, Vec<usize>)>,
        R,
    ),
    ArithErrors,
> {
    let mut sum = R::zero();
    let mut mles = vec![];
    let mut products = Vec::with_capacity(num_products);
    let mut degree = 0;
    let mut current_mle_index = 0;
    for _ in 0..num_products {
        let num_multiplicands = rng.gen_range(num_multiplicands_range.0..num_multiplicands_range.1);
        degree = num_multiplicands.max(degree);
        let (product, product_sum) = random_mle_list(nv, num_multiplicands, rng);

        let coefficient = R::rand(rng);
        mles.extend(product);
        sum += product_sum * coefficient;

        let indices: Vec<usize> =
            (current_mle_index..current_mle_index + num_multiplicands).collect();
        products.push((coefficient, indices));
        current_mle_index += num_multiplicands;
    }

    Ok(((mles, degree), products, sum))
}

impl<R: Ring> DPAuxInfo<R> {
    pub fn new(num_variables: usize, max_degree: usize) -> Self {
        DPAuxInfo {
            max_degree,
            num_variables,
            phantom: PhantomData,
        }
    }
}

/// Evaluate eq polynomial.
pub fn eq_eval<R: Ring>(x: &[R], y: &[R]) -> Result<R, ArithErrors> {
    if x.len() != y.len() {
        return Err(ArithErrors::InvalidParameters(
            "x and y have different length".to_string(),
        ));
    }
    let start = start_timer!(|| "eq_eval");
    let mut res = R::one();
    for (&xi, &yi) in x.iter().zip(y.iter()) {
        let xi_yi = xi * yi;
        res *= xi_yi + xi_yi - xi - yi + R::one();
    }
    end_timer!(start);
    Ok(res)
}

/// This function build the eq(x, r) polynomial for any given r.
///
/// Evaluate
///      eq(x,y) = \prod_i=1^num_var (x_i * y_i + (1-x_i)*(1-y_i))
/// over r, which is
///      eq(x,y) = \prod_i=1^num_var (x_i * r_i + (1-x_i)*(1-r_i))
pub fn build_eq_x_r<R: Ring>(
    r: &[R],
) -> Result<RefCounter<DenseMultilinearExtension<R>>, ArithErrors> {
    let evals = build_eq_x_r_vec(r)?;
    let mle = DenseMultilinearExtension::from_evaluations_vec(r.len(), evals);

    Ok(RefCounter::new(mle))
}
/// This function build the eq(x, r) polynomial for any given r, and output the
/// evaluation of eq(x, r) in its vector form.
///
/// Evaluate
///      eq(x,y) = \prod_i=1^num_var (x_i * y_i + (1-x_i)*(1-y_i))
/// over r, which is
///      eq(x,y) = \prod_i=1^num_var (x_i * r_i + (1-x_i)*(1-r_i))
pub fn build_eq_x_r_vec<R: Ring>(r: &[R]) -> Result<Vec<R>, ArithErrors> {
    // we build eq(x,r) from its evaluations
    // we want to evaluate eq(x,r) over x \in {0, 1}^num_vars
    // for example, with num_vars = 4, x is a binary vector of 4, then
    //  0 0 0 0 -> (1-r0)   * (1-r1)    * (1-r2)    * (1-r3)
    //  1 0 0 0 -> r0       * (1-r1)    * (1-r2)    * (1-r3)
    //  0 1 0 0 -> (1-r0)   * r1        * (1-r2)    * (1-r3)
    //  1 1 0 0 -> r0       * r1        * (1-r2)    * (1-r3)
    //  ....
    //  1 1 1 1 -> r0       * r1        * r2        * r3
    // we will need 2^num_var evaluations

    let mut eval = Vec::new();
    build_eq_x_r_helper(r, &mut eval)?;

    Ok(eval)
}

/// A helper function to build eq(x, r) recursively.
/// This function takes `r.len()` steps, and for each step it requires a maximum
/// `r.len()-1` multiplications.
fn build_eq_x_r_helper<R: Ring>(r: &[R], buf: &mut Vec<R>) -> Result<(), ArithErrors> {
    if r.is_empty() {
        return Err(ArithErrors::InvalidParameters("r length is 0".to_string()));
    } else if r.len() == 1 {
        // initializing the buffer with [1-r_0, r_0]
        buf.push(R::one() - r[0]);
        buf.push(r[0]);
    } else {
        build_eq_x_r_helper(&r[1..], buf)?;

        // suppose at the previous step we received [b_1, ..., b_k]
        // for the current step we will need
        // if x_0 = 0:   (1-r0) * [b_1, ..., b_k]
        // if x_0 = 1:   r0 * [b_1, ..., b_k]
        // let mut res = vec![];
        // for &b_i in buf.iter() {
        //     let tmp = r[0] * b_i;
        //     res.push(b_i - tmp);
        //     res.push(tmp);
        // }
        // *buf = res;

        let mut res = vec![R::zero(); buf.len() << 1];
        cfg_iter_mut!(res).enumerate().for_each(|(i, val)| {
            let bi = buf[i >> 1];
            let tmp = r[0] * bi;
            if (i & 1) == 0 {
                *val = bi - tmp;
            } else {
                *val = tmp;
            }
        });
        *buf = res;
    }

    Ok(())
}

/// Decompose an integer into a binary vector in little endian.
#[cfg(feature = "std")]
pub fn bit_decompose(input: u64, num_var: usize) -> Vec<bool> {
    let mut res = Vec::with_capacity(num_var);
    let mut i = input;
    for _ in 0..num_var {
        res.push((i & 1) == 1);
        i >>= 1;
    }
    res
}